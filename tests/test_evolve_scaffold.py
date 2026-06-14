"""Tests for the Phase 4 evolution scaffolding (seed, adapter, launcher).

These exercise the OpenEvolve integration contract WITHOUT openevolve installed:
the seed generates valid coefficient-bearing candidates, the adapter's
``evaluate(program_path)`` returns a coherent metrics dict (with the MAP-Elites
feature dimensions) using the trusted pipeline, and the launcher parses ``--field``
and degrades gracefully when openevolve is absent.
"""

import importlib.util
import os

import pytest

SEED = os.path.join(os.path.dirname(__file__), "..", "qudit_qec", "evolve",
                    "seed_solution_qudit.py")
SEED = os.path.abspath(SEED)


def _load_seed_with_field(q):
    os.environ["QCODE_FIELD"] = str(q)
    spec = importlib.util.spec_from_file_location(f"seed_q{q}", SEED)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------- #
# seed program
# --------------------------------------------------------------------------- #
def test_seed_generates_triples_over_gf3():
    mod = _load_seed_with_field(3)
    assert mod.FIELD == 3
    cands = mod.generate_candidates(6, 6)
    assert isinstance(cands, list) and len(cands) > 0
    A, B = cands[0]
    assert all(len(t) == 3 for t in A) and all(len(t) == 3 for t in B)
    # coefficients are valid GF(3)* values everywhere
    for a_terms, b_terms in cands:
        for _, _, c in a_terms + b_terms:
            assert 1 <= c <= 2


def test_seed_explores_coefficient_axis_only_for_q_gt_2():
    # over GF(3) the coefficient strategy emits a non-unit coefficient somewhere
    cands3 = _load_seed_with_field(3).generate_candidates(6, 6)
    assert any(c != 1 for a, b in cands3 for _, _, c in a + b)
    # over GF(2) every coefficient is 1 (no coefficient axis)
    cands2 = _load_seed_with_field(2).generate_candidates(6, 6)
    assert all(c == 1 for a, b in cands2 for _, _, c in a + b)


# --------------------------------------------------------------------------- #
# evaluator adapter
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_adapter_evaluate_returns_feature_dims_and_score():
    from qudit_qec.evolve import adapter

    metrics = adapter.evaluate(
        SEED, field=3,
        stage2_lattices=[(6, 6)],          # one small lattice keeps it fast
        max_distance_per_lattice=2,
        milp_total_timeout=30,
    )
    # OpenEvolve contract: combined_score (fitness) + the two MAP-Elites feature dims
    assert "combined_score" in metrics
    assert "lattices_with_high_k" in metrics and "num_high_k" in metrics
    assert metrics["combined_score"] > 0.1          # seed produces valid codes
    assert metrics["num_valid"] >= 1
    assert metrics["field"] == 3.0
    assert metrics["best_code"] is not None


def test_adapter_handles_broken_program(tmp_path):
    from qudit_qec.evolve import adapter

    bad = tmp_path / "bad.py"
    bad.write_text("x = 1\n")  # no generate_candidates
    metrics = adapter.evaluate(str(bad), field=3)
    assert metrics["combined_score"] == 0.0 and "error" in metrics


# --------------------------------------------------------------------------- #
# launcher CLI
# --------------------------------------------------------------------------- #
def test_run_evolution_parser_and_env():
    from qudit_qec.evolve import run_evolution

    args = run_evolution.build_parser().parse_args(["--field", "3", "--iterations", "5"])
    assert args.field == 3 and args.iterations == 5
    run_name = run_evolution.configure_environment(args)
    assert os.environ["QCODE_FIELD"] == "3"
    assert run_name == "qudit_gf3"


def test_run_evolution_rejects_non_prime_field():
    from qudit_qec.evolve import run_evolution

    with pytest.raises(SystemExit):
        run_evolution.main(["--field", "6"])  # composite -> rejected before openevolve


def test_run_evolution_graceful_without_openevolve():
    from qudit_qec.evolve import run_evolution

    # openevolve is not installed in this environment -> main returns 2 (not a crash)
    rc = run_evolution.main(["--field", "3", "--iterations", "1"])
    assert rc == 2
