"""OpenEvolve evaluator adapter for qudit BB code discovery.

Exposes ``evaluate(program_path) -> dict`` -- the fitness function OpenEvolve calls
on each evolved program. It loads the program's ``generate_candidates``, evaluates
its candidates with the trusted qudit pipeline (k-only screen, then prime-q MILP on
the most promising per lattice), and returns a metrics dict whose ``combined_score``
is the fitness and whose ``lattices_with_high_k`` / ``num_high_k`` are the MAP-Elites
behavioral dimensions.

The field order ``q`` comes from the ``QCODE_FIELD`` environment variable (set by
``run_evolution.py``) or the ``field=`` argument. This module depends only on the
``qudit_qec`` pipeline, **not** on ``openevolve`` -- so it is importable and testable
standalone (point ``evaluate`` at any seed file).

Design note (trust): the fitness rewards encoded dimension ``k`` (exact, never
corrupted) and breadth across lattices, plus a bonus for the best **trusted**
(MILP-certified) FOM. The loose GUF bound is never used as a fitness signal -- the
scoping's central lesson (``docs/03``).
"""

from __future__ import annotations

import importlib.util
import math
import os

# Absolute imports (not relative): OpenEvolve loads this file standalone as its
# evaluation_file, so ``from ..evaluator`` would fail -- it has no package context.
# run_evolution.py puts the repo root on sys.path so ``qudit_qec`` is importable.
from qudit_qec.evaluator import evaluate_candidate
from qudit_qec.results import CodeCatalog

# "high-k" threshold for the MAP-Elites feature dimension. NOTE: 8 is the F_2
# convention; the right value over GF(q) should be recalibrated per field (the
# gross-code analogue is k=8 over GF(3), not 12) -- an open campaign-tuning item.
HIGH_K = int(os.environ.get("QCODE_HIGH_K", "8"))

DEFAULT_STAGE2_LATTICES = [(6, 6), (9, 6), (12, 6), (6, 12), (12, 12), (15, 6)]


def _load_generate_candidates(program_path: str):
    spec = importlib.util.spec_from_file_location("evolved_qudit_program", program_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "generate_candidates"):
        raise AttributeError("evolved program is missing generate_candidates(ell, m)")
    return module.generate_candidates


def _field(explicit: int | None) -> int:
    return int(explicit) if explicit is not None else int(os.environ.get("QCODE_FIELD", "2"))


def _error_result(message: str) -> dict:
    return {"combined_score": 0.0, "best_fom": 0.0, "num_valid": 0.0,
            "num_high_k": 0.0, "lattices_with_high_k": 0.0, "total_candidates": 0.0,
            "num_trusted": 0.0, "num_errors": 1.0, "error": message[:500]}


def _pick_top(results: list, limit: int) -> list:
    """Diversify by k value first, then by A polynomial, then fill by k."""
    ranked = sorted(results, key=lambda r: r.k, reverse=True)
    top: list = []
    seen_k: set[int] = set()
    for r in ranked:
        if len(top) >= limit:
            break
        if r.k not in seen_k:
            seen_k.add(r.k)
            top.append(r)
    seen_a = {r.A for r in top}
    for r in ranked:
        if len(top) >= limit:
            break
        if r.A not in seen_a:
            seen_a.add(r.A)
            top.append(r)
    for r in ranked:
        if len(top) >= limit:
            break
        if r not in top:
            top.append(r)
    return top


def _code_repr(code) -> dict | None:
    if code is None:
        return None
    return {"n": code.n, "k": code.k, "d": code.d, "fom": code.fom,
            "trusted": code.trusted, "ell": code.ell, "m": code.m, "q": code.q,
            "A": list(code.A), "B": list(code.B), "decomposable": code.decomposable}


def evaluate(
    program_path: str,
    *,
    field: int | None = None,
    stage2_lattices: list | None = None,
    max_distance_per_lattice: int = 4,
    max_candidates_per_lattice: int = 2000,
    milp_total_timeout: float = 60,
    milp_timeout_per_logical: float = 15,
) -> dict:
    """Evaluate an evolved program; return OpenEvolve metrics incl. ``combined_score``."""
    q = _field(field)
    # Keep the seed's coefficient generation (which reads FIELD = QCODE_FIELD at
    # load time) consistent with the field we evaluate over.
    os.environ["QCODE_FIELD"] = str(q)
    try:
        generate = _load_generate_candidates(program_path)
    except Exception as exc:  # noqa: BLE001
        return _error_result(f"load failed: {exc}")

    lattices = stage2_lattices or DEFAULT_STAGE2_LATTICES
    catalog = CodeCatalog()
    errors: list[str] = []
    total = 0

    for ell, m in lattices:
        try:
            candidates = generate(ell, m)
        except Exception as exc:  # noqa: BLE001
            errors.append(f"({ell},{m}) generate: {type(exc).__name__}: {exc}")
            continue
        if not isinstance(candidates, list):
            errors.append(f"({ell},{m}): generate returned {type(candidates).__name__}, not list")
            continue
        if len(candidates) > max_candidates_per_lattice:
            candidates = candidates[:max_candidates_per_lattice]
        total += len(candidates)

        screened = []
        for pair in candidates:
            try:
                a_terms, b_terms = pair
                r = evaluate_candidate(ell, m, a_terms, b_terms, field=q, distance="none")
            except Exception as exc:  # noqa: BLE001
                errors.append(f"({ell},{m}) eval: {type(exc).__name__}: {exc}")
                continue
            if not r.rejected and r.k > 0:
                screened.append(r)
                catalog.add(r)

        for r in _pick_top(screened, max_distance_per_lattice):
            try:
                full = evaluate_candidate(
                    ell, m, r.A, r.B, field=q, distance="trusted",
                    milp_total_timeout=milp_total_timeout,
                    milp_timeout_per_logical=milp_timeout_per_logical,
                )
                catalog.add(full)
            except Exception as exc:  # noqa: BLE001
                errors.append(f"({ell},{m}) milp: {type(exc).__name__}: {exc}")

    return _aggregate(catalog, q, total, errors)


def _aggregate(catalog: CodeCatalog, q: int, total: int, errors: list[str]) -> dict:
    codes = catalog.codes()
    valid = [c for c in codes if c.k and c.k > 0]
    high_k = [c for c in valid if c.k >= HIGH_K]
    lattices_with_high_k = len({(c.ell, c.m) for c in high_k})
    trusted_foms = [c.fom for c in codes if c.fom and c.trusted]
    best_fom = max(trusted_foms) if trusted_foms else 0.0
    best_rate = max((c.k / c.n for c in valid), default=0.0)
    high_k_quality = math.log1p(len(high_k)) / 10.0

    # Fitness: alive base + encoding rate + high-k breadth + trusted-FOM bonus.
    combined = 0.1 + best_rate + high_k_quality + 0.5 * min(best_fom / 12.0, 2.0)

    best_code = max(codes, key=lambda c: (c.fom or 0.0)) if codes else None
    return {
        "combined_score": float(combined),
        "best_fom": float(best_fom),
        "num_valid": float(len(valid)),
        "num_high_k": float(len(high_k)),
        "lattices_with_high_k": float(lattices_with_high_k),
        "num_trusted": float(len(trusted_foms)),
        "total_candidates": float(total),
        "num_errors": float(len(errors)),
        "field": float(q),
        "best_code": _code_repr(best_code),
    }
