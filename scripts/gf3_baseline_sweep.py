#!/usr/bin/env python3
"""GF(q) baseline sweep: catalog trusted [[n,k,d]]_q bivariate-bicycle codes.

For each requested lattice ``(ell, m)`` this:
  1. generates candidates with the qudit seed's ``generate_candidates`` (over GF(q),
     so the coefficient axis is exercised when q > 2),
  2. k-screens them exactly (GF(q) rank),
  3. runs the **trusted prime-q MILP** distance on the most promising (diversified
     by k and by A polynomial), and
  4. catalogs the results (coefficient-aware dedup, decomposability, FOM).

Output is a JSON shard listing every distance-known code plus the per-lattice and
overall Pareto front. Several invocations over disjoint lattice sets can run in
parallel (one OS process each) and be merged with ``--merge``.

Usage:
  python scripts/gf3_baseline_sweep.py --field 3 --lattices "6,6 9,6" --out shard.json
  python scripts/gf3_baseline_sweep.py --merge a.json b.json --out catalog.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

SEED_PATH = os.path.join(_REPO_ROOT, "qudit_qec", "evolve", "seed_solution_qudit.py")


def parse_lattices(spec: str) -> list[tuple[int, int]]:
    """Parse ``"6,6 9,6 12,6"`` into ``[(6, 6), (9, 6), (12, 6)]``."""
    out = []
    for pair in spec.split():
        ell, m = pair.split(",")
        out.append((int(ell), int(m)))
    return out


def _code_record(r) -> dict:
    return {
        "ell": r.ell, "m": r.m, "q": r.q, "n": r.n, "k": r.k, "d": r.d,
        "d_status": r.d_status, "fom": r.fom, "trusted": r.trusted,
        "decomposable": r.decomposable, "self_dual": r.self_dual,
        "A": [list(t) for t in r.A], "B": [list(t) for t in r.B],
    }


def sweep(lattices, field, max_distance_per_lattice, milp_total_timeout,
          milp_timeout_per_logical, max_candidates) -> dict:
    os.environ["QCODE_FIELD"] = str(field)
    from qudit_qec.evaluator import evaluate_candidate
    from qudit_qec.evolve.adapter import _load_generate_candidates, _pick_top
    from qudit_qec.results import CodeCatalog

    generate = _load_generate_candidates(SEED_PATH)
    catalog = CodeCatalog()
    per_lattice = {}
    errors = []
    t0 = time.monotonic()

    for ell, m in lattices:
        lt0 = time.monotonic()
        try:
            cands = generate(ell, m)
        except Exception as exc:  # noqa: BLE001
            errors.append(f"({ell},{m}) generate: {type(exc).__name__}: {exc}")
            continue
        cands = cands[:max_candidates]
        screened = []
        for a_terms, b_terms in cands:
            try:
                r = evaluate_candidate(ell, m, a_terms, b_terms, field=field, distance="none")
            except Exception as exc:  # noqa: BLE001
                errors.append(f"({ell},{m}) screen: {type(exc).__name__}: {exc}")
                continue
            if not r.rejected and r.k > 0:
                screened.append(r)
                catalog.add(r)  # k-only record (no distance yet)

        # Prefer INDECOMPOSABLE codes for the (costly) distance budget -- direct
        # sums offer no error-correction advantage. Fill remaining slots with
        # decomposable codes only if budget is left.
        indecomp = [r for r in screened if r.decomposable is False]
        decomp = [r for r in screened if r.decomposable]
        selection = _pick_top(indecomp, max_distance_per_lattice)
        if len(selection) < max_distance_per_lattice:
            selection += _pick_top(decomp, max_distance_per_lattice - len(selection))

        n_dist = 0
        for r in selection:
            try:
                full = evaluate_candidate(
                    ell, m, r.A, r.B, field=field, distance="trusted",
                    milp_total_timeout=milp_total_timeout,
                    milp_timeout_per_logical=milp_timeout_per_logical,
                )
                catalog.add(full)
                n_dist += 1
            except Exception as exc:  # noqa: BLE001
                errors.append(f"({ell},{m}) milp: {type(exc).__name__}: {exc}")
        per_lattice[f"{ell},{m}"] = {
            "n": 2 * ell * m, "candidates": len(cands), "valid": len(screened),
            "distance_computed": n_dist, "time_s": round(time.monotonic() - lt0, 1),
        }

    codes = [_code_record(c) for c in catalog.with_distance()]
    pareto = [_code_record(c) for c in catalog.pareto_front()]
    best = catalog.best_by_fom()
    return {
        "field": field,
        "lattices": [f"{e},{m}" for e, m in lattices],
        "num_distinct_codes": len(catalog),
        "num_with_distance": len(codes),
        "per_lattice": per_lattice,
        "codes": codes,
        "pareto_front": pareto,
        "best_by_fom": [_code_record(c) for c in best[:20]],
        "errors": errors,
        "elapsed_s": round(time.monotonic() - t0, 1),
    }


def merge(paths: list[str]) -> dict:
    """Merge shards: union codes (coeff-aware key), recompute Pareto + best-by-FOM."""
    from qudit_qec.evaluator import EvalResult
    from qudit_qec.genotype import pair_key
    from qudit_qec.results import CodeCatalog

    catalog = CodeCatalog()
    fields, all_errors, per_lattice = set(), [], {}
    for path in paths:
        shard = json.load(open(path))
        fields.add(shard.get("field"))
        all_errors += shard.get("errors", [])
        per_lattice.update(shard.get("per_lattice", {}))
        for rec in shard.get("codes", []):
            A = tuple(tuple(t) for t in rec["A"])
            B = tuple(tuple(t) for t in rec["B"])
            key = pair_key(A, B, rec["ell"], rec["m"], rec["q"])
            res = EvalResult(rec["ell"], rec["m"], rec["q"], rec["n"], rec["k"],
                             rec["d"], rec["d_status"], rec["fom"], rec["trusted"],
                             rec["decomposable"], rec["self_dual"], A, B, key,
                             False, "", {})
            catalog.add(res)
    codes = [_code_record(c) for c in catalog.with_distance()]
    return {
        "field": sorted(f for f in fields if f is not None),
        "num_distinct_codes": len(catalog),
        "num_with_distance": len(codes),
        "per_lattice": per_lattice,
        "codes": codes,
        "pareto_front": [_code_record(c) for c in catalog.pareto_front()],
        "best_by_fom": [_code_record(c) for c in catalog.best_by_fom()[:30]],
        "errors": all_errors,
    }


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="GF(q) baseline BB-code sweep.")
    p.add_argument("--field", type=int, default=3)
    p.add_argument("--lattices", default="6,6", help='e.g. "6,6 9,6 12,6"')
    p.add_argument("--max-distance-per-lattice", type=int, default=8)
    p.add_argument("--milp-total-timeout", type=float, default=90)
    p.add_argument("--milp-timeout-per-logical", type=float, default=15)
    p.add_argument("--max-candidates", type=int, default=4000)
    p.add_argument("--out", required=True)
    p.add_argument("--merge", nargs="*", default=None, help="merge these JSON shards instead of sweeping")
    args = p.parse_args(argv)

    if args.merge:
        result = merge(args.merge)
    else:
        result = sweep(parse_lattices(args.lattices), args.field,
                       args.max_distance_per_lattice, args.milp_total_timeout,
                       args.milp_timeout_per_logical, args.max_candidates)
    with open(args.out, "w") as fh:
        json.dump(result, fh, indent=2)
    nd = result.get("num_with_distance", 0)
    print(f"wrote {args.out}: {result.get('num_distinct_codes', '?')} distinct codes, "
          f"{nd} with distance; errors={len(result.get('errors', []))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
