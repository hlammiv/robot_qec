#!/usr/bin/env python3
"""Prime-p magic-state-distillation discovery driver.

Runs the triorthogonal-matrix discovery arm (:mod:`qudit_qec.distill_discovery`):
seeds the validated block-triorthogonal family ``T(p,m,k) -> [[p^2 m - k, k, 2]]_p``,
applies re-validated mutation operators, and catalogs genuine distillation codes by
yield ``gamma = log_d(n/k)``.  First it **validates** the pipeline against the known
arXiv:2403.06228 codes; then it runs a bounded local search and writes
``results/distill_catalog.json`` + ``results/distill_report.md``.

Compute safety (memory rule + docs/08): single process, O(kappa^3) validity gate, NO
distance MILP (distances carried through distance-preserving operators).  Default caps
are deliberately small for the local box; raise ``--iters``/``--max-n`` only on lenore.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from qudit_qec import distill_discovery as dd

# Known codes the arm must reproduce before any search is trusted.
_KNOWN = [
    ("[[20,7,2]]_3", 3, 3, 7, 1.515),
    ("[[14,4,2]]_3", 3, 2, 4, 1.807),
    ("[[13,5,2]]_3", 3, 2, 5, 1.379),
    ("[[23,2,2]]_5", 5, 1, 2, 3.524),
]


def validate_known() -> list[dict]:
    rows = []
    for label, p, m, k, exp in _KNOWN:
        H = dd.triortho_family(p, m, k)
        r = dd.evaluate_distill_candidate(H, p, d=2, d_status="known", op="family")
        ok = bool(r.distills and abs(r.gamma - exp) < 0.01)
        rows.append({"label": label, "n": r.n, "k": r.k, "d": r.d,
                     "gamma": round(r.gamma, 3), "expected_gamma": exp, "ok": ok})
    return rows


def _record(r: dd.DistillResult) -> dict:
    return {"p": r.p, "n": r.n, "k": r.k, "d": r.d, "d_status": r.d_status,
            "gamma": (round(r.gamma, 4) if r.gamma is not None else None),
            "distills": r.distills, "trusted": r.trusted, "sound": r.sound, "op": r.op}


def write_report(path: str, known, cat: dd.DistillCatalog, args) -> None:
    trusted = [r for r in cat.codes() if r.trusted and r.gamma is not None]
    pending = [r for r in cat.codes() if not r.trusted and r.gamma is not None]
    lines = [
        "# Prime-p magic-state-distillation discovery — report",
        "",
        f"Arm: `qudit_qec.distill_discovery` · family `T(p,m,k) -> [[p^2 m - k, k, 2]]_p`.",
        f"Search: primes={args.primes}, m_range={args.m_range}, iters={args.iters}, "
        f"seed={args.seed}, max_n={args.max_n}.  Validity gate: exact cubic triorthogonality. "
        f"Objective: yield gamma=log_d(n/k) (lower better).  No distance MILP (compute-safe).",
        "",
        "## Validation vs known codes (arXiv:2403.06228 / 2408.00436)",
        "",
        "| code | n | k | d | gamma | expected | ok |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in known:
        lines.append(f"| {row['label']} | {row['n']} | {row['k']} | {row['d']} | "
                     f"{row['gamma']} | {row['expected_gamma']} | {'✓' if row['ok'] else '✗'} |")
    # dedup the DISPLAY by (p,n,k,d): distinct genotypes can share parameters.
    seen: dict = {}
    for r in sorted(trusted, key=lambda x: (x.gamma, x.n)):
        seen.setdefault((r.p, r.n, r.k, r.d), r)
    lines += [
        "",
        f"## Catalog: {len(cat)} genuine distillation genotypes "
        f"({len(trusted)} trusted-distance, {len(pending)} distance-pending; "
        f"{len(seen)} distinct [[n,k,d]] parameter sets)",
        "",
        "### Trusted frontier (best by yield gamma; distinct parameter sets)",
        "",
        "| [[n,k,d]]_p | gamma | first-seen op |",
        "|---|---|---|",
    ]
    for r in sorted(seen.values(), key=lambda x: (x.gamma, x.n))[:12]:
        lines.append(f"| [[{r.n},{r.k},{r.d}]]_{r.p} | {r.gamma:.3f} | {r.op} |")
    if pending:
        lines += [
            "",
            "### Distance-pending candidates (gamma is an optimistic lower bound — "
            "verify distance on lenore before trusting)",
            "",
            "| [[n,k,≤d]]_p | gamma(≤) | op |",
            "|---|---|---|",
        ]
        for r in sorted(pending, key=lambda x: (x.gamma, x.n))[:8]:
            lines.append(f"| [[{r.n},{r.k},≤{r.d}]]_{r.p} | {r.gamma:.3f} | {r.op} |")
    lines += [
        "",
        "## Pareto front (n smaller, k larger, d larger; known-distance codes)",
        "",
        "| [[n,k,d]]_p | gamma | op |",
        "|---|---|---|",
    ]
    for r in sorted(cat.pareto_front(), key=lambda x: (x.n, -x.k)):
        g = f"{r.gamma:.3f}" if r.gamma is not None else "—"
        lines.append(f"| [[{r.n},{r.k},{r.d}]]_{r.p} | {g} | {r.op} |")
    fam = [r for r in trusted if r.op == "family"]
    best_fam = min(fam, key=lambda r: r.gamma) if fam else None
    lines += [
        "",
        "## Notes (honest framing)",
        "- Every catalog entry passes the exact cubic triorthogonality gate "
        "(`is_triorthogonal`): nothing is triorthogonal-by-assumption.",
        "- `trusted` distances are carried through distance-preserving operators "
        "(family d=2, column permutation/scaling, direct sum); puncture yields "
        "`d`-upper candidates whose gamma is optimistic.",
        "- **`direct_sum` codes are independent block-stacks of smaller members** — "
        "they have the *same* per-qudit yield gamma as their parts (e.g. `[[26,10,2]]` "
        "= two `[[13,5,2]]`), so they are not genuinely better distillers. The best "
        "per-qudit family member in THIS run is"
        + (f" `[[{best_fam.n},{best_fam.k},{best_fam.d}]]_{best_fam.p}` (gamma={best_fam.gamma:.3f})"
           if best_fam else " (none)")
        + " — but this is a function of the search window: the family's gamma decreases "
        "monotonically toward 1 as `m` grows (optimal `k=pm-1`), so there is no finite "
        "d=2 frontier, only the asymptote.",
        "- **Beating the d=2 regime on *noise suppression* (nu=d=2), not just overhead, "
        "needs distance `d > 2`**: the punctured-Reed–Muller route (`reed_muller_triortho`, "
        "validated on the qubit `[[15,1,3]]`) and the qutrit strange/QR route "
        "(`distill_strange`, the Golay) — see `docs/08`.",
        "- Universality of the transversal T is the separate per-dimension question of "
        "`docs/08` (prime p: any non-Clifford gate is universal — sound here).",
        "",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Prime-p MSD discovery driver.")
    ap.add_argument("--primes", type=int, nargs="+", default=[3])
    ap.add_argument("--m-min", type=int, default=1)
    ap.add_argument("--m-max", type=int, default=2)
    ap.add_argument("--iters", type=int, default=400)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-n", type=int, default=40)
    ap.add_argument("--out-dir", default=os.path.join(_REPO_ROOT, "results"))
    args = ap.parse_args(argv)
    args.primes = tuple(args.primes)
    args.m_range = (args.m_min, args.m_max)

    known = validate_known()
    failed = [r for r in known if not r["ok"]]
    print("validation vs known codes:")
    for row in known:
        print(f"  {row['label']}: gamma={row['gamma']} (exp {row['expected_gamma']}) "
              f"{'OK' if row['ok'] else 'FAIL'}")
    if failed:
        print(f"ABORT: {len(failed)} known-code validation(s) failed; not running search.")
        return 1

    cat = dd.search_distill(primes=args.primes, m_range=args.m_range,
                            iters=args.iters, seed=args.seed, max_n=args.max_n)
    print(f"\nsearch: {len(cat)} genuine distillation codes catalogued")

    os.makedirs(args.out_dir, exist_ok=True)
    catalog_path = os.path.join(args.out_dir, "distill_catalog.json")
    report_path = os.path.join(args.out_dir, "distill_report.md")
    payload = {
        "config": {"primes": list(args.primes), "m_range": list(args.m_range),
                   "iters": args.iters, "seed": args.seed, "max_n": args.max_n},
        "validation": known,
        "codes": sorted((_record(r) for r in cat.codes()),
                        key=lambda d: (d["gamma"] if d["gamma"] is not None else math.inf, d["n"])),
        "pareto_front": [_record(r) for r in cat.pareto_front()],
    }
    with open(catalog_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    write_report(report_path, known, cat, args)
    print(f"wrote {catalog_path}\nwrote {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
