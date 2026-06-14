#!/usr/bin/env python3
"""Qutrit strange-state / QR distillation discovery driver (d > 2).

Runs the strange-state sub-arm (:mod:`qudit_qec.distill_strange`): enumerates
self-orthogonal cyclic codes over F_3, screens them with the cheap Fraction-based
strange-distillation conditions, catalogs the valid distillers, then refines the
exact nu + threshold for the top entries.  Validates by reproducing the 11-qutrit
Golay [[11,1,5]]_3 (cubic nu=3, threshold ~0.387).  Writes
``results/distill_strange_catalog.json`` + ``results/distill_strange_report.md``.

Compute safety (memory rule + docs/08): single process; the weight-enumerator screen
is hard-gated to ``LOCAL_ENUM_DEFAULT`` (~Golay-sized 3^10) so n=13 self-orthogonal
codes that need 3^12 are refused locally and routed to lenore; the sympy threshold
solve runs only for the top ``--refine`` entries.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from qudit_qec import distill_strange as ds


def _record(r: ds.StrangeResult) -> dict:
    return {"p": r.p, "n": r.n, "k": r.k, "nu": r.nu, "distills": r.distills,
            "threshold": (round(r.threshold, 4) if r.threshold is not None else None),
            "min_stabilizer_weight": r.d_lb, "sound": r.sound}


def write_report(path: str, cat: ds.StrangeCatalog, args, golay_ok: bool) -> None:
    ranked = cat.best_by_nu()
    lines = [
        "# Qutrit strange-state / QR distillation discovery — report",
        "",
        "Sub-arm: `qudit_qec.distill_strange` · genotype = self-orthogonal **cyclic "
        "codes over F_3** · CSS distiller `Hx = Hz = C` of the qutrit **strange state**.",
        f"Search: n={args.n_min}..{args.n_max} (gcd(n,3)=1), local enum cap "
        f"{ds.LOCAL_ENUM_DEFAULT}, threshold-refine top {args.refine}.  Objective axis: "
        f"noise-suppression exponent `nu` + threshold `eps_*` (NOT gamma).",
        "",
        f"## Validation: 11-qutrit Golay `[[11,1,5]]_3` reproduced: "
        f"{'✓ (nu=3 cubic, threshold ~0.387)' if golay_ok else '✗ FAILED'}",
        "",
        f"## Catalog: {len(cat)} genuine strange distillers",
        "",
        "| [[n,k]]_3 | nu | threshold eps_* | min stab. weight |",
        "|---|---|---|---|",
    ]
    for r in ranked:
        thr = f"{r.threshold:.3f}" if r.threshold is not None else "— (refine)"
        lines.append(f"| [[{r.n},{r.k}]]_3 | {r.nu} | {thr} | {r.d_lb} |")
    lines += [
        "",
        "## Notes (honest framing)",
        "- `nu` is read from the weight enumerator (lowest power of eps in `3A+B`), "
        "distinct from the code distance `d` (the Golay has `d=5` but `nu=3`); the "
        "`min stab. weight` column is the smallest nonzero stabilizer weight, a "
        "distance hint, not the logical distance.",
        "- **The small-n self-orthogonal *cyclic* F_3 family is sparse** — no distiller "
        "at n=5,7,13(local); the n=11 entries are the two equivalent ternary QR Golay "
        "codes. Broadening to non-cyclic self-orthogonal codes and n>13 is the "
        "**lenore-scale** extension (3^(n-k) exceeds the local enum cap).",
        "- This sub-arm beats the d=2 triorthogonal family on noise suppression "
        "(`nu=3` cubic vs `nu=2`); it is qutrit-only (the strange state is "
        "GF(3)-specific), complementary to the general-prime gamma arm "
        "(`distill_discovery`, docs/08).",
        "",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Qutrit strange/QR distillation driver.")
    ap.add_argument("--n-min", type=int, default=5)
    ap.add_argument("--n-max", type=int, default=13)
    ap.add_argument("--refine", type=int, default=5, help="threshold-refine top-N by nu")
    ap.add_argument("--out-dir", default=os.path.join(_REPO_ROOT, "results"))
    args = ap.parse_args(argv)

    cat = ds.search_strange_cyclic(args.n_min, args.n_max)
    print(f"search n={args.n_min}..{args.n_max}: {len(cat)} strange distillers")
    cat.refine_thresholds(top=args.refine)

    golay = [r for r in cat.codes() if (r.n, r.k) == (11, 1)]
    golay_ok = bool(golay and golay[0].nu == 3 and golay[0].threshold
                    and abs(golay[0].threshold - 0.387) < 0.02)
    print(f"Golay [[11,1,5]]_3 reproduced: {'OK' if golay_ok else 'FAIL'}")
    for r in cat.best_by_nu():
        thr = f"{r.threshold:.3f}" if r.threshold is not None else "(unrefined)"
        print(f"  [[{r.n},{r.k}]]_3 nu={r.nu} threshold={thr}")

    os.makedirs(args.out_dir, exist_ok=True)
    catalog_path = os.path.join(args.out_dir, "distill_strange_catalog.json")
    report_path = os.path.join(args.out_dir, "distill_strange_report.md")
    payload = {
        "config": {"n_min": args.n_min, "n_max": args.n_max, "refine": args.refine,
                   "local_enum_cap": ds.LOCAL_ENUM_DEFAULT},
        "golay_reproduced": golay_ok,
        "codes": [_record(r) for r in cat.best_by_nu()],
    }
    with open(catalog_path, "w") as fh:
        json.dump(payload, fh, indent=2)
    write_report(report_path, cat, args, golay_ok)
    print(f"wrote {catalog_path}\nwrote {report_path}")
    return 0 if golay_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
