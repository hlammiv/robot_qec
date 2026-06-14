#!/usr/bin/env python3
"""Shared discovery-evaluation driver: score candidate BB codes over GF(q).

Reads a JSON list of candidates ``[{"ell":.., "m":.., "A":[[x,y,c]...], "B":[...]}, ...]``
from a file (``--in``) or stdin, evaluates each over GF(``--field``), and writes a JSON
list of results (best available distance + certification status). Used by both the
prime-q and Galois-GF(4) discovery workflows.

Distance / certification, by field type:
  * PRIME q (3,5,7,...):  the full trusted pipeline -- code_distance(certify=True),
    i.e. prime-q MILP + QDistRnd second source + weight-cut certification.
  * PRIME-POWER q (4,8,9,...):  the integer-mod-q MILP is INVALID, so we use the
    QDistRnd both-sector upper bound, plus -- for small codes (q^k <= --exact-cap) --
    qldpc's field-generic exact enumeration run in an ISOLATED `timeout` subprocess
    (fork-free, sidestepping the OpenMP-after-fork abort). Agreement => certified.

Output per code: ell,m,q,n,k,d,certified,method,fom,decomposable,A,B.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)


def _exact_via_subprocess(ell, m, A, B, q, timeout):
    """Run qldpc get_distance_exact in a fresh OS process (no fork) under `timeout`."""
    inner = (
        "import sys; sys.path.insert(0,%r)\n"
        "from qudit_qec import build_bb_code\n"
        "c=build_bb_code(%d,%d,%r,%r,field=%d,validate=False)\n"
        "d=c.get_distance_exact()\n"
        "print(int(d) if d==d else -1)\n"
    ) % (_REPO_ROOT, ell, m, A, B, q)
    try:
        r = subprocess.run(["timeout", str(int(timeout)), sys.executable, "-c", inner],
                           capture_output=True, text=True)
        out = (r.stdout or "").strip().splitlines()
        if out and out[-1].lstrip("-").isdigit():
            v = int(out[-1])
            return v if v >= 0 else None
    except Exception:
        pass
    return None


def evaluate(candidates, q, *, milp_total_timeout, milp_per_logical, exact_cap, exact_timeout):
    from qudit_qec import build_bb_code, code_distance
    from qudit_qec.evaluator import evaluate_candidate
    from qudit_qec.field_utils import is_prime
    from qudit_qec.structure import is_decomposable

    out = []
    for cand in candidates:
        ell, m = int(cand["ell"]), int(cand["m"])
        A = [tuple(t) for t in cand["A"]]
        B = [tuple(t) for t in cand["B"]]
        rec = {"ell": ell, "m": m, "q": q, "A": [list(t) for t in A], "B": [list(t) for t in B]}
        try:
            if is_prime(q):
                r = evaluate_candidate(ell, m, A, B, field=q, distance="trusted", certify=True,
                                       milp_total_timeout=milp_total_timeout,
                                       milp_timeout_per_logical=milp_per_logical)
                if r.rejected:
                    rec.update({"rejected": True, "reason": r.reason}); out.append(rec); continue
                rec.update({"n": r.n, "k": r.k, "d": r.d, "certified": bool(r.trusted),
                            "method": r.d_status, "fom": r.fom, "decomposable": r.decomposable})
            else:
                # prime-power: construct, k, QDistRnd bound, small-code exact certification
                code = build_bb_code(ell, m, A, B, field=q, validate=True)
                n, k = code.num_qudits, code.dimension
                rec.update({"n": n, "k": k, "decomposable": is_decomposable(code)})
                if k == 0:
                    rec.update({"rejected": True, "reason": "k=0"}); out.append(rec); continue
                d_res = code_distance(code, q)  # prime-power -> QDistRnd bound, untrusted
                d, certified, method = d_res.d, False, "qdistrnd_bound"
                if q ** k <= exact_cap:
                    de = _exact_via_subprocess(ell, m, rec["A"], rec["B"], q, exact_timeout)
                    if de is not None:
                        d, certified, method = de, True, "exact_enum"
                rec.update({"d": d, "certified": certified, "method": method,
                            "fom": (k * d * d / n) if d else None})
        except (ValueError, NotImplementedError) as exc:
            rec.update({"rejected": True, "reason": f"{type(exc).__name__}: {exc}"})
        except Exception as exc:  # noqa: BLE001
            rec.update({"rejected": True, "reason": f"error: {type(exc).__name__}: {exc}"})
        out.append(rec)
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Evaluate BB-code candidates over GF(q).")
    p.add_argument("--field", type=int, required=True)
    p.add_argument("--in", dest="infile", default=None, help="candidates JSON file (default stdin)")
    p.add_argument("--out", default=None, help="results JSON file (default stdout)")
    p.add_argument("--milp-total-timeout", type=float, default=90)
    p.add_argument("--milp-per-logical", type=float, default=15)
    p.add_argument("--exact-cap", type=int, default=100000, help="prime-power: try exact if q^k <= this")
    p.add_argument("--exact-timeout", type=float, default=60)
    args = p.parse_args(argv)

    raw = open(args.infile).read() if args.infile else sys.stdin.read()
    candidates = json.loads(raw)
    results = evaluate(candidates, args.field,
                       milp_total_timeout=args.milp_total_timeout,
                       milp_per_logical=args.milp_per_logical,
                       exact_cap=args.exact_cap, exact_timeout=args.exact_timeout)
    payload = json.dumps(results, indent=2)
    if args.out:
        open(args.out, "w").write(payload)
        certs = sum(1 for r in results if r.get("certified"))
        print(f"wrote {args.out}: {len(results)} evaluated, {certs} certified")
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
