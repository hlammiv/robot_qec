"""MILP exact distance for CSS qudit codes over a **prime** field GF(q).

Generalizes the Landahl-Anderson-Rice / Bravyi binary integer-programming
formulation (arXiv:1108.5738, arXiv:2308.07915) to GF(q), solved with HiGHS via
``scipy.optimize.milp`` (thread-safe, no SIGALRM). This is the *trusted* distance
signal of the qudit pipeline; the GUF decoder bound (``distance.py``) is only a
cheap, loose pre-filter.

**Valid only for prime q.** For prime q the integer residue ring Z_q equals the
field GF(q), so the mod-q constraints coincide with field arithmetic. For a
prime-power q = p^m (m > 1) the integer-slack encoding is mathematically
**invalid** (in GF(4) the field product 2*3 = 1, but integer (2*3) % 4 = 2), so
these functions raise ``NotImplementedError``; route prime-power distance through
the GUF bound + exact enumeration instead (see ``distance_qudit``).

Formulation (per logical generator): minimize the number of nonzero qudits
(Hamming weight) of an operator ``x in {0..q-1}^n`` that

* commutes with every check:        ``H[r] . x ≡ 0 (mod q)``   -> ``H[r].x - q*s_r = 0``
* anticommutes with the logical:    ``L . x ≡ 1 (mod q)``       -> ``L.x   - q*t   = 1``

The weight is linearized with a per-qudit binary indicator ``w_j`` and the big-M
constraint ``(q-1)*w_j - x_j >= 0`` (so ``x_j > 0 => w_j = 1``), minimizing
``sum_j w_j``.

The ``L . x ≡ 1`` form is valid for **prime** q: any nontrivial logical pairing to
``r != 0`` can be scaled by the field unit ``r^{-1}`` to pair to 1 with *identical*
support weight, so fixing the residue to 1 loses no minimum-weight operator.
``d = min(d_X, d_Z)``.
"""

from __future__ import annotations

import logging
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

from qldpc.objects import Pauli

from .field_utils import is_prime

logger = logging.getLogger(__name__)


def _require_prime(q: int) -> int:
    q = int(q)
    if not is_prime(q):
        raise NotImplementedError(
            f"the mod-q MILP distance formulation is valid only for prime q; got "
            f"q={q}. For prime-power q = p^m (m>1) the integer-slack encoding is "
            f"invalid (GF(p^m) arithmetic != integer mod q); use the GUF bound + "
            f"exact enumeration (distance_qudit) instead."
        )
    return q


def get_code_matrices(code, q: int):
    """Return ``(hx, hz, lx, lz)`` as integer arrays for a CSS code over GF(q).

    Entries are the galois integer representations (already in ``[0, q)``); for
    prime q these coincide with the field elements, so integer mod-q arithmetic on
    them is exact. ``hx = H_X (m_x, n)``, ``hz = H_Z (m_z, n)``,
    ``lx = X-logicals (k, n)``, ``lz = Z-logicals (k, n)``.
    """
    hx = np.asarray(code.matrix_x, dtype=int) % q
    hz = np.asarray(code.matrix_z, dtype=int) % q
    lx = np.asarray(code.get_logical_ops(Pauli.X), dtype=int) % q
    lz = np.asarray(code.get_logical_ops(Pauli.Z), dtype=int) % q
    return hx, hz, lx, lz


def _build_milp_pieces(check_matrix, logical_op, q: int):
    """Shared (objective-free) constraint system for the GF(q) min-weight-logical MILP.

    Variable layout: ``x(n) in [0,q-1] | w(n) binary | s(m) int slacks | t int slack``.
    Encodes ``H[r].x - q*s_r = 0`` (commute), ``L.x - q*t = 1`` (anticommute), and the
    weight indicator ``(q-1)*w_j - x_j >= 0``. Returns
    ``(nv, ix_w, rows, lb, ub, vlb, vub)``; the caller sets the objective (and may
    append rows, e.g. a weight cap).
    """
    check_matrix = np.asarray(check_matrix, dtype=int)
    logical_op = np.asarray(logical_op, dtype=int)
    m, n = check_matrix.shape
    nv = 2 * n + m + 1
    ix_x, ix_w, ix_t = slice(0, n), slice(n, 2 * n), 2 * n + m

    rows: list[np.ndarray] = []
    lb: list[float] = []
    ub: list[float] = []
    for r in range(m):
        row = np.zeros(nv)
        row[ix_x] = check_matrix[r]
        row[2 * n + r] = -q
        rows.append(row); lb.append(0.0); ub.append(0.0)
    row = np.zeros(nv)
    row[ix_x] = logical_op
    row[ix_t] = -q
    rows.append(row); lb.append(1.0); ub.append(1.0)
    for j in range(n):
        row = np.zeros(nv)
        row[n + j] = q - 1
        row[j] = -1
        rows.append(row); lb.append(0.0); ub.append(np.inf)

    vlb = np.zeros(nv)
    vub = np.zeros(nv)
    vub[ix_x] = q - 1
    vub[ix_w] = 1
    for r in range(m):
        vub[2 * n + r] = int(np.sum(check_matrix[r]) * (q - 1) // q) + 1
    vub[ix_t] = int(np.sum(logical_op) * (q - 1) // q) + 1
    return nv, ix_w, rows, lb, ub, vlb, vub


def ilp_min_weight(check_matrix, logical_op, q: int, timeout: float = 30):
    """Min-weight operator over GF(q) commuting with ``check_matrix``, anticommuting
    with ``logical_op``.

    Args:
        check_matrix: ``(m, n)`` integer matrix of checks (entries in ``[0, q)``).
        logical_op: ``(n,)`` integer logical operator.
        q: a prime field order.
        timeout: solver time limit in seconds (<=0 means unlimited).

    Returns:
        ``(weight, optimal)`` -- ``weight`` is the min weight (int), ``optimal`` is
        True if proven optimal. ``(None, False)`` if no feasible solution was found.
    """
    q = _require_prime(q)
    nv, ix_w, rows, lb, ub, vlb, vub = _build_milp_pieces(check_matrix, logical_op, q)

    c = np.zeros(nv)
    c[ix_w] = 1.0  # minimize the number of nonzero qudits

    opts = {"presolve": True}
    if 0 < timeout < 1e9:
        opts["time_limit"] = float(timeout)

    result = milp(
        c=c,
        constraints=LinearConstraint(np.array(rows), np.array(lb), np.array(ub)),
        integrality=np.ones(nv),
        bounds=Bounds(vlb, vub),
        options=opts,
    )
    if result.x is not None:
        return int(round(result.fun)), bool(result.success)
    return None, False


def ilp_feasible_weight_le(check_matrix, logical_op, q: int, max_weight: int, timeout: float = 20):
    """Decision MILP: does a GF(q) operator exist that commutes with ``check_matrix``,
    anticommutes with ``logical_op``, and has Hamming weight ``<= max_weight``?

    Returns ``'feasible'`` (such an operator exists), ``'infeasible'`` (HiGHS
    *proved* none exists -- a certified lower-bound contribution), or ``'unknown'``
    (the solver hit the time limit). This is the lower-bound counterpart to
    :func:`ilp_min_weight`, used by :func:`certify_distance_geq`.
    """
    q = _require_prime(q)
    nv, ix_w, rows, lb, ub, vlb, vub = _build_milp_pieces(check_matrix, logical_op, q)
    cap = np.zeros(nv)
    cap[ix_w] = 1.0  # sum_j w_j <= max_weight
    rows = rows + [cap]
    lb = lb + [0.0]
    ub = ub + [float(max_weight)]

    opts = {"presolve": True}
    if 0 < timeout < 1e9:
        opts["time_limit"] = float(timeout)

    result = milp(
        c=np.zeros(nv),  # pure feasibility (no objective)
        constraints=LinearConstraint(np.array(rows), np.array(lb), np.array(ub)),
        integrality=np.ones(nv),
        bounds=Bounds(vlb, vub),
        options=opts,
    )
    if result.status == 2:  # scipy.milp: 2 == proven infeasible
        return "infeasible"
    if result.x is not None:
        return "feasible"
    return "unknown"


def certify_distance_geq(
    code,
    q: int,
    w: int,
    *,
    timeout_per_logical: float = 20,
    total_timeout: float = 120,
) -> dict:
    """Attempt to PROVE that the code distance is ``>= w`` over prime GF(q).

    Shows that no nontrivial logical has weight ``<= w-1`` by proving the
    weight-capped decision MILP infeasible for every logical generator in both
    sectors. Paired with an explicit weight-``w`` operator (e.g. from
    :func:`qudit_qec.distance.qdistrnd_bound`), a ``certified=True`` result
    establishes ``d = w`` *exactly* -- a certificate the minimize-weight MILP
    cannot give when it only finds incumbents.

    Returns ``{certified, status, logicals_proven, total_logicals, time_s}`` where
    ``status`` is ``'certified'`` (d >= w proven), ``'refuted'`` (a weight ``<= w-1``
    logical exists, so d < w), or ``'undetermined'`` (a solve timed out).
    """
    q = _require_prime(q)
    k = code.dimension
    if k == 0 or w <= 1:
        return {"certified": True, "status": "certified", "logicals_proven": 0,
                "total_logicals": 0, "time_s": 0.0}

    hx, hz, lx, lz = get_code_matrices(code, q)
    t0 = time.monotonic()
    proven = 0
    undetermined = False
    for check, logicals in ((hx, lx), (hz, lz)):
        for i in range(k):
            remaining = total_timeout - (time.monotonic() - t0)
            if remaining <= 0:
                undetermined = True
                break
            res = ilp_feasible_weight_le(
                check, logicals[i], q, w - 1, timeout=min(timeout_per_logical, remaining)
            )
            if res == "feasible":
                return {"certified": False, "status": "refuted", "logicals_proven": proven,
                        "total_logicals": 2 * k, "time_s": time.monotonic() - t0}
            if res == "infeasible":
                proven += 1
            else:  # unknown -> cannot certify, but keep scanning for a refutation
                undetermined = True
        if undetermined:
            break
    certified = proven == 2 * k
    return {"certified": certified,
            "status": "certified" if certified else "undetermined",
            "logicals_proven": proven, "total_logicals": 2 * k,
            "time_s": time.monotonic() - t0}


def compute_distance_milp(
    code,
    q: int | None = None,
    *,
    timeout_per_logical: float = 30,
    total_timeout: float = 120,
    early_stop: int | None = None,
    verbose: bool = False,
) -> tuple[int, dict]:
    """Compute the exact distance of a CSS code over prime GF(q) via MILP.

    ``d_Z`` (Z-operators commuting with X-checks, anticommuting with X-logicals) is
    computed first, then ``d_X``. With ``early_stop`` set, returns as soon as the
    running distance drops to that threshold (an upper bound, not certified). For a
    *trusted exact* distance pass ``early_stop=None`` (the default).

    Returns ``(d, details)`` where ``details['exact']`` is True iff every logical
    solved to proven optimality within the budget.
    """
    q = _require_prime(q if q is not None else code.field.order)
    if code.field.order != q:
        raise ValueError(f"code is over GF({code.field.order}) but q={q} was given")

    n = code.num_qudits
    k = code.dimension
    if k == 0:
        # No logical qudits: distance is conventionally the block length n (the
        # sentinel "no nontrivial logical exists"); reported as exact.
        return n, {"d_x": n, "d_z": n, "k": 0, "exact": True,
                   "num_logicals_checked": 0, "total_logicals": 0, "time_s": 0.0}

    if timeout_per_logical <= 0:
        timeout_per_logical = float("inf")
    if total_timeout <= 0:
        total_timeout = float("inf")

    hx, hz, lx, lz = get_code_matrices(code, q)
    t_start = time.monotonic()
    checked = optimal = incumbent = 0
    all_solved = True

    def _remaining() -> float:
        return max(0.0, total_timeout - (time.monotonic() - t_start))

    def _sweep(check, logicals, label: str, best: int) -> int:
        nonlocal checked, optimal, incumbent, all_solved
        for i in range(k):
            rem = _remaining()
            if rem <= 0:
                all_solved = False
                break
            w, opt = ilp_min_weight(check, logicals[i], q, timeout=min(timeout_per_logical, rem))
            checked += 1
            if w is not None:
                best = min(best, w)
                if opt:
                    optimal += 1
                else:
                    incumbent += 1
                    all_solved = False
                if verbose:
                    logger.info("%s[%d]: w=%d%s", label, i, w, "" if opt else " (incumbent)")
            else:
                all_solved = False
            if early_stop is not None and best <= early_stop:
                break
        return best

    d_z = _sweep(hx, lx, "Z", n)
    if early_stop is not None and d_z <= early_stop:
        return d_z, {"d_x": n, "d_z": d_z, "k": k, "exact": False,
                     "d_x_computed": False, "num_logicals_checked": checked,
                     "logicals_optimal": optimal, "logicals_incumbent": incumbent,
                     "total_logicals": 2 * k, "time_s": time.monotonic() - t_start}

    d_x = _sweep(hz, lz, "X", n)
    d = min(d_x, d_z)
    return d, {"d_x": d_x, "d_z": d_z, "k": k, "exact": all_solved,
               "d_x_computed": True, "num_logicals_checked": checked,
               "logicals_optimal": optimal, "logicals_incumbent": incumbent,
               "total_logicals": 2 * k, "time_s": time.monotonic() - t_start}
