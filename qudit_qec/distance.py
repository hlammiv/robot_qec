"""Decoder-bound (upper-bound) and exact distance for CSS qudit codes.

Two field-aware, *non-trusted-or-slow* estimators that complement the trusted MILP
(``distance_milp``):

* :func:`decoder_bound` -- a fast stochastic **upper bound** via qldpc's decoder.
  For ``q == 2`` it uses the binary BP-OSD path (``bp_method='product_sum'``); for
  ``q > 2`` qldpc routes to the generalized Union-Find decoder (``GUFDecoder``),
  which rejects the binary kwargs -- so we call it **kwarg-free** and gate the
  BP-OSD kwargs on ``field.order == 2``. The per-sector GUF bound is cheap but
  **loose** (≈3x on GF(3)); use it only as a pre-filter, never as a trusted value.
  ``tight=True`` uses the slower both-sector ``get_distance_bound`` (≈true d).

* :func:`compute_distance_exact` -- qldpc's brute-force exact distance, run in a
  forked subprocess with an **OS-level timeout** (in-process ``signal.alarm`` around
  qldpc C calls crashes the interpreter). Gated by qudit count because the cost
  grows ~``q^k``; returns ``None`` on timeout / when skipped.
"""

from __future__ import annotations

import multiprocessing as mp

from qldpc.objects import Pauli


def _safe_int(value, default: int) -> int:
    """Convert a (possibly NaN/inf) distance to int, falling back to ``default``."""
    try:
        if value != value:  # NaN
            return default
        v = int(value)
        return v if v >= 0 else default
    except (TypeError, ValueError, OverflowError):
        return default


def decoder_bound(code, num_trials: int = 50, *, tight: bool = False) -> int:
    """Stochastic upper bound on the distance of a CSS code (field-aware).

    Args:
        code: a CSS code (``BBCode``/``CSSCode``) over GF(q).
        num_trials: decoder trials (split across X and Z sectors).
        tight: if True, use the slower both-sector ``get_distance_bound`` (tighter);
            otherwise use the cheap per-sector bound (default; loose for q > 2).

    Returns:
        An integer upper bound on ``d`` (``num_qudits`` if the decoder returns no
        finite value or the code is trivial).
    """
    n = code.num_qudits
    if code.dimension == 0:
        return 0

    if tight:
        return _safe_int(code.get_distance_bound(num_trials), n)

    q = code.field.order
    trials_x = num_trials // 2
    trials_z = num_trials - trials_x
    if q == 2:
        d_x = code.get_distance_bound_with_decoder(Pauli.X, trials_x, bp_method="product_sum")
        d_z = code.get_distance_bound_with_decoder(Pauli.Z, trials_z, bp_method="product_sum")
    else:
        # q > 2 -> GUFDecoder; the binary bp_method/osd kwargs raise TypeError here.
        d_x = code.get_distance_bound_with_decoder(Pauli.X, trials_x)
        d_z = code.get_distance_bound_with_decoder(Pauli.Z, trials_z)
    return min(_safe_int(d_x, n), _safe_int(d_z, n))


def qdistrnd_bound(code, num_trials: int = 200) -> int:
    """Independent both-sector UPPER bound on distance via qldpc's ``get_distance_bound``.

    This path is backed by GAP/QDistRnd -- a peer-reviewed randomized GF(q) distance
    tool that constructs *explicit* low-weight codewords -- so it is a genuinely
    independent second source from the MILP, and is typically much tighter than a
    loose MILP incumbent (it caught the GF(3) ``[[108,6,<=15]] -> d<=12`` overestimate).
    It is still an UPPER bound; pair it with
    :func:`qudit_qec.distance_milp.certify_distance_geq` to certify exactness.

    Portability: qldpc's ``get_distance_bound`` requires the GAP ``GUAVA`` package and
    will *prompt to install it* if missing (which raises ``EOFError`` with no TTY, or
    hangs on a TTY). If GAP/GUAVA is unavailable we fall back to the GAP-free
    per-sector GUF bound (:func:`decoder_bound`) -- looser, but never crashes or
    blocks. Install GUAVA in GAP to get the tight bound on such machines.
    """
    if code.dimension == 0:
        return 0
    try:
        return _safe_int(code.get_distance_bound(num_trials), code.num_qudits)
    except Exception:  # noqa: BLE001 - GAP/GUAVA missing (EOFError on its prompt) or any failure
        return decoder_bound(code, max(int(num_trials), 50), tight=False)


def _exact_worker(code, queue) -> None:
    try:
        d = code.get_distance_exact()
        queue.put(_safe_int(d, None) if d == d else None)
    except Exception:  # noqa: BLE001 - any failure -> "unavailable"
        queue.put(None)


def compute_distance_exact(
    code,
    timeout: float = 60,
    *,
    max_num_qudits: int = 48,
) -> int | None:
    """Exact distance via qldpc brute force, with an OS-level subprocess timeout.

    Runs in a forked child process so a hung/slow enumeration can be terminated
    cleanly (no SIGALRM, which crashes qldpc C calls). Gated by ``max_num_qudits``
    because exact enumeration cost grows roughly as ``q^k``.

    Returns the exact distance, or ``None`` if skipped (too large), timed out, or
    the worker failed.

    Known limitation (prime-power q): when the *parent* has already exercised
    prime-power ``galois`` arithmetic (e.g. via ``decoder_bound`` on a GF(p^m)
    code), it has initialized GNU OpenMP, and ``fork()`` then aborts the child
    ("fork() called from a process already using GNU OpenMP"). The child exits
    nonzero with an empty queue and this returns ``None`` -- so exact corroboration
    is currently unavailable for prime-power q. Prime q never triggers OpenMP, so
    the trusted prime-q path is unaffected. A forkserver/spawn worker that
    reconstructs the code from its (picklable) symplectic matrix is the planned
    Phase 7 fix; prime-power exact distance is also slow regardless.
    """
    if code.num_qudits > max_num_qudits:
        return None
    ctx = mp.get_context("fork")
    queue: mp.Queue = ctx.Queue()
    proc = ctx.Process(target=_exact_worker, args=(code, queue))
    proc.start()
    proc.join(timeout)
    if proc.is_alive():  # timed out -> kill
        proc.terminate()
        proc.join(5)
        return None
    # Note: an aborted child (e.g. OpenMP-after-fork on prime-power q) leaves the
    # queue empty with a nonzero exitcode -- indistinguishable here from "skipped",
    # and likewise returned as None. That is acceptable for a best-effort corroborator.
    try:
        return queue.get_nowait()
    except Exception:  # noqa: BLE001
        return None
