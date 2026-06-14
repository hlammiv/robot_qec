"""Prime-``p`` magic-state-distillation (MSD) discovery arm.

A genotype -> evaluate -> catalog search whose *family* is **triorthogonal
matrices** over ``GF(p)`` (not BB pairs) and whose *objective* is the distillation
**yield** ``gamma = log_d(n/k)`` (lower is better) at fixed distance, with ``(n, k, d)``
as the Pareto axes.  It reuses the validated primitives in
:mod:`qudit_qec.distillation` (``is_triorthogonal`` as the hard validity gate,
``magic_state_yield`` for ``gamma``) and mirrors the structure of the BB arm
(:mod:`qudit_qec.genotype` / :mod:`~qudit_qec.evaluator` / :mod:`~qudit_qec.results`).
See ``docs/08`` for the plan and the dimension-dependent universality layer.

What is searched
----------------
The seed family is the **general prime-p block-triorthogonal family**

    T(p, m, k)  ->  [[ p^2 m - k ,  k ,  2 ]]_p      (k = 1 .. p*m - 1)

a length-``p^2 m`` triorthogonal *matrix* with ``k`` magic (``H_1``) rows and a
transversal third-level ``T = sum_j exp(2*pi*i*j/p^2)|j><j|`` gate.  For ``p=3`` this
is exactly the validated arXiv:2403.06228 family (``[[20,7,2]]_3`` ``gamma=1.51``,
``[[14,4,2]]_3`` ``gamma=1.81``); the generalization to every prime is verified by
``is_triorthogonal`` (see ``tests/test_distill_discovery.py``).  On top of the family
the arm applies mutation operators (column permutation, direct sum, puncturing,
column scaling) -- **every candidate is re-validated by the exact cubic gate**, so
nothing in the catalog is triorthogonal-by-assumption.

Distance / yield honesty (mirrors the trusted/untrusted split of the BB arm)
----------------------------------------------------------------------------
``gamma`` needs the distance ``d``.  Distance provenance is tracked per candidate:

* ``"known"``  -- distance carried exactly through a distance-preserving operator
  (the family seed has ``d=2``; column permutation and column scaling are monomial
  transforms that preserve ``d``; direct sum gives ``d = min(d1, d2)``).  ``gamma``
  is exact and ``trusted=True``.
* ``"upper"``  -- puncturing can only *lower* the distance, so ``d`` is recorded as
  ``<= d_parent``; ``gamma`` is reported against that upper bound (an optimistic
  *lower* bound on the true ``gamma``) and flagged ``trusted=False``.
* ``"trusted"`` -- opt-in only: build the CSS code and run the full prime-``p``
  distance pipeline.  HEAVY (MILP); off by default and meant for **lenore_remote**
  batch runs, never the local box (see ``docs/08`` compute-safety + the memory rule).

Compute safety
--------------
The validity gate is ``O(kappa^3)`` and the search is single-process with hard caps
(``max_n``, ``max_iters``); no distance MILP runs unless explicitly requested.  Keep
the local box light -- route ``distance="trusted"`` sweeps to lenore.

Public API: :func:`triortho_family`, :func:`qutrit_Tm`, the operators
(:func:`permute_columns`, :func:`direct_sum`, :func:`puncture_columns`,
:func:`scale_columns`, :func:`mutate`), :func:`evaluate_distill_candidate`,
:class:`DistillResult`, :class:`DistillCatalog`, :func:`search_distill`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field as _dc_field
from typing import Optional, Sequence

import numpy as np

from .distillation import is_triorthogonal, magic_state_yield
from .field_utils import is_prime

# Hard local cap on the code length the search will materialize (compute safety).
DEFAULT_MAX_N = 64


# --------------------------------------------------------------------------- #
# Genotype: a triorthogonal matrix over GF(p)
# --------------------------------------------------------------------------- #
def _as_matrix(H, p: int) -> np.ndarray:
    return np.asarray(H, dtype=int) % p


def genotype_key(H, p: int) -> tuple:
    """Column-multiset dedup key for a triorthogonal matrix over GF(p).

    Invariant under **column permutation** (the chief mutation) and distinguishes
    codes with different column multisets.  Conservative: two row-basis-different
    generator matrices of the same row space get *distinct* keys (never wrongly
    merged), since the magic/``H_0`` split is basis-dependent and is the genotype.
    """
    M = _as_matrix(H, p)
    cols = sorted(tuple(int(x) for x in M[:, j]) for j in range(M.shape[1]))
    return (int(p), M.shape[0], tuple(cols))


# --------------------------------------------------------------------------- #
# Constructive seed family (validated; general prime p)
# --------------------------------------------------------------------------- #
def triortho_family(p: int, m: int, k: int) -> np.ndarray:
    """The block-triorthogonal matrix ``T(p,m,k)`` -> ``[[p^2 m - k, k, 2]]_p``.

    Length ``p^2 m`` over ``p*m`` blocks of ``p`` columns: a "counting" row
    ``w=(0,1,..,p-1)`` repeated, plus block-indicator rows (block ``a`` = ones, last
    block = ``p-1``), then punctured at the first coordinate of the first ``k`` blocks
    to create ``k`` magic (nonzero-cube) rows.  Valid (``is_triorthogonal`` matrix
    with self-orthogonal ``H_0``) for ``1 <= k <= p*m - 1``.  ``p=3`` reproduces the
    arXiv:2403.06228 qutrit family exactly.
    """
    if not is_prime(p):
        raise ValueError(f"triortho_family is defined for prime p; got p={p}")
    if m < 1:
        raise ValueError(f"need m >= 1; got m={m}")
    nblocks = p * m
    if not (1 <= k <= nblocks - 1):
        raise ValueError(f"need 1 <= k <= p*m - 1 = {nblocks - 1}; got k={k}")
    n = p * p * m
    w = np.array([i % p for i in range(n)], dtype=int)
    rows = [w]
    for a in range(1, nblocks):                       # a = 1 .. p*m - 1
        v = np.zeros(n, dtype=int)
        v[(a - 1) * p:a * p] = 1                       # block a -> ones
        v[(nblocks - 1) * p:nblocks * p] = p - 1       # last block -> p-1
        rows.append(v)
    T = np.array(rows, dtype=int) % p
    del_cols = {p * j for j in range(k)}               # first coord of first k blocks
    keep = [c for c in range(n) if c not in del_cols]
    return T[:, keep] % p


def qutrit_Tm(m: int, k: int) -> np.ndarray:
    """Qutrit ``[[9m - k, k, 2]]_3`` triorthogonal matrix (``triortho_family(3,m,k)``)."""
    return triortho_family(3, m, k)


# --------------------------------------------------------------------------- #
# Mutation operators (each candidate is re-validated downstream)
# --------------------------------------------------------------------------- #
def permute_columns(H, perm: Sequence[int], p: int) -> np.ndarray:
    """Relabel columns by ``perm`` -- a code equivalence: preserves triorthogonality
    AND distance exactly."""
    M = _as_matrix(H, p)
    return M[:, list(perm)] % p


def scale_columns(H, scales: Sequence[int], p: int) -> np.ndarray:
    """Multiply column ``j`` by nonzero ``scales[j]`` -- a monomial transform that
    preserves the code distance, but NOT necessarily triorthogonality (re-validate)."""
    M = _as_matrix(H, p)
    s = np.asarray(scales, dtype=int) % p
    if np.any(s == 0):
        raise ValueError("column scales must be nonzero in GF(p)")
    return (M * s[None, :]) % p


def puncture_columns(H, cols: Sequence[int], p: int) -> np.ndarray:
    """Delete columns ``cols`` -- can only *lower* the distance; re-validate."""
    M = _as_matrix(H, p)
    drop = set(int(c) for c in cols)
    keep = [c for c in range(M.shape[1]) if c not in drop]
    return M[:, keep] % p


def direct_sum(H1, H2, p: int) -> np.ndarray:
    """Block-diagonal sum -- provably triorthogonal if both inputs are (disjoint
    supports kill every cross triple), with ``d = min(d1, d2)`` and magic rows added."""
    A = _as_matrix(H1, p)
    B = _as_matrix(H2, p)
    out = np.zeros((A.shape[0] + B.shape[0], A.shape[1] + B.shape[1]), dtype=int)
    out[:A.shape[0], :A.shape[1]] = A
    out[A.shape[0]:, A.shape[1]:] = B
    return out % p


def mutate(H, p: int, rng, *, parent_d: Optional[int] = None) -> tuple[np.ndarray, str, str, Optional[int]]:
    """Apply one random operator to ``H``.  Returns ``(H2, op_name, d_status, d_hint)``
    where ``d_status`` is ``"known"`` (permute/scale preserve d) or ``"upper"``
    (puncture: ``d <= parent_d``)."""
    M = _as_matrix(H, p)
    n = M.shape[1]
    op = rng.choice(["permute", "scale", "puncture"])
    if op == "permute":
        return permute_columns(M, rng.permutation(n), p), op, "known", parent_d
    if op == "scale":
        scales = rng.integers(1, p, size=n)            # nonzero scalars
        return scale_columns(M, scales, p), op, "known", parent_d
    # puncture 1..2 columns
    npunc = int(rng.integers(1, 3))
    cols = rng.choice(n, size=min(npunc, n - 1), replace=False)
    return puncture_columns(M, cols, p), op, "upper", parent_d


# --------------------------------------------------------------------------- #
# Evaluation
# --------------------------------------------------------------------------- #
@dataclass
class DistillResult:
    """Outcome of evaluating one triorthogonal-matrix distillation candidate."""

    p: int
    n: int
    kappa: int                 # number of generator rows
    k: int                     # number of magic (H_1) logicals = output qudits
    d: Optional[int]           # distance (None if unknown)
    d_status: str              # 'known' | 'upper' | 'trusted' | 'unknown'
    gamma: Optional[float]     # yield log_d(n/k); lower is better
    distills: bool             # passes the triorthogonal-matrix validity gate
    trusted: bool              # gamma rests on a known/certified distance
    sound: bool                # prime p (the fully sound regime)
    op: str                    # operator that produced it ('family' for seeds)
    key: tuple
    rejected: bool
    reason: str
    caveats: list = _dc_field(default_factory=list)


def build_distill_css(H, p: int):
    """Build the qldpc ``CSSCode`` of the triorthogonal matrix ``H`` (for distance).

    ``S_x = span(H_0)`` (the zero-cube rows), ``S_z = H^perp`` (right null space),
    matching the construction validated in ``tests/test_distillation.py``.  Heavy
    callers only -- importing/constructing qldpc + galois is not free.
    """
    import galois
    from qldpc.codes import CSSCode

    GF = galois.GF(p)
    M = _as_matrix(H, p)
    cubes = np.array([int(np.sum((M[r] * M[r] * M[r]) % p)) % p for r in range(M.shape[0])])
    H0 = M[cubes == 0]
    Hperp = np.asarray(GF(M).null_space()) % p
    return CSSCode(np.asarray(H0, dtype=int), Hperp.astype(int), field=p)


def evaluate_distill_candidate(
    H,
    p: int,
    *,
    d: Optional[int] = None,
    d_status: str = "known",
    distance: str = "given",
    op: str = "family",
    max_n: int = DEFAULT_MAX_N,
) -> DistillResult:
    """Validate + score one triorthogonal-matrix candidate over ``GF(p)``.

    Validity gate: :func:`~qudit_qec.distillation.is_triorthogonal` must report a
    triorthogonal *matrix* with ``>= 1`` magic row and self-orthogonal ``H_0``.
    Objective: ``gamma = log_d(n/k)`` via
    :func:`~qudit_qec.distillation.magic_state_yield`.

    ``distance``:
      * ``"given"``   -- use the passed ``d`` / ``d_status`` (default; light).
      * ``"trusted"`` -- build the CSS code and run the prime-``p`` distance pipeline
        (HEAVY MILP; route to lenore, see module docstring).
    """
    M = _as_matrix(H, p)
    n = M.shape[1]
    kappa = M.shape[0]
    key = genotype_key(M, p)

    if n > max_n:
        return DistillResult(p, n, kappa, 0, None, "unknown", None, False, False,
                             is_prime(p), op, key, True, f"n>{max_n} (compute cap)")

    rep = is_triorthogonal(None, p=p, generators=M)
    if not (rep.is_matrix and rep.n_magic_rows > 0 and rep.h0_self_orthogonal):
        reason = ("not a triorthogonal matrix" if not rep.is_matrix else
                  "no magic (H_1) row" if rep.n_magic_rows == 0 else
                  "H_0 not self-orthogonal")
        return DistillResult(p, n, kappa, rep.n_magic_rows, None, "unknown", None,
                             False, False, rep.sound, op, key, True, reason,
                             caveats=list(rep.caveats))
    k = rep.n_magic_rows

    if distance == "trusted":
        from .distance_qudit import code_distance

        code = build_distill_css(M, p)
        res = code_distance(code, p)
        d, d_status = res.d, ("trusted" if res.trusted else "upper")
    elif distance != "given":
        raise ValueError(f"distance must be 'given' or 'trusted', got {distance!r}")

    if d is None:
        return DistillResult(p, n, kappa, k, None, "unknown", None, True, False,
                             rep.sound, op, key, False, "distance unknown",
                             caveats=list(rep.caveats) +
                             ["distance pending -- run distance='trusted' on lenore"])

    my = magic_state_yield(None, p=p, mode="triorthogonal", d=d, generators=M)
    trusted = d_status in ("known", "trusted")
    caveats = list(my.caveats)
    if d_status == "upper":
        caveats.append(f"d recorded as <= {d} (puncture); gamma is an optimistic lower bound.")
    return DistillResult(
        p=p, n=n, kappa=kappa, k=k, d=d, d_status=d_status, gamma=my.gamma,
        distills=my.distills, trusted=trusted, sound=rep.sound, op=op, key=key,
        rejected=False, reason="", caveats=caveats,
    )


# --------------------------------------------------------------------------- #
# Catalog (dedup + Pareto front + gamma ranking)
# --------------------------------------------------------------------------- #
def _better(a: DistillResult, b: DistillResult) -> bool:
    """Keep the better record for the same genotype key: trusted over untrusted,
    then lower gamma, then larger known distance."""
    if a.trusted != b.trusted:
        return a.trusted
    ga = a.gamma if a.gamma is not None else math.inf
    gb = b.gamma if b.gamma is not None else math.inf
    if ga != gb:
        return ga < gb
    da = a.d if a.d is not None else -1
    db = b.d if b.d is not None else -1
    return da > db


def _dominates(a: DistillResult, b: DistillResult) -> bool:
    """``a`` dominates ``b``: no worse on (n smaller, k larger, d larger), better on
    one.  Only over candidates with a known distance."""
    if a.d is None or b.d is None:
        return False
    no_worse = a.n <= b.n and a.k >= b.k and a.d >= b.d
    strictly = a.n < b.n or a.k > b.k or a.d > b.d
    return no_worse and strictly


@dataclass
class DistillCatalog:
    """Deduplicated catalog of distillation candidates, keyed by genotype."""

    by_key: dict = _dc_field(default_factory=dict)

    def add(self, result: DistillResult) -> bool:
        """Add ``result``; return True if a new genotype, False if duplicate/rejected.
        On a duplicate key the better record (see :func:`_better`) is retained."""
        if result.rejected or result.key is None:
            return False
        existing = self.by_key.get(result.key)
        if existing is None:
            self.by_key[result.key] = result
            return True
        if _better(result, existing):
            self.by_key[result.key] = result
        return False

    def __len__(self) -> int:
        return len(self.by_key)

    def codes(self) -> list:
        return list(self.by_key.values())

    def with_gamma(self) -> list:
        return [r for r in self.by_key.values() if r.gamma is not None]

    def pareto_front(self) -> list:
        """Codes not dominated by any other on (n, k, d). Known-distance codes only."""
        codes = [r for r in self.by_key.values() if r.d is not None]
        return [a for a in codes if not any(_dominates(b, a) for b in codes if b is not a)]

    def best_by_gamma(self, top: Optional[int] = None) -> list:
        """Known-gamma codes sorted by yield ascending (lower gamma first)."""
        ranked = sorted(self.with_gamma(), key=lambda r: (r.gamma, r.n, -r.k))
        return ranked[:top] if top is not None else ranked


# --------------------------------------------------------------------------- #
# Search driver (bounded, single-process, compute-safe)
# --------------------------------------------------------------------------- #
def seed_family(primes=(3,), m_range=(1, 2), max_n: int = DEFAULT_MAX_N) -> list[DistillResult]:
    """Evaluate every valid family member ``T(p,m,k)`` over the given grid (d=2 known)."""
    out: list[DistillResult] = []
    for p in primes:
        for m in range(m_range[0], m_range[1] + 1):
            for k in range(1, p * m):
                H = triortho_family(p, m, k)
                if H.shape[1] > max_n:
                    continue
                out.append(evaluate_distill_candidate(
                    H, p, d=2, d_status="known", op="family", max_n=max_n))
    return out


def search_distill(
    primes=(3,),
    m_range=(1, 2),
    *,
    iters: int = 200,
    seed: int = 0,
    max_n: int = DEFAULT_MAX_N,
    allow_direct_sum: bool = True,
) -> DistillCatalog:
    """Bounded local search over triorthogonal matrices, seeded from the family.

    Single process, ``O(kappa^3)`` validity gate, no distance MILP (distances are
    carried through distance-preserving operators; puncture yields ``d``-upper
    candidates).  Hard caps: ``max_n`` on length, ``iters`` on mutation steps.
    Returns a :class:`DistillCatalog`.
    """
    rng = np.random.default_rng(seed)
    cat = DistillCatalog()

    # pool entries: (H, p, d, d_status) -- distance carried so gamma stays computable.
    pool: list[tuple[np.ndarray, int, Optional[int], str]] = []
    for p in primes:
        for m in range(m_range[0], m_range[1] + 1):
            for k in range(1, p * m):
                H = triortho_family(p, m, k)
                if H.shape[1] > max_n:
                    continue
                r = evaluate_distill_candidate(H, p, d=2, d_status="known",
                                               op="family", max_n=max_n)
                if not r.rejected:
                    cat.add(r)
                    pool.append((H, p, 2, "known"))

    for _ in range(iters):
        if not pool:
            break
        idx = int(rng.integers(len(pool)))
        H, p, pd, _pstat = pool[idx]
        # occasionally recombine two pool members of the same prime via direct sum
        if allow_direct_sum and rng.random() < 0.15:
            same = [e for e in pool if e[1] == p]
            H2, _, pd2, _ = same[int(rng.integers(len(same)))]
            if H.shape[1] + H2.shape[1] <= max_n:
                cand = direct_sum(H, H2, p)
                cd = min(pd, pd2) if (pd is not None and pd2 is not None) else None
                cstat, dhint = "known", cd
                op = "direct_sum"
            else:
                continue
        else:
            cand, op, cstat, dhint = mutate(H, p, rng, parent_d=pd)

        res = evaluate_distill_candidate(cand, p, d=dhint, d_status=cstat,
                                         op=op, max_n=max_n)
        if res.rejected:
            continue
        is_new = cat.add(res)
        if is_new and res.n <= max_n:
            pool.append((cand, p, res.d, res.d_status))

    return cat
