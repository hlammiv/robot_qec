"""Qutrit strange-state / QR magic-state-distillation discovery sub-arm (d > 2).

The second distillation family of ``docs/08``, complementary to the d=2 triorthogonal
``T``-gate arm in :mod:`qudit_qec.distill_discovery`.  Where that arm optimizes the
yield ``gamma`` of triorthogonal codes, this one searches **self-orthogonal codes over
``F_3``** whose CSS distiller ``H_x = H_z = C`` distills the qutrit **strange state**
``|S> = (|1> - |2>)/sqrt 2`` -- the route of Prakash-Singhal (arXiv:2408.00436) whose
canonical member is the **11-qutrit Golay** ``[[11,1,5]]_3`` (distance ``d=5``, *cubic*
noise suppression ``nu=3``, threshold ``eps_* ~ 0.387``).  This beats the d=2 family
on noise suppression (``nu = 3`` vs ``nu = 2``); the objective axis here is the
**noise-suppression exponent / threshold**, NOT ``gamma`` (``gamma = log_d(n/k)`` does
not apply).  Note ``nu != d`` for the strange route: ``nu`` is read from the weight
enumerator (the lowest power of ``eps`` in ``3A + B``), distinct from the code
distance ``d`` (the Golay has ``d=5`` but ``nu=3``).

Genotype: a **cyclic code over ``F_3``**, ``C = <g(x)>`` for a monic divisor
``g(x) | x^n - 1`` (``n`` coprime to 3).  A self-orthogonal ``C`` (``C subset
C^perp``) gives a valid strange CSS distiller.  Validity + ``nu`` come from the
**validated** weight-enumerator pipeline of :mod:`qudit_qec.distillation`
(``weight_enumerator`` + the strange conditions ``B(-1/2) != 0`` and
``3 A'(-1/2) + B'(-1/2) = 0``).

Compute safety (memory rule + docs/08)
--------------------------------------
The strange weight enumerator enumerates ``3^{2k}`` (``= 3^{n-k_q}``) stabilizer
elements; this arm **hard-gates** that to ``MAX_ENUM`` (``~1e6``, i.e. classical
``k <= 6``, ``n <= 13`` locally) and refuses larger codes, routing them to
**lenore_remote**.  The threshold solve (``sympy.solve``) is the expensive step and is
**opt-in** (``compute_threshold=True`` / :meth:`StrangeCatalog.refine_thresholds`):
the cheap screen is ``nu`` + the two algebraic conditions.

Honest scope: the small-``n`` self-orthogonal *cyclic* family over ``F_3`` is
**sparse** -- it is essentially Golay-dominated (no self-orthogonal cyclic distiller at
``n=5,7``; the two equivalent QR codes at ``n=11`` reproduce the Golay).  Broadening to
non-cyclic self-orthogonal codes and ``n > 13`` is the lenore-scale extension.

Public API: :func:`cyclic_codes`, :func:`self_orthogonal_cyclic_codes`,
:func:`build_strange_css`, :func:`evaluate_strange_candidate`, :class:`StrangeResult`,
:class:`StrangeCatalog`, :func:`search_strange_cyclic`.
"""

from __future__ import annotations

from dataclasses import dataclass, field as _dc_field
from fractions import Fraction
from itertools import combinations
from typing import Optional

import numpy as np

from .distillation import MAX_ENUM, _strange_conditions, weight_enumerator

_P = 3  # qutrit-only: the strange state, its Wigner function and z(eps) map are GF(3)-specific

# Local-search enum default: the per-code weight enumerator iterates 3^{n-k} elements
# in Python, slow well below the absolute MAX_ENUM refuse cap.  Keep the *search* light
# (~Golay-sized, 3^10) locally; codes above this are refused with "route to lenore".
LOCAL_ENUM_DEFAULT = 60_000


# --------------------------------------------------------------------------- #
# Cyclic-code genotype over F_3
# --------------------------------------------------------------------------- #
def cyclic_codes(n: int, p: int = _P):
    """All cyclic codes ``C = <g(x)>`` of length ``n`` over ``F_p`` as generator
    matrices, one per monic divisor ``g(x) | x^n - 1``.  Requires ``gcd(n, p) = 1``."""
    import galois

    if n % p == 0:
        raise ValueError(f"need gcd(n,p)=1 for cyclic codes; got n={n}, p={p}")
    GF = galois.GF(p)
    xn_minus_1 = galois.Poly([1] + [0] * (n - 1) + [(-1) % p], field=GF)  # x^n - 1
    facs, mults = galois.factors(xn_minus_1)
    irr = []
    for f, mlt in zip(facs, mults):
        irr.extend([f] * mlt)
    out = []
    seen = set()
    for r in range(len(irr) + 1):
        for combo in combinations(range(len(irr)), r):
            g = galois.Poly([1], field=GF)
            for idx in combo:
                g = g * irr[idx]
            key = str(g)
            if key in seen:
                continue
            seen.add(key)
            G = _cyclic_generator_matrix(g, n, GF, p)
            if G is not None:
                out.append(G)
    return out


def _cyclic_generator_matrix(g, n: int, GF, p: int):
    """``k x n`` generator matrix of ``<g(x)>`` (``k = n - deg g``), or None if empty."""
    deg = g.degree
    k = n - deg
    if k <= 0:
        return None
    gcoef = np.array(g.coeffs, dtype=int)[::-1]  # low -> high degree
    G = np.zeros((k, n), dtype=int)
    for i in range(k):
        G[i, i:i + deg + 1] = gcoef
    return np.asarray(G) % p


def _is_self_orthogonal(G, p: int = _P) -> bool:
    """``C subset C^perp``: every pair of rows (and each row with itself) is orthogonal."""
    M = np.asarray(G, dtype=int) % p
    return bool(np.all((M @ M.T) % p == 0))


def self_orthogonal_cyclic_codes(n: int, p: int = _P):
    """Cyclic codes of length ``n`` over ``F_p`` that are self-orthogonal (valid strange
    CSS distillers ``H_x = H_z = C``)."""
    return [G for G in cyclic_codes(n, p) if _is_self_orthogonal(G, p)]


def build_strange_css(G, p: int = _P):
    """Build the strange CSS distiller ``CSSCode(H_x = H_z = C)`` from a self-orthogonal
    generator matrix ``C`` over ``F_p``."""
    from qldpc.codes import CSSCode

    M = np.asarray(G, dtype=int) % p
    if not _is_self_orthogonal(M, p):
        raise ValueError("generator matrix is not self-orthogonal (C subset C^perp required)")
    return CSSCode(M, M, field=p)


def _cheap_strange_conditions(A: list, B: list, n: int) -> tuple:
    """Screen the strange-distillation conditions with exact ``Fraction`` arithmetic on
    the enumerator coefficients -- no sympy, no ``z(eps)`` substitution.

    Returns ``(distills, nu_hint, info)``.  ``nu_hint`` is the *guaranteed* suppression
    lower bound (2 from cond1+cond2; 3 when ``n`` is odd and ``3A+B(-1/2)=0`` gives the
    automatic cubic of arXiv:2408.00436).  The *exact* ``nu`` (possibly higher) and the
    threshold are the opt-in sympy refine path.
    """
    half = Fraction(-1, 2)
    Bm = sum(Fraction(int(c)) * half ** i for i, c in enumerate(B))
    Am = sum(Fraction(int(c)) * half ** i for i, c in enumerate(A))
    Ad = sum(i * Fraction(int(c)) * half ** (i - 1) for i, c in enumerate(A) if i >= 1)
    Bd = sum(i * Fraction(int(c)) * half ** (i - 1) for i, c in enumerate(B) if i >= 1)
    cond1 = Bm != 0                       # B(-1/2) != 0
    cond2 = (3 * Ad + Bd) == 0            # 3A'+B' = 0  -> >= quadratic
    cond0 = (3 * Am + Bm) == 0            # 3A+B = 0    (P vanishes at the fixed point)
    # nu = order of vanishing of P=3A+B at z0=-1/2 needs BOTH P(z0)=0 (cond0) AND
    # P'(z0)=0 (cond2) to reach >=2; cond2 alone does NOT force cond0 (it can fail for
    # k>=2 codes). Require cond0 so this matches the sympy verdict exactly and nu_hint
    # stays a true LOWER bound, not an overstatement.
    distills = bool(cond1 and cond2 and cond0)
    nu_hint = (3 if (distills and n % 2 == 1) else (2 if distills else None))
    info = {"B(-1/2)!=0": cond1, "3A'+B'(-1/2)=0": cond2, "3A+B(-1/2)=0": cond0,
            "nu_lower_bound": nu_hint}
    return distills, nu_hint, info


def _genotype_key(G, p: int) -> tuple:
    """Row-space canonical key (RREF) -- dedups row-equivalent generators of the same
    code.  Permutation-equivalent codes (e.g. the two length-11 QR Golay forms) keep
    distinct keys (conservative)."""
    import galois

    GF = galois.GF(p)
    R = np.asarray(GF(np.asarray(G, dtype=int) % p).row_reduce())
    rows = tuple(tuple(int(x) for x in R[i]) for i in range(R.shape[0]) if R[i].any())
    return (int(p), R.shape[1], rows)


# --------------------------------------------------------------------------- #
# Evaluation
# --------------------------------------------------------------------------- #
@dataclass
class StrangeResult:
    """Outcome of evaluating one strange-state qutrit distillation candidate."""

    p: int
    n: int
    k: int                       # logical (output) qudits of the CSS distiller
    nu: Optional[int]            # noise-suppression exponent (cubic=3, quintic=5, ...)
    distills: bool               # passes the strange-distillation conditions
    threshold: Optional[float]   # eps_* (only if compute_threshold=True)
    d_lb: Optional[int]          # min nonzero stabilizer weight (distance hint)
    conditions: dict
    sound: bool                  # p == 3 (the only sound regime for the strange route)
    key: tuple
    rejected: bool
    reason: str
    caveats: list = _dc_field(default_factory=list)


def evaluate_strange_candidate(
    code_or_G,
    p: int = _P,
    *,
    compute_threshold: bool = False,
    max_enum: int = LOCAL_ENUM_DEFAULT,
) -> StrangeResult:
    """Validate + score a strange-state distillation candidate (qutrit only).

    ``code_or_G`` is either a self-orthogonal generator matrix ``C`` over ``F_3`` (the
    CSS distiller ``H_x = H_z = C`` is built) or a ready CSS code.  By default a cheap
    exact-``Fraction`` screen of the strange conditions gives ``distills`` + a ``nu``
    lower bound; ``compute_threshold=True`` instead runs the validated sympy pipeline
    (``_strange_conditions``) for the exact ``nu`` + threshold (expensive).  Refuses
    (``rejected``) when the weight-enumerator enumeration ``3^{n-k}`` exceeds
    ``max_enum``.
    """
    if p != _P:
        return StrangeResult(p, 0, 0, None, False, None, None, {}, False, (int(p),),
                             True, "strange route is qutrit-only (p=3)")
    # build / accept the CSS code
    if hasattr(code_or_G, "matrix_x") or hasattr(code_or_G, "num_qudits"):
        # content key (NOT id(): addresses are reused after GC -> false catalog
        # collisions). For a strange distiller Hx=Hz=C, key on the X-stabilizer C so
        # the key matches the generator-matrix path and stays rebuildable.
        code = code_or_G
        key = _genotype_key(np.asarray(code.matrix_x, dtype=int), p)
    else:
        G = np.asarray(code_or_G, dtype=int) % p
        if not _is_self_orthogonal(G, p):
            return StrangeResult(p, G.shape[1], 0, None, False, None, None, {}, True,
                                 _genotype_key(G, p), True, "C not self-orthogonal")
        key = _genotype_key(G, p)
        code = build_strange_css(G, p)

    n, kq = code.num_qudits, code.dimension
    if kq < 1:
        return StrangeResult(p, n, kq, None, False, None, None, {}, True, key, True,
                             "k=0 (no logical qudit)")
    enum = p ** (n - kq)
    if enum > max_enum:
        return StrangeResult(p, n, kq, None, False, None, None, {}, True, key, True,
                             f"weight-enum {enum} > cap {max_enum}; route to lenore",
                             caveats=[f"strange enum 3^(n-k)={enum} exceeds local cap"])

    we = weight_enumerator(code, p=p)
    d_lb = next((i for i, c in enumerate(we.A) if i > 0 and c), None)
    if compute_threshold:
        # exact nu + threshold via the validated sympy pipeline (expensive)
        distills, nu, threshold, info = _strange_conditions(we, compute_threshold=True)
    else:
        # cheap Fraction screen: conditions + nu lower bound, no sympy
        distills, nu, info = _cheap_strange_conditions(we.A, we.B, n)
        threshold = None
    return StrangeResult(
        p=p, n=n, k=kq, nu=nu, distills=distills, threshold=threshold, d_lb=d_lb,
        conditions=info, sound=True, key=key, rejected=False, reason="",
        caveats=["strange-state route is qutrit-only (p=3); assumes transversal "
                 "Clifford + trivial-syndrome projection."]
        + (["nu is a guaranteed lower bound (cheap screen); exact nu via "
            "compute_threshold=True"] if not compute_threshold else []),
    )


# --------------------------------------------------------------------------- #
# Catalog
# --------------------------------------------------------------------------- #
def _better(a: StrangeResult, b: StrangeResult) -> bool:
    """Keep the better record per genotype: higher nu, then higher threshold, then
    larger distance hint."""
    na = a.nu if a.nu is not None else -1
    nb = b.nu if b.nu is not None else -1
    if na != nb:
        return na > nb
    ta = a.threshold if a.threshold is not None else -1.0
    tb = b.threshold if b.threshold is not None else -1.0
    if ta != tb:
        return ta > tb
    return (a.d_lb or -1) > (b.d_lb or -1)


@dataclass
class StrangeCatalog:
    """Deduplicated catalog of strange-state distillation candidates."""

    by_key: dict = _dc_field(default_factory=dict)

    def add(self, result: StrangeResult) -> bool:
        if result.rejected or not result.distills:
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

    def best_by_nu(self, top: Optional[int] = None) -> list:
        """Distillers ranked by noise suppression (nu desc), then threshold, then n."""
        ranked = sorted(self.by_key.values(),
                        key=lambda r: (-(r.nu or 0), -(r.threshold or 0.0), r.n))
        return ranked[:top] if top is not None else ranked

    def refine_thresholds(self, top: int = 5) -> None:
        """Compute the (expensive) exact ``nu`` + threshold for the top-``top``
        distillers in place (sympy ``z(eps)`` series + threshold solve)."""
        for r in self.best_by_nu(top):
            if r.threshold is None:
                refined = evaluate_strange_candidate(
                    build_strange_css_from_result(r), compute_threshold=True)
                r.threshold = refined.threshold
                if refined.nu is not None:
                    r.nu = refined.nu


def build_strange_css_from_result(r: StrangeResult):
    """Rebuild the CSS code of a catalog entry from its row-space key (for threshold
    refinement)."""
    rows = r.key[2]
    G = np.array(rows, dtype=int)
    return build_strange_css(G, r.p)


# --------------------------------------------------------------------------- #
# Search driver (bounded, single-process, compute-safe)
# --------------------------------------------------------------------------- #
def search_strange_cyclic(
    n_min: int = 5,
    n_max: int = 13,
    *,
    p: int = _P,
    compute_threshold: bool = False,
    max_enum: int = LOCAL_ENUM_DEFAULT,
) -> StrangeCatalog:
    """Enumerate self-orthogonal cyclic codes over ``F_3`` for ``n_min <= n <= n_max``
    (``gcd(n,3)=1``) and catalog the valid strange distillers.  Single process,
    compute-gated by ``max_enum`` (default :data:`LOCAL_ENUM_DEFAULT`, Golay-sized);
    thresholds are opt-in.  Reproduces the Golay ``[[11,1,5]]_3`` at ``n=11``."""
    cat = StrangeCatalog()
    for n in range(n_min, n_max + 1):
        if n % p == 0:
            continue
        for G in self_orthogonal_cyclic_codes(n, p):
            res = evaluate_strange_candidate(
                G, p=p, compute_threshold=compute_threshold, max_enum=max_enum)
            cat.add(res)
    return cat
