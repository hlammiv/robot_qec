"""Magic-state-distillation (MSD) suitability checks for CSS qudit codes over GF(p).

This module decides, for a CSS code, **(a)** whether it admits a transversal
non-Clifford gate (via *triorthogonality*) and **(b)** how good a magic-state
distiller it is (via *weight enumerators* -> yield ``gamma`` / threshold).  The two
criteria come from the Prakash et al. cluster:

* **Triorthogonality** (Krishna-Tillich qudit form of Bravyi-Haah), arXiv:2403.06228
  Sec. II.  A length-``n`` dimension-``kappa`` subspace ``T`` of ``GF(p)^n`` is a
  *triorthogonal space* iff for every triple of vectors ``h^a, h^b, h^c in T``

      sum_i  h^a_i h^b_i h^c_i  ==  0   (mod p).                         (cubic)

  A *triorthogonal matrix* ``H = [H_1 ; H_0]`` (the object actually used to build a
  distillation code) is weaker: ``H_0`` is fully self-orthogonal and every triple
  product over the rows of ``H`` vanishes **except** the diagonal self-cube of each
  ``H_1`` ("magic") row, ``sum_i (h^a_i)^3 != 0``.  That nonzero cube is exactly what
  forces the transversal gate to leave the third level of the Clifford hierarchy.
  ``S_x = span(H_0)``, ``S_z = span(H)^perp``, logical Z's are the ``H_1`` rows.

* **Weight-enumerator => performance** (qutrit strange state), arXiv:2408.00436.
  For the qutrit "strange" state ``|S> = (|1>-|2>)/sqrt2`` the success probability,
  output error and threshold of an ``[[n,k]]_3`` distiller are fixed by the *simple
  weight enumerator* ``A(z)=w(S;1,z)`` of the stabilizer group and its MacWilliams
  dual ``B(z)=w(S^perp;1,z)``.  Distillation requires ``B(-1/2) != 0`` and a
  better-than-linear noise-suppression exponent (``3A'(-1/2)+B'(-1/2)=0`` gives
  >= quadratic; for odd ``n`` that automatically gives cubic).

----------------------------------------------------------------------------------
DIMENSION SOUNDNESS  (read this before trusting any "universal magic" claim)
----------------------------------------------------------------------------------
Only the **prime** ``p`` case is fully sound here, and even there with caveats:

* ``is_triorthogonal`` / ``transversal_gate_level`` work over any prime field
  ``GF(p)`` (the cubic form is genuinely mod ``p``).  They are a *necessary*
  structural predicate for the standard MSD route.

* ``magic_state_yield`` has two modes:
    - ``mode="triorthogonal"``: ``gamma = log_d(n/k)`` with noise-suppression
      exponent ``nu = d`` (distance).  Valid for any **prime** ``p`` *given* a
      certified triorthogonal matrix.  Verified vs the paper: ``[[20,7,2]]_3`` ->
      1.51, ``[[14,4,2]]_3`` -> 1.81, ``[[15,1,3]]_2`` -> 2.46.
    - ``mode="strange"``: the full weight-enumerator A/B/threshold pipeline.
      **Qutrit-only (p=3).**  The strange state, its discrete Wigner function and the
      ``z(eps)`` map are specific to ``GF(3)``; do NOT call this for ``p != 3``.

* **Prime-power GF(p^m), m>1, and composite Z_d: NOT covered, and NOT merely a
  coding detail.**  Triorthogonality guarantees a *transversal third-level-Clifford
  gate*, but whether that yields a *universal* magic resource is dimension-dependent.
  The decisive fact (Borda-Rincon-Galindo arXiv:2512.20787; classical basis
  Nebe-Rains-Sloane math/0001038) is that the single-qudit Clifford group is a
  MAXIMAL finite subgroup of U(d)/phase **iff d is PRIME**:
    - For prime ``p``, Clifford acts irreducibly on ``sl(p)``, so ANY non-Clifford
      gate (e.g. the transversal T from triorthogonality) forces universality.
      Sound.
    - For ``GF(p^m)``, m>1 (the native/modular convention), Clifford is
      NON-maximal: its adjoint action on ``sl(d)`` is REDUCIBLE, splitting into
      p-adic "sectors" ``W_k = span{P_{u,v}: gcd(u,v,p^m)=p^k}``.  A non-Clifford
      gate gives universality only if it COUPLES the sectors.  Whether
      ``{Clifford + CCZ}`` does so for ``p^m`` is **OPEN**
      (``qudit_qec.universality`` has a decidable sector-coverage checker; CCZ's
      natural single-qudit reductions are Clifford and do NOT couple sectors).
      Triorthogonality is necessary, **not sufficient**, here.
    - Galois/Kronecker route: a ``GF(p^m)`` qudit is m prime-p qudits; per-factor
      prime-p magic is universal per factor, but full universality ALSO needs a
      cross-factor entangling 2-qudit gate (Brylinski).  Per-factor magic alone is
      NOT enough.
    - Composite ``Z_d`` is a *trichotomy*, NOT "magic on every factor": for
      pairwise-COPRIME prime-power factors, inter-factor generalized-CNOTs give
      universality with NO diagonal magic at all (arXiv:2512.20787 Thm 48); only a
      ``p^m`` block (m>1) still needs its own sector-coupling/Galois resource.
  Accordingly the functions below **refuse** (raise / flag ``sound=False``) for
  non-prime fields rather than silently returning a number, and the returned
  dataclasses carry an explicit ``caveats`` list.  See ``docs/08`` for the full
  per-dimension universality layer.  Never read "transversal gate" as "universal"
  for ``m>1`` or composite ``d``.

Public API:

* :func:`is_triorthogonal`        -- cubic mod-p test (+ which conditions hold).
* :func:`transversal_gate_level`  -- which Clifford-hierarchy level is transversal.
* :func:`weight_enumerator`       -- simple weight enumerator A(z) (SMALL-n only).
* :func:`magic_state_yield`       -- yield gamma (+ threshold, qutrit strange mode).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from itertools import combinations_with_replacement, product
from typing import Optional

import numpy as np

from .field_utils import is_prime

# enumeration safety cap: refuse to materialize a weight enumerator that would
# require iterating more than this many stabilizer elements (local machine is tiny).
MAX_ENUM = 1_000_000


# --------------------------------------------------------------------------- #
# Code-matrix extraction (prime-p, integer representation == field element)
# --------------------------------------------------------------------------- #
def _matrices(code, p: Optional[int]):
    """Return ``(hx, hz, lx, lz, p)`` as int arrays in ``[0, p)`` for a prime CSS code.

    Accepts either a qldpc CSS code object (uses ``matrix_x``/``matrix_z`` and
    ``get_logical_ops``) or a plain mapping/object exposing those attributes.  For
    prime ``p`` the galois integer representation coincides with the field element,
    so integer ``% p`` arithmetic is exact (same guarantee ``get_code_matrices``
    relies on).  Logical ops are best-effort: returned as ``None`` if unavailable.
    """
    if p is None:
        p = int(getattr(getattr(code, "field", None), "order", 0)) or None
    if p is None:
        raise ValueError("could not infer field order p; pass p= explicitly")
    p = int(p)

    hx = np.asarray(code.matrix_x, dtype=int) % p
    hz = np.asarray(code.matrix_z, dtype=int) % p
    lx = lz = None
    try:  # qldpc CSSCode exposes logicals; not all callers will
        from qldpc.objects import Pauli

        lx = np.asarray(code.get_logical_ops(Pauli.X), dtype=int) % p
        lz = np.asarray(code.get_logical_ops(Pauli.Z), dtype=int) % p
    except Exception:
        pass
    return hx, hz, lx, lz, p


def _require_prime(p: int) -> int:
    p = int(p)
    if not is_prime(p):
        raise NotImplementedError(
            f"distillation suitability is implemented only for PRIME p (where "
            f"Z_p == GF(p)); got p={p}. For prime-power GF(p^m) (m>1) or composite "
            f"Z_d the Clifford hierarchy / resource theory of magic differ, "
            f"triorthogonality is necessary but NOT sufficient for a universal magic "
            f"state, and a CRT-factored Z_d needs magic on EACH prime-power factor. "
            f"See the module docstring."
        )
    return p


# --------------------------------------------------------------------------- #
# (1) Triorthogonality
# --------------------------------------------------------------------------- #
@dataclass
class TriorthogonalReport:
    """Result of :func:`is_triorthogonal`.

    Attributes:
        is_space: True iff the *strict* triorthogonal-**space** condition holds --
            every triple product over the X-generator rows is ``0 mod p`` (no magic
            rows).  This means the rows span a self-orthogonal triorthogonal space.
        is_matrix: True iff the *triorthogonal-**matrix*** condition holds -- the
            rows split into ``H_0`` (fully self-orthogonal, all triples zero) and
            ``H_1`` (every triple zero except its own diagonal cube), so the code
            carries a transversal non-Clifford gate with ``len(H_1)`` magic
            logicals.  A nonempty ``H_1`` is what a *distillation* code needs.
        triorthogonal: ``is_space or is_matrix`` -- admits a transversal non-Clifford
            gate from triorthogonal structure.
        n_magic_rows: number of ``H_1`` (nonzero-cube) rows = number of transversal
            magic logicals.
        h0_self_orthogonal: True iff every double product among the inferred ``H_0``
            rows vanishes (required for ``S_x subset S_z``).
        conditions: human-readable pass/fail of each sub-condition.
        sound: True iff ``p`` is prime (the only fully sound regime).
        caveats: dimension caveats (non-prime, etc.).
    """

    is_space: bool
    is_matrix: bool
    triorthogonal: bool
    n_magic_rows: int
    h0_self_orthogonal: bool
    conditions: dict
    sound: bool
    caveats: list = field(default_factory=list)

    def __bool__(self) -> bool:  # truthy iff it admits a transversal non-Clifford gate
        return self.triorthogonal


def _rows_for_triortho(code, p, generators):
    """Pick the GF(p) row set whose cubic form defines triorthogonality.

    For a triorthogonal CSS code the magic structure lives in the X-generator space
    ``H = [H_1 ; H_0]`` (``S_x = span(H_0)``, logicals = ``H_1``).  We therefore use
    the X-stabilizer rows together with the X-logical operators when those are
    available, since the ``H_1`` (magic) rows are X-logicals of the constructed code.
    Callers may override with an explicit ``generators`` matrix (e.g. the raw
    triorthogonal matrix ``H`` straight from the paper).
    """
    if generators is not None:
        return np.asarray(generators, dtype=int) % p
    hx, hz, lx, lz, _ = _matrices(code, p)
    rows = [hx]
    if lx is not None and lx.size:
        rows.append(lx)
    H = np.vstack(rows) % p
    return H


def is_triorthogonal(code, p: Optional[int] = None, generators=None) -> TriorthogonalReport:
    """Test the (cubic, mod-p) triorthogonality of a CSS code's X-generator space.

    Computes, over the rows ``h^a`` of the X-generator space ``H`` (= X-stabilizers
    plus X-logicals, unless ``generators`` is given explicitly):

        T[a,b,c] = sum_i h^a_i h^b_i h^c_i  (mod p)     for all a<=b<=c.

    and classifies:

    * **strict space** (``is_space``): ALL ``T[a,b,c] == 0`` -- the rows span a
      triorthogonal *space* (Eq. (cubic) of arXiv:2403.06228).
    * **triorthogonal matrix** (``is_matrix``): every ``T[a,b,c] == 0`` *except*
      possibly the diagonal cubes ``T[a,a,a]`` of some rows; those nonzero-cube rows
      form ``H_1`` (magic logicals) and the rest form a self-orthogonal ``H_0``.
      This is the object Bravyi-Haah / Krishna-Tillich actually use to build a code
      with a transversal non-Clifford gate.

    Args:
        code: a CSS code (qldpc ``CSSCode`` or any object with ``matrix_x`` /
            ``get_logical_ops``).  Ignored if ``generators`` is supplied.
        p: prime field order; inferred from ``code.field.order`` if omitted.
        generators: optional explicit ``(kappa, n)`` GF(p) matrix to test directly
            (e.g. a published triorthogonal matrix ``H``); bypasses ``code``.

    Returns:
        :class:`TriorthogonalReport`.  Truthy iff it admits a transversal
        non-Clifford gate.

    Caveat:
        Sound only for **prime** ``p``.  For ``GF(p^m)`` (m>1) the cube must be the
        *field* cube and -- more importantly -- a nonzero field cube does NOT by
        itself certify a universal magic state (see module docstring); this function
        flags ``sound=False`` there but still reports the mod-``p`` cubic structure.
    """
    if p is None:
        p = int(getattr(getattr(code, "field", None), "order", 0)) or None
    sound = bool(p is not None and is_prime(int(p)))
    if p is None:
        raise ValueError("could not infer field order p; pass p= explicitly")
    p = int(p)

    H = _rows_for_triortho(code, p, generators)
    kappa = H.shape[0]

    magic_rows: list[int] = []      # rows with nonzero self-cube
    offdiag_violations = 0          # any non-cube triple that is nonzero -> not triortho
    cube_vals: dict[int, int] = {}
    for a, b, c in combinations_with_replacement(range(kappa), 3):
        s = int(np.sum((H[a] * H[b] * H[c]) % p)) % p
        if a == b == c:
            cube_vals[a] = s
            if s != 0:
                magic_rows.append(a)
        elif s != 0:
            offdiag_violations += 1

    is_matrix = offdiag_violations == 0          # all non-cube triples vanish
    is_space = is_matrix and not magic_rows      # ... and no nonzero cubes either

    # H_0 = rows with zero cube; verify they are pairwise (double-product) orthogonal,
    # which is required for S_x subset S_z (so S_x sets the distance).
    h0 = [r for r in range(kappa) if r not in magic_rows]
    h0_self_orthogonal = all(
        int(np.sum((H[i] * H[j]) % p)) % p == 0
        for i, j in combinations_with_replacement(h0, 2)
    )

    conditions = {
        "all_offdiagonal_triples_zero": offdiag_violations == 0,
        "n_nonzero_cube_rows(H1/magic)": len(magic_rows),
        "H0_double_products_zero": h0_self_orthogonal,
        "strict_triorthogonal_space": is_space,
        "triorthogonal_matrix": is_matrix,
    }
    caveats = []
    if not sound:
        caveats.append(
            f"p={p} is not prime: the cube is taken mod p, not in GF(p^m); a "
            f"transversal third-level gate here does NOT certify a universal magic "
            f"state (needs per-prime-power magic). Treat result as structural only."
        )
    if is_matrix and magic_rows:
        caveats.append(
            f"{len(magic_rows)} magic (H1) rows -> transversal non-Clifford gate; "
            f"this is a distillation code only if the H1 rows are genuine logicals."
        )

    return TriorthogonalReport(
        is_space=is_space,
        is_matrix=is_matrix,
        triorthogonal=is_matrix,
        n_magic_rows=len(magic_rows),
        h0_self_orthogonal=h0_self_orthogonal,
        conditions=conditions,
        sound=sound,
        caveats=caveats,
    )


# --------------------------------------------------------------------------- #
# (2) Transversal gate level
# --------------------------------------------------------------------------- #
@dataclass
class TransversalGate:
    """Which diagonal gate is transversal for a (tri)orthogonal CSS code.

    Attributes:
        level: Clifford-hierarchy level (1=Pauli, 2=Clifford, 3=non-Clifford "T").
        gate: short label, e.g. ``"T (3rd level)"`` or ``"Clifford/none"``.
        p: field order.
        diagonal_phase: the diagonal-gate phase angle factor (``2*pi/p^2`` for the
            qudit T of arXiv:2403.06228, ``T = sum_k exp(2*pi*i*k/p^2)|k><k|``;
            for p=3 this is the ``exp(2*pi*i*k/9)`` gate).
        sound: True iff p is prime.
        caveats: dimension caveats.
    """

    level: int
    gate: str
    p: int
    diagonal_phase: Optional[str]
    sound: bool
    caveats: list = field(default_factory=list)


def transversal_gate_level(code, p: Optional[int] = None, generators=None) -> TransversalGate:
    """Return which Clifford-hierarchy gate is transversal, from triorthogonal structure.

    A *triorthogonal matrix* (some ``H_1`` magic rows present) gives a transversal
    **third-level** (non-Clifford) diagonal gate -- for qudits the
    ``T = sum_k exp(2*pi*i*k/p^2)|k><k|`` of arXiv:2403.06228 -- which is the gate
    that powers magic-state distillation.  A strict triorthogonal *space* (no magic
    rows) carries no non-Clifford transversal gate from this structure (level <= 2).

    Sound only for prime ``p``; for ``p^m``/composite the existence of the mod-p
    diagonal gate does not certify a universal non-Clifford gate of the qudit (see
    module docstring), so ``sound=False`` and the level is reported as a structural
    upper hint, not a universality claim.
    """
    rep = is_triorthogonal(code, p=p, generators=generators)
    # recover p the same way is_triorthogonal did
    if p is None:
        p = int(getattr(getattr(code, "field", None), "order", 0))
    p = int(p)

    if rep.is_matrix and rep.n_magic_rows > 0:
        level = 3
        gate = "T (3rd-level non-Clifford diagonal)"
        phase = f"2*pi/{p}^2  (T = sum_k exp(2*pi*i*k/{p**2})|k><k|)"
    elif rep.is_space:
        # self-orthogonal triorthogonal space => transversal gate is (at most) Clifford
        level = 2
        gate = "Clifford-or-lower (no non-Clifford gate from this structure)"
        phase = None
    else:
        level = 0
        gate = "none (not triorthogonal)"
        phase = None

    caveats = list(rep.caveats)
    return TransversalGate(
        level=level,
        gate=gate,
        p=p,
        diagonal_phase=phase,
        sound=rep.sound,
        caveats=caveats,
    )


# --------------------------------------------------------------------------- #
# (3) Weight enumerator (SMALL-n only)
# --------------------------------------------------------------------------- #
@dataclass
class WeightEnumerator:
    """Simple weight enumerators of a CSS stabilizer code (small-n).

    Attributes:
        A: ``A(z) = w(S; 1, z)`` -- coefficient list, ``A[w]`` = #stabilizer
            elements of symplectic Hamming weight ``w`` (qudit nontrivial on ``w``
            sites).  ``A[0]=1``, ``sum A = |S| = p^(n-k)``.
        B: ``B(z) = w(S^perp; 1, z)`` via the MacWilliams identity
            ``B(z) = (1+(p^2-1)z)^n / p^(n-k) * A((1-z)/(1+(p^2-1)z))``.
            Coefficient list; ``sum B = |S^perp| = p^(n+k)``.
        n, k, p: code parameters.
        method: how A was obtained (``"enumerate"``).
    """

    A: list
    B: list
    n: int
    k: int
    p: int
    method: str

    def A_poly(self):
        import sympy as sp

        z = sp.symbols("z")
        return sum(int(c) * z**i for i, c in enumerate(self.A))

    def B_poly(self):
        import sympy as sp

        z = sp.symbols("z")
        return sum(int(c) * z**i for i, c in enumerate(self.B))


def _stabilizer_simple_weight_enum(hx, hz, p, n, k):
    """Enumerate the simple weight enumerator A(z) of the stabilizer group S.

    ``S`` is generated symplectically by ``(hx_row | 0)`` and ``(0 | hz_row)``.  A
    group element ``(x | z)`` has *symplectic* Hamming weight = number of sites ``i``
    with ``(x_i, z_i) != (0, 0)``.  We iterate the ``p^(rank Hx) * p^(rank Hz)``
    elements (capped by :data:`MAX_ENUM`).  Exact, prime-p only.
    """
    import galois

    GF = galois.GF(p)
    Hx = GF(hx % p)
    Hz = GF(hz % p)
    # reduce to row bases so the count is p^rank
    Hx_b = Hx[~np.all(np.asarray(Hx) == 0, axis=1)] if Hx.shape[0] else Hx
    Hz_b = Hz[~np.all(np.asarray(Hz) == 0, axis=1)] if Hz.shape[0] else Hz
    rx = int(np.linalg.matrix_rank(Hx_b)) if Hx_b.shape[0] else 0
    rz = int(np.linalg.matrix_rank(Hz_b)) if Hz_b.shape[0] else 0
    total = p**rx * p**rz
    if total > MAX_ENUM:
        raise ValueError(
            f"weight-enumerator enumeration needs {total} elements > cap {MAX_ENUM}; "
            f"refuse on this small machine (n={n}). Raise MAX_ENUM or use a bigger box."
        )

    # all X-parts and Z-parts (row span)
    x_space = [np.asarray(GF(list(a)) @ Hx_b) % p for a in product(range(p), repeat=Hx_b.shape[0])] if Hx_b.shape[0] else [np.zeros(n, dtype=int)]
    z_space = [np.asarray(GF(list(b)) @ Hz_b) % p for b in product(range(p), repeat=Hz_b.shape[0])] if Hz_b.shape[0] else [np.zeros(n, dtype=int)]

    counts = np.zeros(n + 1, dtype=object)
    for xv in x_space:
        xnz = xv != 0
        for zv in z_space:
            w = int(np.count_nonzero(xnz | (zv != 0)))
            counts[w] += 1
    return [int(c) for c in counts]


def _macwilliams(A, p, n, k):
    """Compute B(z) from A(z) via the GF(p^2) MacWilliams identity (2408.00436 eq.)."""
    import sympy as sp

    z = sp.symbols("z")
    A_poly = sum(int(c) * z**i for i, c in enumerate(A))
    pf = p**2 - 1
    sub = (1 - z) / (1 + pf * z)
    B_poly = sp.expand((1 + pf * z) ** n / p ** (n - k) * A_poly.subs(z, sub))
    B_poly = sp.Poly(sp.simplify(B_poly), z)
    coeffs = [0] * (n + 1)
    for (deg,), c in B_poly.terms():
        coeffs[deg] = int(c)
    return coeffs


def weight_enumerator(code, p: Optional[int] = None) -> WeightEnumerator:
    """Compute the simple weight enumerators ``A(z)``, ``B(z)`` of a CSS code (small-n).

    Uses the classical ``GF(p^2)`` <-> stabilizer correspondence of arXiv:2408.00436:
    ``A(z) = w(S;1,z)`` enumerates the symplectic Hamming weights of the stabilizer
    group ``S``; ``B(z) = w(S^perp;1,z)`` follows from the MacWilliams identity.

    **Small-n only.**  Enumerating ``S`` is exponential; this refuses (raises) when
    ``|S| = p^(n-k)`` exceeds :data:`MAX_ENUM` (~1e6), matching the way the pipeline
    gates exact distance.  The papers' MSD searches stay at ``n <= 23`` qutrits.

    Sound only for prime ``p`` (integer == field element).  Raises for non-prime ``p``.
    """
    hx, hz, lx, lz, p = _matrices(code, p)
    _require_prime(p)
    n = hx.shape[1] if hx.shape[0] else hz.shape[1]
    # k = n - rank(Hx) - rank(Hz)
    import galois

    GF = galois.GF(p)
    rx = int(np.linalg.matrix_rank(GF(hx % p))) if hx.shape[0] else 0
    rz = int(np.linalg.matrix_rank(GF(hz % p))) if hz.shape[0] else 0
    k = n - rx - rz

    A = _stabilizer_simple_weight_enum(hx, hz, p, n, k)
    B = _macwilliams(A, p, n, k)
    return WeightEnumerator(A=A, B=B, n=n, k=k, p=p, method="enumerate")


# --------------------------------------------------------------------------- #
# (4) Magic-state yield
# --------------------------------------------------------------------------- #
@dataclass
class MagicYield:
    """Yield / threshold of a magic-state distillation routine.

    Attributes:
        gamma: yield parameter; overhead is ``O(log^gamma eps^-1)``.  Lower is better;
            ``gamma -> 1`` is the low-overhead frontier (``gamma < 1`` sublog).
        noise_suppression_exponent: ``nu`` in ``eps_out ~ eps_in^nu``.  For a
            triorthogonal code ``nu = d`` (distance).
        distills: True iff the code qualifies as a distillation routine
            (triorthogonal: has a magic logical; strange: conditions (1)&(2) hold).
        threshold: depolarizing/strange-state threshold ``eps_*`` if cheaply
            computable (``mode="strange"`` only), else ``None``.
        mode: ``"triorthogonal"`` or ``"strange"``.
        p, n, k, d: code parameters used.
        sound: True iff the regime is fully sound (prime p; strange => p==3).
        caveats: honest dimension/validity caveats.
    """

    gamma: Optional[float]
    noise_suppression_exponent: Optional[float]
    distills: bool
    threshold: Optional[float]
    mode: str
    p: int
    n: int
    k: int
    d: Optional[int]
    sound: bool
    caveats: list = field(default_factory=list)


def _strange_conditions(we: WeightEnumerator):
    """Evaluate arXiv:2408.00436 conditions (1),(2),(cubic) and threshold for p=3.

    Returns ``(distills, nu, threshold, info)``.  ``nu`` is the noise-suppression
    exponent read from the lowest power of eps dividing ``3A(z(eps))+B(z(eps))``.
    """
    import sympy as sp

    z, eps = sp.symbols("z eps")
    A = sum(int(c) * z**i for i, c in enumerate(we.A))
    B = sum(int(c) * z**i for i, c in enumerate(we.B))
    half = sp.Rational(-1, 2)

    Bm = B.subs(z, half)
    cond1 = (Bm != 0)                                   # B(-1/2) != 0
    Ad = sp.diff(A, z).subs(z, half)
    Bd = sp.diff(B, z).subs(z, half)
    cond2 = sp.simplify(3 * Ad + Bd) == 0               # >= quadratic suppression
    Am = A.subs(z, half)
    cond0 = sp.simplify(3 * Am + Bm) == 0               # auto for odd n

    # noise-suppression exponent nu: lowest power of eps in 3A(z(eps))+B(z(eps)),
    # with z(eps) = (3-eps)/(8 eps - 6).
    zeps = (3 - eps) / (8 * eps - 6)
    num = sp.expand((3 * A + B).subs(z, zeps).rewrite(sp.Pow))
    num = sp.together(num)
    n_num, n_den = sp.fraction(num)
    # divide out the denominator's eps-order; nu = ord_eps(numerator) - ord_eps(den)
    def ord_eps(expr):
        expr = sp.expand(expr)
        if expr == 0:
            return sp.oo
        poly = sp.Poly(expr, eps)
        return min(m for m in range(poly.degree() + 1) if poly.nth(m) != 0)

    nu = None
    try:
        nu = int(ord_eps(n_num) - ord_eps(n_den))
    except Exception:
        nu = None

    distills = bool(cond1 and cond2 and (nu is not None and nu >= 2))

    # threshold eps_*: smallest positive real root of eps'(eps) = eps with eps'<eps below.
    threshold = None
    try:
        # eps'(eps) = 3*(3A+B)/(4B), per msd-formula
        epsp = sp.simplify(3 * (3 * A + B).subs(z, zeps) / (4 * B.subs(z, zeps)))
        roots = sp.solve(sp.Eq(epsp, eps), eps)
        cand = [sp.re(r) for r in roots if abs(sp.im(sp.N(r))) < 1e-9 and 0 < sp.re(sp.N(r)) < 1]
        if cand:
            threshold = float(min(sp.N(c) for c in cand))
    except Exception:
        threshold = None

    info = {
        "cond1_B(-1/2)!=0": bool(cond1),
        "cond0_3A+B(-1/2)==0": bool(cond0),
        "cond2_3A'+B'(-1/2)==0": bool(cond2),
        "nu(noise_suppression)": nu,
    }
    return distills, nu, threshold, info


def magic_state_yield(
    code,
    p: Optional[int] = None,
    mode: str = "auto",
    d: Optional[int] = None,
    generators=None,
) -> MagicYield:
    """Yield parameter ``gamma`` (and threshold where cheap) for an MSD routine.

    Two modes:

    * ``mode="triorthogonal"`` -- for a triorthogonal code (noise-suppression
      exponent ``nu = d``), the yield is ``gamma = log_d(n / k)`` with ``k`` =
      number of magic (output) qudits.  Valid for any **prime** ``p``.  Pass ``d``
      (distance) explicitly to avoid heavy distance computation; ``generators`` may
      be the raw triorthogonal matrix ``H`` so ``k`` is read as the number of magic
      ``H_1`` rows.  Verified vs arXiv:2403.06228: ``[[20,7,2]]_3``->1.51,
      ``[[14,4,2]]_3``->1.81, ``[[15,1,3]]_2``->2.46.

    * ``mode="strange"`` -- **qutrit-only (p=3)** weight-enumerator pipeline of
      arXiv:2408.00436 for the strange state: compute A/B, check distillation
      conditions, read ``nu`` and (cheaply) the threshold ``eps_*``.  ``gamma`` for
      the strange-state routine is ``log_nu(n/k)`` (here ``k=1``) when ``nu`` is
      finite, else ``None`` (e.g. linear suppression -> infinite yield).

    ``mode="auto"`` uses ``"triorthogonal"`` if the code is a triorthogonal matrix,
    else falls back to ``"strange"`` for ``p==3``.

    Returns :class:`MagicYield` with an explicit ``sound`` flag and ``caveats``.
    Never interprets ``gamma`` as a *universality* claim for non-prime ``p``.
    """
    # infer p
    if p is None:
        p = int(getattr(getattr(code, "field", None), "order", 0)) or None
    if p is None and generators is None:
        raise ValueError("could not infer field order p; pass p= explicitly")

    rep = is_triorthogonal(code, p=p, generators=generators)
    p = int(p)
    _require_prime(p)

    if mode == "auto":
        mode = "triorthogonal" if (rep.is_matrix and rep.n_magic_rows > 0) else "strange"

    caveats: list = []

    if mode == "triorthogonal":
        # n, k from the triorthogonal matrix / code
        if generators is not None:
            H = np.asarray(generators, dtype=int) % p
            n = H.shape[1]
            k = rep.n_magic_rows
        else:
            hx, hz, lx, lz, _ = _matrices(code, p)
            n = hx.shape[1] if hx.shape[0] else hz.shape[1]
            k = int(getattr(code, "dimension", 0)) or rep.n_magic_rows
        if d is None:
            d = int(getattr(code, "_distance", 0)) or None
        nu = d
        if not (rep.is_matrix and rep.n_magic_rows > 0):
            caveats.append("code is NOT a triorthogonal matrix; gamma below is not valid.")
        if k == 0 or n == 0:
            gamma = None
        elif d is None:
            gamma = None
            caveats.append("distance d unknown; gamma=log_d(n/k) needs d (pass d=).")
        else:
            gamma = math.log(n / k) / math.log(d)
        if not rep.sound:
            caveats.append(
                f"p={p} non-prime: a transversal 3rd-level gate does NOT certify a "
                f"universal magic state (per-prime-power magic / CRT needed)."
            )
        return MagicYield(
            gamma=gamma,
            noise_suppression_exponent=nu,
            distills=bool(rep.is_matrix and rep.n_magic_rows > 0),
            threshold=None,
            mode="triorthogonal",
            p=p, n=n, k=k, d=d,
            sound=rep.sound,
            caveats=caveats + rep.caveats,
        )

    # strange-state mode: qutrit only
    if p != 3:
        raise NotImplementedError(
            f"strange-state weight-enumerator yield is defined only for qutrits "
            f"(p=3); got p={p}. The strange state |S>, its Wigner function and the "
            f"z(eps) map are GF(3)-specific. Use mode='triorthogonal' for general "
            f"prime p, or extend the formalism for other dimensions first."
        )
    we = weight_enumerator(code, p=p)
    distills, nu, threshold, info = _strange_conditions(we)
    n, k = we.n, we.k
    if nu and nu >= 1 and k > 0 and n > 0:
        gamma = math.log(n / k) / math.log(nu) if nu > 1 else math.inf
    else:
        gamma = None
    caveats.append(
        "strange-state mode is qutrit-only (p=3) and assumes a complete set of "
        "transversal Clifford gates (trivial-syndrome projection)."
    )
    caveats.append(f"weight-enumerator conditions: {info}")
    return MagicYield(
        gamma=gamma,
        noise_suppression_exponent=nu,
        distills=distills,
        threshold=threshold,
        mode="strange",
        p=p, n=n, k=k, d=None,
        sound=True,
        caveats=caveats,
    )
