"""CRT layer: arbitrary **square-free** qudit dimension ``d`` via prime factors.

A modular-qudit stabilizer code over ``Z_d`` factors *exactly* through the Chinese
Remainder Theorem into its prime-power reductions (verified in
``docs/05-arbitrary-dimension-crt.md``): a ``Z_d`` check matrix ``H`` corresponds to
the tuple ``(H mod q_i)`` over the coprime prime-power moduli ``q_i``, and
stabilizers, logicals, ``k`` and distance all decompose per factor. In particular:

* the code's distance is ``d = min_i d_i`` (a nontrivial logical can live entirely
  in the weakest factor -- a composite-d code is only as strong as its weakest
  prime-power factor), and
* the encoded dimension is ``k_i`` per factor (in general a heterogeneous product,
  not a single uniform ``k``).

When ``d`` is **square-free** (a product of *distinct* primes), every CRT factor is
a prime *field* ``GF(p_i) = Z_{p_i}``, so we simply reduce the ``Z_d`` genotype's
coefficients mod each prime and reuse the prime-q field pipeline (construction +
trusted MILP distance) per factor. This module delivers exactly that.

Non-square-free ``d`` (a repeated prime factor, e.g. 12 = 2^2 * 3) has a ``Z_{p^a}``
*ring* factor that ``galois``/``qldpc`` cannot represent as a field; it is rejected
here and deferred to the Phase 7b Smith-normal-form ring backend.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field as _dc_field

from .evaluator import EvalResult, evaluate_candidate
from .field_utils import prime_factorization
from .genotype import normalize_terms


# --------------------------------------------------------------------------- #
# dimension classification
# --------------------------------------------------------------------------- #
def crt_moduli(d: int) -> list[int]:
    """Coprime prime-power moduli ``[p_i^{a_i}]`` of ``d`` (CRT factors), sorted."""
    return sorted(p**a for p, a in prime_factorization(d).items())


def is_squarefree(d: int) -> bool:
    """True iff ``d`` is a product of *distinct* primes (every exponent is 1)."""
    return all(a == 1 for a in prime_factorization(d).values())


def classify(d: int) -> str:
    """Classify ``d``: ``'prime'`` | ``'prime_power'`` | ``'squarefree'`` | ``'composite'``.

    * ``prime``: a single prime (3, 5, 7, ...).
    * ``prime_power``: a single prime to a power > 1 (4, 8, 9, ...).
    * ``squarefree``: >= 2 *distinct* primes, all exponent 1 (6, 10, 15, 30, ...).
    * ``composite``: >= 2 distinct primes with a repeated factor (12, 18, 20, ...).
    """
    factors = prime_factorization(d)
    if len(factors) == 1:
        return "prime" if next(iter(factors.values())) == 1 else "prime_power"
    return "squarefree" if is_squarefree(d) else "composite"


def canonicalize_zd(terms, ell: int, m: int, d: int) -> tuple:
    """Canonical ``Z_d`` (ring) genotype: reduce exponents mod (ell, m), combine like
    monomials by **integer mod-d** addition (the ring arithmetic, valid for any d),
    drop zero coefficients, sort.

    Distinct from ``genotype.canonicalize``, which combines over the *field* GF(q)
    and so requires prime-power q; here the modulus ``d`` may be composite.
    """
    acc: dict[tuple[int, int], int] = {}
    for xe, ye, c in normalize_terms(terms):
        key = (xe % ell, ye % m)
        acc[key] = (acc.get(key, 0) + int(c)) % d
    return tuple(sorted((xe, ye, c) for (xe, ye), c in acc.items() if c != 0))


def split_terms(terms, p: int) -> list[tuple[int, int, int]]:
    """Reduce a ``Z_d`` genotype's coefficients mod prime ``p`` -> a GF(p) genotype.

    Monomials whose coefficient reduces to 0 mod ``p`` vanish in that factor.
    """
    out: list[tuple[int, int, int]] = []
    for xe, ye, c in normalize_terms(terms):
        cp = int(c) % p
        if cp != 0:
            out.append((xe, ye, cp))
    return out


# --------------------------------------------------------------------------- #
# CRT-combined result
# --------------------------------------------------------------------------- #
@dataclass
class CRTResult:
    """A modular-qudit ``Z_d`` code assembled by CRT from its prime-factor codes."""

    dim: int                              # the qudit dimension d
    ell: int
    m: int
    n: int
    A: tuple                              # canonical Z_d genotype
    B: tuple
    moduli: list[int]                     # the prime factors p_i
    factors: dict                         # p_i -> EvalResult over GF(p_i)
    distance: int | None                  # min_i d_i (the code distance)
    k_per_factor: dict                    # p_i -> k_i
    d_per_factor: dict                    # p_i -> d_i (or None)
    trusted: bool                         # all factors have a trusted distance
    rejected: bool
    reason: str = ""

    @property
    def min_k(self) -> int:
        """Number of *uniform* d-dimensional logical qudits (min k over factors)."""
        return min(self.k_per_factor.values()) if self.k_per_factor else 0

    @property
    def logical_qubit_equivalents(self) -> float:
        """Total logical content in qubit-equivalents: ``sum_i k_i * log2(p_i)``.

        A field-agnostic cross-factor measure (a ``Z_d`` code's per-factor logical
        qudits weighted by their bit content). FOM normalization across fields is an
        open product decision (docs/05); this is the provisional convention.
        """
        return sum(k * math.log2(p) for p, k in self.k_per_factor.items())

    @property
    def fom(self) -> float | None:
        """Provisional FOM = logical_qubit_equivalents * d^2 / n (see docs/05 open Q)."""
        if self.distance is None:
            return None
        return self.logical_qubit_equivalents * self.distance**2 / self.n

    def summary(self) -> str:
        parts = ", ".join(
            f"GF({p}):[k={self.k_per_factor.get(p)},d={self.d_per_factor.get(p)}]"
            for p in self.moduli
        )
        tag = "trusted" if self.trusted else "untrusted"
        return f"[[{self.n}, k_per_factor, d={self.distance}]]_{self.dim} ({tag}) {{{parts}}}"


# --------------------------------------------------------------------------- #
# build / evaluate
# --------------------------------------------------------------------------- #
def _require_squarefree(d: int) -> list[int]:
    if d < 2:
        raise ValueError(f"qudit dimension d must be >= 2, got {d}")
    cls = classify(d)
    if cls in ("prime", "squarefree"):
        return crt_moduli(d)
    raise NotImplementedError(
        f"dimension d={d} is '{cls}': it has a prime-power factor Z_(p^a) with a>1, "
        f"which is a ring (not a field) that the CRT field pipeline cannot handle. "
        f"That needs the Phase 7b Smith-normal-form ring backend; see "
        f"docs/05-arbitrary-dimension-crt.md."
    )


def evaluate_crt_candidate(
    ell: int,
    m: int,
    A_terms,
    B_terms,
    dim: int,
    **eval_kwargs,
) -> CRTResult:
    """Evaluate a ``Z_dim`` CSS BB candidate by CRT-splitting over its prime factors.

    Each prime factor is built and evaluated with the field pipeline
    (:func:`~qudit_qec.evaluator.evaluate_candidate`); results are combined with
    ``distance = min_i d_i`` and ``k`` reported per factor. ``dim`` must be prime or
    square-free; extra keyword arguments are forwarded to ``evaluate_candidate``.
    """
    primes = _require_squarefree(dim)
    n = 2 * ell * m
    canon_A = canonicalize_zd(A_terms, ell, m, dim)
    canon_B = canonicalize_zd(B_terms, ell, m, dim)

    factors: dict[int, EvalResult] = {}
    for p in primes:
        factors[p] = evaluate_candidate(
            ell, m, split_terms(A_terms, p), split_terms(B_terms, p), field=p,
            **eval_kwargs,
        )

    rejected = [p for p, r in factors.items() if r.rejected]
    if rejected:
        reasons = "; ".join(f"GF({p}): {factors[p].reason}" for p in rejected)
        return CRTResult(dim, ell, m, n, canon_A, canon_B, primes, factors, None,
                         {}, {}, False, True, f"factor(s) rejected -> {reasons}")

    k_per_factor = {p: r.k for p, r in factors.items()}
    d_per_factor = {p: r.d for p, r in factors.items()}
    known = [d for d in d_per_factor.values() if d is not None]
    distance = min(known) if len(known) == len(primes) else None
    trusted = all(r.trusted for r in factors.values())

    return CRTResult(dim, ell, m, n, canon_A, canon_B, primes, factors, distance,
                     k_per_factor, d_per_factor, trusted, False, "")
