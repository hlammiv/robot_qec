"""Field substrate and guards for qudit code construction over GF(q).

This module centralizes **all** finite-field arithmetic for the qudit pipeline so
that the single most dangerous qudit pitfall lives in exactly one place:

    Applying a raw integer ``% q`` to a ``galois.FieldArray`` over a *prime-power*
    field GF(p^m) (m > 1) silently corrupts the arithmetic without raising.
    e.g. in GF(4): ``3 + 3 = 0`` and ``2 * 3 = 1`` in-field, but ``(3+3) % 4 = 2``
    and ``(2*3) % 4 = 2``.

For **prime** q the integer residue ring Z_q coincides with the field GF(q), so
``% q`` happens to be correct there — which is exactly why the prime-q MVP is the
safe starting point. Everywhere else, route field arithmetic through this module.

Dimension policy (see ``docs/05-arbitrary-dimension-crt.md``):

* ``q`` prime power  -> :func:`get_field` returns the Galois field GF(q).
* ``q`` composite (not a prime power, e.g. 6, 10, 12) -> :func:`get_field` raises
  ``NotImplementedError``.  Such dimensions are handled by CRT-factoring into
  prime-power factors (``qudit_qec.crt``, Phase 4.5) and, for the physical
  modular-qudit interpretation, a ``Z_{p^a}`` Smith-normal-form ring backend
  (Phase 7b) -- not by treating GF(d) as a field (it is not one).

Public API:

* :func:`prime_factorization`, :func:`is_prime_power` -- pure integer helpers.
* :func:`get_field` -- validated ``galois.GF(q)`` with a clear composite-d guard.
* :func:`to_field_element` -- map an integer coefficient into GF(q) correctly.
* :func:`combine_like_terms` -- merge duplicate monomials by **field** addition.
* :func:`terms_to_poly` -- exponent/coefficient terms -> sympy polynomial.
* :func:`symplectic_conjugate` -- ``[X|Z] -> [-Z|X]`` over the array's field.
* :func:`assert_is_stabilizer_code` -- in-field commutativity + non-subsystem check.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import galois
import sympy
from sympy.abc import x, y

# qldpc's own symplectic conjugation uses GF(q) arithmetic internally; reuse it
# so our commutativity test matches qldpc's stabilizer-matrix block convention.
from qldpc.math import symplectic_conjugate as _qldpc_symplectic_conjugate


# A genotype term is (x_exp, y_exp) [legacy, coeff=1] or (x_exp, y_exp, coeff).
Term = tuple


# --------------------------------------------------------------------------- #
# Integer helpers (dimension classification)
# --------------------------------------------------------------------------- #
def prime_factorization(n: int) -> dict[int, int]:
    """Return the prime factorization ``{p: exponent}`` of ``n >= 2``."""
    if n < 2:
        raise ValueError(f"prime_factorization requires n >= 2, got {n}")
    factors: dict[int, int] = {}
    m = int(n)
    d = 2
    while d * d <= m:
        while m % d == 0:
            factors[d] = factors.get(d, 0) + 1
            m //= d
        d += 1
    if m > 1:
        factors[m] = factors.get(m, 0) + 1
    return factors


def is_prime_power(n: int) -> bool:
    """True iff ``n = p^a`` for a single prime ``p`` and ``a >= 1`` (n >= 2)."""
    if not isinstance(n, (int, np.integer)) or n < 2:
        return False
    return len(prime_factorization(int(n))) == 1


def is_prime(n: int) -> bool:
    """True iff ``n`` is a prime (n >= 2).

    Distinguished from :func:`is_prime_power` because the mod-q MILP distance
    formulation is valid only for *prime* q (where Z_q == GF(q)), not for
    prime-power q = p^m with m > 1.
    """
    if not isinstance(n, (int, np.integer)) or n < 2:
        return False
    factors = prime_factorization(int(n))
    return len(factors) == 1 and next(iter(factors.values())) == 1


# --------------------------------------------------------------------------- #
# Field construction
# --------------------------------------------------------------------------- #
def get_field(q: int) -> type[galois.FieldArray]:
    """Return the Galois field ``GF(q)`` as a ``galois.FieldArray`` subclass.

    Args:
        q: A prime power (2, 3, 4, 5, 7, 8, 9, ...). The qudit/Galois dimension.

    Raises:
        ValueError: if ``q`` is not an integer >= 2.
        NotImplementedError: if ``q`` is composite (not a prime power). GF(q) is
            not a field for such q; these dimensions are reached by CRT factoring
            (``qudit_qec.crt``) and a Z_{p^a} ring backend, not by this function.
    """
    if not isinstance(q, (int, np.integer)) or q < 2:
        raise ValueError(f"field order q must be an integer >= 2, got {q!r}")
    q = int(q)
    if not is_prime_power(q):
        fac = " * ".join(
            f"{p}^{a}" if a > 1 else f"{p}"
            for p, a in sorted(prime_factorization(q).items())
        )
        raise NotImplementedError(
            f"dimension {q} = {fac} is not a prime power, so GF({q}) is not a "
            f"field. Composite dimensions are handled by CRT-factoring into "
            f"prime-power factors (qudit_qec.crt, Phase 4.5) and a Z_(p^a) ring "
            f"backend for physical modular qudits (Phase 7b); see "
            f"docs/05-arbitrary-dimension-crt.md. Do not treat GF({q}) as a field."
        )
    return galois.GF(q)


def to_field_element(c: int, field: type[galois.FieldArray]) -> galois.FieldArray:
    """Map an integer coefficient ``c`` to an element of ``field`` correctly.

    For a prime field (degree 1) the integer is reduced modulo the order. For an
    extension field GF(p^m) the integer must already index a valid element
    ``0 <= c < q`` -- there is no meaningful ``% q`` reduction in that case.
    """
    c = int(c)
    if field.degree == 1:  # prime field: Z_p == GF(p), reduce the residue
        return field(c % field.order)
    if not (0 <= c < field.order):  # extension field: c is an element index
        raise ValueError(
            f"coefficient {c} is out of range for GF({field.order}); extension-"
            f"field coefficients must be element indices in [0, {field.order})."
        )
    return field(c)


# --------------------------------------------------------------------------- #
# Polynomial / term helpers
# --------------------------------------------------------------------------- #
def _as_triple(term: Term) -> tuple[int, int, int]:
    """Normalize a 2-tuple ``(x_exp, y_exp)`` to ``(x_exp, y_exp, 1)``."""
    if len(term) == 2:
        xe, ye = term
        return int(xe), int(ye), 1
    if len(term) == 3:
        xe, ye, c = term
        return int(xe), int(ye), int(c)
    raise ValueError(f"term must be (x_exp, y_exp[, coeff]), got {term!r}")


def combine_like_terms(terms: Iterable[Term], q: int) -> list[tuple[int, int, int]]:
    """Merge duplicate monomials by **field** addition of coefficients over GF(q).

    Coefficients are added in GF(q) (never via integer ``% q``), and monomials
    whose coefficients sum to the field zero are dropped. Exponents are compared
    as-is -- this function does *not* reduce them modulo a lattice (see
    :func:`qudit_qec.genotype.canonicalize`, which reduces first).

    Returns a list of ``(x_exp, y_exp, coeff)`` triples with ``coeff`` a nonzero
    GF(q) element (its integer representation), sorted by ``(x_exp, y_exp)``.
    """
    field = get_field(q)
    acc: dict[tuple[int, int], galois.FieldArray] = {}
    for term in terms:
        xe, ye, c = _as_triple(term)
        key = (xe, ye)
        val = to_field_element(c, field)
        acc[key] = (acc[key] + val) if key in acc else val
    out: list[tuple[int, int, int]] = []
    for (xe, ye), val in acc.items():
        if int(val) != 0:  # drop field-zero coefficients
            out.append((xe, ye, int(val)))
    out.sort(key=lambda t: (t[0], t[1]))
    return out


def terms_to_poly(terms: Iterable[Term], q: int | None = None) -> sympy.Expr:
    """Convert ``(x_exp, y_exp[, coeff])`` terms to a sympy polynomial in x, y.

    2-tuples are treated as coefficient 1. Uses the ``sympy.abc`` x, y symbols
    that qldpc's BBCode expects.

    If ``q`` is given, each integer coefficient is first normalized into a valid
    GF(q) element index via :func:`to_field_element` (prime-field coefficients are
    reduced modulo q; extension-field coefficients are range-checked), keeping this
    in lock-step with :func:`combine_like_terms` / ``genotype.canonicalize``.

    If ``q`` is ``None``, coefficients are emitted as-is and the caller is
    responsible for keeping them in range: ``qldpc.codes.BBCode`` interprets each
    integer as a GF(q) element index and **raises** on out-of-range values -- it
    does *not* auto-reduce, not even for prime fields. Prefer passing ``q`` (or
    routing terms through ``genotype.canonicalize`` first) to avoid that.
    """
    field = get_field(q) if q is not None else None
    expr: sympy.Expr = sympy.Integer(0)
    for term in terms:
        xe, ye, c = _as_triple(term)
        if field is not None:
            c = int(to_field_element(c, field))
        expr += int(c) * x**xe * y**ye
    return expr


# --------------------------------------------------------------------------- #
# Symplectic / stabilizer checks (field-correct)
# --------------------------------------------------------------------------- #
def symplectic_conjugate(matrix: galois.FieldArray) -> galois.FieldArray:
    """Return ``[X|Z] -> [-Z|X]`` over the array's field (qldpc convention).

    Thin pass-through to ``qldpc.math.symplectic_conjugate`` so our commutativity
    test uses exactly qldpc's stabilizer-matrix block layout and field arithmetic.
    """
    return _qldpc_symplectic_conjugate(matrix)


def assert_is_stabilizer_code(code, *, dimension: int | None = None) -> bool:
    """Assert that ``code`` is a genuine stabilizer code over its field.

    Checks, using GF(q) arithmetic:

    1. all stabilizer generators pairwise commute under the symplectic form
       (``M @ symplectic_conjugate(M).T == 0``), and
    2. ``code`` is not a subsystem (gauge) code, and
    3. optionally, that ``code.field.order == dimension``.

    Returns ``True`` on success; raises ``AssertionError`` otherwise. This is the
    guard that catches a code accidentally built over the wrong field (e.g. a
    ``QuditCode`` constructed without ``field=q``), which otherwise reports a
    plausible-but-wrong ``k``.
    """
    field = code.field
    if dimension is not None and field.order != dimension:
        raise AssertionError(
            f"code is over GF({field.order}) but GF({dimension}) was expected"
        )
    matrix = code.matrix  # FieldArray, shape (num_checks, 2n)
    comm = np.asarray(matrix @ symplectic_conjugate(matrix).T)
    nnz = int(np.count_nonzero(comm))
    if nnz != 0:
        raise AssertionError(
            f"stabilizers do not commute over GF({field.order}): "
            f"{nnz} nonzero symplectic products"
        )
    if getattr(code, "is_subsystem_code", False):
        raise AssertionError(
            "code is a subsystem code (gauge degrees of freedom present), not a "
            "pure stabilizer code -- usually a sign/field error in construction"
        )
    return True
