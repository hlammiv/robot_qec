"""CSS bivariate-bicycle code construction over GF(q).

This is the qudit, coefficient-aware analogue of the reference pipeline's
``evaluation/bb_code.py``. A bivariate bicycle (BB) code is defined by two
polynomials ``A, B`` over ``GF(q)[x, y] / (x^ell - 1, y^m - 1)``; over GF(q) the
qldpc library constructs the code field-generically and auto-inserts the
antipode + sign (``H_Z = [B^T, -A^T]``) that makes the stabilizers commute, so we
only thread ``field=q`` and a coefficient-carrying genotype through to it.

Two things make the construction correct over GF(q) where the F_2 path was free:

1. **Coefficients are combined in-field before construction.** Two genotype terms
   on the same monomial (e.g. ``1*y`` and ``2*y``) must be merged by GF(q)
   addition *first*: over GF(3) they sum to ``3*y == 0`` and the term vanishes.
   Handing the un-combined ``3*y`` to ``qldpc`` would call ``galois.GF(3)(3)``,
   which **raises** (galois never auto-reduces, not even for prime fields). We
   therefore route every polynomial through ``genotype.canonicalize`` first.
2. **Coefficients are load-bearing.** ``x^3 + y + y^2`` gives ``k = 8`` over GF(3),
   but ``x^3 + 2y + y^2`` gives ``k = 0`` -- a single coefficient flip. The
   genotype carries ``(x_exp, y_exp, coeff)`` triples accordingly.

Public API:

* :func:`validate_terms` -- canonicalize + range/term-count checks (raises).
* :func:`build_bb_code` -- construct a ``qldpc.codes.BBCode`` over GF(q).
* :func:`get_code_params_fast` -- ``(n, k)`` without distance.
* :func:`code_params` -- build and return ``(n, k)`` in one call.
"""

from __future__ import annotations

from typing import Iterable

from qldpc import codes
from sympy.abc import x, y

from .field_utils import assert_is_stabilizer_code, get_field, terms_to_poly
from .genotype import Triple, canonicalize


def _check_lattice(ell: int, m: int) -> None:
    if not (isinstance(ell, int) and isinstance(m, int)) or ell < 1 or m < 1:
        raise ValueError(f"lattice orders (ell, m) must be integers >= 1, got ({ell!r}, {m!r})")


def validate_terms(
    ell: int,
    m: int,
    terms: Iterable[Triple],
    q: int,
    name: str = "polynomial",
    min_terms: int = 1,
    max_terms: int = 6,
) -> tuple[tuple[int, int, int], ...]:
    """Validate and canonicalize a polynomial genotype over GF(q) on (ell, m).

    Canonicalizes (reduce exponents mod the lattice, combine like terms by GF(q)
    addition, drop field-zero coefficients, sort), then enforces:

    * the polynomial is not the zero polynomial (all coefficients cancelled), and
    * the number of *distinct* monomials lies in ``[min_terms, max_terms]``.

    Unlike the F_2 reference, duplicate monomials are **combined**, not rejected
    (over GF(q) ``c1*x^a + c2*x^a`` is a legitimate single term of coefficient
    ``c1 + c2``). Coefficient range is validated as part of canonicalization
    (``get_field``/``to_field_element``), so composite ``q`` raises
    ``NotImplementedError`` here too.

    Returns the canonical ``(x_exp, y_exp, coeff)`` tuple; raises ``ValueError``
    on a violation.
    """
    _check_lattice(ell, m)
    canon = canonicalize(terms, ell, m, q)
    if not canon:
        raise ValueError(
            f"{name} is the zero polynomial over GF({q}) (all coefficients cancelled)"
        )
    if not (min_terms <= len(canon) <= max_terms):
        raise ValueError(
            f"{name} has {len(canon)} distinct monomials over GF({q}); "
            f"expected between {min_terms} and {max_terms}"
        )
    return canon


def build_bb_code(
    ell: int,
    m: int,
    A_terms: Iterable[Triple],
    B_terms: Iterable[Triple],
    field: int = 2,
    *,
    validate: bool = True,
    min_terms: int = 1,
    max_terms: int = 6,
    verify_stabilizer: bool = False,
) -> codes.BBCode:
    """Construct a CSS bivariate-bicycle code over GF(``field``).

    Args:
        ell, m: cyclic group orders for x and y (``n = 2*ell*m`` qudits).
        A_terms, B_terms: genotype polynomials as ``(x_exp, y_exp[, coeff])``
            terms; 2-tuples mean coefficient 1.
        field: the qudit dimension / field order (prime power). Default 2 keeps
            legacy qubit behavior.
        validate: enforce term-count bounds (``min_terms..max_terms``) per
            polynomial. Coefficients are always canonicalized regardless.
        verify_stabilizer: if True, additionally assert the constructed code's
            stabilizers commute over GF(q) and it is not a subsystem code.

    Returns:
        A ``qldpc.codes.BBCode`` over GF(field).

    Raises:
        NotImplementedError: if ``field`` is composite (not a prime power).
        ValueError: on a malformed/zero polynomial or out-of-bounds term count.
    """
    _check_lattice(ell, m)
    q = int(field)
    get_field(q)  # fail fast & clearly on composite q before touching qldpc

    if validate:
        a_canon = validate_terms(ell, m, A_terms, q, "A", min_terms, max_terms)
        b_canon = validate_terms(ell, m, B_terms, q, "B", min_terms, max_terms)
    else:
        a_canon = canonicalize(A_terms, ell, m, q)
        b_canon = canonicalize(B_terms, ell, m, q)
        if not a_canon or not b_canon:
            raise ValueError("A or B is the zero polynomial over GF(q) after reduction")

    # a_canon/b_canon coefficients are already in-field range and duplicate-free,
    # so terms_to_poly cannot produce an out-of-range coefficient for qldpc.
    code = codes.BBCode(
        {x: ell, y: m}, terms_to_poly(a_canon), terms_to_poly(b_canon), field=q
    )
    if code.field.order != q:  # post-construction field guard
        raise AssertionError(
            f"constructed code is over GF({code.field.order}), expected GF({q})"
        )
    if verify_stabilizer:
        assert_is_stabilizer_code(code, dimension=q)
    return code


def get_code_params_fast(code: codes.BBCode) -> tuple[int, int]:
    """Return ``(n, k)`` without computing distance (field-aware via qldpc)."""
    return code.num_qudits, code.dimension


def code_params(
    ell: int,
    m: int,
    A_terms: Iterable[Triple],
    B_terms: Iterable[Triple],
    field: int = 2,
    **kwargs,
) -> tuple[int, int]:
    """Build a BB code over GF(field) and return ``(n, k)`` in one call."""
    return get_code_params_fast(build_bb_code(ell, m, A_terms, B_terms, field, **kwargs))
