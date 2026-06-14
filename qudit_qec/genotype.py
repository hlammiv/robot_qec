"""Normalized, coefficient-carrying genotype for qudit BB codes.

The reference (qubit) pipeline represents a polynomial as a list of exponent
tuples ``(x_exp, y_exp)`` with an implicit coefficient of 1 -- which loses nothing
over F_2, where 1 is the only nonzero scalar. Over GF(q) coefficients are
load-bearing: a single monomial coefficient flip can change ``k`` (verified:
``x^3 + y + y^2`` gives ``k=8`` over GF(3) but ``x^3 + 2y + y^2`` gives ``k=0``).

This module defines the GF(q) genotype as ``(x_exp, y_exp, coeff)`` triples, with
a back-compat shim that promotes legacy 2-tuples to coefficient 1, and a single
:func:`canonicalize` routine that every cache/dedup site must route through so
that two coefficient-differing codes never silently collide.

``canonicalize`` reduces exponents modulo the lattice ``(ell, m)``, combines like
monomials by **field** addition over GF(q) (via
:func:`qudit_qec.field_utils.combine_like_terms`, never integer ``% q``), drops
field-zero coefficients, and sorts -- producing a hashable canonical key.

Public API:

* :func:`as_triple`, :func:`normalize_terms` -- the 2-tuple -> 3-tuple shim.
* :func:`canonicalize` -- lattice-reduced, field-combined, sorted tuple of triples.
* :func:`poly_key`, :func:`pair_key`, :func:`tuple_key` -- hashable dedup/cache keys.
"""

from __future__ import annotations

from typing import Iterable, Sequence

from .field_utils import _as_triple, combine_like_terms

# A canonical term is an (x_exp, y_exp, coeff) triple of ints.
Triple = tuple


def as_triple(term: Triple) -> tuple[int, int, int]:
    """Promote a legacy ``(x_exp, y_exp)`` 2-tuple to ``(x_exp, y_exp, 1)``.

    Pass-through for ``(x_exp, y_exp, coeff)`` 3-tuples.
    """
    return _as_triple(term)


def normalize_terms(terms: Iterable[Triple]) -> list[tuple[int, int, int]]:
    """Normalize every term in ``terms`` to a ``(x_exp, y_exp, coeff)`` triple."""
    return [as_triple(t) for t in terms]


def canonicalize(
    terms: Iterable[Triple], ell: int, m: int, q: int
) -> tuple[tuple[int, int, int], ...]:
    """Return the canonical, hashable form of a polynomial over GF(q) on (ell, m).

    Steps: reduce exponents modulo ``(ell, m)`` (since ``x^ell = y^m = 1`` in the
    quotient ring), combine like monomials by GF(q) field addition, drop
    field-zero coefficients, and sort by ``(x_exp, y_exp)``.

    Two polynomials are equal as ring elements on the ``(ell, m)`` lattice iff
    their canonical forms are identical, so the returned tuple is a sound dedup
    and cache key.
    """
    reduced = [(xe % ell, ye % m, c) for xe, ye, c in normalize_terms(terms)]
    combined = combine_like_terms(reduced, q)
    return tuple(combined)


def poly_key(
    terms: Iterable[Triple], ell: int, m: int, q: int
) -> tuple[tuple[int, int, int], ...]:
    """Hashable canonical key for a single polynomial (alias of canonicalize)."""
    return canonicalize(terms, ell, m, q)


def pair_key(
    a_terms: Iterable[Triple],
    b_terms: Iterable[Triple],
    ell: int,
    m: int,
    q: int,
) -> tuple:
    """Hashable canonical key for a CSS candidate ``(A, B)`` pair over GF(q)."""
    return (
        canonicalize(a_terms, ell, m, q),
        canonicalize(b_terms, ell, m, q),
    )


def tuple_key(
    terms_list: Sequence[Iterable[Triple]],
    ell: int,
    m: int,
    q: int,
) -> tuple:
    """Hashable canonical key for an n-tuple of polynomials (e.g. PBB ``(A,B,C,D)``)."""
    return tuple(canonicalize(terms, ell, m, q) for terms in terms_list)
