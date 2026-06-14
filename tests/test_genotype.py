"""Unit tests for qudit_qec.genotype — the coefficient-carrying canonical key."""

import pytest

from qudit_qec import genotype as gt


def test_as_triple_shim():
    assert gt.as_triple((1, 2)) == (1, 2, 1)
    assert gt.as_triple((1, 2, 4)) == (1, 2, 4)


def test_normalize_terms():
    assert gt.normalize_terms([(0, 1), (2, 3, 2)]) == [(0, 1, 1), (2, 3, 2)]


def test_canonicalize_reduces_exponents_mod_lattice():
    # x^6 == 1 on a (6, 6) lattice; (7, 8) -> (1, 2)
    assert gt.canonicalize([(6, 0)], 6, 6, 3) == ((0, 0, 1),)
    assert gt.canonicalize([(7, 8)], 6, 6, 3) == ((1, 2, 1),)


def test_canonicalize_combines_after_reduction():
    # (0,1) and (6,1) collapse to the same monomial on (6,6): 1 + 1 = 2 over GF(3)
    assert gt.canonicalize([(0, 1), (6, 1)], 6, 6, 3) == ((0, 1, 2),)


def test_canonicalize_order_independent():
    a = gt.canonicalize([(2, 0), (0, 1), (1, 1)], 6, 6, 5)
    b = gt.canonicalize([(1, 1), (2, 0), (0, 1)], 6, 6, 5)
    assert a == b


def test_canonicalize_drops_field_zero():
    # 1 + 2 = 0 over GF(3) -> empty polynomial
    assert gt.canonicalize([(0, 1, 1), (0, 1, 2)], 6, 6, 3) == ()


def test_coefficients_are_not_collapsed_anti_collision():
    # THE load-bearing property: coeff-differing polynomials get distinct keys.
    assert gt.poly_key([(0, 1, 1)], 6, 6, 3) != gt.poly_key([(0, 1, 2)], 6, 6, 3)


def test_two_tuple_shim_equals_coeff_one():
    assert gt.canonicalize([(0, 1)], 6, 6, 3) == gt.canonicalize([(0, 1, 1)], 6, 6, 3)


def test_pair_key_and_tuple_key():
    pk = gt.pair_key([(3, 0)], [(0, 3)], 6, 6, 3)
    assert pk == (((3, 0, 1),), ((0, 3, 1),))
    tk = gt.tuple_key([[(3, 0)], [(0, 3)], [(1, 1, 2)]], 6, 6, 3)
    assert tk == (((3, 0, 1),), ((0, 3, 1),), ((1, 1, 2),))
