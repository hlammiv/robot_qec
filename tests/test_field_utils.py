"""Unit tests for qudit_qec.field_utils — the field substrate and its guards."""

import numpy as np
import pytest
import sympy
from sympy.abc import x, y

from qudit_qec import field_utils as fu


# --------------------------------------------------------------------------- #
# integer helpers
# --------------------------------------------------------------------------- #
def test_prime_factorization():
    assert fu.prime_factorization(12) == {2: 2, 3: 1}
    assert fu.prime_factorization(7) == {7: 1}
    assert fu.prime_factorization(72) == {2: 3, 3: 2}
    assert fu.prime_factorization(2) == {2: 1}


def test_prime_factorization_rejects_small():
    with pytest.raises(ValueError):
        fu.prime_factorization(1)


@pytest.mark.parametrize("n", [2, 3, 4, 5, 7, 8, 9, 16, 25, 27])
def test_is_prime_power_true(n):
    assert fu.is_prime_power(n) is True


@pytest.mark.parametrize("n", [1, 0, -3, 6, 10, 12, 15, 100])
def test_is_prime_power_false(n):
    assert fu.is_prime_power(n) is False


# --------------------------------------------------------------------------- #
# get_field
# --------------------------------------------------------------------------- #
def test_get_field_prime():
    GF = fu.get_field(3)
    assert GF.order == 3 and GF.degree == 1 and GF.characteristic == 3


def test_get_field_prime_power():
    GF = fu.get_field(4)
    assert GF.order == 4 and GF.degree == 2 and GF.characteristic == 2


def test_get_field_composite_raises_with_crt_pointer():
    with pytest.raises(NotImplementedError) as exc:
        fu.get_field(6)
    msg = str(exc.value)
    assert "not a prime power" in msg
    assert "CRT" in msg  # the guard points at the CRT path


@pytest.mark.parametrize("bad", [1, 0, -2, 2.5, "3"])
def test_get_field_bad_input_raises_value_error(bad):
    with pytest.raises(ValueError):
        fu.get_field(bad)


# --------------------------------------------------------------------------- #
# to_field_element
# --------------------------------------------------------------------------- #
def test_to_field_element_prime_reduces():
    GF = fu.get_field(3)
    assert int(fu.to_field_element(5, GF)) == 2  # 5 % 3
    assert int(fu.to_field_element(-1, GF)) == 2  # -1 mod 3


def test_to_field_element_extension_in_range():
    GF = fu.get_field(4)
    assert int(fu.to_field_element(3, GF)) == 3


def test_to_field_element_extension_out_of_range_raises():
    GF = fu.get_field(4)
    with pytest.raises(ValueError):
        fu.to_field_element(5, GF)  # no meaningful % 4 over GF(4)


# --------------------------------------------------------------------------- #
# combine_like_terms — the field-arithmetic core
# --------------------------------------------------------------------------- #
def test_combine_drops_field_zero_prime():
    # 1 + 2 = 0 in GF(3) -> the monomial vanishes
    assert fu.combine_like_terms([(0, 1, 1), (0, 1, 2)], 3) == []


def test_combine_sums_prime():
    assert fu.combine_like_terms([(0, 1, 1), (0, 1, 1)], 3) == [(0, 1, 2)]


def test_combine_uses_field_arithmetic_not_integer_mod_gf4():
    # 3 + 3 = 0 in GF(4) (char 2); integer (3+3) % 4 = 2 would wrongly keep it.
    assert fu.combine_like_terms([(0, 1, 3), (0, 1, 3)], 4) == []
    # 2 + 3 = 1 in GF(4) (XOR); integer (2+3) % 4 = 1 coincides here, value must be 1.
    assert fu.combine_like_terms([(1, 0, 2), (1, 0, 3)], 4) == [(1, 0, 1)]


def test_combine_two_tuple_shim():
    # legacy 2-tuples are coefficient 1: 1 + 1 = 2 over GF(3)
    assert fu.combine_like_terms([(0, 1), (0, 1)], 3) == [(0, 1, 2)]


def test_combine_sorts_output():
    out = fu.combine_like_terms([(2, 0, 1), (0, 0, 1), (1, 1, 1)], 5)
    assert out == sorted(out, key=lambda t: (t[0], t[1]))


# --------------------------------------------------------------------------- #
# terms_to_poly
# --------------------------------------------------------------------------- #
def test_terms_to_poly_matches_sympy():
    expr = fu.terms_to_poly([(3, 0), (0, 1), (0, 2)])
    assert sympy.expand(expr) == sympy.expand(x**3 + y + y**2)


def test_terms_to_poly_carries_coefficients():
    expr = fu.terms_to_poly([(0, 1, 2)])
    assert sympy.expand(expr) == sympy.expand(2 * y)


def test_terms_to_poly_normalizes_when_field_given():
    # prime field: 5 -> 2 (mod 3), -1 -> 2 (mod 3); keeps genotype/construction in sync
    assert sympy.expand(fu.terms_to_poly([(0, 1, 5)], q=3)) == sympy.expand(2 * y)
    assert sympy.expand(fu.terms_to_poly([(0, 1, -1)], q=3)) == sympy.expand(2 * y)
    # without q, out-of-range coeffs pass through unchanged (caller's responsibility)
    assert sympy.expand(fu.terms_to_poly([(0, 1, 5)])) == sympy.expand(5 * y)


def test_terms_to_poly_extension_field_range_checks():
    # GF(4) has no meaningful '% 4' reduction; 5 is not a valid element index
    with pytest.raises(ValueError):
        fu.terms_to_poly([(0, 1, 5)], q=4)


# --------------------------------------------------------------------------- #
# stabilizer-code guard
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_assert_is_stabilizer_code_passes_gf3():
    from qldpc import codes

    code = codes.BBCode({x: 6, y: 6}, x**3 + y + y**2, y**3 + x + x**2, field=3)
    assert fu.assert_is_stabilizer_code(code) is True
    assert fu.assert_is_stabilizer_code(code, dimension=3) is True


@pytest.mark.slow
def test_assert_is_stabilizer_code_dimension_mismatch_raises():
    from qldpc import codes

    code = codes.BBCode({x: 6, y: 6}, x**3 + y + y**2, y**3 + x + x**2, field=3)
    with pytest.raises(AssertionError):
        fu.assert_is_stabilizer_code(code, dimension=5)
