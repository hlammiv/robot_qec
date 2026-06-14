"""Tests for qudit_qec.construct — CSS BB construction over GF(q).

Includes the Phase 1 milestone (docs/04): build_bb_code reproduces [[72,8]] over
GF(3), is not a subsystem code, and a single-monomial coefficient flip changes k.
"""

import pytest

from qudit_qec import construct as cn
from qudit_qec import field_utils as fu

# Gross-code generator polynomials as coefficient-carrying genotype terms.
A = [(3, 0, 1), (0, 1, 1), (0, 2, 1)]  # x^3 + y + y^2
B = [(0, 3, 1), (1, 0, 1), (2, 0, 1)]  # y^3 + x + x^2


# --------------------------------------------------------------------------- #
# validate_terms (fast — no qldpc construction)
# --------------------------------------------------------------------------- #
def test_validate_terms_canonicalizes():
    assert cn.validate_terms(6, 6, A, 3, "A") == ((0, 1, 1), (0, 2, 1), (3, 0, 1))


def test_validate_terms_combines_duplicates_over_field():
    # (0,1,1) + (0,1,2) = 0 over GF(3) -> that monomial vanishes, leaving x^3 + y^2
    out = cn.validate_terms(6, 6, [(3, 0, 1), (0, 1, 1), (0, 1, 2), (0, 2, 1)], 3, "A")
    assert out == ((0, 2, 1), (3, 0, 1))


def test_validate_terms_rejects_zero_polynomial():
    with pytest.raises(ValueError, match="zero polynomial"):
        cn.validate_terms(6, 6, [(0, 1, 1), (0, 1, 2)], 3, "A")  # 1+2=0 mod 3


def test_validate_terms_enforces_term_count():
    with pytest.raises(ValueError, match="distinct monomials"):
        cn.validate_terms(6, 6, A, 3, "A", min_terms=1, max_terms=2)


def test_validate_terms_composite_field_raises():
    with pytest.raises(NotImplementedError):
        cn.validate_terms(6, 6, A, 6, "A")


def test_validate_terms_bad_lattice_raises():
    with pytest.raises(ValueError):
        cn.validate_terms(0, 6, A, 3, "A")


# --------------------------------------------------------------------------- #
# build_bb_code  (slow — qldpc construction)
# --------------------------------------------------------------------------- #
@pytest.mark.slow
@pytest.mark.parametrize("q,nk", [(2, (72, 12)), (3, (72, 8)), (4, (72, 12))])
def test_build_reproduces_known_params(q, nk):
    code = cn.build_bb_code(6, 6, A, B, field=q, verify_stabilizer=True)
    assert cn.get_code_params_fast(code) == nk
    assert code.is_subsystem_code is False


@pytest.mark.slow
def test_milestone_single_coeff_flip_changes_k():
    # docs/04 Phase 1 milestone regression: a one-monomial coefficient flip
    # collapses k from 8 to 0 over GF(3).
    base_k = cn.build_bb_code(6, 6, A, B, field=3).dimension
    flipped_A = [(3, 0, 1), (0, 1, 2), (0, 2, 1)]  # x^3 + 2y + y^2
    flip_k = cn.build_bb_code(6, 6, flipped_A, B, field=3).dimension
    assert base_k == 8
    assert flip_k == 0
    assert base_k != flip_k


@pytest.mark.slow
def test_build_combines_cancelling_duplicate_without_crash():
    # The construction-time hazard: an un-combined duplicate would become 3*y and
    # make galois.GF(3)(3) raise. canonicalize-in-build must prevent that.
    dup_A = [(3, 0, 1), (0, 1, 1), (0, 1, 2), (0, 2, 1)]  # y-term cancels -> x^3 + y^2
    code = cn.build_bb_code(6, 6, dup_A, B, field=3, validate=False, verify_stabilizer=True)
    # equals the code built directly from the canonical x^3 + y^2
    ref = cn.build_bb_code(6, 6, [(3, 0, 1), (0, 2, 1)], B, field=3)
    assert cn.get_code_params_fast(code) == cn.get_code_params_fast(ref)


@pytest.mark.slow
def test_build_reduces_out_of_range_exponent():
    # x^7 == x on a (6,6) lattice; building from either must give the same code.
    oor = cn.build_bb_code(6, 6, [(7, 0, 1), (0, 1, 1), (0, 2, 1)], B, field=3)
    red = cn.build_bb_code(6, 6, [(1, 0, 1), (0, 1, 1), (0, 2, 1)], B, field=3)
    assert cn.get_code_params_fast(oor) == cn.get_code_params_fast(red)


@pytest.mark.slow
def test_build_default_field_is_qubit():
    code = cn.build_bb_code(6, 6, A, B)  # field defaults to 2
    assert code.field.order == 2
    assert cn.get_code_params_fast(code) == (72, 12)


@pytest.mark.slow
def test_code_params_convenience():
    assert cn.code_params(6, 6, A, B, field=3) == (72, 8)


def test_build_composite_field_raises_fast():
    # composite q is rejected before any qldpc work
    with pytest.raises(NotImplementedError):
        cn.build_bb_code(6, 6, A, B, field=6)
