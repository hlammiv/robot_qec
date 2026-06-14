"""Tests for qudit_qec.crt — arbitrary square-free dimension via CRT factoring.

Includes the Phase 4.5 milestone: a Z_6 (qubit (x) qutrit) BB code built and
MILP-verified per prime factor and reported as one code with distance = min(d_2, d_3).
"""

import pytest

from qudit_qec import crt
from qudit_qec.crt import classify, crt_moduli, evaluate_crt_candidate, is_squarefree, split_terms

ONE_X = [(0, 0, 1), (1, 0, 1)]   # 1 + x
ONE_Y = [(0, 0, 1), (0, 1, 1)]   # 1 + y
GROSS_A = [(3, 0, 1), (0, 1, 1), (0, 2, 1)]
GROSS_B = [(0, 3, 1), (1, 0, 1), (2, 0, 1)]


# --------------------------------------------------------------------------- #
# classification
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("d,cls", [
    (3, "prime"), (5, "prime"), (7, "prime"),
    (4, "prime_power"), (8, "prime_power"), (9, "prime_power"),
    (6, "squarefree"), (10, "squarefree"), (15, "squarefree"), (30, "squarefree"),
    (12, "composite"), (18, "composite"), (20, "composite"),
])
def test_classify(d, cls):
    assert classify(d) == cls


@pytest.mark.parametrize("d,mods", [
    (6, [2, 3]), (10, [2, 5]), (30, [2, 3, 5]), (12, [3, 4]), (5, [5]),
])
def test_crt_moduli(d, mods):
    assert crt_moduli(d) == mods


def test_is_squarefree():
    assert is_squarefree(6) and is_squarefree(30) and is_squarefree(5)
    assert not is_squarefree(4) and not is_squarefree(12)


def test_split_terms_reduces_and_drops():
    # over the prime 2: coeff 2 -> 0 (dropped), coeff 3 -> 1
    assert split_terms([(1, 0, 2), (0, 1, 3), (2, 0, 1)], 2) == [(0, 1, 1), (2, 0, 1)]
    # over the prime 3: coeff 3 -> 0 (dropped), coeff 4 -> 1
    assert split_terms([(1, 0, 3), (0, 1, 4)], 3) == [(0, 1, 1)]


# --------------------------------------------------------------------------- #
# square-free guard
# --------------------------------------------------------------------------- #
def test_rejects_non_squarefree_dimension():
    with pytest.raises(NotImplementedError):
        evaluate_crt_candidate(2, 2, ONE_X, ONE_Y, 12)  # 12 = 2^2 * 3 has a ring factor
    with pytest.raises(NotImplementedError):
        evaluate_crt_candidate(2, 2, ONE_X, ONE_Y, 4)


# --------------------------------------------------------------------------- #
# Phase 4.5 milestone + structure
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_milestone_z6_code_distance_is_min_over_factors():
    res = evaluate_crt_candidate(2, 2, ONE_X, ONE_Y, 6)  # qubit (x) qutrit, [[8,...]]_6
    assert res.rejected is False
    assert res.moduli == [2, 3]
    assert res.n == 8
    # each factor MILP-verified
    assert res.factors[2].trusted and res.factors[3].trusted
    assert res.trusted is True
    # distance is the minimum over factors
    assert res.distance == min(res.d_per_factor.values())
    assert res.distance == res.d_per_factor[2] == res.d_per_factor[3] == 2
    assert res.k_per_factor == {2: 2, 3: 2}


@pytest.mark.slow
def test_z6_distance_is_min_when_factors_differ():
    # gross (6,6): GF(2) gives d=6, GF(3) gives a (smaller-or-equal) bound; the code
    # distance is the minimum. Use the cheap bound stage to keep it fast.
    res = evaluate_crt_candidate(6, 6, GROSS_A, GROSS_B, 6, distance="bound")
    assert res.rejected is False and res.moduli == [2, 3]
    assert res.n == 72
    assert res.k_per_factor == {2: 12, 3: 8}  # heterogeneous logical structure
    assert res.distance == min(res.d_per_factor.values())


@pytest.mark.slow
def test_crt_rejects_when_a_factor_degenerates():
    # coeff-2 polynomial vanishes mod 2 -> the GF(2) factor is the zero polynomial.
    res = evaluate_crt_candidate(2, 2, [(0, 0, 2), (1, 0, 2)], ONE_Y, 6, distance="none")
    assert res.rejected is True and "GF(2)" in res.reason


def test_crt_result_fom_and_summary():
    res = evaluate_crt_candidate(2, 2, ONE_X, ONE_Y, 6, distance="none")
    # distance None under k-only, so fom is None but k_per_factor is populated
    assert res.k_per_factor == {2: 2, 3: 2}
    assert res.distance is None and res.fom is None
    assert "_6" in res.summary()
