"""Tests for qudit_qec.universality -- the p-adic sector-coverage check.

Validates the decidable instantiation (arXiv:2512.20787) of "does {Clifford + G}
restore irreducibility at d=p^m": Clifford gates must NOT couple the p-adic sectors,
a generic non-Clifford gate CAN, and -- the in-repo result -- CCZ's single-qudit
reductions are Pauli and do NOT couple sectors.
"""

import pytest

from qudit_qec import universality as U


@pytest.mark.parametrize("d", [4, 8, 16, 9, 27])
def test_clifford_gates_do_not_couple_sectors(d):
    # identity, Pauli Z, and the (correct) Clifford S all preserve the stratification
    for phases in (U.identity_phases(d), U.pauli_Z_power_phases(d, 1), U.clifford_S_phases(d)):
        sc = U.sector_coverage_diagonal(phases, d)
        assert sc.is_clifford_like is True
        assert sc.couples_all_sectors is False


@pytest.mark.parametrize("d", [4, 8, 16])
def test_positive_control_couples_all_sectors(d):
    sc = U.sector_coverage_diagonal(U.level_bump_phases(d, 1, 1, 4), d)
    assert sc.couples_all_sectors is True
    assert sc.is_clifford_like is False


@pytest.mark.parametrize("d", [4, 8, 16])
@pytest.mark.parametrize("bc", [(1, 1), (2, 3), (3, 1), (2, 2)])
def test_ccz_single_qudit_reductions_are_clifford(d, bc):
    # every computational-basis reduction of CCZ is a Pauli Z^{b*c} -> no sector coupling
    b, c = bc
    sc = U.sector_coverage_diagonal(U.ccz_single_qudit_reduction_phases(d, b, c), d)
    assert sc.couples_all_sectors is False
    assert sc.is_clifford_like is True


def test_prime_dimension_is_trivially_irreducible():
    # d prime -> m=1: sector test trivially satisfied (any non-Clifford gate is universal)
    sc = U.sector_coverage_diagonal(U.level_bump_phases(3, 1, 1, 3), 3)
    assert sc.m == 1 and sc.couples_all_sectors is True


def test_composite_dimension_refused():
    with pytest.raises(NotImplementedError):
        U.sector_coverage_diagonal(U.identity_phases(6), 6)


def test_sector_of_basic():
    # d=4=2^2: gcd(2,0,4)=2 -> valuation 1 (k=1); gcd(1,2,4)=1 -> k=0
    assert U.sector_of(2, 0, 4, 2, 2) == 1
    assert U.sector_of(1, 2, 4, 2, 2) == 0
