"""Tests for qudit_qec.evaluator — the field-threaded evaluation cascade.

Includes the Phase 3 milestone: an end-to-end ``evaluate_candidate(..., field=3)``
returning a coherent ``{n, k, d, fom}`` with a trusted (MILP) distance.
"""

import pytest

from qudit_qec import evaluate_candidate

SWAP_A = [(2, 0, 1), (0, 1, 1), (0, 3, 1)]   # x^2 + y + y^3
SWAP_B = [(0, 2, 1), (1, 0, 1), (3, 0, 1)]   # y^2 + x + x^3
GROSS_A = [(3, 0, 1), (0, 1, 1), (0, 2, 1)]  # x^3 + y + y^2
GROSS_B = [(0, 3, 1), (1, 0, 1), (2, 0, 1)]  # y^3 + x + x^2


@pytest.mark.slow
def test_milestone_end_to_end_trusted_gf3():
    r = evaluate_candidate(4, 4, SWAP_A, SWAP_B, field=3)  # [[32,2,4]]_3
    assert (r.n, r.k, r.d) == (32, 2, 4)
    assert r.trusted is True and r.d_status == "milp"
    assert r.fom == pytest.approx(2 * 16 / 32)  # k*d^2/n = 1.0
    assert r.decomposable is False
    assert r.rejected is False
    assert r.q == 3


def test_k_only_stage():
    r = evaluate_candidate(6, 6, GROSS_A, GROSS_B, field=3, distance="none")
    assert r.k == 8 and r.d is None and r.fom is None
    assert not r.rejected and r.key is not None


@pytest.mark.slow
def test_bound_stage_is_untrusted_upper_bound():
    r = evaluate_candidate(4, 4, SWAP_A, SWAP_B, field=3, distance="bound")
    assert r.d_status == "bound" and r.trusted is False
    assert r.d >= 4  # GUF pre-filter is an upper bound on the true d=4


def test_rejects_zero_polynomial():
    r = evaluate_candidate(6, 6, [(0, 1, 1), (0, 1, 2)], GROSS_B, field=3, distance="none")
    assert r.rejected and "invalid" in r.reason


def test_rejects_composite_field():
    r = evaluate_candidate(6, 6, GROSS_A, GROSS_B, field=6, distance="none")
    assert r.rejected and "invalid" in r.reason


@pytest.mark.slow
def test_rejects_low_k():
    r = evaluate_candidate(4, 4, SWAP_A, SWAP_B, field=2, distance="none")  # k=0 over GF(2)
    assert r.k == 0 and r.rejected and "k<" in r.reason


def test_self_dual_marker():
    r = evaluate_candidate(6, 6, GROSS_A, GROSS_A, field=3, distance="none")  # A == B
    assert r.self_dual is True


def test_invalid_distance_mode_raises():
    with pytest.raises(ValueError):
        evaluate_candidate(6, 6, GROSS_A, GROSS_B, field=3, distance="bogus")
