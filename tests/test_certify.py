"""Tests for the QDistRnd second-source + weight-cut certification (Phase 2.5).

Verifies the lower-bound proof (ilp_feasible_weight_le / certify_distance_geq), the
QDistRnd upper bound, and that the dispatcher pins an exact distance via
QDistRnd-upper-bound + weight-cut-lower-bound when the minimize-MILP only finds an
incumbent.
"""

import numpy as np
import pytest

from qudit_qec import construct as cn
from qudit_qec import distance as dist
from qudit_qec import distance_milp as dm
from qudit_qec.distance_qudit import code_distance

ONE_X = [(0, 0, 1), (1, 0, 1)]   # 1 + x
ONE_Y = [(0, 0, 1), (0, 1, 1)]   # 1 + y
SWAP_A = [(2, 0, 1), (0, 1, 1), (0, 3, 1)]
SWAP_B = [(0, 2, 1), (1, 0, 1), (3, 0, 1)]


# --------------------------------------------------------------------------- #
# weight-cut feasibility
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_feasible_weight_le_decision():
    code = cn.build_bb_code(2, 2, ONE_X, ONE_Y, field=3)  # [[8,2,2]]_3
    hx, hz, lx, lz = dm.get_code_matrices(code, 3)
    # a weight-2 logical exists -> feasible at max_weight>=2, infeasible at max_weight=1
    assert dm.ilp_feasible_weight_le(hz, lz[0], 3, 2) == "feasible"
    assert dm.ilp_feasible_weight_le(hz, lz[0], 3, 1) == "infeasible"


@pytest.mark.slow
def test_certify_distance_geq_brackets_d2():
    code = cn.build_bb_code(2, 2, ONE_X, ONE_Y, field=3)  # d = 2
    assert dm.certify_distance_geq(code, 3, 2)["certified"] is True   # d >= 2 proven
    refute = dm.certify_distance_geq(code, 3, 3)                       # d >= 3 is false
    assert refute["certified"] is False and refute["status"] == "refuted"


@pytest.mark.slow
def test_certify_distance_geq_for_d4_code():
    code = cn.build_bb_code(4, 4, SWAP_A, SWAP_B, field=3)  # [[32,2,4]]_3
    assert dm.certify_distance_geq(code, 3, 4)["certified"] is True    # d >= 4
    assert dm.certify_distance_geq(code, 3, 5)["certified"] is False   # d >= 5 false (d=4)


# --------------------------------------------------------------------------- #
# QDistRnd second source
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_qdistrnd_bound_is_upper_bound():
    code = cn.build_bb_code(4, 4, SWAP_A, SWAP_B, field=3)  # true d = 4
    assert dist.qdistrnd_bound(code, num_trials=100) >= 4


# --------------------------------------------------------------------------- #
# dispatcher: certification pins an exact distance
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_dispatcher_certifies_d4_via_cut_when_milp_capped():
    code = cn.build_bb_code(4, 4, SWAP_A, SWAP_B, field=3)
    # tiny MILP budget forces an incumbent, so the QDistRnd+cut path must certify
    res = code_distance(
        code, 3,
        milp_total_timeout=0.5, milp_timeout_per_logical=0.5,  # force incumbent
        certify=True, certify_total_timeout=90,
    )
    assert res.d == 4
    assert res.trusted is True
    assert res.method in ("qdistrnd+milp_cut", "milp")  # certified either way


@pytest.mark.slow
def test_dispatcher_reports_tightened_upper_bound_without_certify():
    code = cn.build_bb_code(4, 4, SWAP_A, SWAP_B, field=3)
    res = code_distance(
        code, 3,
        milp_total_timeout=0.5, milp_timeout_per_logical=0.5,
        certify=False,
    )
    # no certification -> untrusted, but the reported d is a real upper bound (>= true 4)
    assert res.trusted is False
    assert res.d >= 4
    assert "qdistrnd_bound" in res.details
