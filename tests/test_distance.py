"""Tests for the Phase 2 distance backend over GF(q).

The trusted signal is the prime-q mod-q MILP; here it is cross-checked against an
*independent* brute-force enumeration and qldpc's exact distance on small codes,
and the prime-power guard / dispatcher trust gate are exercised. Includes the
Phase 2 milestone: a trusted (MILP-certified) exact distance for a GF(3) CSS code.
"""

import itertools

import numpy as np
import pytest
from qldpc.objects import Pauli

from qudit_qec import construct as cn
from qudit_qec import distance as dist
from qudit_qec import distance_milp as dm
from qudit_qec.distance_qudit import code_distance

# --- small CSS BB codes as coefficient-carrying genotypes -------------------- #
ONE_X = [(0, 0, 1), (1, 0, 1)]          # 1 + x
ONE_Y = [(0, 0, 1), (0, 1, 1)]          # 1 + y
TRI_X = [(0, 0, 1), (1, 0, 1), (2, 0, 1)]  # 1 + x + x^2
TRI_Y = [(0, 0, 1), (0, 1, 1), (0, 2, 1)]  # 1 + y + y^2
SWAP_A = [(2, 0, 1), (0, 1, 1), (0, 3, 1)]  # x^2 + y + y^3
SWAP_B = [(0, 2, 1), (1, 0, 1), (3, 0, 1)]  # y^2 + x + x^3


def _brute_distance(code, q, max_w=4):
    """Independent ground truth: min weight of a nontrivial logical, by enumeration."""
    hx = np.asarray(code.matrix_x, int) % q
    hz = np.asarray(code.matrix_z, int) % q
    lx = np.asarray(code.get_logical_ops(Pauli.X), int) % q
    lz = np.asarray(code.get_logical_ops(Pauli.Z), int) % q
    n = code.num_qudits

    def sector(check, logs):
        for w in range(1, max_w + 1):
            for cols in itertools.combinations(range(n), w):
                for vals in itertools.product(range(1, q), repeat=w):
                    v = np.zeros(n, int)
                    for col, val in zip(cols, vals):
                        v[col] = val
                    if np.any(check @ v % q):
                        continue
                    if np.any(logs @ v % q):
                        return w
        return None

    cands = [d for d in (sector(hx, lx), sector(hz, lz)) if d is not None]
    return min(cands) if cands else None


# --------------------------------------------------------------------------- #
# MILP vs independent ground truth
# --------------------------------------------------------------------------- #
@pytest.mark.slow
@pytest.mark.parametrize(
    "ell,m,A,B,q,expected",
    [
        (2, 2, ONE_X, ONE_Y, 3, 2),   # [[8,2,2]]_3
        (2, 2, ONE_X, ONE_Y, 5, 2),   # [[8,2,2]]_5
        (3, 3, TRI_X, TRI_Y, 3, 2),   # [[18,8,2]]_3
        (4, 4, SWAP_A, SWAP_B, 3, 4),  # [[32,2,4]]_3  (the d>2 case)
    ],
)
def test_milp_matches_bruteforce(ell, m, A, B, q, expected):
    code = cn.build_bb_code(ell, m, A, B, field=q)
    d, details = dm.compute_distance_milp(code, q, early_stop=None)
    brute = _brute_distance(code, q, max_w=5)
    assert d == expected
    assert d == brute
    assert details["exact"] is True


@pytest.mark.slow
def test_milp_matches_qldpc_exact_tiny():
    code = cn.build_bb_code(2, 2, ONE_X, ONE_Y, field=3)  # [[8,2,2]]_3
    d, _ = dm.compute_distance_milp(code, 3, early_stop=None)
    assert d == int(code.get_distance_exact())


@pytest.mark.slow
def test_milp_q2_reproduces_binary_distance():
    # the generalized MILP at q=2 must agree with brute force on a qubit code
    code = cn.build_bb_code(2, 2, ONE_X, ONE_Y, field=2)  # [[8,2,2]]_2, k>0
    assert code.dimension > 0
    d, _ = dm.compute_distance_milp(code, 2, early_stop=None)
    assert d == _brute_distance(code, 2, max_w=5)


# --------------------------------------------------------------------------- #
# prime-power guard
# --------------------------------------------------------------------------- #
def test_milp_rejects_prime_power():
    with pytest.raises(NotImplementedError):
        dm.ilp_min_weight(np.zeros((1, 4), int), np.ones(4, int), 4)


@pytest.mark.slow
def test_compute_distance_milp_rejects_prime_power_code():
    code = cn.build_bb_code(2, 2, ONE_X, ONE_Y, field=4)
    with pytest.raises(NotImplementedError):
        dm.compute_distance_milp(code, 4)


# --------------------------------------------------------------------------- #
# decoder bound (pre-filter)
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_decoder_bound_works_over_gf3_no_typeerror():
    # the binary bp_method path would raise TypeError on GF(3); gating must avoid it
    code = cn.build_bb_code(6, 6, [(3, 0, 1), (0, 1, 1), (0, 2, 1)],
                            [(0, 3, 1), (1, 0, 1), (2, 0, 1)], field=3)  # [[72,8]]_3
    ub = dist.decoder_bound(code, num_trials=20)
    assert isinstance(ub, int) and ub > 0


@pytest.mark.slow
def test_decoder_bound_is_upper_bound():
    code = cn.build_bb_code(4, 4, SWAP_A, SWAP_B, field=3)  # true d = 4
    assert dist.decoder_bound(code, num_trials=40) >= 4


# --------------------------------------------------------------------------- #
# dispatcher trust gate + MILESTONE
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_milestone_trusted_exact_distance_gf3():
    # Phase 2 milestone: a TRUSTED (MILP-certified) exact distance for a GF(3) CSS code.
    code = cn.build_bb_code(4, 4, SWAP_A, SWAP_B, field=3)  # [[32,2,4]]_3
    res = code_distance(code, milp_total_timeout=120)
    assert res.d == 4
    assert res.trusted is True and res.exact is True
    assert res.method == "milp"
    assert res.upper_bound >= res.d  # GUF pre-filter is an upper bound


@pytest.mark.slow
def test_dispatcher_prime_power_is_untrusted():
    code = cn.build_bb_code(2, 2, ONE_X, ONE_Y, field=4)  # [[8,?]]_4
    res = code_distance(code)
    assert res.trusted is False
    assert res.method == "guf_bound"


@pytest.mark.slow
def test_dispatcher_exact_corroboration():
    code = cn.build_bb_code(2, 2, ONE_X, ONE_Y, field=3)  # [[8,2,2]]_3
    res = code_distance(code, want_exact=True, exact_max_num_qudits=16)
    assert res.d == 2 and res.trusted is True
    assert res.details.get("exact_enum") == 2
