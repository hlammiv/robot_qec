"""Tests for qudit_qec.structure — Tanner-graph decomposability (direct-sum)."""

import pytest

from qudit_qec import build_bb_code
from qudit_qec.structure import (
    connected_components,
    is_decomposable,
    num_connected_components,
)

# all x-exponents even on a (4,4) lattice -> even/odd-x sublattices decouple
EVEN_A = [(2, 0, 1), (0, 1, 1), (0, 2, 1)]
EVEN_B = [(0, 2, 1), (2, 0, 1), (0, 1, 1)]
# x/y-swap trinomials -> indecomposable
SWAP_A = [(2, 0, 1), (0, 1, 1), (0, 3, 1)]
SWAP_B = [(0, 2, 1), (1, 0, 1), (3, 0, 1)]


@pytest.mark.slow
def test_even_exponent_code_is_decomposable():
    code = build_bb_code(4, 4, EVEN_A, EVEN_B, field=3, validate=False)
    assert num_connected_components(code) == 2
    assert is_decomposable(code) is True


@pytest.mark.slow
def test_swap_code_is_indecomposable():
    code = build_bb_code(4, 4, SWAP_A, SWAP_B, field=3, validate=False)
    assert num_connected_components(code) == 1
    assert is_decomposable(code) is False


@pytest.mark.slow
def test_components_partition_all_qudits():
    code = build_bb_code(4, 4, SWAP_A, SWAP_B, field=3, validate=False)
    comps = connected_components(code)
    flat = sorted(j for comp in comps for j in comp)
    assert flat == list(range(code.num_qudits))
    assert sum(len(c) for c in comps) == code.num_qudits  # disjoint
