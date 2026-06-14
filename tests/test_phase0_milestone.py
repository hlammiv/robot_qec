"""Phase 0 acceptance milestone (docs/04-implementation-roadmap.md):

    "field_utils round-trips GF(3)/GF(4) codes; composite-d raises cleanly."

We build the gross-code polynomials over GF(2)/GF(3)/GF(4) *through our own*
``terms_to_poly`` (the coefficient-carrying genotype path) and qldpc's field-aware
BBCode, then assert the result is a genuine commuting stabilizer code over the
intended field with the expected (n, k). This is the foundation the Phase 1
construction wrapper and the Phase 2 distance backend build on.
"""

import pytest
from sympy.abc import x, y

from qudit_qec import field_utils as fu

# Gross-code generator polynomials, as coefficient-carrying genotype terms.
A_TERMS = [(3, 0, 1), (0, 1, 1), (0, 2, 1)]  # x^3 + y + y^2
B_TERMS = [(0, 3, 1), (1, 0, 1), (2, 0, 1)]  # y^3 + x + x^2

# (lattice 6x6) field -> expected (n, k), verified live against qldpc 0.2.1.
EXPECTED = {2: (72, 12), 3: (72, 8), 4: (72, 12)}


@pytest.mark.slow
@pytest.mark.parametrize("q", [2, 3, 4])
def test_round_trip_gf_q(q):
    from qldpc import codes

    poly_a = fu.terms_to_poly(A_TERMS)
    poly_b = fu.terms_to_poly(B_TERMS)
    code = codes.BBCode({x: 6, y: 6}, poly_a, poly_b, field=q)

    n_exp, k_exp = EXPECTED[q]
    assert (code.num_qudits, code.dimension) == (n_exp, k_exp)
    # commuting stabilizers over GF(q), not a subsystem code, right field
    assert fu.assert_is_stabilizer_code(code, dimension=q) is True


def test_composite_dimension_raises_cleanly():
    # The physical d=6 qudit is reached by CRT (Phase 4.5), not by GF(6).
    with pytest.raises(NotImplementedError):
        fu.get_field(6)
    with pytest.raises(NotImplementedError):
        fu.get_field(12)
