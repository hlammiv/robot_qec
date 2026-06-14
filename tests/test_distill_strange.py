"""Tests for qudit_qec.distill_strange -- the qutrit strange/QR distillation sub-arm.

Validates that the self-orthogonal cyclic-F_3 genotype reproduces the 11-qutrit Golay
[[11,1,5]]_3 strange distiller (cubic, threshold ~0.387), that the cheap Fraction-based
condition screen agrees with the validated sympy pipeline, that the compute gate refuses
oversize codes, and that the (honestly sparse) cyclic search recovers exactly the Golay.
"""

import numpy as np
import pytest

from qudit_qec import distill_strange as ds


@pytest.fixture(scope="module")
def golay_gen():
    """A self-orthogonal cyclic [11,5] generator over F_3 (the ternary QR Golay)."""
    codes = ds.self_orthogonal_cyclic_codes(11)
    assert codes, "expected self-orthogonal cyclic codes at n=11"
    return codes[0]


# --------------------------------------------------------------------------- #
# Cyclic-code genotype
# --------------------------------------------------------------------------- #
def test_cyclic_codes_basic():
    codes = ds.cyclic_codes(7)
    assert all(np.asarray(G).shape[1] == 7 for G in codes)
    assert len(codes) >= 2                       # at least the trivial divisors


def test_cyclic_requires_coprime_length():
    with pytest.raises(ValueError):
        ds.cyclic_codes(9)                       # gcd(9,3) != 1


def test_self_orthogonal_golay_exists(golay_gen):
    # the two equivalent ternary QR Golay codes are self-orthogonal [11,5]
    so = ds.self_orthogonal_cyclic_codes(11)
    assert len(so) == 2
    assert np.asarray(golay_gen).shape == (5, 11)
    assert ds._is_self_orthogonal(golay_gen)


def test_sparse_family_no_distiller_at_5_7():
    # honest sparsity: no self-orthogonal cyclic distiller at small n
    assert ds.self_orthogonal_cyclic_codes(5) == []
    assert ds.self_orthogonal_cyclic_codes(7) == []


# --------------------------------------------------------------------------- #
# CSS build + validity
# --------------------------------------------------------------------------- #
def test_build_strange_css_golay(golay_gen):
    code = ds.build_strange_css(golay_gen)
    assert code.num_qudits == 11 and code.dimension == 1


def test_build_rejects_non_self_orthogonal():
    G = np.array([[1, 0, 0], [0, 1, 0]])         # not self-orthogonal over F_3
    with pytest.raises(ValueError):
        ds.build_strange_css(G)


# --------------------------------------------------------------------------- #
# Evaluation (cheap screen)
# --------------------------------------------------------------------------- #
def test_evaluate_golay_cheap(golay_gen):
    r = ds.evaluate_strange_candidate(golay_gen)    # no threshold
    assert r.distills is True and r.sound is True
    assert r.nu == 3                                # cubic lower bound (odd n)
    assert r.k == 1 and r.n == 11
    assert r.threshold is None                      # cheap screen defers threshold
    assert r.d_lb == 6                              # min nonzero stabilizer weight


def test_evaluate_rejects_non_self_orthogonal():
    G = np.array([[1, 0, 0], [0, 1, 0]])
    r = ds.evaluate_strange_candidate(G)
    assert r.rejected is True and "self-orthogonal" in r.reason


def test_evaluate_qutrit_only():
    r = ds.evaluate_strange_candidate(np.zeros((1, 4), dtype=int), p=4)
    assert r.rejected is True and "qutrit-only" in r.reason


def test_compute_gate_refuses_oversize(golay_gen):
    r = ds.evaluate_strange_candidate(golay_gen, max_enum=100)   # 3^10 >> 100
    assert r.rejected is True and "lenore" in r.reason


def test_cheap_conditions_match_known_golay(golay_gen):
    # the cheap Fraction screen must agree with the distillation verdict
    we_distills, nu_hint, info = ds._cheap_strange_conditions(
        *_golay_AB(golay_gen), n=11)
    assert we_distills is True and nu_hint == 3
    assert info["B(-1/2)!=0"] and info["3A'+B'(-1/2)=0"]


def _golay_AB(G):
    from qudit_qec.distillation import weight_enumerator
    code = ds.build_strange_css(G)
    we = weight_enumerator(code, p=3)
    return we.A, we.B


def test_cheap_conditions_require_cond0():
    # cond2 (3A'+B'=0) alone does NOT imply distillation: nu>=2 also needs cond0
    # (3A+B(-1/2)=0). A case with cond2 True but cond0 False must screen as NON-distilling
    # and report no nu (regression for the cheap-screen false positive / nu overstatement).
    A = [1, 0, 4]
    B = [1, 12, 0]          # 3A+B has P'(-1/2)=0 but P(-1/2)=1 != 0  (cond0 fails)
    distills, nu_hint, info = ds._cheap_strange_conditions(A, B, n=4)
    assert info["3A'+B'(-1/2)=0"] is True and info["3A+B(-1/2)=0"] is False
    assert distills is False and nu_hint is None


# --------------------------------------------------------------------------- #
# Catalog + search
# --------------------------------------------------------------------------- #
def test_search_recovers_golay():
    cat = ds.search_strange_cyclic(5, 13)
    assert len(cat) >= 1
    params = {(r.n, r.k) for r in cat.codes()}
    assert (11, 1) in params
    assert all(r.distills for r in cat.codes())
    assert all(r.nu and r.nu >= 2 for r in cat.codes())


def test_catalog_best_by_nu_ranking():
    cat = ds.search_strange_cyclic(11, 11)
    ranked = cat.best_by_nu()
    nus = [r.nu for r in ranked]
    assert nus == sorted(nus, reverse=True)         # nu descending


# --------------------------------------------------------------------------- #
# Opt-in exact path (slower: one sympy threshold solve)
# --------------------------------------------------------------------------- #
def test_golay_threshold_refine():
    cat = ds.search_strange_cyclic(11, 11)
    cat.refine_thresholds(top=1)
    best = cat.best_by_nu(1)[0]
    assert best.nu == 3                              # exact cubic
    assert best.threshold == pytest.approx(0.387, abs=0.02)
