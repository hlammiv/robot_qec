"""Tests for qudit_qec.distill_discovery -- the prime-p MSD discovery arm.

Validates that (1) the constructive family reproduces the arXiv:2403.06228 codes and
generalizes to every prime, (2) the validity gate and operators never admit a
non-triorthogonal matrix to the catalog, (3) distance provenance is tracked honestly
(known vs upper-bound vs trusted), and (4) the bounded search is deterministic,
compute-capped, and emits only genuine distillation codes.
"""

import math

import numpy as np
import pytest

from qudit_qec import distill_discovery as dd


# --------------------------------------------------------------------------- #
# Constructive family
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("p,m,k,n,exp_gamma", [
    (3, 3, 7, 20, 1.515),   # paper headline [[20,7,2]]_3
    (3, 2, 4, 14, 1.807),   # paper [[14,4,2]]_3
    (3, 2, 5, 13, 1.379),   # extended frontier member of the same family
])
def test_family_reproduces_paper(p, m, k, n, exp_gamma):
    H = dd.triortho_family(p, m, k)
    r = dd.evaluate_distill_candidate(H, p, d=2, d_status="known", op="family")
    assert r.n == n and r.k == k and r.d == 2
    assert r.distills is True and r.trusted is True and r.sound is True
    assert r.gamma == pytest.approx(exp_gamma, abs=0.01)


@pytest.mark.parametrize("p", [3, 5, 7])
def test_family_valid_every_prime(p):
    # k = 1 .. p*m-1 all give genuine triorthogonal matrices with k magic rows
    for k in range(1, p):           # m=1 -> valid k up to p-1
        H = dd.triortho_family(p, 1, k)
        r = dd.evaluate_distill_candidate(H, p, d=2)
        assert r.distills is True
        assert r.k == k and r.n == p * p - k


def test_family_rejects_bad_params():
    with pytest.raises(ValueError):
        dd.triortho_family(4, 1, 1)          # non-prime
    with pytest.raises(ValueError):
        dd.triortho_family(3, 1, 3)          # k > p*m - 1


# --------------------------------------------------------------------------- #
# Validity gate
# --------------------------------------------------------------------------- #
def test_gate_rejects_random_matrix():
    rng = np.random.default_rng(0)
    M = rng.integers(0, 3, size=(4, 14))
    r = dd.evaluate_distill_candidate(M, 3, d=2)
    assert r.rejected is True and r.distills is False


def test_gate_respects_max_n_cap():
    H = dd.triortho_family(3, 3, 7)          # n=20
    r = dd.evaluate_distill_candidate(H, 3, d=2, max_n=10)
    assert r.rejected is True and "compute cap" in r.reason


def test_unknown_distance_is_not_trusted():
    H = dd.triortho_family(3, 2, 4)
    r = dd.evaluate_distill_candidate(H, 3, d=None)
    assert r.distills is True            # it IS a valid triorthogonal matrix
    assert r.gamma is None and r.trusted is False
    assert any("distance pending" in c for c in r.caveats)


def test_non_prime_flagged_not_sound():
    M = np.array([[0, 1, 2, 3], [1, 1, 0, 0]])   # GF(4) matrix
    r = dd.evaluate_distill_candidate(M, 4, d=2)
    assert r.sound is False


# --------------------------------------------------------------------------- #
# Operators
# --------------------------------------------------------------------------- #
def test_permute_preserves_validity_and_key():
    H = dd.triortho_family(3, 2, 4)
    rng = np.random.default_rng(1)
    Hp = dd.permute_columns(H, rng.permutation(H.shape[1]), 3)
    assert dd.evaluate_distill_candidate(Hp, 3, d=2).distills is True
    # column-multiset key is permutation-invariant
    assert dd.genotype_key(H, 3) == dd.genotype_key(Hp, 3)


def test_direct_sum_provably_triorthogonal():
    H = dd.triortho_family(3, 1, 2)          # [[7,2,2]]
    Hds = dd.direct_sum(H, H, 3)
    r = dd.evaluate_distill_candidate(Hds, 3, d=2)
    assert r.distills is True
    assert r.n == 14 and r.k == 4            # n and magic rows add


def test_scale_columns_distance_preserving():
    H = dd.triortho_family(3, 2, 4)
    Hs = dd.scale_columns(H, [2] * H.shape[1], 3)   # uniform nonzero scale
    r = dd.evaluate_distill_candidate(Hs, 3, d=2, d_status="known")
    # uniform scaling keeps triorthogonality here; if valid, distance is preserved
    if r.distills:
        assert r.trusted is True


def test_scale_rejects_zero():
    with pytest.raises(ValueError):
        dd.scale_columns(dd.triortho_family(3, 1, 1), [0, 1, 1, 1, 1, 1, 1, 1], 3)


def test_puncture_marks_distance_upper():
    H = dd.triortho_family(3, 3, 7)          # [[20,7,2]]
    # find a puncture that stays triorthogonal, then check it is flagged upper-bound
    rng = np.random.default_rng(3)
    found = False
    for _ in range(50):
        Hp = dd.puncture_columns(H, [int(rng.integers(H.shape[1]))], 3)
        r = dd.evaluate_distill_candidate(Hp, 3, d=2, d_status="upper")
        if r.distills:
            assert r.d_status == "upper" and r.trusted is False
            assert any("optimistic lower bound" in c for c in r.caveats)
            found = True
            break
    assert found, "expected at least one triorthogonality-preserving single puncture"


# --------------------------------------------------------------------------- #
# Catalog
# --------------------------------------------------------------------------- #
def test_catalog_dedup_and_better():
    cat = dd.DistillCatalog()
    H = dd.triortho_family(3, 2, 4)
    untrusted = dd.evaluate_distill_candidate(H, 3, d=2, d_status="upper")
    trusted = dd.evaluate_distill_candidate(H, 3, d=2, d_status="known")
    assert cat.add(untrusted) is True        # new genotype
    assert cat.add(trusted) is False         # same key -> dedup
    assert cat.codes()[0].trusted is True    # trusted record retained


def test_catalog_pareto_and_gamma_ranking():
    cat = dd.DistillCatalog()
    for k in range(1, 5):
        H = dd.triortho_family(3, 2, k)
        cat.add(dd.evaluate_distill_candidate(H, 3, d=2, d_status="known"))
    ranked = cat.best_by_gamma()
    gammas = [r.gamma for r in ranked]
    assert gammas == sorted(gammas)          # ascending (lower gamma better)
    assert all(r.d is not None for r in cat.pareto_front())


# --------------------------------------------------------------------------- #
# Search
# --------------------------------------------------------------------------- #
def test_search_emits_only_valid_codes():
    cat = dd.search_distill(primes=(3,), m_range=(1, 2), iters=150, seed=1, max_n=30)
    assert len(cat) > 0
    assert all(r.distills for r in cat.codes())       # every catalog entry is genuine
    assert all(r.n <= 30 for r in cat.codes())        # compute cap respected


def test_search_deterministic():
    a = dd.search_distill(primes=(3,), m_range=(1, 2), iters=100, seed=7, max_n=30)
    b = dd.search_distill(primes=(3,), m_range=(1, 2), iters=100, seed=7, max_n=30)
    assert {r.key for r in a.codes()} == {r.key for r in b.codes()}


def test_search_finds_family_frontier():
    # the search must at least catalog the validated [[20,7,2]]_3 family member
    cat = dd.search_distill(primes=(3,), m_range=(1, 3), iters=50, seed=0, max_n=24)
    params = {(r.n, r.k, r.d) for r in cat.codes()}
    assert (20, 7, 2) in params


# --------------------------------------------------------------------------- #
# Trusted-distance route (opt-in; tiny case only)
# --------------------------------------------------------------------------- #
def test_trusted_distance_small_member():
    H = dd.triortho_family(3, 1, 1)          # [[8,1,2]]_3
    r = dd.evaluate_distill_candidate(H, 3, distance="trusted")
    assert r.d == 2 and r.d_status == "trusted" and r.trusted is True
    assert r.gamma == pytest.approx(math.log(8) / math.log(2), abs=1e-6)
