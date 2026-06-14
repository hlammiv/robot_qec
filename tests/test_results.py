"""Tests for qudit_qec.results — catalog + Pareto with coefficient-aware dedup.

Includes the Phase 3 milestone: two GF(3) codes differing only in a coefficient are
NOT merged (the silent-merge hazard the reference's exponent-only keys would hit).
"""

from qudit_qec import CodeCatalog, EvalResult, evaluate_candidate


def _mk(n, k, d, *, q=3, key=None, trusted=True, fom=None):
    fom = fom if fom is not None else (k * d * d / n if d is not None else None)
    return EvalResult(1, 1, q, n, k, d, "milp", fom, trusted, False, False,
                      (), (), key or (n, k, d), False, "", {})


def test_milestone_coeff_differing_codes_not_merged():
    A1 = [(3, 0, 1), (0, 1, 1), (0, 2, 1)]          # x^3 + y + y^2
    A2 = [(3, 0, 2), (0, 1, 2), (0, 2, 2)]          # 2*(x^3 + y + y^2), distinct over GF(3)
    B = [(0, 3, 1), (1, 0, 1), (2, 0, 1)]
    r1 = evaluate_candidate(6, 6, A1, B, field=3, distance="none")
    r2 = evaluate_candidate(6, 6, A2, B, field=3, distance="none")
    assert not r1.rejected and not r2.rejected
    assert r1.key != r2.key  # coefficients distinguish the keys
    cat = CodeCatalog()
    assert cat.add(r1) is True
    assert cat.add(r2) is True
    assert len(cat) == 2


def test_duplicate_key_keeps_better():
    cat = CodeCatalog()
    untrusted = _mk(32, 2, 4, key="X", trusted=False)
    trusted = _mk(32, 2, 4, key="X", trusted=True)
    assert cat.add(untrusted) is True
    assert cat.add(trusted) is False  # same key -> not a new code
    assert len(cat) == 1
    assert cat.codes()[0].trusted is True  # the trusted record is retained


def test_rejected_results_are_not_added():
    cat = CodeCatalog()
    r = EvalResult(1, 1, 3, 0, 0, None, "none", None, False, None, False,
                   (), (), None, True, "bad", {})
    assert cat.add(r) is False
    assert len(cat) == 0


def test_pareto_front():
    cat = CodeCatalog()
    for r in (_mk(72, 12, 12, key="dom"),   # dominates "sub"
              _mk(72, 8, 8, key="sub"),     # dominated (same n, smaller k and d)
              _mk(144, 20, 6, key="oth")):  # not dominated (bigger n, bigger k)
        cat.add(r)
    front = {r.key for r in cat.pareto_front()}
    assert front == {"dom", "oth"}


def test_best_by_fom():
    cat = CodeCatalog()
    cat.add(_mk(72, 12, 12, key="a"))  # FOM = 24
    cat.add(_mk(32, 2, 4, key="b"))    # FOM = 1
    best = cat.best_by_fom(1)
    assert len(best) == 1 and best[0].key == "a"
