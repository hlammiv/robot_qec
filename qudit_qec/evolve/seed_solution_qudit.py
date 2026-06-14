"""Seed program for LLM-guided evolutionary search of **qudit** BB codes.

This is the file OpenEvolve loads and mutates. It exposes:

``TARGET_LATTICES`` : list[tuple[int, int]]
    Candidate ``(ell, m)`` lattice dimensions to search.

``STRUCTURAL_SEEDS`` : list[dict]
    Field-agnostic ``(ell, m, A_terms, B_terms)`` polynomial *structures* (all
    coefficient 1) drawn from known good qubit BB codes. Over GF(q) their exact
    ``[[n,k,d]]`` differs and is computed by the evaluator; here they serve as
    starting points and a k>0 safety net (coefficient-1 trinomials are non-degenerate
    over the small prime fields we target).

``generate_candidates(ell, m)`` : the function the LLM evolves (between the
    ``# EVOLVE-BLOCK-START`` / ``# EVOLVE-BLOCK-END`` markers). It returns a list of
    ``(A_terms, B_terms)`` candidate pairs, where each term is a
    ``(x_exp, y_exp, coeff)`` triple with ``coeff in 1..q-1`` (a 2-tuple
    ``(x_exp, y_exp)`` is also accepted and means coeff 1).

The qudit field order ``q`` is read from the ``QCODE_FIELD`` environment variable
(set by ``run_evolution.py``); this is what lets the evolved function explore the
**coefficient axis** that does not exist over F_2. Note ``q`` is a *prime* for the
direct search; composite (square-free) dimensions are reached via CRT at evaluation
time. See ``docs/`` for the algebra.
"""

from __future__ import annotations

import os

# GF(q) field order for this campaign (prime). Set by run_evolution.py via env.
FIELD = int(os.environ.get("QCODE_FIELD", "2"))


# Field-agnostic structural seeds (coefficient 1), from known qubit BB codes.
STRUCTURAL_SEEDS = [
    # x/y-swap "gross" structure: A = x^3 + y + y^2, B = y^3 + x + x^2
    {"ell": 6, "m": 6,
     "A_terms": [(3, 0, 1), (0, 1, 1), (0, 2, 1)],
     "B_terms": [(0, 3, 1), (1, 0, 1), (2, 0, 1)]},
    {"ell": 12, "m": 6,
     "A_terms": [(3, 0, 1), (0, 1, 1), (0, 2, 1)],
     "B_terms": [(0, 3, 1), (1, 0, 1), (2, 0, 1)]},
    {"ell": 9, "m": 6,
     "A_terms": [(3, 0, 1), (0, 1, 1), (0, 2, 1)],
     "B_terms": [(0, 3, 1), (1, 0, 1), (2, 0, 1)]},
    # constant-monomial (univariate/HGP) structure: A = 1 + y + y^2, B = 1 + x^c + x^2c
    {"ell": 15, "m": 12,
     "A_terms": [(0, 0, 1), (0, 1, 1), (0, 2, 1)],
     "B_terms": [(0, 0, 1), (5, 0, 1), (10, 0, 1)]},
    # GF(3)-verified discovery: [[72,6,8]]_3 (mixed-monomial x^5*y term), the
    # highest CERTIFIED-EXACT FOM (5.33) from the GF(3) baseline sweep + independent
    # verification (see results/gf3_report.md). A good qutrit-specific seed.
    {"ell": 6, "m": 6,
     "A_terms": [(0, 2, 1), (3, 0, 1), (5, 1, 1)],
     "B_terms": [(0, 3, 1), (1, 0, 1), (2, 0, 1)]},
]


# Candidate lattice dimensions (kept small-to-moderate; the evaluator can subset).
TARGET_LATTICES = [
    (6, 6), (9, 6), (12, 6), (6, 12), (12, 12), (15, 6),
    (15, 12), (18, 6), (24, 6), (30, 6),
]


# EVOLVE-BLOCK-START
def generate_candidates(ell, m):
    """Generate candidate ``(A_terms, B_terms)`` polynomial pairs over GF(q).

    Each term is ``(x_exp, y_exp, coeff)`` with ``coeff in 1..q-1``. The LLM
    evolves this function to discover algebraic patterns -- in both the exponents
    *and* the GF(q) coefficients -- that yield high k (logical qudits) and high d
    (distance), maximizing FOM = k*d^2/n.

    Args:
        ell, m: cyclic group orders for x and y (n = 2*ell*m).

    Returns:
        list of (A_terms, B_terms) pairs.
    """
    q = FIELD
    nonzero_coeffs = list(range(1, q))  # GF(q)* = {1, ..., q-1}
    candidates = []
    seen = set()

    def _add(a_terms, b_terms):
        key = (tuple(sorted(map(tuple, a_terms))), tuple(sorted(map(tuple, b_terms))))
        if key not in seen:
            seen.add(key)
            candidates.append(([tuple(t) for t in a_terms], [tuple(t) for t in b_terms]))

    # Strategy 1: x/y-swap symmetric, all coefficients 1 (the proven structure).
    # A = x^a + y^b + y^(2b),  B = y^d + x^e + x^(2e).
    max_a = ell // 2 + 1
    max_b = m // 2 + 1
    for a in range(1, max_a):
        for b in range(1, max_b):
            c = (2 * b) % m
            A = [(a, 0, 1), (0, b, 1), (0, c, 1)]
            if len({(t[0], t[1]) for t in A}) != 3:
                continue
            for d in range(1, max_b):
                for e in range(1, max_a):
                    f = (2 * e) % ell
                    B = [(0, d, 1), (e, 0, 1), (f, 0, 1)]
                    if len({(t[0], t[1]) for t in B}) == 3:
                        _add(A, B)

    # Strategy 2: the QUDIT coefficient axis -- vary one coefficient over GF(q)*.
    # For q = 2 this is a no-op (only coeff 1 exists); for q > 2 it explores codes
    # unreachable by the qubit search. We perturb the middle monomial of A.
    if q > 2:
        for a in range(1, max_a):
            for b in range(1, max_b):
                c = (2 * b) % m
                if len({(a, 0), (0, b), (0, c)}) != 3:
                    continue
                for coeff in nonzero_coeffs[1:]:  # skip coeff 1 (Strategy 1)
                    A = [(a, 0, 1), (0, b, coeff), (0, c, 1)]
                    B = [(0, 1, 1), (1, 0, 1), (2, 0, 1)]
                    if 1 < ell and 2 < ell:
                        _add(A, B)

    # Strategy 3: structural seeds at this lattice + small exponent perturbations.
    for spec in STRUCTURAL_SEEDS:
        if spec["ell"] != ell or spec["m"] != m:
            continue
        base_A, base_B = spec["A_terms"], spec["B_terms"]
        _add(base_A, base_B)
        for delta in (-1, 1):
            for i in range(len(base_A)):
                pa = [list(t) for t in base_A]
                pa[i][0] = (pa[i][0] + delta) % ell
                pa_t = [tuple(t) for t in pa]
                if len({(t[0], t[1]) for t in pa_t}) == len(pa_t):
                    _add(pa_t, base_B)

    # Strategy 4: small-exponent x/y-swap search (coeff 1), independent exponents.
    mx = min(ell, 5)
    my = min(m, 5)
    for a in range(1, mx):
        for b in range(0, my):
            for c in range(b + 1, my):
                A = [(a, 0, 1), (0, b, 1), (0, c, 1)]
                if len({(t[0], t[1]) for t in A}) != 3:
                    continue
                for d in range(1, my):
                    for e in range(0, mx):
                        for f in range(e + 1, mx):
                            B = [(0, d, 1), (e, 0, 1), (f, 0, 1)]
                            if len({(t[0], t[1]) for t in B}) == 3:
                                _add(A, B)

    return candidates
# EVOLVE-BLOCK-END
