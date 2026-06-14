"""Structural analysis of a stabilizer code: Tanner-graph decomposability.

A code is *decomposable* (a direct sum) when its qudits partition into groups with
no stabilizer spanning two groups -- the Tanner graph disconnects. Such a code
offers no error-correction advantage over its independent pieces (the paper's
``[[288,24,12]] = gross + gross`` is the canonical example), so the search should
detect and discount them.

The check is purely structural -- it depends only on the *support* (nonzero
positions) of the stabilizers, not on the GF(q) coefficient values -- so it is
field-agnostic and works unchanged for qudits.
"""

from __future__ import annotations

import numpy as np


def qudit_check_incidence(code) -> np.ndarray:
    """Boolean ``(num_checks, n)`` matrix: does check ``r`` touch qudit ``j``?

    A qudit is touched if the stabilizer has X- or Z-support on it. Uses the
    symplectic stabilizer matrix ``code.matrix = [X | Z]``; the X/Z block order is
    irrelevant to support.
    """
    matrix = np.asarray(code.matrix, dtype=int)
    n = code.num_qudits
    x_part = matrix[:, :n]
    z_part = matrix[:, n : 2 * n]
    return (x_part != 0) | (z_part != 0)


def connected_components(code) -> list[list[int]]:
    """Partition qudit indices by Tanner-graph connectivity.

    Two qudits are connected if some stabilizer touches both; components are the
    transitive closure. Returns a list of sorted qudit-index lists.
    """
    incidence = qudit_check_incidence(code)
    n = incidence.shape[1]
    parent = list(range(n))

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for r in range(incidence.shape[0]):
        cols = np.nonzero(incidence[r])[0]
        for j in cols[1:]:
            union(int(cols[0]), int(j))

    groups: dict[int, list[int]] = {}
    for j in range(n):
        groups.setdefault(find(j), []).append(j)
    return [sorted(v) for v in groups.values()]


def num_connected_components(code) -> int:
    """Number of Tanner-graph connected components among the qudits."""
    return len(connected_components(code))


def is_decomposable(code) -> bool:
    """True if the code is a direct sum (its Tanner graph has > 1 component)."""
    return num_connected_components(code) > 1
