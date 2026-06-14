"""Catalog and Pareto front for discovered qudit codes, with coefficient-aware dedup.

``CodeCatalog`` accumulates :class:`~qudit_qec.evaluator.EvalResult` records keyed by
the **coefficient-aware** canonical genotype key (``genotype.pair_key``, which also
carries the lattice and field). Two codes that share exponents but differ in a GF(q)
coefficient therefore get distinct keys and are **never merged** -- the silent-merge
hazard the reference pipeline's exponent-only keys would hit over GF(q).

The catalog records the field order ``q`` per entry and exposes the Pareto front
(non-dominated in qudit count ``n``, encoded dimension ``k``, and distance ``d``) and
a best-by-FOM ranking.
"""

from __future__ import annotations

from dataclasses import dataclass, field as _dc_field


def _better(a, b) -> bool:
    """True if result ``a`` should replace ``b`` for the same genotype key.

    Prefer a trusted distance over an untrusted one, then a larger FOM, then a
    larger (definitely-known) distance.
    """
    if a.trusted != b.trusted:
        return a.trusted
    fa = a.fom if a.fom is not None else -1.0
    fb = b.fom if b.fom is not None else -1.0
    if fa != fb:
        return fa > fb
    da = a.d if a.d is not None else -1
    db = b.d if b.d is not None else -1
    return da > db


def _dominates(a, b) -> bool:
    """True if code ``a`` dominates code ``b``: no worse on (n, k, d), better on one.

    Smaller ``n`` (fewer physical qudits), larger ``k``, larger ``d`` are better.
    Only defined for codes with a known distance.
    """
    if a.d is None or b.d is None:
        return False
    no_worse = a.n <= b.n and a.k >= b.k and a.d >= b.d
    strictly = a.n < b.n or a.k > b.k or a.d > b.d
    return no_worse and strictly


@dataclass
class CodeCatalog:
    """A deduplicated collection of evaluated codes, keyed by coefficient-aware genotype."""

    by_key: dict = _dc_field(default_factory=dict)

    def add(self, result) -> bool:
        """Add ``result``; return True if it was a new code, False if a duplicate key.

        Rejected results and results without a key are ignored (returns False). On a
        duplicate key the better record (see :func:`_better`) is retained.
        """
        if result.rejected or result.key is None:
            return False
        existing = self.by_key.get(result.key)
        if existing is None:
            self.by_key[result.key] = result
            return True
        if _better(result, existing):
            self.by_key[result.key] = result
        return False

    def __len__(self) -> int:
        return len(self.by_key)

    def codes(self) -> list:
        """All catalogued (deduplicated) results."""
        return list(self.by_key.values())

    def with_distance(self) -> list:
        """Catalogued results that have a known distance."""
        return [r for r in self.by_key.values() if r.d is not None]

    def pareto_front(self) -> list:
        """Codes not dominated by any other (on n, k, d). Distance-known codes only."""
        codes = self.with_distance()
        return [a for a in codes if not any(_dominates(b, a) for b in codes if b is not a)]

    def best_by_fom(self, top: int | None = None) -> list:
        """Distance-known codes sorted by FOM descending (optionally top-``top``)."""
        ranked = sorted(self.with_distance(), key=lambda r: r.fom, reverse=True)
        return ranked[:top] if top is not None else ranked
