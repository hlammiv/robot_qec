"""Single-qudit gate-set universality: the p-adic *sector-coverage* check over Z_d.

Decidable instantiation (arXiv:2512.20787, Borda-Rincon-Galindo) of *what extra gate
makes ``{Clifford + G}`` universal* on a single d-level qudit, beyond the prime case.

Background.  The single-qudit Clifford group is a maximal finite subgroup of
``U(d)`` modulo phase **iff d is prime** (Nebe-Rains-Sloane math/0001038;
Borda-Rincon-Galindo arXiv:2512.20787).  For prime ``d`` the adjoint representation
on ``sl(d)`` is irreducible, so *any* non-Clifford gate forces universality.  For
``d = p^m`` (m>1) it is **reducible**: ``sl(d) = (+)_k W_k`` with

    W_k = span{ P_{a,b} : (a,b) != (0,0),  v_p(gcd(a,b,d)) = k },   k = 0..m-1,

and Clifford acts transitively *within* each ``W_k`` but cannot mix them.  A gate
``G`` can restore irreducibility (a *necessary* condition for ``{Clifford+G}`` to be
universal) only if its adjoint action ``Ad_G`` **couples** the sectors so the
sector graph on ``{0..m-1}`` is connected.  (The full Sawicki-Karnas criterion adds
a density/threshold condition on top; this module checks the decidable
sector-coupling part -- failure => definitely not universal; success => candidate.)

For a single-qudit **diagonal** gate ``G = diag(exp(2*pi*i*phase_j))`` (the relevant
case for T-type and CCZ-induced gates), ``G X^a G^\\dagger = D_a X^a`` with
``D_a = diag(exp(2*pi*i*(phase_j - phase_{j-a})))``; expanding
``D_a = sum_c hatD_a(c) Z^c`` gives ``Ad_G : (a,b) -> { (a, b+c) : hatD_a(c) != 0 }``.
So sector-coupling is read off the Fourier supports of the ``D_a`` -- an O(d^2)
computation.  Pure ``Z_d`` arithmetic, tiny compute.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Sequence

import numpy as np

from .field_utils import prime_factorization

_TOL = 1e-9


def _prime_power(d: int) -> tuple[int, int]:
    pf = prime_factorization(d)
    if len(pf) != 1:
        raise NotImplementedError(
            f"sector analysis is defined for prime-power d=p^m; got d={d} "
            f"= {dict(pf)} (composite). Use the CRT/inter-factor route (docs/08)."
        )
    p, m = next(iter(pf.items()))
    return p, m


def _valuation(x: int, p: int) -> int:
    """p-adic valuation v_p(x) for x>0; large sentinel for x==0."""
    if x == 0:
        return 1 << 30
    k = 0
    while x % p == 0:
        x //= p
        k += 1
    return k


def sector_of(a: int, b: int, d: int, p: int, m: int) -> int:
    """Sector index k = v_p(gcd(a,b,d)) in 0..m-1 for (a,b) != (0,0)."""
    return min(_valuation(math.gcd(math.gcd(a % d, b % d), d), p), m - 1)


def _diagonal_adjoint_supports(phases: Sequence[float], d: int) -> dict[int, set[int]]:
    """For diagonal G=diag(exp(2*pi*i*phases[j])), return a -> Fourier support C_a of
    ``D_a = diag(exp(2*pi*i*(phases[j]-phases[(j-a)%d])))`` (the Z^c content of
    ``G X^a G^\\dagger``)."""
    phases = np.asarray(phases, dtype=float)
    j = np.arange(d)
    supports: dict[int, set[int]] = {}
    for a in range(d):
        diag = np.exp(2j * np.pi * (phases - phases[(j - a) % d]))
        # Fourier coefficients over Z_d: hatD_a(c) = (1/d) sum_j diag_j omega^{-c j}
        fft = np.fft.fft(diag) / d  # numpy fft uses exp(-2pi i c j / d): exactly hatD_a(c)
        supports[a] = {int(c) for c in range(d) if abs(fft[c]) > _TOL}
    return supports


@dataclass
class SectorCoverage:
    """Result of the p-adic sector-coupling test for a single-qudit gate over Z_d."""

    d: int
    p: int
    m: int
    couples_all_sectors: bool   # sector graph on {0..m-1} connected (necessary for univ.)
    sector_edges: list          # coupled sector pairs (k, k')
    is_clifford_like: bool       # True if Ad_G never leaves a sector (preserves stratification)
    note: str = ""
    caveats: list = field(default_factory=list)


def sector_coverage_diagonal(phases: Sequence[float], d: int) -> SectorCoverage:
    """Does a single-qudit DIAGONAL gate over Z_d couple all p-adic Pauli sectors?

    ``phases`` are the diagonal phases in *turns* (G_jj = exp(2*pi*i*phases[j])).
    Returns a :class:`SectorCoverage`. ``couples_all_sectors`` is a NECESSARY
    condition for ``{Clifford + G}`` to be universal on a native ``d=p^m`` qudit;
    if False, ``{Clifford+G}`` is definitely not universal (G stays within the
    reducible sector structure).
    """
    p, m = _prime_power(d)
    supports = _diagonal_adjoint_supports(phases, d)

    edges: set[tuple[int, int]] = set()
    mixes = False
    for a in range(d):
        for b in range(d):
            if a == 0 and b == 0:
                continue
            k0 = sector_of(a, b, d, p, m)
            for c in supports[a]:
                b2 = (b + c) % d
                if a == 0 and b2 == 0:
                    continue
                k1 = sector_of(a, b2, d, p, m)
                if k0 != k1:
                    mixes = True
                    edges.add((min(k0, k1), max(k0, k1)))

    # connectivity of {0..m-1} under the coupling edges
    parent = list(range(m))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for k0, k1 in edges:
        parent[find(k0)] = find(k1)
    connected = len({find(k) for k in range(m)}) == 1 if m > 1 else True

    caveats = [
        "Sector-coupling is a NECESSARY condition for universality (the Sawicki-Karnas "
        "criterion adds a density/threshold check on top); it is not by itself sufficient.",
    ]
    if m == 1:
        caveats.append("d is prime: Clifford is already irreducible, so any non-Clifford "
                       "gate is universal regardless of this sector test.")
    return SectorCoverage(d=d, p=p, m=m, couples_all_sectors=bool(connected and (m == 1 or mixes)),
                          sector_edges=sorted(edges), is_clifford_like=not mixes,
                          note=f"d={p}^{m}, {m} sectors, {'mixes' if mixes else 'no mixing'}",
                          caveats=caveats)


# --------------------------------------------------------------------------- #
# Reference single-qudit diagonal gates (phases in turns)
# --------------------------------------------------------------------------- #
def identity_phases(d: int) -> list[float]:
    return [0.0] * d


def pauli_Z_power_phases(d: int, t: int) -> list[float]:
    """Z^t = diag(exp(2*pi*i*t*j/d)) -- a Pauli (Clifford); should NOT couple sectors."""
    return [(t * j % d) / d for j in range(d)]


def clifford_S_phases(d: int) -> list[float]:
    """The Clifford phase gate S (quadratic phase), parity-correct so it is periodic
    mod d: for odd d, ``S = diag(omega^{(2^{-1}) j(j-1)})`` with ``2^{-1}=(d+1)/2``;
    for even d, ``S = diag(omega_{2d}^{j^2}) = diag(exp(2*pi*i*j^2/(2d)))``. Its adjoint
    sends ``X^a -> (phase) Z^a`` (a single Fourier frequency), so it PRESERVES the
    p-adic stratification and must NOT couple sectors -- the key negative control."""
    if d % 2 == 1:
        inv2 = (d + 1) // 2  # 2^{-1} mod d for odd d
        return [((inv2 * j * (j - 1)) % d) / d for j in range(d)]
    return [(j * j % (2 * d)) / (2.0 * d) for j in range(d)]


def level_bump_phases(d: int, j0: int = 1, num: int = 1, den: int = 4) -> list[float]:
    """A single-level phase bump (phase num/den on |j0>, 0 elsewhere) -- a manifestly
    non-quadratic, non-Clifford diagonal that SHOULD couple sectors. Positive control
    validating that the checker can return couples_all_sectors=True at a prime power."""
    return [(num / den) if j == j0 else 0.0 for j in range(d)]


def ccz_single_qudit_reduction_phases(d: int, b: int, c: int) -> list[float]:
    """The single-qudit gate induced by CCZ (CCZ|j,k,l> = omega^{jkl}|j,k,l>) when the
    two control qudits are fixed to computational |b>,|c>: this is Z^{b*c}, a PAULI.
    So every computational-basis single-qudit reduction of CCZ is Clifford -- evidence
    that CCZ does NOT supply a single-qudit sector-coupling resource by such reductions
    (the magic-state-injection effective gate is a separate, open question)."""
    return pauli_Z_power_phases(d, (b * c) % d)
