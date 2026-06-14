"""Field-aware distance dispatcher for qudit CSS codes.

Combines the three backends into one trust-aware result:

1. a **cheap GUF/BP-OSD pre-filter** upper bound (``distance.decoder_bound``) --
   fast but loose (~3x for q>2), used only to screen, never trusted;
2. the **prime-q mod-q MILP** (``distance_milp.compute_distance_milp``) -- the
   *trusted* distance signal when q is prime and all logicals certify optimal;
3. optional **exact enumeration** (``distance.compute_distance_exact``) -- for
   small codes, to corroborate (or supply, for prime-power q where MILP is invalid).

The trust gate is the key design point from the scoping: a result is ``trusted``
only when MILP-certified or exact-enumerated. The loose GUF bound alone is never
trusted (a 3x-loose d inflates ``FOM = k*d^2/n`` ~9x).
"""

from __future__ import annotations

from dataclasses import dataclass, field as _dc_field

from .distance import compute_distance_exact, decoder_bound
from .distance_milp import compute_distance_milp
from .field_utils import is_prime


@dataclass
class DistanceResult:
    """Outcome of a distance computation.

    Attributes:
        d: the best distance value found (exact if ``trusted``, else an upper bound).
        method: which backend produced ``d`` (``'milp'``, ``'milp_incumbent'``,
            ``'exact'``, or ``'guf_bound'``).
        trusted: True only if ``d`` is MILP-certified or exact-enumerated.
        exact: alias of ``trusted`` (``d`` is the proven minimum distance).
        upper_bound: the cheap GUF/BP-OSD pre-filter bound (always an upper bound).
        q: the field order.
        details: backend-specific diagnostics.
    """

    d: int
    method: str
    trusted: bool
    exact: bool
    upper_bound: int
    q: int
    details: dict = _dc_field(default_factory=dict)


def code_distance(
    code,
    q: int | None = None,
    *,
    prefilter_trials: int = 50,
    milp_timeout_per_logical: float = 30,
    milp_total_timeout: float = 180,
    want_exact: bool = False,
    exact_timeout: float = 60,
    exact_max_num_qudits: int = 48,
) -> DistanceResult:
    """Compute a trust-aware distance for a CSS code over GF(q).

    For **prime** q the MILP provides the trusted signal; for prime-power q (where
    the integer MILP is invalid) the result falls back to the loose GUF bound and
    is left untrusted unless exact enumeration corroborates it.
    """
    q = int(q if q is not None else code.field.order)
    if code.dimension == 0:
        return DistanceResult(code.num_qudits, "trivial", True, True,
                              code.num_qudits, q, {"k": 0})

    upper = decoder_bound(code, prefilter_trials)  # cheap, loose pre-filter

    if is_prime(q):
        d_milp, details = compute_distance_milp(
            code, q,
            timeout_per_logical=milp_timeout_per_logical,
            total_timeout=milp_total_timeout,
            early_stop=None,  # certify exact
        )
        trusted = bool(details.get("exact"))
        method = "milp" if trusted else "milp_incumbent"
        d = d_milp
    else:
        # prime-power: integer mod-q MILP is invalid; only the bound is available.
        d = upper
        trusted = False
        method = "guf_bound"
        details = {"note": "prime-power q: MILP invalid; bound only (see docs/05)"}

    if want_exact:
        d_exact = compute_distance_exact(
            code, timeout=exact_timeout, max_num_qudits=exact_max_num_qudits
        )
        if d_exact is not None:
            details["exact_enum"] = d_exact
            if not trusted:
                d, trusted, method = d_exact, True, "exact"
            elif d_exact != d:
                # MILP-certified and exact disagree -> surface loudly; trust exact
                # (exact enumeration is itself a trusted method, so the existing
                # trusted/exact=True flags remain correct for the adopted value).
                details["milp_exact_mismatch"] = {"milp": d, "exact": d_exact}
                d, method = d_exact, "exact"

    return DistanceResult(
        d=int(d), method=method, trusted=trusted, exact=trusted,
        upper_bound=int(upper), q=q, details=details,
    )
