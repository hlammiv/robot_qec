"""Field-threaded evaluation cascade for qudit CSS BB codes.

``evaluate_candidate`` turns a genotype ``(ell, m, A, B, field=q)`` into a coherent
scored result through the staged cascade that mirrors the reference pipeline, but
field- and coefficient-aware throughout:

* **Stage 0 -- validate / canonicalize.** Reject malformed or zero polynomials;
  combine like terms over GF(q). Composite ``q`` is rejected here (CRT/ring path
  not yet wired -- see ``docs/05``).
* **Stage 1 -- build + read (n, k).** ``k`` is exact (field-aware qldpc dimension);
  candidates with ``k < min_k`` are rejected cheaply before any distance work.
* **Structural -- decomposability.** A cheap Tanner-graph check flags direct-sum
  codes (no error-correction advantage).
* **Stage 2 -- distance.** ``"none"`` (k-only), ``"bound"`` (cheap loose GUF/BP-OSD
  upper bound, untrusted), or ``"trusted"`` (prime-q MILP; the default).
* **FOM** = ``k * d^2 / n``, reported only when a distance is available; ``trusted``
  is True only when the distance is MILP-certified or exact-enumerated.

The result carries the coefficient-aware canonical key (so distinct-coefficient
codes never collide downstream) and a ``self_dual`` marker (``A == B``), which the
reference identifies as a distance trap.
"""

from __future__ import annotations

from dataclasses import dataclass, field as _dc_field

from .construct import build_bb_code, validate_terms
from .distance import decoder_bound
from .distance_qudit import code_distance
from .genotype import pair_key
from .structure import is_decomposable


@dataclass
class EvalResult:
    """Outcome of evaluating one CSS BB candidate over GF(q)."""

    ell: int
    m: int
    q: int
    n: int
    k: int
    d: int | None
    d_status: str           # 'none'|'bound'|'milp'|'milp_incumbent'|'exact'|'trivial'
    fom: float | None
    trusted: bool
    decomposable: bool | None
    self_dual: bool
    A: tuple
    B: tuple
    key: tuple | None
    rejected: bool
    reason: str
    details: dict = _dc_field(default_factory=dict)


def _rejected(ell, m, q, reason, *, n=0, k=0, A=(), B=(), key=None,
              self_dual=False, decomposable=None) -> EvalResult:
    return EvalResult(ell, m, q, n, k, None, "none", None, False, decomposable,
                      self_dual, A, B, key, True, reason, {})


def evaluate_candidate(
    ell: int,
    m: int,
    A_terms,
    B_terms,
    field: int = 2,
    *,
    distance: str = "trusted",
    want_exact: bool = False,
    milp_timeout_per_logical: float = 20,
    milp_total_timeout: float = 120,
    prefilter_trials: int = 50,
    min_k: int = 1,
    check_decomposable: bool = True,
) -> EvalResult:
    """Evaluate a CSS BB candidate over GF(``field``); see module docstring.

    ``distance`` selects the cascade depth: ``"none"`` (Stage 1, k only),
    ``"bound"`` (cheap untrusted upper bound), or ``"trusted"`` (MILP, default).
    """
    q = int(field)

    # Stage 0 -- validate + canonicalize (rejects composite q, zero/oversize polys)
    try:
        A = validate_terms(ell, m, A_terms, q, "A")
        B = validate_terms(ell, m, B_terms, q, "B")
    except (ValueError, NotImplementedError) as exc:
        return _rejected(ell, m, q, f"invalid_genotype: {exc}")

    key = pair_key(A, B, ell, m, q)
    self_dual = A == B  # reference distance trap (A == B)

    # Stage 1 -- build + (n, k)
    code = build_bb_code(ell, m, A, B, field=q, validate=False)
    n, k = code.num_qudits, code.dimension
    if k < min_k:
        return EvalResult(ell, m, q, n, k, None, "none", 0.0, False, None,
                          self_dual, A, B, key, True, f"k<{min_k}", {})

    decomposable = is_decomposable(code) if check_decomposable else None

    if distance == "none":
        return EvalResult(ell, m, q, n, k, None, "none", None, False,
                          decomposable, self_dual, A, B, key, False, "", {})

    if distance == "bound":
        ub = decoder_bound(code, prefilter_trials)
        return EvalResult(ell, m, q, n, k, ub, "bound", k * ub * ub / n, False,
                          decomposable, self_dual, A, B, key, False, "",
                          {"upper_bound": ub})

    if distance != "trusted":
        raise ValueError(f"distance must be 'none'|'bound'|'trusted', got {distance!r}")

    # Stage 2 -- trusted distance (MILP, with optional exact corroboration)
    res = code_distance(
        code, q,
        prefilter_trials=prefilter_trials,
        milp_timeout_per_logical=milp_timeout_per_logical,
        milp_total_timeout=milp_total_timeout,
        want_exact=want_exact,
    )
    d = res.d
    fom = (k * d * d / n) if d is not None else None
    return EvalResult(ell, m, q, n, k, d, res.method, fom, res.trusted,
                      decomposable, self_dual, A, B, key, False, "", res.details)
