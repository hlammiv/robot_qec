"""robot_qec — LLM-guided discovery of qudit (GF(q)) quantum error-correcting codes.

This package houses the qudit extension of the IBM ``qcode-discovery`` pipeline.

Design principle: lean on the field-generic ``qldpc`` library (which already
constructs qudit BB/CSS/non-CSS codes over GF(q) via a ``field=`` argument)
rather than reimplementing code construction. Our code concentrates on the parts
that are GF(2)-specific in the reference pipeline: the evolutionary genotype
(GF(q) coefficients), GF(q) distance/decoding, and field-aware dedup/equivalence.

Phase 0 (this release) provides the field substrate and genotype foundations; see
``docs/04-implementation-roadmap.md``.
"""

from __future__ import annotations

from . import (
    construct,
    crt,
    distance,
    distance_milp,
    distance_qudit,
    evaluator,
    field_utils,
    genotype,
    results,
    structure,
)
from .crt import (
    CRTResult,
    canonicalize_zd,
    classify,
    crt_moduli,
    evaluate_crt_candidate,
    is_squarefree,
    split_terms,
)
from .construct import (
    build_bb_code,
    code_params,
    get_code_params_fast,
    validate_terms,
)
from .distance import compute_distance_exact, decoder_bound, qdistrnd_bound
from .distance_milp import (
    certify_distance_geq,
    compute_distance_milp,
    ilp_feasible_weight_le,
    ilp_min_weight,
)
from .distance_qudit import DistanceResult, code_distance
from .evaluator import EvalResult, evaluate_candidate
from .results import CodeCatalog
from .structure import connected_components, is_decomposable, num_connected_components
from .field_utils import (
    assert_is_stabilizer_code,
    combine_like_terms,
    get_field,
    is_prime,
    is_prime_power,
    prime_factorization,
    symplectic_conjugate,
    terms_to_poly,
    to_field_element,
)
from .genotype import (
    as_triple,
    canonicalize,
    normalize_terms,
    pair_key,
    poly_key,
    tuple_key,
)

__version__ = "0.1.0"

__all__ = [
    "field_utils",
    "genotype",
    "construct",
    "distance",
    "distance_milp",
    "distance_qudit",
    "evaluator",
    "results",
    "structure",
    "crt",
    # construct
    "validate_terms",
    "build_bb_code",
    "get_code_params_fast",
    "code_params",
    # distance
    "decoder_bound",
    "qdistrnd_bound",
    "compute_distance_exact",
    "ilp_min_weight",
    "ilp_feasible_weight_le",
    "compute_distance_milp",
    "certify_distance_geq",
    "code_distance",
    "DistanceResult",
    # field_utils
    "prime_factorization",
    "is_prime_power",
    "is_prime",
    "get_field",
    "to_field_element",
    "combine_like_terms",
    "terms_to_poly",
    "symplectic_conjugate",
    "assert_is_stabilizer_code",
    # genotype
    "as_triple",
    "normalize_terms",
    "canonicalize",
    "poly_key",
    "pair_key",
    "tuple_key",
    # evaluator / results / structure
    "evaluate_candidate",
    "EvalResult",
    "CodeCatalog",
    "is_decomposable",
    "connected_components",
    "num_connected_components",
    # crt (arbitrary square-free dimension)
    "classify",
    "crt_moduli",
    "is_squarefree",
    "split_terms",
    "evaluate_crt_candidate",
    "CRTResult",
]
