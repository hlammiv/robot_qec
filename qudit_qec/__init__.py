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

from . import construct, field_utils, genotype
from .construct import (
    build_bb_code,
    code_params,
    get_code_params_fast,
    validate_terms,
)
from .field_utils import (
    assert_is_stabilizer_code,
    combine_like_terms,
    get_field,
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
    # construct
    "validate_terms",
    "build_bb_code",
    "get_code_params_fast",
    "code_params",
    # field_utils
    "prime_factorization",
    "is_prime_power",
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
]
