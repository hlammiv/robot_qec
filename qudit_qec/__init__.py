"""robot_qec — LLM-guided discovery of qudit (GF(q)) quantum error-correcting codes.

This package will house the qudit extension of the IBM `qcode-discovery`
pipeline. It is currently a scaffold; modules are added as the implementation
roadmap (`docs/04-implementation-roadmap.md`) is executed.

Design principle: lean on the field-generic `qldpc` library (which already
constructs qudit BB/CSS/non-CSS codes over GF(q) via a ``field=`` argument)
rather than reimplementing code construction. Our code concentrates on the
parts that are GF(2)-specific in the reference pipeline: the evolutionary
genotype (GF(q) coefficients), GF(q) distance/decoding, and field-aware
dedup/equivalence.
"""

__version__ = "0.0.0"
