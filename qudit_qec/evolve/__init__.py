"""Evolutionary-search scaffolding for qudit BB code discovery (Phase 4).

This subpackage wires the trusted ``qudit_qec`` evaluation pipeline into an
OpenEvolve / MAP-Elites loop:

* ``seed_solution_qudit.py`` -- the seed program OpenEvolve mutates (a
  coefficient-bearing ``generate_candidates`` over GF(q)).
* ``adapter.py`` -- ``evaluate(program_path)``, the fitness function (uses only the
  ``qudit_qec`` pipeline; importable/testable without ``openevolve``).
* ``run_evolution.py`` -- the CLI that configures and launches a campaign
  (``--field`` -> ``QCODE_FIELD``; imports ``openevolve`` lazily).
* ``config_qudit.yaml`` / ``prompt_context_qudit.md`` -- OpenEvolve config and the
  GF(q) domain knowledge fed to the LLM.

``adapter.evaluate`` is the only piece imported here, to keep this subpackage
importable without ``openevolve`` installed.
"""

from .adapter import evaluate

__all__ = ["evaluate"]
