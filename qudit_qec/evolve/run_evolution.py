"""Launch an OpenEvolve campaign to discover qudit (GF(q)) BB codes.

Usage (requires ``openevolve`` + ``litellm`` installed and an LLM proxy running):

    python -m qudit_qec.evolve.run_evolution --field 3 --iterations 200 \
        --model anthropic/claude-opus-4-8 --api-base http://localhost:4000/v1

The ``--field q`` value is threaded to the evaluator and seed via the
``QCODE_FIELD`` environment variable (OpenEvolve calls ``evaluate(program_path)``
with no way to pass extra args). ``q`` must be **prime** for the direct search (the
trusted MILP distance is prime-only); square-free composite dimensions are reached
via the CRT layer at evaluation time.

``openevolve`` is imported lazily inside :func:`main`, so this module (and its
argument parser) import fine without it -- which is what the scaffolding tests rely
on. When ``openevolve`` is absent, ``main`` prints install instructions and exits.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_REPO_ROOT = _HERE.parents[1]  # .../robot_qec
DEFAULT_SEED = _HERE / "seed_solution_qudit.py"
DEFAULT_CONFIG = _HERE / "config_qudit.yaml"
DEFAULT_ADAPTER = _HERE / "adapter.py"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Run an OpenEvolve campaign for qudit (GF(q)) BB code discovery."
    )
    p.add_argument("--field", type=int, default=2,
                   help="GF(q) field order (prime) for the search. Default 2 (qubits).")
    p.add_argument("--iterations", type=int, default=100,
                   help="Number of evolutionary iterations.")
    p.add_argument("--model", action="append", default=None,
                   help="LLM model id (repeat for an ensemble). Overrides config models.")
    p.add_argument("--api-base", default=None,
                   help="LLM API base URL (e.g. a LiteLLM proxy at http://localhost:4000/v1).")
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG,
                   help=f"OpenEvolve config YAML (default: {DEFAULT_CONFIG.name}).")
    p.add_argument("--seed", type=Path, default=DEFAULT_SEED,
                   help=f"Seed program to evolve (default: {DEFAULT_SEED.name}).")
    p.add_argument("--run-name", default=None,
                   help="Campaign name (defaults to qudit_gf<q>).")
    p.add_argument("--high-k", type=int, default=None,
                   help="Override the high-k feature threshold (default 8; recalibrate per q).")
    return p


def _validate_prime_field(q: int) -> None:
    # Local import keeps the parser usable even if the package isn't on sys.path yet.
    from qudit_qec.field_utils import is_prime

    if not is_prime(q):
        raise SystemExit(
            f"--field must be a prime for the direct search, got {q}. "
            f"Prime-power and composite dimensions are not yet wired into the loop "
            f"(prime powers: Phase 7; square-free composites: use the CRT layer)."
        )


def configure_environment(args: argparse.Namespace) -> str:
    """Set the QCODE_* env channel and ensure the repo root is importable."""
    if str(_REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(_REPO_ROOT))
    run_name = args.run_name or f"qudit_gf{args.field}"
    os.environ["QCODE_FIELD"] = str(args.field)
    os.environ["QCODE_RUN_NAME"] = run_name
    if args.high_k is not None:
        os.environ["QCODE_HIGH_K"] = str(args.high_k)
    return run_name


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    _validate_prime_field(args.field)
    run_name = configure_environment(args)

    try:
        import openevolve  # noqa: F401
    except ImportError:
        print(
            "openevolve is not installed. This launcher is ready, but a campaign "
            "needs:\n"
            "  pip install openevolve litellm\n"
            "and a running LLM proxy (e.g. LiteLLM at --api-base). "
            "The seed, evaluator (adapter.py), config, and prompt are all in place; "
            "see qudit_qec/evolve/README.md.",
            file=sys.stderr,
        )
        return 2

    # --- OpenEvolve wiring (mirrors the reference run_evolution.py) ---------- #
    from openevolve import Config, run_evolution
    from openevolve.config import LLMModelConfig

    config = Config.from_yaml(str(args.config))
    if args.api_base:
        config.llm.api_base = args.api_base
    if args.model:
        config.llm.models = [LLMModelConfig(name=name, weight=1.0) for name in args.model]

    print(f"[run_evolution] campaign '{run_name}' over GF({args.field}), "
          f"{args.iterations} iters, seed={args.seed.name}, evaluator={DEFAULT_ADAPTER.name}")

    run_evolution(
        initial_program=str(args.seed),
        evaluation_file=str(DEFAULT_ADAPTER),
        config=config,
        iterations=args.iterations,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
