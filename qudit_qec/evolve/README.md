# `qudit_qec.evolve` — LLM-guided discovery loop (Phase 4 scaffolding)

Wires the trusted `qudit_qec` evaluation pipeline into an OpenEvolve / MAP-Elites
search so an LLM *discovers* new qudit (GF(q)) BB codes.

| File | Role |
|------|------|
| `seed_solution_qudit.py` | The seed program OpenEvolve mutates — a coefficient-bearing `generate_candidates(ell, m)` over GF(q) (EVOLVE-BLOCK). Reads `QCODE_FIELD`. |
| `adapter.py` | `evaluate(program_path) → dict` — the fitness function. k-screens all candidates, runs the **trusted prime-q MILP** on the most promising, returns `combined_score` + the MAP-Elites feature dims (`lattices_with_high_k`, `num_high_k`). Uses only `qudit_qec`; no `openevolve` dependency. |
| `run_evolution.py` | CLI launcher. `--field q` → `QCODE_FIELD`; imports `openevolve` lazily. |
| `config_qudit.yaml` | OpenEvolve config (GF(q) system message, models, feature dims). |
| `prompt_context_qudit.md` | GF(q) domain knowledge for the LLM (antipode/sign, the coefficient axis, structural families, traps). |

## Status

Scaffolding is **complete and tested** (the seed generates valid GF(q) candidates;
`adapter.evaluate` runs end-to-end and returns the OpenEvolve contract). What it
does **not** yet have is the external runtime: `openevolve`, `litellm`, and an LLM
endpoint. The launcher degrades gracefully (prints install instructions) until
those are present.

## Running a campaign (once the runtime is available)

```bash
pip install openevolve litellm        # not currently installed
# start a LiteLLM proxy exposing your model(s) at, e.g., http://localhost:4000/v1

python -m qudit_qec.evolve.run_evolution \
    --field 3 \
    --iterations 200 \
    --model anthropic/claude-opus-4-8 \
    --api-base http://localhost:4000/v1
```

- `--field` must be **prime** (the trusted MILP distance is prime-only). Square-free
  composite dimensions are reached by evaluating through the CRT layer; prime powers
  are Phase 7.
- `--model` may be repeated for an ensemble.
- The fitness rewards exact encoded dimension `k`, breadth across lattices, and the
  best **trusted (MILP-certified)** FOM — never a loose decoder bound.

## Testing without the runtime

`tests/test_evolve_scaffold.py` exercises the whole contract with no `openevolve`:
seed generation over GF(2)/GF(3), `adapter.evaluate` returning the feature dims, and
the launcher's `--field` plumbing + graceful no-openevolve exit.
