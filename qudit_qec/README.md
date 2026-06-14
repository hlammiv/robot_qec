# `qudit_qec` (scaffold)

Our code for the qudit extension. Empty for now beyond this map — modules land
here as we execute `docs/04-implementation-roadmap.md`.

## Planned module map (ours ⟶ what it extends in `qcode-discovery`)

| Planned module | Extends / replaces | Responsibility for GF(q) |
|----------------|--------------------|--------------------------|
| `construct.py` | `evaluation/bb_code.py`, `pbb_code.py` | Build qudit BB/PBB codes via `qldpc ...BBCode(..., field=q)` / `QuditCode(..., field=q)`; carry GF(q) coefficients. |
| `genotype.py` | the exponent-tuple convention in `bb_code.py` + seeds | Term = `(x_exp, y_exp, coeff∈GF(q)*)`; `terms_to_poly`, `validate_terms`. |
| `distance.py` | `evaluation/distance.py` | Field-aware distance: exact-enum (small n), qldpc qudit decoder bound, fallthrough. |
| `distance_milp.py` | `evaluation/distance_milp.py` | GF(q) minimum-weight-logical MILP (mod-q constraints, symplectic weight). |
| `evaluator.py` | `evaluation/evaluator.py` | Field-threaded cascade: k → quick distance → exact. |
| `dedup.py` | `evaluation/tanner_equivalence.py` | BLISS coloring that encodes GF(q) edge coefficients. |
| `equivalence.py` | `evaluation/clifford_equivalence.py` | Qudit Clifford / LC equivalence (likely post-MVP / research). |
| `evolve/` | `evolve/seed_solution*.py`, `openevolve_evaluator*.py`, `prompt_context*.md`, `config*.yaml` | Qudit-aware seeds, fitness, MAP-Elites features, and LLM prompt domain knowledge. |

> This table is provisional and is finalized by `docs/03-qudit-extension-scope.md`
> once the scoping workflow completes.
