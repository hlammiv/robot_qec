# `qudit_qec`

Our code for the qudit extension — a trusted GF(q) code-discovery pipeline, two
magic-state-distillation discovery arms, and a gate-set-universality layer. See the
[top-level README](../README.md) for the capability overview and
[`docs/08`](../docs/08-distillation-arm-and-universality.md) for the distillation +
universality design.

## Module map (ours ⟶ what it extends in `qcode-discovery`)

The discovery-pipeline modules were finalized by
[`docs/03`](../docs/03-qudit-extension-scope.md) /
[`docs/04`](../docs/04-implementation-roadmap.md); the distillation + universality
modules by [`docs/07`](../docs/07-magic-state-distillation-scope.md) /
[`docs/08`](../docs/08-distillation-arm-and-universality.md).

Status: **✅ implemented** · 🟡 built, needs an external endpoint to run · — deferred.

### Trusted code-discovery pipeline

| Planned module | Extends / replaces | Responsibility for GF(q) | Phase |
|---|---|---|---|
| ✅ `field_utils.py` | — | `get_field(q)` prime-power guard (composite→CRT), in-field `assert_is_stabilizer_code`, `terms_to_poly`, `combine_like_terms`, `to_field_element`; **forbids raw `%q` on FieldArrays** | 0 ✓ |
| ✅ `genotype.py` | exponent-tuple convention in `bb_code.py` + seeds | Term = `(x_exp, y_exp, coeff∈GF(q))`; `canonicalize`/`poly_key`/`pair_key`/`tuple_key`; reused by all cache/dedup keys | 0 ✓ |
| ✅ `construct.py` | `evaluation/bb_code.py` | Build qudit CSS BB via `BBCode(field=q)` (qldpc auto-handles antipode/sign); `build_bb_code`/`validate_terms`/`code_params`; canonicalizes coeffs in-field before construction | 1 ✓ |
| ✅ `distance.py` | `evaluation/distance.py` | `decoder_bound` gates BP-OSD on `field.order==2`, kwarg-free GUF for q>2 (cheap, loose pre-filter); `compute_distance_exact` (forked, OS-timeout) | 2 ✓ |
| ✅ `distance_milp.py` | `evaluation/distance_milp.py` | **Prime-q** mod-q MILP (`−q` slack, big-M weight indicator, `L·x=1` unit-scaling, prime guard); `ilp_feasible_weight_le` + `certify_distance_geq` (weight-cut **lower-bound proof**) | 2 ✓ |
| ✅ `distance_qudit.py` *(new)* | — | `code_distance` + `DistanceResult`: GUF pre-filter → MILP + **QDistRnd** (2nd independent source) → tightest bound; **certifies exact** via QDistRnd-upper + weight-cut-lower; trust gate | 2 ✓ |
| ✅ `evaluator.py` | `evaluation/evaluator.py` | `evaluate_candidate(field=q)` cascade (validate→k→bound/MILP→FOM); `EvalResult` with trust gate, `self_dual` marker, coeff-aware key | 3 ✓ |
| ✅ `structure.py` *(new)* | (part of `tanner_equivalence`) | Tanner-graph decomposability (union-find): `is_decomposable`/`connected_components` — direct-sum detection, field-agnostic | 3 ✓ |
| ✅ `results.py` | `evaluation/results.py` | `CodeCatalog`: coeff-aware dedup (records q), Pareto front (n/k/d), best-by-FOM | 3 ✓ |
| ✅ `crt.py` *(new)* | — | Arbitrary **square-free** `d` via CRT: `classify`/`crt_moduli`/`split_terms`, `evaluate_crt_candidate` → per-prime-factor field pipeline, distance = minᵢ dᵢ, `CRTResult` | 4.5 ✓ |
| 🟡 `evolve/` | `seed_solution*.py`, `openevolve_evaluator*.py`, `prompt_context*.md`, `config*.yaml` | **Scaffolding built+tested:** coeff-bearing `seed_solution_qudit.py`, `adapter.evaluate`, `run_evolution --field`/`QCODE_FIELD`, `config_qudit.yaml`, `prompt_context_qudit.md`. Live campaign needs `openevolve`+LLM endpoint | 4 🟡 |

### Distillation + universality (is the code *useful*?)  — [`docs/08`](../docs/08-distillation-arm-and-universality.md)

| Module | Responsibility | Validated against |
|---|---|---|
| ✅ `distillation.py` | v0 MSD suitability: `is_triorthogonal` (cubic mod-p), `transversal_gate_level`, `weight_enumerator` (+MacWilliams over GF(p²)), `magic_state_yield` (γ / strange threshold). Prime-only, soundness-gated. | `[[20,7,2]]₃` γ=1.51, `[[14,4,2]]₃` γ=1.81, 11-qutrit Golay A(z)/B(z), `[[15,1,3]]₂` γ=2.46 |
| ✅ `distill_discovery.py` | d=2 triorthogonal **T-gate arm** (odd prime): `triortho_family`, `reed_muller_triortho`, re-validated operators, `DistillCatalog`, `search_distill`. Objective = yield γ. | family vs arXiv:2403.06228; `reed_muller_triortho(2,4,1)` = qubit `[[15,1,3]]` (d=3 trusted MILP) |
| ✅ `distill_strange.py` | d>2 **qutrit strange/QR** sub-arm: self-orthogonal cyclic-F₃ genotype → strange-state CSS; cheap Fraction screen + opt-in sympy threshold; `StrangeCatalog`, `search_strange_cyclic`. | 11-qutrit Golay `[[11,1,5]]₃` (cubic ν=3, threshold ≈0.387) |
| ✅ `universality.py` | p-adic **sector-coverage checker**: does a single-qudit diagonal gate couple all Pauli sectors `W_k` at `d=p^m`? (necessary for `{Clifford+G}` universality). | Clifford gates don't couple (d=4,8,16,9,27); CCZ reductions are Pauli → don't couple |

### Deferred research tracks

| Planned module | Extends / replaces | Responsibility for GF(q) | Phase |
|---|---|---|---|
| `pbb_construct.py` | `evaluation/pbb_code.py` | One-line `[Bᵀ,−Aᵀ]%q` sign fix + `field=q`; in-field commutativity | 5 — |
| `dedup.py` | `evaluation/tanner_equivalence.py` | BLISS **structural value-gadget** encoding GF(q) edge coefficients | 6 — |
| `qudit_clifford.py` *(new)* | `evaluation/clifford_equivalence.py` | `Sp(2,q)` Clifford / LC↔CSS equivalence (research) | 7 — |
