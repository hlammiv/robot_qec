# `qudit_qec` (scaffold)

Our code for the qudit extension. Empty for now beyond this map — modules land
here as we execute `docs/04-implementation-roadmap.md`.

## Module map (ours ⟶ what it extends in `qcode-discovery`)

Finalized by [`docs/03`](../docs/03-qudit-extension-scope.md) /
[`docs/04`](../docs/04-implementation-roadmap.md). ✓ = MVP (CSS, prime q); — = deferred.

Status: **✅ implemented** · 🔜 next · — deferred.

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
| `pbb_construct.py` | `evaluation/pbb_code.py` | One-line `[Bᵀ,−Aᵀ]%q` sign fix + `field=q`; in-field commutativity | 5 — |
| `dedup.py` | `evaluation/tanner_equivalence.py` | BLISS **structural value-gadget** encoding GF(q) edge coefficients | 6 — |
| `qudit_clifford.py` *(new)* | `evaluation/clifford_equivalence.py` | `Sp(2,q)` Clifford / LC↔CSS equivalence (research) | 7 — |
