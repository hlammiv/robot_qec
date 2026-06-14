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
| 🔜 `construct.py` | `evaluation/bb_code.py` | Build qudit CSS BB via `BBCode(field=q)` (qldpc auto-handles antipode/sign) | 1 ✓ |
| `distance.py` | `evaluation/distance.py` | Gate BP-OSD on `field.order==2`; kwarg-free GUF bound for q>2; OS-level timeout; `q^k`-gated exact | 2 ✓ |
| `distance_milp.py` | `evaluation/distance_milp.py` | **Prime-q** mod-q MILP (`−q` slack, big-M weight indicator, prime guard) — the trusted signal | 2 ✓ |
| `distance_qudit.py` *(new)* | — | Dispatcher: GUF pre-filter → MILP trusted → exact corroboration | 2 ✓ |
| `evaluator.py` | `evaluation/evaluator.py` | Field-threaded cascade; MILP-corroboration trust gate; coeff-aware keys; `A==B⇒d=2` gated on q==2 | 3 ✓ |
| `evolve/` | `seed_solution*.py`, `openevolve_evaluator*.py`, `prompt_context*.md`, `config*.yaml` | Coeff-bearing seeds + mutation, `QCODE_FIELD` plumbing, GF(q) prompts/baselines, retuned thresholds | 4 ✓ |
| `pbb_construct.py` | `evaluation/pbb_code.py` | One-line `[Bᵀ,−Aᵀ]%q` sign fix + `field=q`; in-field commutativity | 5 — |
| `dedup.py` | `evaluation/tanner_equivalence.py` | BLISS **structural value-gadget** encoding GF(q) edge coefficients | 6 — |
| `qudit_clifford.py` *(new)* | `evaluation/clifford_equivalence.py` | `Sp(2,q)` Clifford / LC↔CSS equivalence (research) | 7 — |
