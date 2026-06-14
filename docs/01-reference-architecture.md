# 01 — Reference architecture: how `qcode-discovery` works

> Source: IBM Research `qcode-discovery`
> (https://github.com/qiskit-community/qcode-discovery), the code behind
> arXiv:2606.02418. This document is the verified output of the `Map` phase of
> [`workflows/qudit-qec-scope.js`](../workflows/qudit-qec-scope.js) (8 parallel
> readers over the codebase), reconciled with live experiments
> ([02](02-grounding-experiments.md)).

## What the system does

An LLM (OpenEvolve / MAP-Elites, FunSearch-style) mutates a Python program
`generate_candidates(ell, m)` that emits candidate polynomial tuples. Each
candidate is turned into a concrete `qldpc` code object and scored through a
multi-stage cascade. Data flows top-down; **the field-specific math lives in the
`evaluation/*` layer, while `evolve/*` is orchestration + genotype.**

```
evolve/  (LLM mutates generate_candidates → exponent-tuple polynomials)
   │   genotype: CSS = (A,B) pairs;  non-CSS PBB = (A,B,C,D) 4-tuples
   ▼
evaluation/bb_code.py | pbb_code.py        Stage 0  build qldpc code (FIELD FIXED HERE)
   ▼
evaluation/evaluator.py                    Stage 1  validate → read (n,k) via qldpc
   ▼                                        Stage 2  BP-OSD distance UPPER bound (ldpc)
   ▼                                        Stage 3  MILP EXACT distance (HiGHS/scipy)
   ▼
evaluation/results.py, tracking.py         Pareto front, JSONL logging
   ▼
post-campaign: tanner_equivalence.py (BLISS dedup) · clifford_equivalence.py (LC↔CSS)
```

FOM = `k·d²/n`. Five campaigns discovered 465 distinct codes (97 CSS, 368
non-CSS PBB) over `F₂[x,y]/(x^ℓ-1, y^m-1)`.

## Cluster by cluster

### 1. construction — `evaluation/bb_code.py`, `pbb_code.py`, `mirror_code.py`
Entry point. A genotype (exponent-tuple polynomials) becomes a concrete qldpc
code object.

- **CSS path** (`bb_code.py`): `terms_to_poly` builds a sympy poly,
  `build_bb_code` calls `codes.BBCode({x:ell, y:m}, poly_a, poly_b)` — currently
  with **no `field` arg**, defaulting to GF(2). *This is where the field is
  fixed.* Internally qldpc routes BB construction through its two-block code
  `TBCode`: `matrix_x = [A, B]`, `matrix_z = [Bᵀ, −Aᵀ]`.
- **Non-CSS PBB path** (`pbb_code.py`): a 4-tuple `(A,B,C,D)` → per-poly
  circulants via `bb.eval(poly).lift()` → tiled into a `2n`-wide symplectic
  stabilizer matrix → `QuditCode`. Contains hand-rolled GF(2) linear algebra
  (`_gf2_rref`, `_gf2_nullspace`, `_compute_logicals_gf2`). The current
  `block2_z = [mat_B.T % 2, mat_A.T % 2]` (≈line 159) is the naive GF(2) form
  that is **wrong over GF(q)** (it produces a subsystem code).
- **mirror_code.py**: Khesin–Lu baseline/test helper only; low priority.

### 2. evaluator-cascade — `evaluation/evaluator.py`, `results.py`, `tracking.py`
Orchestration spine. `evaluate_candidate` runs: Stage 1 validate → Stage 2 build
+ read `(n,k)` → Stage 3/4 BP-OSD distance → Stage 5 exact (qldpc brute force)
or MILP. Most field math is delegated to qldpc; the residual GF(2) commitments
are: coefficient-free cache keys (`_milp_cache_key`), the self-dual `A==B ⇒ d=2`
gate, and the `d/√n` trust thresholds. `tracking.py` is field-agnostic;
`results.py` Pareto math (`k/n`, `d`, `1/n`) is field-independent — only its
`_code_key` dedup is coefficient-blind.

### 3. distance-bposd — `evaluation/distance.py`, `distance_bposd_noncss.py`
Stage-2 distance estimator returning an **upper bound**.

- `distance.py` (CSS): delegates to `code.get_distance_bound_with_decoder(
  Pauli.X/Z, …, bp_method='product_sum', osd_method=…, osd_order=…)`. **Indirect
  break over GF(q):** qldpc routes non-binary matrices to `GUFDecoder`, which
  rejects those kwargs (`TypeError`). Note `get_distance_bound_with_decoder`
  exists **only on `CSSCode`/`BBCode`, not on `QuditCode`**.
- `distance_bposd_noncss.py` (non-CSS): re-implements a Bravyi-style randomized
  BP-OSD bound calling `ldpc.BpOsdDecoder` (**binary-only**) directly, with
  uint8 mod-2 XOR syndrome hashing and a 3-Pauli-per-qubit (X/Z/Y) model.
  **Direct, deep break — research-grade rewrite.**

### 4. distance-milp — `evaluation/distance_milp.py`
Stage-3 exact distance via HiGHS MILP (no decoder, thread-safe). A CSS
formulation (`ilp_min_weight`, solving `d_X` and `d_Z` separately) and a non-CSS
symplectic formulation (`ilp_min_weight_symplectic`). **Hard-wired to GF(2):**
every matrix `% 2`, binary `{0,1}` vars, mod-2 commutation via
`sum − 2·slack = 0`, Hamming-weight objective, anticommutation `= 1`. Portable
for **prime q** via concrete deltas; **invalid for prime-power q** (a MILP
cannot do GF(pᵐ) arithmetic). This is the true critical path — see
[03](03-qudit-extension-scope.md).

### 5. equivalence-dedup — `evaluation/tanner_equivalence.py`, `clifford_equivalence.py`
Post-campaign. `tanner_equivalence.py`: BLISS canonical-hash permutation dedup
via a colored Tanner graph (python-igraph). `clifford_equivalence.py`: LC-to-CSS
equivalence via the 6-element single-qubit Clifford group `{I,S,H,HS,SH,HSH}` and
union-find 2-coloring. **Tanner dedup is salvageable for GF(q)** (a structural
value-gadget encodes coefficients, verified); **the Clifford machinery is
research-grade** (the single-qudit Clifford group is `Sp(2,q)` of order
`q(q²−1)`, not 6).

### 6. evolve — `evolve/seed_solution*.py`, `openevolve_evaluator*.py`, `config*.yaml`
OpenEvolve/MAP-Elites mutates `generate_candidates`. **No direct qldpc/galois
calls.** GF(2)-specificity lives in: the coefficient-free genotype, GF(2)-tuned
scoring constants (`k ≥ 8`, the `d/√n` trust band, Bravyi FOMs), and the GF(2)
`KNOWN_CODES` catalog. The field is threaded via the established `QCODE_*`
env-var channel (OpenEvolve's `evaluate(program_path)` takes no extra args).

### 7. prompts-entry-scripts — `evolve/prompt_context*.md`, `main.py`, `scripts/verify_*.py`
LLM-facing prompts (the F₂ genotype contract, FOM, F₂ algebra) + the CLI driver
+ post-campaign verification harnesses (which carry **drifted inlined copies** of
the BLISS hash, each with its own `% 2`).

## The load-bearing takeaway

The construction library **`qldpc` is field-generic**: `BBCode(orders, A, B,
field=q)` builds a *qudit* BB code over GF(q) today, and auto-inserts the
antipode+sign that makes the checks commute for odd `q`. So the qudit extension
is a focused set of edits — concentrated in the genotype (coefficients) and the
distance layer (GF(q) decoding + MILP) — not a rebuild. The evidence is in
[02](02-grounding-experiments.md); the plan is in
[03](03-qudit-extension-scope.md) and [04](04-implementation-roadmap.md).
