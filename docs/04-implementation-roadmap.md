# 04 — Implementation roadmap

> The phased plan to reach an MVP that **evolves AND MILP/exact-verifies at least
> one new qudit code**, then expands toward the full pipeline. Derived from
> [03](03-qudit-extension-scope.md). **Start with q=3, CSS arm only.**

## MVP definition

**Minimal believable MVP: prime q (start q=3), CSS arm only, that evolves AND
trustworthily MILP/exact-verifies at least one NEW qudit code.**

- **Construction:** `field=q` threaded through `build_bb_code` +
  coefficient-carrying genotype `(x_exp, y_exp, coeff)` with a 2-tuple→coeff-1
  back-compat shim + coefficient-aware dedup/cache keys. qldpc auto-handles the
  CSS antipode/sign, so no hand-derived sign is needed on the CSS path.
- **Distance:** kwarg-free GUF per-sector bound as a **cheap pre-filter only**
  (~3× loose); prime-q CSS mod-q MILP (`−q` slack, big-M weight indicator,
  residue-disjunction anticommutation, hard prime guard) as the **trusted**
  ranking/save signal; `q^k`-gated `get_distance_exact` with an OS-level
  subprocess timeout for final corroboration. Keep a corroboration trust gate.
- **Evolution:** `--field`/`QCODE_FIELD` plumbed through `run_evolution.py` and
  both evaluators; a GF(3)-verified seed catalog; a GF(q)-rewritten prompt
  teaching q-ary coefficients and the antipode; recalibrated k/d/trust thresholds;
  a one-off GF(3) lattice-viability sweep.

**Explicitly OUT of the MVP (deferred to research):** the entire non-CSS/PBB arm,
all prime-power (GF(4)/GF(8)/GF(9)) support, composite-`d` (Z_d) qudits, the full
`distance_bposd_noncss.py` rewrite, and all Clifford/LC-to-CSS equivalence. The
MVP delivers a real, trusted, novel GF(3) CSS code with exact/MILP distance
certification plus a Pareto front + permutation-dedup + decomposability accounting.

**Realistic effort: ~1.5–2.5 weeks** for the cross-concern MVP (not ~1 week). The
prime-q MILP is the ~2–3 day critical path.

## Phases

### Phase 0 — Field substrate + guards (~0.5 day)
- New `evaluation/field_utils.py`: `get_field(q)` (galois, prime-power
  validation, reject composite-`d` with a `NotImplementedError` pointing at the
  Howell-normal-form research path), `assert_is_stabilizer_code(code)` (in-field
  `symplectic_conjugate(M)@M.T == 0` + `not is_subsystem_code`), `terms_to_poly_q`,
  `combine_like_terms`. Centralizes the prime-power `%q`-on-`FieldArray` pitfall.
- New `evaluation/genotype.py`: normalized `(x,y,coeff)` term type + `canonicalize`
  (sorted triples), reused by every cache/dedup key.
- **Milestone:** `field_utils` round-trips GF(3)/GF(4) codes; composite-`d` raises
  cleanly.

### Phase 1 — CSS construction over GF(q) (~0.5 day, de-risked)
- `bb_code.py`: coeff-carrying `terms_to_poly` (2-tuple→coeff=1 shim), `field=q`
  in `build_bb_code`, coeff-aware `validate_terms`, post-construction
  `assert code.field.order == q`.
- **Milestone:** `build_bb_code(6,6,A,B,field=3)` reproduces `[[72,8]]`,
  `is_subsystem_code=False`; a single-monomial coeff flip changes `k` (regression
  test).

### Phase 2 — Prime-q distance backend (~2–3 days, **critical path**)
- `distance.py`: gate BP-OSD kwargs on `field.order==2`; kwarg-free GUF CSS bound
  for q>2; OS-level timeout replacing SIGALRM; `q^k`-gated exact.
- `distance_milp.py`: prime-q `ilp_min_weight` (CSS) — big-M weight indicator,
  `−q` slack, residue-loop/big-M anticommutation, **prime guard**. Thread q.
- New `evaluation/distance_qudit.py` dispatcher: prime vs prime-power branch, GUF
  pre-filter, MILP trusted ranking, exact corroboration.
- `tests/test_distance_gf.py`: GUF returns a finite bound over GF(3) (no
  TypeError); MILP min-weight matches hand-checked GF(3) examples; GF(4) MILP
  raises/skips; GUF+exact agree on a tiny GF(3) code.
- **Milestone:** MILP returns a TRUSTED exact distance for a GF(3) CSS BB code;
  GUF only screens. **First MILP-verified qudit code.**

### Phase 3 — Cascade + dedup/cache keys (~1 day)
- `evaluator.py`: add q; thread to result dict; coeff-aware self-dual gate (gated
  on q==2); coeff in `_milp_cache_key`; trust-filter gated to MILP/exact
  corroboration for q>2.
- `results.py`: coeff in `_code_key`; record q.
- `clifford_equivalence.py:247-278`: galois GF(q) rank → decomposability works for
  free (verified `rank 32+32=64` over GF(3)).
- **Milestone:** end-to-end `evaluate_candidate(…, field=3)` returns a coherent
  `{n,k,d,fom}` with `d` from MILP; two coeff-differing GF(3) codes are not merged.

### Phase 4 — Evolution loop + prompts (~2–3 days)
- `run_evolution.py`: `--field` → `QCODE_FIELD` env. Evaluators read q, thread
  `field=q`, coeff-aware dedup/classify/feedback, distance-path gating,
  recalibrated k/d/trust thresholds.
- New `evolve/seed_solution_qudit.py` with coeff-bearing genotype + coeff-mutation
  + GF(3)-verified `KNOWN_CODES` (re-verify the `_safety_net_codes` `k>0`
  guarantee over GF(3) or the cascade starves).
- New `evolve/prompt_context_qudit.md` + `config_qudit.yaml`: F_q ring, q-ary
  coeffs, antipode/sign, qudit baselines; remove F₂-only lore (`1+y+y²`
  irreducibility, `A=B⇒d=2`, doubling `c=2b`).
- One-off GF(3) lattice-viability sweep regenerating `STAGE*_LATTICES`.
- **🎯 MVP MILESTONE:** a short GF(3) CSS campaign evolves candidates,
  MILP-verifies distances, and surfaces at least one **new** `[[n,k,d]]₃` code not
  in the GF(2) catalog, with permutation-dedup + decomposability accounting.

### Phase 4.5 — CRT layer for arbitrary square-free `d` ✅ DONE (`qudit_qec/crt.py`)
> Added by the 2026-06-14 decision to target arbitrary `d` via CRT factoring.
> See [05](05-arbitrary-dimension-crt.md). Depends only on the prime-`q` CSS MVP.
- New `qudit_qec/crt.py`: `factor_dimension`, `crt_moduli`, `classify(d)`,
  `split_genotype` (reduce `Z_d` coeffs mod each distinct prime), `build_crt_code`
  (per-factor dispatch), `combine` (`k` per factor, `d = minᵢ dᵢ`, FOM).
- Genotype is `Z_d`; evaluation CRT-splits to prime fields and reuses the Phase-2
  prime-field distance backend per factor (MILP valid: each factor is a prime
  field). New search signal: prefer balanced `kᵢ`/`dᵢ` across factors.
- **Milestone:** an end-to-end `Z_6` (qubit⊗qutrit) BB code is built by CRT,
  each factor MILP-verified, and reported as a single `d=6` code with
  `distance = min(d₂, d₃)`. Delivers **arbitrary square-free dimension** with no
  new field/ring math.

### Phase 5 — Non-CSS PBB arm, prime q (~3–5 days, research-grade)
- `pbb_code.py`: one-line `[Bᵀ | −Aᵀ] % q` sign fix + `field=q`; in-field
  commutativity check; kwarg-free `get_logical_ops()` for logicals (port the
  galois fallback only if it fails on a given code).
- Symplectic MILP (`ilp_min_weight_symplectic`) as the non-CSS trusted distance
  (no cheap GUF non-CSS bound exists).
- Non-CSS evaluator/seed/prompt arm; verify scripts threaded with `field=q`.
- **Milestone:** a non-CSS PBB GF(3) code is constructed (`subsystem=False`) and
  MILP-verified.

### Phase 6 — Post-campaign GF(q) dedup (~2–3 days)
- `tanner_equivalence.py`: structural value-gadget BLISS encoding (verified live on
  `[[72,8]]₃`); de-dup the two drifted inlined copies in verify scripts into the
  shared module; optional scalar canonicalization for per-qudit GF(q)* equivalence.
- **Milestone:** a coeff flip changes the canonical hash; a qudit permutation
  preserves it; GF(2) behavior unchanged.

### Phase 7+ — Research expansion (~1–2 weeks each)
- Full non-CSS q-ary randomized estimator (`distance_bposd_noncss.py` rewrite
  around a GUFDecoder wrapper).
- **Prime-power dimensions — both tracks, ordered cheapest-first** (decided
  2026-06-14; see [05](05-arbitrary-dimension-crt.md)):
  - **Phase 7a — Galois-qudit `GF(p^a)`** (medium, *do first*): construction + `k`
    already free via qldpc/galois; add only the `GF(p^a)`-valid distance path
    (audit every `%q` for field arithmetic; prime-subfield MILP, or GUF + exact).
    Delivers `GF(4)/GF(8)/GF(9)` field qudits.
  - **Phase 7b — Modular-qudit `Z_{p^a}`** (large, research, *after 7a*): a
    ring-linear-algebra backend — `k`/logicals via Smith/Howell normal form over
    `Z_{p^a}` (galois cannot represent it); the distance MILP is valid here (integer
    mod-`p^a` *is* the ring arithmetic). Yields the *physical* clock-mod-`d` qudit
    and, with Phase 4.5's CRT layer, **completes arbitrary integer dimension `d`**.
- Qudit Clifford/LC-to-CSS equivalence (new `qudit_clifford.py`: `Sp(2,q)`
  generators + orbit search) — needed to reproduce non-CSS publication claims over
  GF(q).

## Critical path & sequencing notes

- **Phase 2 (prime-q MILP) is the bottleneck** and gates "trusted distance," which
  gates the whole MVP claim. Everything else is comparatively cheap.
- Phases 0–1 are low-risk and unlock testing for Phase 2.
- Phases 5–7 are independent research tracks; none blocks the CSS MVP.
- Decide **open question #1 (FOM normalization)** before Phase 4 — it touches both
  evaluators, the Pareto front, prompts, and thresholds simultaneously.
