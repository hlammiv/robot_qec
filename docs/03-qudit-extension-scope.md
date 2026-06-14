# 03 — Qudit (GF(q)) extension: scope

> The verified scoping synthesis from
> [`workflows/qudit-qec-scope.js`](../workflows/qudit-qec-scope.js), reconciled
> with the live evidence in [02](02-grounding-experiments.md) and the literature.
> Read [01](01-reference-architecture.md) first for the pipeline map; see
> [04](04-implementation-roadmap.md) for the phased plan.

## Executive summary

Extending `qcode-discovery` from qubit (F₂) to qudit (GF(q)) is **feasible and
front-loaded.** The construction layer is the entry point and the **CSS path is
nearly free**, because qldpc's `BBCode` is field-generic and auto-inserts the
antipode + sign (`matrix_z = [Bᵀ, −Aᵀ]`) needed for odd-`q` CSS commutativity —
confirmed live that `BBCode(field=3)` yields `[[72,8]]` with
`is_subsystem_code=False` and zero symplectic residual.

The single most **pervasive** edit is the genotype change from coefficient-free
`(x_exp, y_exp)` tuples to coefficient-carrying `(x_exp, y_exp, coeff)` triples —
small at each site but load-bearing and touching ~6 clusters. The **real work**
is the distance layer: the CSS BP-OSD bound breaks indirectly (hardcoded
`bp_method` kwargs reach `GUFDecoder` over GF(q) → `TypeError`; fix by gating on
`field.order == 2`), and the **prime-q mod-q MILP is the trustworthy distance
signal and the true critical path** (~2–3 days hardened).

**MVP = prime q (start q=3), CSS arm only, evolving + MILP/exact-verifying at
least one new qudit code.** The full pipeline (non-CSS PBB, prime-power, qudit
Clifford/LC equivalence) is research-grade and reasonably deferred. A believable
cross-concern MVP is **~1.5–2.5 weeks**, not the ~1 week implied by summing
per-concern minimums.

---

## Math framing: qudit BB = the GF(q) case of 2BGA / GB codes

Qudit bivariate bicycle (BB) codes are the abelian-cyclic special case of
**two-block group-algebra (2BGA)** codes over the group algebra `F_q[G]` with
`G = Z_ℓ × Z_m`, i.e. the ring `R = F_q[x,y]/(x^ℓ−1, y^m−1)`. Grounded in
Lin & Pryadko, *Quantum two-block group algebra codes* (arXiv:2306.16400, Phys.
Rev. A 109, 022407), and the qudit twisted-torus BB construction
(arXiv:2602.04443).

**Field requirement — prime-power `q` only.** `q = pᵐ` must be a **prime power**
so that GF(q) is a field (both qldpc and `galois` require this). **Composite-`d`
(the ring `Z_d`, `d` not a prime power — e.g. 4, 6) is a distinct modular-qudit
family that qldpc cannot reach** and that needs ring linear algebra (Smith /
Howell normal form). Composite-`d` is **out of scope** — flag and defer.

**CSS stabilizer form (verified).** For circulants `A = L(a)`, `B = R(b)` from
`a, b ∈ R`:
- `H_X = (A, B)`
- `H_Z = (Bᵀ, −Aᵀ)` — equivalently `(b̂, −â)` where `^` is the **antipode**
  `g → g⁻¹` (`x → x⁻¹`, `y → y⁻¹`). Lin–Pryadko Eq. (41): `L(a)ᵀ = L(â)`, so the
  transpose of a circulant *is* the antipode; Eq. (16): the minus sign sits on
  the second block.

Load-bearing facts verified live on GF(3) `[[72,8]]`:
- The **relative opposite sign** between the two Z-blocks (plus antipode on both)
  is essential; a global sign is convention. Over F₂, `−1 = +1`, so it collapses
  to `(Bᵀ, Aᵀ)` — the minus is invisible only in characteristic 2.
- qldpc **auto-inserts this** on the CSS/`BBCode` path (`matrix_z` carries
  `2 = −1 mod 3` on half its nonzero entries). **No hand-derived sign needed.**
- For the hand-rolled **PBB symplectic**, `[Bᵀ | −Aᵀ] % q` (and `[−Bᵀ | Aᵀ] % q`)
  both give `k=8`, `is_subsystem_code=False`, zero residual; the current
  `[Bᵀ | Aᵀ]` (naive GF(2)) gives a wrong `k=36` subsystem code. **One-line fix.**

**Genotype must carry GF(q) coefficients.** 2BGA elements live in `F_q[G]` with
coefficients in `F_q`, not just `{0,1}`. A term is `c · x^a · y^b` with
`c ∈ GF(q)*`. Coefficients are load-bearing: over GF(3), `x³+y+y²` gives `k=8`,
but flipping a **single** monomial coefficient (`x³+2y+y²`) drives `k → 0`.
(Scaling the **whole** polynomial by a unit preserves `k` — a distinction an
implementer must respect when writing tests.)

**Commutativity / dimension over GF(q).** Abelian `G ⇒ L(a)=R(a) ⇒` commutativity
is automatic (Lin–Pryadko Sec. IV.2). Field-aware `k = n − rank_q(H_X) −
rank_q(H_Z)` is exactly `code.dimension` (verified `72−32−32=8`).

**Prime vs prime-power arithmetic (a correctness landmine, not a knob).** For
`q = pᵐ` you must use field arithmetic, **never integer `% q`**: in GF(4),
`2·3 = 1` and `3+3 = 0` in-field but `2` by integer mod — any residual `% q` on a
`galois.FieldArray` silently corrupts arithmetic without crashing. For **prime
`q` (2,3,5,7), integer mod-p coincides with field arithmetic**, which is why
prime `q` is the safe MVP.

---

## Change inventory

Difficulty: **T**=trivial, **E**=easy, **M**=moderate, **H**=hard, **R**=research.
MVP: ✓ = in the CSS prime-q MVP; — = deferred.

| File / location | Change | Diff | MVP |
|---|---|---|---|
| `bb_code.py:35-47` `terms_to_poly` | Accept `(x,y,coeff)` triples; build `coeff*x**xe*y**ye`; 2-tuple→coeff=1 shim | E | ✓ |
| `bb_code.py:84-103` `build_bb_code` | Add `field:int=2`; pass `field=` to BBCode; assert `code.field.order==field` | T | ✓ |
| `bb_code.py:50-81` `validate_terms` | Coefficient-aware: combine like monomials mod q, drop zero-sums; decouple term-count bound | E | ✓ |
| `bb_code.py:106-112` `get_code_params_fast` | No math change (already field-aware); docstring qubits→qudits | T | ✓ |
| `pbb_code.py:159` `block2_z` | **One-line sign fix** `[mat_B.T, -mat_A.T] % q`; pass `field=q` to QuditCode | E | — |
| `pbb_code.py:47-65` poly/matrix | Accept coeff triples; build over `field=q` so `lift()` returns GF(q); drop `% 2` | E | — |
| `pbb_code.py:68-73` `check_commutativity` | In-field `symplectic_conjugate(M)@M.T==0` + assert `not is_subsystem_code` | M | — |
| `pbb_code.py:173-263` `_gf2_*` logicals | Replace with `galois.GF(q)` row_reduce/null_space, or prefer kwarg-free `get_logical_ops()` | M | — |
| `distance.py:103-107` `estimate_distance` | Gate `bp_method` on `field.order==2`; for q>2 call kwarg-free `get_distance_bound_with_decoder(Pauli.X/Z, trials)` (CSS only) | E | ✓ |
| `distance.py:127-133` `estimate_distance_osd_cs` | Gate on q==2; alias to `estimate_distance` for q>2 (no q-ary OSD-CS in ldpc) | E | ✓ |
| `distance.py:61-83,139-188` exact/timeout | Replace SIGALRM with OS-level subprocess timeout; gate exact to small `q^k` | E | ✓ |
| `distance_milp.py:38-65` `get_code_matrices` | `% 2`→`% q`; weight via `np.count_nonzero` not `np.sum` | E | ✓ |
| `distance_milp.py:68-138` `ilp_min_weight` (CSS) | **PRIME q only**: ub=q−1, big-M indicator `(q-1)w_j ≥ x_j`, slack `-2`→`-q`, residue-loop/big-M anticommutation; **assert q prime** | H | ✓ |
| `distance_milp.py:337-471` symplectic MILP | Same deltas + **signed** `s_x·z − s_z·x mod q`; PRIME q only | H | — |
| `distance_milp.py:141-329,474-598` drivers | Thread q; hard guard refusing prime-power q | M | ✓ |
| `distance_bposd_noncss.py` (whole) | Re-architect: `q²−1` per-site Paulis, mod-q hash as linear-combination search, GUFDecoder wrapper | R | — |
| `evaluator.py:86-87,218-262` trust/distance | Add q; for q>2 set `distance_trusted` only on MILP/exact corroboration; branch BP-OSD off q==2 | M | ✓ |
| `evaluator.py:162-169` self-dual gate | Compare `(exp,coeff)` tuples; gate `A==B⇒d=2` on q==2 | M | ✓ |
| `evaluator.py:545-549` `_milp_cache_key` | Include coeff in canonical key | E | ✓ |
| `results.py:158-181` `_code_key` | Include coeff per term; record q | E | ✓ |
| `tanner_equivalence.py:132-183,297-336` | Structural value-gadget for coeff edges; `==1`→`!=0`; `%2`→`%q`; optional scalar canonicalization | H | — |
| `clifford_equivalence.py:247-278` rank | galois GF(q) rank (decomposability works for free) | E | ✓ (light) |
| `clifford_equivalence.py:89-787` LC machinery | Replace 6-Clifford/2-coloring with `Sp(2,q)` + multiplier orbit search | R | — |
| `evolve/seed_solution*.py` | `(x,y)→(x,y,coeff)`; read q from `QCODE_FIELD`; coeff-mutation; recompute KNOWN_CODES per q | M–H | ✓ |
| `evolve/openevolve_evaluator*.py` | Read `QCODE_FIELD`; thread `field=q`; coeff-aware dedup/classify; retune k/d/trust; distance-path gating | M | ✓ |
| `evolve/run_evolution.py` | `--field` arg; set `QCODE_FIELD` env | E | ✓ |
| `evolve/prompt_context*.md`, config `system_message` | Rewrite F₂→F_q: q-ary coeffs, antipode/sign, qudit baselines; remove F₂-only lore | H (prose) | ✓ |
| `scripts/verify_*.py` BLISS copies | De-dup into shared field-aware `tanner_equivalence`; thread q | M | — |
| `main.py:48,128,135` | Add `--field`; thread into `evaluate_batch` | E | ✓ |

**New files:** `evaluation/field_utils.py` (centralize `get_field(q)` prime-power
guard, in-field `assert_is_stabilizer_code`, `terms_to_poly_q`,
`combine_like_terms` — forbids raw `%q` on FieldArrays); `evaluation/genotype.py`
(normalized coeff-carrying term type + `canonicalize`, reused by all cache/dedup
keys); `evaluation/distance_qudit.py` (the field-aware dispatcher);
`tests/test_distance_gf.py`, `tests/test_tanner_gf_q.py`; GF(3)/GF(5) reference
fixtures replacing the GF(2) `KNOWN_CODES`.

---

## Distance / decoder strategy over GF(q)

The distance layer breaks in two distinct ways and needs a field-and-structure
dispatcher. The four backends:

**(A) CSS BP-OSD bound (`distance.py`) — indirect break, easy fix.** Hardcoded
`bp_method='product_sum'` reaches `GUFDecoder` over GF(q) and raises `TypeError`.
Gate the kwargs on `code.field.order == 2`; for `q>2` call the **kwarg-free** CSS
GUF bound. ~10 lines, verified working.

**(B) GUF per-sector bound — cheap pre-filter only.** Fast (~0.1 s, mild
`q`-scaling to q=7) but **loose by ~3×** (returns 18 vs true ≈6 on GF(3)
`[[72,8]]`). Since `FOM = k·d²/n`, a 3×-loose `d` inflates FOM ~9×. **Never pass
`max_weight` for q>2** (exponential null-vector enumeration). Use only as a screen;
require MILP/exact corroboration before trusting/saving.

**(C) Prime-q mod-q MILP (`distance_milp.py`) — the trustworthy signal, true
critical path.** Two validated deltas over scipy/HiGHS (no new solver):
1. replace the `−2` commutation/anticommutation slack with `−q`
   (`check·x − q·s == 0`);
2. replace the binary-OR weight indicator with a **big-M indicator**: q-ary vars
   `ub = q−1`, binary `w_j`, constraint `(q−1)·w_j − x_j ≥ 0` (and `≥ z_j` for
   symplectic), minimize `Σ w_j`.

Use the **signed** symplectic form `s_x·z − s_z·x mod q` (sign matters for q>2).
Anticommutation `≠ 0 mod q` is **not** one linear constraint: either loop the ILP
once per residue `r ∈ {1..q−1}` (q−1 exact solves) or add a big-M indicator
(1 looser-presolve solve). **Hard guard: assert q is prime** — the integer-slack
encoding is mathematically **invalid for GF(pᵐ)**.

**(D) Exact enumeration (`get_distance_exact`).** Field-generic but ~`q^k`
(brute-forces over GF(q)). Gate to tiny `(q,n,k)`; **replace SIGALRM with an
OS-level subprocess timeout** (in-process `signal.alarm` around qldpc C calls
crashes the interpreter, exit 144).

**MVP distance recipe (CSS, prime q):** GUF kwarg-free per-sector bound as a cheap
pre-filter → prime-q CSS mod-q MILP as the trusted ranking/save signal → `q^k`-gated
exact (OS timeout) for top candidates. **Keep a corroboration trust gate** (require
MILP confirmation before saving) rather than blindly disabling the `d/√n` filter —
the loose GUF bound with no guard is the single most under-rated risk.

**Non-CSS distance (deferred):** no cheap field-generic non-CSS decoder bound
exists (`QuditCode` lacks `with_decoder`; `get_distance_bound` is ~15 s/code and
loose). The only trustworthy non-CSS prime-q distance is the symplectic MILP
(research-grade). The full `distance_bposd_noncss.py` rewrite is the long pole.

---

## Top risks and mitigations

- **Loose GUF bound with no trust gate (the true under-rated risk).** ~3× loose ⇒
  FOM inflated ~9×; combined with disabling the `d/√n` filter, the search chases
  phantom-high FOMs. → GUF as pre-filter only; require MILP corroboration before
  any FOM is trusted/saved.
- **Prime vs prime-power MILP invalidity (correctness landmine).** Integer-slack
  mod-q MILP is invalid for GF(pᵐ). → hard `assert q prime` in the MILP path;
  route prime powers to GUF + exact only.
- **Prime-power galois-vs-integer pitfall everywhere.** Raw `% q` on a
  `FieldArray` over GF(pᵐ) silently corrupts. → centralize all field arithmetic
  in `field_utils`; for MVP stay on prime q.
- **Coefficient-blind canonicalization (silent merges).** `_milp_cache_key`,
  `_code_key`, openevolve dedup keys, `_classify_pattern`, and the `A==B⇒d=2`
  gate canonicalize on exponents only → over GF(q) distinct codes collide,
  collapsing MAP-Elites diversity. → one `genotype.canonicalize` helper including
  coeff; gate the `d=2` shortcut on q==2 (the `A=B⇒d=2` theorem is GF(2)-specific).
- **API-name errors (would crash before the field issue).**
  `get_distance_bound_with_decoder` is CSS-only; non-CSS distance uses the
  symplectic MILP; logicals via kwarg-free `get_logical_ops()`.
- **SIGALRM crashes the interpreter (exit 144).** → OS-level subprocess timeouts
  for all q>2 exact/bound calls.
- **Seed/baseline mismatch starving the cascade.** GF(2) `KNOWN_CODES` (n,k,d)
  differ over GF(q) (the gross-code analogue is `k=8`, not 12, over GF(3)); the
  `_safety_net_codes` `k>0` guarantee may not hold. → re-verify the safety net and
  recompute baselines on the target q *before* launching.
- **Effort optimism.** A believable cross-concern MVP is ~1.5–2.5 weeks; the
  prime-q MILP is the ~2–3 day critical path. → plan to the realistic estimate;
  gate scope to CSS prime-q first.

---

## Open questions (need a decision before/early in implementation)

1. **FOM normalization across fields:** keep `FOM = k·d²/n`, or weight by
   `log₂(q)` so a GF(q) code's `k` logical qudits count as `k·log₂(q)`
   logical-qubit-equivalents? Propagates to both evaluators, `results.py` Pareto,
   prompts, baselines, trust/early-stop thresholds. **Product call before Phase 4.**
2. **Coefficient search-space sizing:** open the full GF(q)* coefficient axis
   immediately, or de-risk with a coeff=1-only intermediate milestone first?
3. **MILP anticommutation encoding:** residue-loop (q−1 exact solves) vs big-M
   indicator (1 solve)? Benchmark at q=3,5 to pick the default.
4. **q fixed per campaign vs q in the genotype/MAP-Elites feature space?** Affects
   whether q lives in env or in each candidate.
5. **Genotype coefficient type:** `int ∈ 1..q−1` (prime q only) vs a true galois
   field element (required for GF(pᵐ))? Decide now to avoid a second migration.
6. **Which qudit best-known BB/PBB codes seed the prompts and references?**
   arXiv:2606.02418, Lin–Pryadko (2306.16400), and a GF(q) baseline sweep must
   supply verified `[[n,k,d]]_q` tuples before prompts can be authored truthfully.
7. **Stage-1 lattice viability is field-dependent** (rank determines k): a one-time
   GF(q) sweep is needed to rebuild `STAGE*_LATTICES` (current exclusions are GF(2)
   artifacts).
8. **Is qldpc's kwarg-free non-CSS `get_logical_ops()` trustworthy for all PBB
   codes over GF(q)**, or is a galois-ported fallback still needed?
9. **Target equivalence relation for dedup:** bare permutation, permutation +
   per-qudit GF(q)* scaling, or full local-Clifford? Decides whether the research
   `qudit_clifford` module is on the critical path for the publication headline.
10. **Scope confirmation:** is prime-power q (GF(4)/GF(8)/GF(9)) wanted near-term,
    and is composite-`d` (Z_d ring) explicitly out for v1?
