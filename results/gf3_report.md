# GF(3) qutrit BB-code baseline sweep — report

> Produced by the trusted `qudit_qec` pipeline (`scripts/gf3_baseline_sweep.py`)
> plus an independent adversarial verification workflow (14 verifier agents +
> reviewer). Raw catalog: [`gf3_catalog.json`](gf3_catalog.json); machine-readable
> trustworthy baselines: [`gf3_known_codes.json`](gf3_known_codes.json).

## What was run

A sweep over 7 lattices (`n = 36 … 144`) generated candidates with the qudit seed
(GF(3), so the coefficient axis is exercised), exact-screened **~7,700 valid
qutrit codes** by GF(3) rank, and ran the trusted prime-q MILP on the most
promising **indecomposable** candidates per lattice — **39 distinct codes with
computed distances**. The 14 highest-FOM / certified codes were then independently
re-verified by three orthogonal methods: qldpc's GUF + GAP/QDistRnd decoder bound,
a **from-scratch GF(3) information-set-decoding (ISD)** search using only the raw
check matrices, and MILP weight-cut certification.

## Headline (honest)

**The sweep's top-FOM entries were inflated and are now corrected.** The three
highest-FOM catalog codes were uncertified MILP *incumbents*; verification found
explicit weight-12 logicals (validated `H·v ≡ 0 mod 3`, anticommuting with a
logical) by two independent methods, refuting their distances:

| Code (claimed) | claimed FOM | verified | corrected FOM |
|---|---|---|---|
| `[[108,6,≤15]]₃` | 12.5 | `d ≤ 12` | ≤ 8.0 |
| `[[144,6,≤16]]₃` | 10.7 | `d ≤ 12` | ≤ 6.0 |
| `[[108,6,≤13]]₃` | 9.4 | `d ≤ 12` | ≤ 8.0 |

This is the paper's distance-overestimation lesson reappearing as loose MILP
incumbents — and the pipeline's **trust gate flagged all three as `trusted=false`,
never reporting the inflated FOMs as fact.** For BB codes `d = O(√n)` (√108 ≈ 10),
so the claimed d=15/16 were correctly suspect.

- **Highest CERTIFIED-EXACT FOM = 5.33** — the `[[72,6,8]]₃` codes at lattice (6,6).
- **Highest TRUSTWORTHY FOM = 8.0**, but only as an *upper bound* (`d ≤ 12` at
  n=108/144) — a ceiling, not an achievement.
- All certified-exact codes are small (`n ∈ {36, 54, 72}`); **no n ≥ 108 distance
  could be certified** (MILP proved 0 logicals optimal there; brute force / exact
  enumeration intractable). Treat every `d_status = milp_incumbent` as suspect.

## Three trust tiers

### Tier 1 — CERTIFIED EXACT (proven distance)
| Code | lattice | A | B | FOM | indecomp. |
|---|---|---|---|---|---|
| `[[72,6,8]]₃` | (6,6) | `x³+x⁵y+y²` | `y³+x+x²` | **5.33** | yes |
| `[[72,6,8]]₃` | (6,6) | `x³+x⁵y²+y` | `y³+x+x²` | 5.33 | — |
| `[[54,6,6]]₃` | (9,3) | `1+y+x⁴` | `1+y+x` | 4.0 | yes |
| `[[36,4,6]]₃` | (3,6) | `1+y+x` | `1+y²+x²` | 4.0 | yes |
| `[[36,4,6]]₃` | (6,3) | `1+y+x³` | `y+x+x²` | 4.0 | yes |
| `[[36,4,6]]₃` | (6,3) | `1+y²+x³` | `y+x+x²` | 4.0 | — |

(`[[72,6,8]]a` proven by two-sided weight-cut MILP infeasibility — no nontrivial
logical of weight ≤ 7 in either sector — plus a weight-8 witness, GUF-corroborated.)

### Tier 2 — CONFIRMED UPPER BOUND (multi-method agreement, not proven)
| Code | lattice | FOM (≤) | notes |
|---|---|---|---|
| `[[144,8,≤12]]₃` | (12,6) | ≤ 8.0 | GUF + MILP + ISD all converge on 12 |
| `[[144,8,≤12]]₃` | (6,12) | ≤ 8.0 | |
| `[[108,6,≤12]]₃` | (9,6) | ≤ 8.0 | |
| `[[108,6,≤11]]₃` | (9,6) | ≤ 6.72 | **d bracketed in [7,11]** — `d ≥ 7` *certified* |
| `[[144,6,≤10]]₃` | (6,12) | ≤ 4.17 | |

### Tier 3 — REFUTED INCUMBENTS
The three codes in the headline table. Their corrected `d ≤ 12` is itself an upper
bound (not certified) — the true distance could be lower; only the *claim* is
disproven.

## Notable

- **`[[72,6,8]]₃` is a genuine certified discovery** of our pipeline (mixed-monomial
  `x⁵y` term) — added to the campaign seed as a structural baseline.
- One verifier *certified a lower bound* `d ≥ 7` for `[[108,6,11]]₃` (all 12 logicals
  proven infeasible at weight ≤ 6) — the only sub-exact code with a nontrivial
  certified lower bound.
- **`GAP`/`QDistRnd` is installed** in the environment — a peer-reviewed independent
  GF(q) distance tool (the engine behind qldpc's `get_distance_bound`). Worth adopting
  as a second trusted distance source alongside the MILP.

## Cautions / limitations

- **Do not seed or report the catalog's `best_by_fom` leaders** (FOM 12.5/10.7/9.4):
  they are refuted incumbents. `gf3_catalog.json` ranks by inflated incumbents.
- All `n ≥ 108` distances are **upper bounds**, not certified — a hidden lower-weight
  logical cannot be rigorously excluded.
- MILP incumbents are **systematically loose** for GF(3) BB codes at `n ≥ 108`
  (HiGHS couldn't certify any logical even at 90–500 s); always cross-check with an
  independent method.
- Coverage is partial: 14 of 39 distinct codes verified (the high-FOM tail);
  lower-FOM entries inherit their original MILP status.
- A few large-code A/B triples were reconciled against `gf3_catalog.json`; the two
  `(6,3) [[36,4,6]]` variants may coincide under relabeling.

## Methodology note

This sweep validates the **end-to-end discovery+verification toolchain** over GF(3),
not a maximal search: only the trusted-MILP-confirmable regime (small n) yields
certified results, and the most valuable outcome is the demonstrated ability to
*discover candidates and then correctly distrust the loose ones*. Reproduce with:

```bash
python scripts/gf3_baseline_sweep.py --field 3 --lattices "6,6 9,6 12,6" --out shard.json
```
