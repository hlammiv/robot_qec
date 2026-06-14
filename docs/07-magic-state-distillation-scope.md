# 07 — Scoping: a "can this code distill magic states?" capability

> Investigation (no heavy compute) of how hard it would be to add a magic-state-
> distillation (MSD) capability to the discovery pipeline. Grounded in five papers
> in `literature/` (Prakash–Saha–Singhal–Zurel–de Silva cluster, 2024–2026).

## The computable property that decides "can it distill"

A code distills magic states (via the standard route) iff it admits a **transversal
non-Clifford gate**, which for a CSS code over GF(p) is exactly the condition that
its X-type generator matrix `G` (rows `h⁽ᵃ⁾ ∈ GF(p)ⁿ`) spans a **triorthogonal
space** (Bravyi–Haah for qubits; Krishna–Tillich qudit generalization, used in
arXiv:2403.06228):

```
  ∑_i  h⁽ᵃ⁾_i · h⁽ᵇ⁾_i · h⁽ᶜ⁾_i  ≡  0   (mod p)      for all triples (a,b,c)
```
(plus the per-row weight / lower-order conditions splitting "magic" rows from
stabilizer rows). Triorthogonality ⇒ the diagonal gate `T = Σ_k e^{2πik/p²}|k⟩⟨k|`
is transversal ⇒ the code distills the magic state `|M₀⟩ = Σ_j e^{2πij/p²}|j⟩`, with
noise suppression `ε_out ~ ε_in^d`.

**Crucially this is a *cheap, exact* check** — an `O(κ³·n)` cubic-form test in GF(p)
arithmetic (κ = #generators), the same `galois` linear algebra we already use. It is
a yes/no structural predicate, on par with our existing `is_decomposable`.

### How *well* it distills (the performance/quality criterion)

arXiv:2408.00436 proves the MSD performance (output error, threshold, yield) is
captured by the code's **complete weight enumerator** — and for the qutrit "strange"
state `|S⟩`, by its **simple weight enumerator** (a 2-variable polynomial via the
GF(p²) classical-code correspondence; the relevant object is `B_S(z) − A_S(z)`, the
logical-coset enumerator). From it you read the **yield parameter** `γ` (overhead
`O(log^γ ε⁻¹)`) and the threshold. Computing weight enumerators is exponential in
general but routine at the sizes MSD actually uses (their searches go to **n ≤ 23
qutrits**). So "how good a distiller" is a moderate, **small-n** computation.

## What the five papers give us

| Paper | Contributes | Use to us |
|---|---|---|
| **2403.06228** Low-Overhead Qutrit MSD (Prakash–Saha) | `[[9m−k,k,2]]₃` **triorthogonal** family; `[[20,7,2]]₃` has γ=1.51 (beats all qubit triorthogonal) | the triorthogonality construction + a concrete target to beat |
| **2408.00436** High-Threshold Qutrit MSD Search (Prakash–Singhal) | **weight-enumerator ⇒ performance** theorem; extensive search to n=23 | the exact *objective function* for an MSD search |
| **2510.10852** Sublog Distillation, all prime dims (Saha–Prakash) | punctured **Reed–Muller** codes give γ<1 for qudits | an asymptotic family / alternative genotype |
| **2603.18560** High-threshold MSD with **quantum QR codes** (Zurel et al.) | unifies 5-qubit/Steane/Golay + new high-threshold QR codes | a second family + "is it a QR code" route |
| **2605.30108** Asymptotic MSD, almost-linear rate (Ehara–Takagi) | overhead-exponent vs asymptotic-rate theory | context; not a per-code test |

## Integration difficulty — two very different answers

### Option A — add an MSD-suitability *check* to the evaluator: **EASY (~0.5–1 day)**
A new `qudit_qec/distillation.py`:
- `is_triorthogonal(code, p)` — the cubic mod-p condition on the generators
  (reuses `get_code_matrices` + `galois`); cheap, exact.
- `transversal_T_level(code)` — which diagonal gate level is transversal.
- `magic_state_yield(code)` — weight-enumerator → `γ` and threshold (small-n only;
  gate by `p^k`/`n` like we gate exact distance).
Hook it into `EvalResult` as a flag/score, exactly like `decomposable`. **Low risk,
fully in our wheelhouse.**

**Honest caveat (the important part):** our **BB/qLDPC codes are essentially never
triorthogonal.** BB codes are high-distance, low-weight LDPC *memories*;
triorthogonal/RM/QR distillation codes are a *different, small-distance, specially-
structured* family (the triorthogonal family above has `d=2`!). So bolting an MSD
flag onto the current BB search would return "no" almost always — cheap but
uninformative. Distillation suitability and qLDPC-memory quality are **largely
disjoint design regimes.**

### Option B — a distillation-code *discovery arm*: **MODERATE (~1–2 weeks)**
The valuable version: point the *same* discovery framework at the MSD family.
Reuses the whole skeleton (genotype → construct → evaluate-a-property → catalog),
swapping:
- **genotype** — triorthogonal generator matrices / punctured-Reed–Muller puncturing
  patterns / QR-code parameters (instead of BB polynomial pairs);
- **objective** — the yield `γ` / threshold from the weight enumerator (instead of
  `FOM = kd²/n`), with `is_triorthogonal` as the validity gate;
- **prompt context** — the triorthogonality + weight-enumerator criteria above.
This is a genuinely novel, in-scope target: arXiv:2408.00436 and 2603.18560
**literally do manual/automated searches** for such codes, so an **LLM-guided search
for triorthogonal/QR qutrit distillation codes** — e.g. beating `[[20,7,2]]₃`'s
γ=1.51, or finding an n<23 strange-state distiller better than their search found —
directly parallels (and could extend) the published work. The certification machinery
(exact distance, GF(p) linear algebra) transfers; the weight-enumerator yield is the
new evaluator piece.

## Effort & risk summary

| Piece | Effort | Risk | Notes |
|---|---|---|---|
| `is_triorthogonal` + transversal-gate check | ~0.5 day | low | cubic mod-p; reuses galois |
| weight-enumerator → γ/threshold (small n) | ~2–4 days | low–med | exponential, so `n`-gated; validate vs the papers' tabulated γ |
| evaluator/`EvalResult` hook (Option A) | ~0.5 day | low | mirror `decomposable` |
| MSD genotype + seed + prompt (Option B arm) | ~1 wk | med | new family (triorthogonal/RM/QR), not BB |
| run an MSD discovery campaign | — | (compute) | small-n, so *lighter* than the qLDPC MILP sweeps |

## Recommendation

- **The check itself is easy and low-risk** — and it's correct, prime-and-prime-
  power-friendly (the criterion is mod-p; for prime-power Galois it's the field
  version), and reuses our GF(q) stack.
- But **don't expect our BB codes to distill** — that's a category difference. The
  high-value move is **Option B**: a *distillation-focused discovery arm* that
  searches the triorthogonal/QR/punctured-RM families for high-yield qutrit (and
  general prime-p) distillation codes, validated by triorthogonality + weight-
  enumerator γ. It reuses ~80% of the pipeline and targets exactly what the
  literature searches by hand — a credible novel contribution, and notably **lighter
  compute** than the qLDPC distance sweeps (n ≤ ~23).

## A v0 that needs no heavy compute

`is_triorthogonal` + `magic_state_yield` for **small, known** codes (verify against
the papers: `[[20,7,2]]₃` γ=1.51; the 11-qutrit Golay; 5-qubit/Steane) is a
self-contained, low-compute first step that proves the evaluator piece before any
search. That can be built and unit-tested locally without touching lenore or running
MILP.
