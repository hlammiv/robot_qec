# 08 — Distillation discovery arm + universality layer

> The plan for a magic-state-distillation (MSD) discovery arm, and the
> dimension-dependent **universality layer** it requires. Built on the v0 in
> `qudit_qec/distillation.py` + `qudit_qec/universality.py`, and on two adversarially
> verified, web-grounded surveys (see references). Scopes [`docs/07`](07-magic-state-distillation-scope.md).

## v0 status — built and validated

`qudit_qec/distillation.py` (prime-`p` only, by design):
- `is_triorthogonal` — cubic mod-`p` test `Σᵢ hᵃᵢ hᵇᵢ hᶜᵢ ≡ 0` (Krishna–Tillich), distinguishing the strict triorthogonal *space* from the triorthogonal *matrix* (with `H₁` "magic" rows).
- `transversal_gate_level` — detects the level-3 non-Clifford `T`.
- `weight_enumerator` — simple `A(z)` by capped stabilizer enumeration + MacWilliams `B(z)` over GF(`p²`); refuses above `MAX_ENUM=1e6`.
- `magic_state_yield` — `γ = log_d(n/k)` (triorthogonal mode, any prime `p`) or the qutrit-only strange-state pipeline.

**Validated exactly against the papers**: `[[20,7,2]]₃` γ=1.51, `[[14,4,2]]₃` γ=1.81, the 11-qutrit Golay's exact `A(z)`/`B(z)` + threshold ≈0.387, qubit `[[15,1,3]]₂` γ=2.46. Independently re-derived by an adversarial verifier (no bug found). **Prime-only soundness is enforced** (flags `sound=False` / raises for prime-power & composite).

## The universality question — `{Clifford + CCZ}` is dimension-dependent

**Decisive fact:** the single-qudit Clifford group is a maximal finite subgroup of `U(d)/phase` **iff `d` is prime** (Nebe–Rains–Sloane; Borda–Rincon–Galindo arXiv:2512.20787). This — not the conjectural Wigner-negativity statement — is what makes any non-Clifford gate universal at prime `d`.

| Dimension class | Examples | Clifford maximal? | `{Clifford+CCZ}` universal? | Why |
|---|---|---|---|---|
| **Prime `d`** | 2, 3, 5, 7 | **Yes** (iff prime) | **Yes** | adjoint rep on `sl(d)` irreducible ⇒ any non-Clifford gate universal; `q`-ary CCZ is level-3 non-Clifford |
| **Prime-power, powers of 2** | **4, 8, 16** | **No** (`m≥2`) | **OPEN** (leaning *insufficient*) | adjoint rep reducible: `sl(d)=⊕ₖWₖ` by `p`-adic valuation; need a gate that **couples** the sectors. CCZ never analyzed in the literature; CCZ's single-qudit reductions are Clifford (below) |
| **Prime-power, odd** | 9, 25, 27 | **No** (`m≥2`) | **OPEN** | same reducibility obstruction |
| **Composite, coprime factors** | 6, 15, 35 | n/a (CRT register) | **Universal *without* magic** | inter-factor generalized-CNOTs alone restore irreducibility (2512.20787 Thm 48) — no diagonal magic needed |
| **Composite with a `p^m` block** | 12, 24 | mixed | CRT-CNOTs + per-block resource | coprime CNOTs stitch distinct primes; each `p^m` block (`m≥2`) still needs its own sector-coupling/Galois resource |

**Two inequivalent `p^m` conventions (must be declared per campaign):**
- **Native / modular** (Clifford `= Sp₂(Z_d)`): the reducibility obstruction applies; `{modular-Clifford + CCZ}` is **OPEN**.
- **Galois / Kronecker** (Clifford `= Sp₂(F_d)`, a 2-design for all prime powers; arXiv:2603.02659): a `GF(p^m)` qudit *is* `m` prime-`p` qudits — per-factor prime-`p` magic is universal per factor, **but full universality additionally needs a genuine cross-factor entangling 2-qudit gate** (Brylinski). Per-factor magic alone is **not** sufficient.

**Load-bearing caveat:** *transversal-gate-exists ≠ gate-set-universal* (Eastin–Knill / Bravyi–König). Golowich–Guruswami (arXiv:2408.09254, transversal CCZ over GF(`q`) incl. `q=2`) is a **code/MSD** result — do **not** cite it as a `p^m`-universality theorem.

### The sector-coverage checker — built, and applied to CCZ
`qudit_qec/universality.py` implements the **decidable** part of the 2512.20787 criterion for a single-qudit *diagonal* gate `G` over `Z_d`: does `Ad_G` **couple** all `p`-adic sectors `Wₖ`? (`O(d²)`, via the Fourier supports of `G X^a G†`.) Necessary for `{Clifford+G}` to be universal at `d=p^m`; failure ⇒ definitely not.

Validated (`d=4,8,16,9,27`, 23 tests):
- **Clifford gates** (identity, Pauli `Z`, correct `S`): **never couple** sectors ✓ (they preserve the stratification).
- **Non-Clifford positive control** (a level-bump diagonal): **couples all sectors** ✓.
- **CCZ's single-qudit reductions** are `Z^{bc}` (Pauli) for any fixed controls ⇒ **do NOT couple sectors** for `d=4,8,16`.

**In-repo finding:** CCZ's natural single-qudit reductions are Clifford and supply *no* sector-coupling resource — concrete (necessary-condition) evidence that `{Clifford + CCZ}` does **not** give native `p^m` universality by these reductions. *Still open:* whether a magic-state-injection *effective* gate from CCZ (level 3) couples the sectors — not captured by simple reductions.

## What the discovery pipeline should TARGET (per dimension)

- **Prime `p` (qutrit, GF(5), GF(7)):** any level-3 non-Clifford (qudit-`T` / CCZ) **+ an entangling 2-qudit gate** (our BB-code CNOT structure supplies this). The **clean, lowest-risk** case — our GF(3) catalog + triorthogonality check feed it directly. **Scope the MSD arm here first.**
- **Prime-power `p^m`:** prefer the **Galois-decomposition route** — per-factor prime-`p` magic **+ a genuine cross-factor entangler** + an explicit Galois-convention declaration. (Native route alternative: a *sector-coupling* single-qudit gate — a `T_s` with `s∤d`, or a transposition — found via the sector checker.)
- **Composite, coprime factors:** **inter-factor CNOTs, no diagonal magic** — do *not* spend search budget on triorthogonal/CCZ codes there.

## The MSD discovery arm — plan

Reuse the genotype → evaluate → catalog skeleton; swap the family and objective:
- **genotype**: triorthogonal generator matrices / punctured-Reed–Muller patterns / quantum-QR parameters (not BB pairs).
- **objective (Pareto, *not* a scalar)**: yield `γ` (triorthogonal family) **and** strange/QR **threshold** as **separate** axes — `ν=d` holds *only* for the triorthogonal family (Golay has `d=5` but `ν=3`), so `γ=log_d(n/k)` must not be applied to strange/QR candidates.
- **validity gate**: `is_triorthogonal`; **universality** handled by the layer above (separately — a transversal gate is MSD *input*, not a universality certificate).

**Effort (regime-scoped, honest):**
- Prime-`p` triorthogonal sub-arm + universality layer (prime + corrected `p^m`/composite framing): **~3.5–4.5 days MVP** (high reuse of validated primitives).
- + punctured-RM family: +1–1.5 d (needs the `3r < m(p−1)` bound + one small `d>2` `(n,k,d)` validation first).
- + QR-over-`F_{d²}` family: +1.5–2 d (**new** code — the `ι⁻¹` → prime-`d` stabilizer map; `weight_enumerator` currently raises for prime-power).
- Composite `Z_d` wrapper: +1 d, **new** (the CRT path is BB-hardwired — `crt.evaluate_crt_candidate` consumes `(ell,m,A,B)` BB pairs, *not* raw `H`; cannot be billed as reuse).
- **Total ~6.5–9 days** full.

## Compute safety (after the 2026-06-14 OOM crash)

- Strange-mode weight enumeration is `p^{n−k}` — `3^{n−1}` is **4.8M at n=14**, **31 billion at n=23**. The v0 already caps at `1e6` and refuses; the *plan* must **hard-gate local strange-mode to `n≤13`** and route larger weight-enumerator/distance work to **lenore_remote**. Tier-1 (`γ` with distance passed as an argument) is the only local objective.
- The sector checker is `O(d²)` (`d≤16` → ≤256) — trivially local.

## Open questions

- Native `{modular-Clifford + CCZ}` universality for `p^m` (incl. 4,8,16): **open**; 2512.20787 gives the criterion but never instantiates CCZ. The magic-injection *effective* gate (vs the Clifford simple reductions we checked) is the remaining sub-question.
- The weight-enumerator ⇒ yield theorem: the *complete* enumerator is general odd-prime, but the *simple*-enumerator collapse is **qutrit-only**; `p>3` needs the full `p²`-variable enumerator.
- 2512.20787 is a Dec-2025 preprint — re-verify the `p^m` sector theorems against any published version (its "maximal iff prime" core rests on solid classical results).
- The canonical universality-completing magic target per odd prime `p` (strange-state analogue) is not pinned down.

## References

- Borda–Rincon–Galindo, *Quantum universality in composite dimensions* — arXiv:2512.20787 (2025)
- Nebe–Rains–Sloane, *The invariants of the Clifford groups* — arXiv:math/0001038 (2001)
- Golowich–Guruswami, *Asymptotically good codes with transversal CCZ* — arXiv:2408.09254 (2024)
- Howard–Vala, *Qudit versions of the π/8 gate* — arXiv:1206.1598 (2012)
- Sawicki–Karnas, *Universality of single-qudit gates* — arXiv:1609.05780 (2017)
- Brylinski–Brylinski, *Universal quantum gates* — Phys. Rev. Lett. 89, 247902 (2002)
- *Qudit designs* (Galois-Clifford 2-design) — arXiv:2603.02659 (2026)
- Prakash–Saha, *Low-overhead qutrit MSD* — arXiv:2403.06228 · Prakash–Singhal, *High-threshold qutrit MSD search* — arXiv:2408.00436
