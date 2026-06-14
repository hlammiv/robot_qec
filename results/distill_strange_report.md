# Qutrit strange-state / QR distillation discovery — report

Sub-arm: `qudit_qec.distill_strange` · genotype = self-orthogonal **cyclic codes over F_3** · CSS distiller `Hx = Hz = C` of the qutrit **strange state**.
Search: n=5..13 (gcd(n,3)=1), local enum cap 60000, threshold-refine top 5.  Objective axis: noise-suppression exponent `nu` + threshold `eps_*` (NOT gamma).

## Validation: 11-qutrit Golay `[[11,1,5]]_3` reproduced: ✓ (nu=3 cubic, threshold ~0.387)

## Catalog: 2 genuine strange distillers

| [[n,k]]_3 | nu | threshold eps_* | min stab. weight |
|---|---|---|---|
| [[11,1]]_3 | 3 | 0.387 | 6 |
| [[11,1]]_3 | 3 | 0.387 | 6 |

## Notes (honest framing)
- `nu` is read from the weight enumerator (lowest power of eps in `3A+B`), distinct from the code distance `d` (the Golay has `d=5` but `nu=3`); the `min stab. weight` column is the smallest nonzero stabilizer weight, a distance hint, not the logical distance.
- **The small-n self-orthogonal *cyclic* F_3 family is sparse** — no distiller at n=5,7,13(local); the n=11 entries are the two equivalent ternary QR Golay codes. Broadening to non-cyclic self-orthogonal codes and n>13 is the **lenore-scale** extension (3^(n-k) exceeds the local enum cap).
- This sub-arm beats the d=2 triorthogonal family on noise suppression (`nu=3` cubic vs `nu=2`); it is qutrit-only (the strange state is GF(3)-specific), complementary to the general-prime gamma arm (`distill_discovery`, docs/08).
