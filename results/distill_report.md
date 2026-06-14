# Prime-p magic-state-distillation discovery — report

Arm: `qudit_qec.distill_discovery` · family `T(p,m,k) -> [[p^2 m - k, k, 2]]_p`.
Search: primes=(3, 5), m_range=(1, 2), iters=500, seed=0, max_n=40.  Validity gate: exact cubic triorthogonality. Objective: yield gamma=log_d(n/k) (lower better).  No distance MILP (compute-safe).

## Validation vs known codes (arXiv:2403.06228 / 2408.00436)

| code | n | k | d | gamma | expected | ok |
|---|---|---|---|---|---|---|
| [[20,7,2]]_3 | 20 | 7 | 2 | 1.515 | 1.515 | ✓ |
| [[14,4,2]]_3 | 14 | 4 | 2 | 1.807 | 1.807 | ✓ |
| [[13,5,2]]_3 | 13 | 5 | 2 | 1.379 | 1.379 | ✓ |
| [[23,2,2]]_5 | 23 | 2 | 2 | 3.524 | 3.524 | ✓ |

## Catalog: 31 genuine distillation genotypes (29 trusted-distance, 2 distance-pending; 23 distinct [[n,k,d]] parameter sets)

### Trusted frontier (best by yield gamma; distinct parameter sets)

| [[n,k,d]]_p | gamma | first-seen op |
|---|---|---|
| [[13,5,2]]_3 | 1.379 | family |
| [[26,10,2]]_3 | 1.379 | direct_sum |
| [[39,15,2]]_3 | 1.379 | direct_sum |
| [[33,12,2]]_3 | 1.459 | direct_sum |
| [[40,14,2]]_3 | 1.515 | direct_sum |
| [[27,9,2]]_3 | 1.585 | direct_sum |
| [[34,11,2]]_3 | 1.628 | direct_sum |
| [[7,2,2]]_3 | 1.807 | family |
| [[14,4,2]]_3 | 1.807 | family |
| [[21,6,2]]_3 | 1.807 | direct_sum |
| [[28,8,2]]_3 | 1.807 | direct_sum |
| [[29,7,2]]_3 | 2.051 | direct_sum |

### Distance-pending candidates (gamma is an optimistic lower bound — verify distance on lenore before trusting)

| [[n,k,≤d]]_p | gamma(≤) | op |
|---|---|---|
| [[29,7,≤2]]_3 | 2.051 | direct_sum |
| [[22,5,≤2]]_3 | 2.138 | direct_sum |

## Pareto front (n smaller, k larger, d larger; known-distance codes)

| [[n,k,d]]_p | gamma | op |
|---|---|---|
| [[7,2,2]]_3 | 1.807 | family |
| [[13,5,2]]_3 | 1.379 | family |
| [[21,6,2]]_3 | 1.807 | direct_sum |
| [[26,10,2]]_3 | 1.379 | direct_sum |
| [[33,12,2]]_3 | 1.459 | direct_sum |
| [[39,15,2]]_3 | 1.379 | direct_sum |

## Notes (honest framing)
- Every catalog entry passes the exact cubic triorthogonality gate (`is_triorthogonal`): nothing is triorthogonal-by-assumption.
- `trusted` distances are carried through distance-preserving operators (family d=2, column permutation/scaling, direct sum); puncture yields `d`-upper candidates whose gamma is optimistic.
- **`direct_sum` codes are independent block-stacks of smaller members** — they have the *same* per-qudit yield gamma as their parts (e.g. `[[26,10,2]]` = two `[[13,5,2]]`), so they are not genuinely better distillers. The best per-qudit family member in THIS run is `[[13,5,2]]_3` (gamma=1.379) — but this is a function of the search window: the family's gamma decreases monotonically toward 1 as `m` grows (optimal `k=pm-1`), so there is no finite d=2 frontier, only the asymptote.
- **Beating the d=2 regime on *noise suppression* (nu=d=2), not just overhead, needs distance `d > 2`**: the punctured-Reed–Muller route (`reed_muller_triortho`, validated on the qubit `[[15,1,3]]`) and the qutrit strange/QR route (`distill_strange`, the Golay) — see `docs/08`.
- Universality of the transversal T is the separate per-dimension question of `docs/08` (prime p: any non-Clifford gate is universal — sound here).
