# 02 — Grounding experiments (live evidence)

> Verified output of the `Probe` + `Verify` phases of
> [`workflows/qudit-qec-scope.js`](../workflows/qudit-qec-scope.js). Every claim
> below was produced by running `qldpc 0.2.1`, `ldpc 2.4.0`, `galois 0.4.6` in
> the base environment (`/home/hlamm/miniforge3/bin/python3`), not from memory.
> These facts anchor the scoping in [03](03-qudit-extension-scope.md).

## E1 — qldpc builds **commuting** qudit BB/CSS/non-CSS codes over GF(q)

`qldpc.codes.BBCode(orders, A, B, field=q)` is field-generic. With
`A = x³+y+y²`, `B = y³+x+x²` at `(ℓ,m)=(6,6)`:

| q | n | k | checks |
|---|---|---|---|
| 2 | 72 | 12 | 72 |
| 3 | 72 | 8 | 72 |
| 4 | 72 | 12 | 72 |
| 5 | 72 | 8 | 72 |

- **The stabilizer checks genuinely commute over GF(q)** for all tested
  `q ∈ {2,3,4,5}`. Internally qldpc uses its two-block code `TBCode`:
  `matrix_x = [A, B]`, `matrix_z = [Bᵀ, −Aᵀ]` (`quantum.py:145–146,170–171`).
  The `−Aᵀ` term is qldpc **automatically inserting the antipode (transpose) +
  sign (negation)** required for odd-`q` CSS commutativity — visible as the value
  `2 (= −1 mod 3)` filling exactly half of `matrix_z`'s nonzero entries over
  GF(3), while `matrix_x` stays all-`1`s. **No hand-derived sign is needed on the
  CSS path.**
- **Non-CSS** `QuditCode(matrix, field=q)` also works: a hand-built symplectic
  `[X|Z]` for commuting generators yields a valid stabilizer code over GF(3)/GF(5).
- **Pitfall (confirmed):** computing the symplectic inner product with integer
  `numpy % q` on a GF(4)=GF(2²) `FieldArray` *falsely* reports non-commuting.
  This is a field-vs-integer arithmetic artifact — qldpc's internal
  `symplectic_conjugate` uses true GF(q) arithmetic (the `−1` is the field
  additive inverse). **Lesson: never apply raw `% q` to a `galois.FieldArray`
  over a prime-power field.**

## E2 — Field-aware `k` is free; `galois` does the GF(q) linear algebra

- `k = n − rank_q(H_X) − rank_q(H_Z)`, computed by ordinary Gaussian elimination
  over the field — and this is exactly what `code.dimension` returns. Verified
  live over GF(3): `72 − 32 − 32 = 8 = code.dimension`.
- `galois.GF(q)` supplies everything the pipeline needs: `matrix_rank`,
  `row_reduce` (RREF), `null_space`, `left_null_space`, element inverse — for
  **prime and prime-power** `q` (e.g. in GF(4): `2·3 = 1`, `2⁻¹ = 3`).
- `code.get_logical_ops()` over GF(3) returns a `(2k × 2n)` GF(3) `FieldArray`
  (first half = X-support, second half = Z-support); symplectic weight = number
  of qudits with nonzero X- or Z-support.

## E3 — Distance over GF(q): what qldpc actually exposes

Measured on the GF(3) `[[72,8]]` BB code:

| Method | CSS / `BBCode` | non-CSS `QuditCode` | Field-generic | Cost / result (GF(3)) |
|---|---|---|---|---|
| `get_distance_bound_with_decoder(Pauli.X/Z, trials)` **kwarg-free** | **yes** | **no** (AttributeError) | yes (GUF) | ~0.1–0.6 s, returns **18 (loose)** |
| `get_distance_bound(trials)` (both sectors) | yes | yes | yes (GUF) | ~2.5–5.7 s, returns **6 (≈ true d)** |
| `get_distance_exact()` | yes | yes | yes (brute force) | **did not return in 45–60 s**; ~`q^k` |
| `get_logical_ops()` **kwarg-free** | yes | yes | yes → `(2k,2n)` | fast |

Key facts:
- **The binary signature fails over GF(q):** `get_distance_bound_with_decoder(
  Pauli.X, …, bp_method='product_sum', osd_method='osd_cs', osd_order=7)` raises
  `TypeError: GUFDecoder.__init__() got an unexpected keyword argument bp_method`.
  qldpc auto-routes any `FieldArray` with `order ≠ 2` to `GUFDecoder`
  (generalized Union-Find, arXiv:2103.08049), whose only kwargs are `max_weight`
  and `symplectic`. **Fix:** gate the BP-OSD kwargs on `code.field.order == 2`.
- **The cheap GUF per-sector bound is ~3× loose** (returns 18 vs true ≈6) and
  scales mildly with `q` (`Pauli.X, 3 trials`: q=2 → 16 in 0.025 s; q=3 → 21;
  q=5 → 24; q=7 → 26 in 0.26 s). Because `FOM = k·d²/n`, a 3×-loose `d` inflates
  FOM ~9×, so this must be a **pre-filter only**, never the trusted signal.
- **Never pass `max_weight` for `q>2`** — it triggers GUFDecoder's
  worst-case-exponential null-vector enumeration (`field.elements ** len(null)`),
  which timed out at 30 s.
- **Exact distance is the expensive part** (qldpc itself warns "computing the
  exact distance of a non-binary code may take a (very) long time"). Gate it to
  tiny `q^k`, and use an **OS-level subprocess timeout** — in-process
  `signal.alarm` around qldpc C calls crashes the interpreter (observed exit 144).

## E4 — Dedup and equivalence over GF(q)

- **BLISS Tanner dedup is salvageable.** The current code binarizes (`% 2`,
  `== 1`) and would *wrongly merge* GF(3) codes differing only by coefficients
  (verified: binary `hash([[1,1]]) == hash([[1,2]])`). The working fix is a
  **structural value-gadget**: subdivide each `(qudit, check)` edge with a
  dedicated-color coefficient vertex carrying a length-`c` marker path, then hash
  the plain canonical edge list (verified on the real `[[72,8]]₃`). *Caveat:* a
  naive "extra colored vertex per coeff with colors in a separate tuple" encoding
  is broken with igraph 1.0.0 — use the structural gadget.
- **The structural hash does not quotient by per-qudit GF(q)\* scaling** (a qudit
  multiplication-gate Clifford). Whether to canonicalize that out is an
  equivalence-relation choice — see open questions in [03](03-qudit-extension-scope.md).
- **`clifford_equivalence.py` is entirely qubit-specific:** it hardcodes
  `(x,z) ∈ F₂²` with `I/X/Z/Y`, `% 2`, XOR/union-find 2-coloring, and Hadamard
  swapping `(x,z)→(z,x)`. None of this holds over GF(q), where the single-qudit
  Clifford group is `Sp(2,q)` (order `q(q²−1)`). **Research-grade rewrite.**

## What the verifiers confirmed (and corrected)

Three adversarial verifiers (literature, implementation, feasibility) reviewed
the scoping. Net result:

- ✅ **Math framing confirmed** against Lin–Pryadko *Quantum two-block group
  algebra codes* (arXiv:2306.16400): 2BGA/GB codes over `F_q[G]`, antipode
  `[L(a)]ᵀ = L(â)` (Eq. 41), `H_Zᵀ = (B, −A)ᵀ` with the minus on the second block
  (Eq. 16), abelian ⇒ commutativity automatic (Sec. IV.2). Over F₂ the minus is
  invisible (`−1 = +1`). Galois-qudit CSS requires **prime-power q (a field)**;
  composite-`d` is a *distinct* modular-qudit family over the ring `Z_d`.
- ✅ **qldpc reality confirmed live** (commutativity nnz = 0 over GF(3);
  `n − rank_X − rank_Z = 8 = dimension`; non-CSS `QuditCode` commutes
  symplectically; `get_distance_exact` field-generic but slow).
- ⚠️ **Corrections folded into the scoping:** `get_distance_bound_with_decoder`
  is CSS-only; the cheap-non-CSS-GUF idea has no foundation (use the symplectic
  MILP). Coefficients are load-bearing but the *example* matters: scaling the
  **whole** polynomial by a unit preserves `k`, whereas flipping a **single**
  monomial's coefficient (`x³+2y+y²`) drives `k: 8 → 0` over GF(3).
