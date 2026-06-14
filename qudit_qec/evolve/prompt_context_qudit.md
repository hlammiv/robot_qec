# Domain context: qudit (GF(q)) bivariate bicycle codes

Reference knowledge for the LLM mutating `generate_candidates(ell, m)`. This is the
GF(q) analogue of the qubit BB code domain; the new degree of freedom is the
**coefficients**.

## The object

A bivariate bicycle (BB) code over the field `GF(q)` (q prime here) is defined by
two polynomials `A, B` in the quotient ring `R = GF(q)[x, y] / (x^ℓ − 1, y^m − 1)`:

```
H_X = [ A | B ],   H_Z = [ B^T | −A^T ]
```

with `n = 2·ℓ·m` physical qudits. The transpose is the **antipode** (`x → x⁻¹`,
`y → y⁻¹`), and over GF(q) the **minus sign matters** (over F₂ it vanishes, since
−1 = 1). The construction library inserts the antipode and sign automatically — you
only choose `A` and `B`. Because the group `Z_ℓ × Z_m` is abelian, the stabilizers
commute automatically, so every `(A, B)` is a valid CSS code.

This is the abelian case of **two-block group-algebra (2BGA) / generalized bicycle
codes** (Lin–Pryadko, arXiv:2306.16400).

## Genotype

A term is `(x_exp, y_exp, coeff)` with `coeff ∈ {1, …, q−1}`; a 2-tuple
`(x_exp, y_exp)` means `coeff = 1`. The field order `q` is the module constant
`FIELD`. Example: over GF(3), `A = x³ + 2y + y²` is
`[(3, 0, 1), (0, 1, 2), (0, 2, 1)]`.

## Parameters and figure of merit

- `k` = logical qudits, exact via GF(q) rank: `k = n − rank_q(H_X) − rank_q(H_Z)`.
- `d` = code distance, computed **exactly** by a trusted mod-q MILP (not a loose
  decoder bound). `d = 2` is a trap; aim for `d ≥ 6`.
- `FOM = k·d²/n`. Reference qubit codes reach FOM ≈ 12 (the gross code) up to ≈ 19.

## What is genuinely new over GF(q): the coefficient axis

Coefficients are **load-bearing**, unlike F₂ where 1 is the only nonzero scalar:

- `x³ + y + y²` over GF(3) → `k = 8`; flip one coefficient to `x³ + 2y + y²` → `k = 0`.
- Multiplying a *whole* polynomial by a unit `c ∈ GF(q)*` gives the *same* code
  (a trivial relabeling) — do not waste the search on global rescalings.
- Mixing coefficients **within** a polynomial, or **differently between** `A` and
  `B`, accesses codes with no F₂ analogue. This is the region most worth exploring.

## Useful structural families (starting points)

- **x/y-swap**: `A = x^a + y^b + y^c`, `B = y^d + x^e + x^f` (mixed x and y in each).
  The only qubit trinomial family reaching `d ≥ 6`; a strong starting structure.
- **Doubling**: `c = 2b`, `f = 2e` (e.g. the gross code) — often high distance.
- **Constant-monomial / HGP**: `A = f(y)`, `B = g(x)` with constant term 1 — reaches
  high `k` but distance collapses to `d ≤ 4`. High-rate, low-distance.

## Traps the evaluator flags (avoid)

- `A == B` ⇒ low distance (a structural distance trap).
- **Decomposable** codes (Tanner graph splits into a direct sum) offer no advantage
  over their independent pieces — the evaluator reports `decomposable`; prefer
  **indecomposable** codes.
- High-`k` codes are usually low-`d`; the rate–distance tradeoff is real. The prize
  is an *indecomposable* code with both high `k` and high `d`.

## Practical guidance for `generate_candidates`

- Return a list of `(A_terms, B_terms)` pairs; keep counts per lattice modest
  (hundreds to low thousands) — the evaluator k-screens all of them and runs the
  exact MILP on the most promising.
- Seed from the structural families above, then **mutate exponents and
  coefficients**. A single mutation that expresses an algebraic pattern (e.g. "use
  `x^{ℓ/3}`" or "set the middle coefficient to 2") generalizes across lattices.
- Respect the lattice: `0 ≤ x_exp < ℓ`, `0 ≤ y_exp < m` (the library reduces
  out-of-range exponents, but staying in range keeps the genotype clean).
