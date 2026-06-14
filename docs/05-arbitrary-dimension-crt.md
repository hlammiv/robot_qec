# 05 — Arbitrary dimension via CRT factoring

> Decision (2026-06-14): the pipeline should ultimately target **arbitrary qudit
> dimension `d`**, reached by **CRT-factoring** `Z_d` into prime-power pieces and
> reusing the per-factor machinery. This document makes that route precise,
> grounds it in a verified computation, and amends the scope/roadmap
> ([03](03-qudit-extension-scope.md), [04](04-implementation-roadmap.md)).

## The CRT decomposition (verified, exact, complete)

A physical qudit of dimension `d` is a `d`-level system with clock/shift Paulis
mod `d` — i.e. a **modular qudit** over the ring `Z_d`. For `d = ∏ᵢ pᵢ^{aᵢ}`,
the Chinese Remainder Theorem gives a ring isomorphism

```
Z_d  ≅  Z_{q₁} × Z_{q₂} × ···        with  qᵢ = pᵢ^{aᵢ}  (pairwise coprime)
```

and a matching tensor factorization of the Hilbert space `C^d ≅ ⊗ᵢ C^{qᵢ}` and of
the generalized Pauli group `P_d ≅ ⊗ᵢ P_{qᵢ}`. (A `d=6` qudit *is* a qubit ⊗
qutrit.)

For any stabilizer code defined by a check matrix `H` over `Z_d` — which is
exactly the BB/CSS setting — **the whole code factors through CRT**:

```
ker_{Z_d}(H)  ≅  ∏ᵢ ker_{Z_{qᵢ}}(H mod qᵢ)
```

so stabilizers, logicals, `k`, and distance all decompose per prime-power factor.
This is **not a heuristic that covers "most" codes — it is exact and complete for
check-matrix-defined modular-qudit codes.** Verified live (brute-force kernel
counts over `Z_d^4`):

| d | factors | `|ker_{Z_d}|` | `∏ᵢ|ker_{Z_qᵢ}|` | match |
|---|---|---|---|---|
| 6 | 2·3 | 36 | 36 | ✅ |
| 10 | 2·5 | 100 | 100 | ✅ |
| 15 | 3·5 | 225 | 225 | ✅ |
| 30 | 2·3·5 | 3600 | 3600 | ✅ |

**Consequences we get for free from CRT:**
- **`k` per factor:** the code encodes `kᵢ` logical `Z_{qᵢ}`-qudits in factor `i`.
  A *uniform* `d`-dimensional logical qudit exists only where the factor logicals
  line up; in general the logical space is a product of possibly-different pieces.
- **Distance = min over factors:** `d(Z_d code) = minᵢ d(Z_{qᵢ} code)`. A
  nontrivial logical can live entirely in the weakest factor, so a composite-`d`
  code is **only as strong as its weakest prime-power factor.** (This is the same
  phenomenon as the paper's direct-sum `[[288,24,12]] = gross ⊕ gross`
  decomposability — CRT is a graded version of it.)

## The one subtlety that decides the backend: square-free vs not

CRT splits `d` into **prime-power moduli** `Z_{pᵢ^{aᵢ}}` — and
`Z_{p^a}` for `a>1` is **a ring, not the field `GF(p^a)`** (`Z_4 ≠ GF(4)`: `Z_4`
has zero divisors; `GF(4)` does not). `galois`/`qldpc` represent the **field**
`GF(p^a)`, never the **ring** `Z_{p^a}`.

So CRT reduces composite `d` to *prime fields* **iff `d` is square-free**:

| `d` | CRT factors | every factor a field? | Covered by prime-field machinery? |
|---|---|---|---|
| 6, 10, 15, 21, 30, 35, … (**square-free**) | distinct primes `Z_{pᵢ}=GF(pᵢ)` | ✅ yes | ✅ **fully** (just add a CRT layer) |
| 4, 8, 9, 16, 25, … (**prime power `p^a`**) | single `Z_{p^a}`, CRT-irreducible | ❌ no (`Z_{p^a}` ring) | ❌ needs a decision (below) |
| 12, 18, 20, 24, … (mixed) | e.g. `Z_4 × Z_3` | ❌ partial | ❌ ring factor remains |

Verified: `galois.GF(d)` and `qldpc BBCode(field=d)` succeed for prime-power `d`
(2,3,4,5,7,8,9,16) and fail for composite non-prime-power `d`
(6,10,12,15) — *"order must be a prime power"*. The prime-power successes are the
**field** `GF(p^a)`, which is a *different code family* from the **modular ring**
`Z_{p^a}`.

## What "arbitrary d via CRT" therefore delivers, in tiers

- **Tier 1 — square-free `d` (6, 10, 15, 21, 30, …): FULLY in scope, cheap.**
  Once the prime-`q` CSS MVP exists, a thin CRT layer (`crt.py`) delivers arbitrary
  square-free dimensions: reduce the `Z_d` genotype mod each distinct prime `pᵢ`,
  run the existing **prime-field** pipeline per factor, recombine (`k` per factor,
  `d = minᵢ dᵢ`, distance MILP is valid because each factor is a prime field).
  **This is the bulk of "arbitrary dimension" at low marginal cost.**

- **Tier 2 — prime-power factors `p^a` (a>1): a genuine fork** (the remaining
  open decision). For a factor like `Z_4`/`Z_8`/`Z_9` you must choose the object:
  - **(2a) Galois-qudit `GF(p^a)`** — a `p^a`-level system with *field* Paulis.
    `qldpc`/`galois` build it and give `k` **for free**; the added work is only the
    `GF(p^a)`-valid distance path (prime-subfield MILP or GUF+exact — the existing
    "Phase 7 prime-power" item). Cheaper, but it is **not** the clock-mod-`d`
    physical qudit.
  - **(2b) Modular-qudit `Z_{p^a}`** — the *physical* clock/shift-mod-`d` qudit.
    Requires a **new ring backend**: `k`/logicals via **Smith/Howell normal form**
    over `Z_{p^a}` (galois/qldpc cannot help); the distance MILP is, conversely,
    *valid* here because integer mod-`p^a` **is** the native ring arithmetic.
    This is the larger research track.

## Proposed CRT architecture

```
qudit_qec/crt.py
  factor_dimension(d)            -> {p: a}            (trial division)
  crt_moduli(d)                  -> [q_i = p_i^a_i]   (coprime prime-power moduli)
  classify(d)                    -> 'prime' | 'square_free' | 'prime_power' | 'mixed'
  split_genotype(terms_Zd, q_i)  -> terms over each factor (coeffs reduced mod q_i)
  build_crt_code(d, A, B, ...)   -> [per-factor code objects]   # dispatch field/ring
  combine(results)               -> {n, k_per_factor, k_uniform, d = min(d_i), fom}
```

The **genotype is `Z_d`** (coefficients in `Z_d`); evaluation CRT-splits, scores
each factor with the appropriate (field or ring) backend, and recombines. A
natural new search objective falls out: favor `Z_d` codes whose factor reductions
have **balanced `k` and distance**, so every logical is a full `d`-dimensional
qudit rather than a lopsided product.

## Roadmap impact (amends [04](04-implementation-roadmap.md))

- **New Phase 4.5 — CRT layer for square-free `d`** (~1–2 days *after* the prime-q
  CSS MVP): `crt.py` + per-factor dispatch + recombination + a `[[·]]_6`
  (qubit⊗qutrit) end-to-end demo. Delivers **arbitrary square-free dimension** with
  no new field/ring math — pure reuse of the prime-field pipeline. High value,
  low cost.
- **Phase 7 (was "prime-power")** splits into the **(2a) `GF(p^a)` field** track
  (medium — distance path only) and the **(2b) `Z_{p^a}` ring** track (large — SNF
  backend). The choice is the open decision below.
- Composite `d` moves from *out of scope* to **in scope via CRT**: square-free now
  (Tier 1), prime-power factors per the Tier-2 fork.

## Resolution (2026-06-14): both tracks, ordered

Prime-power dimensions will be supported **both ways, sequenced cheapest-first**:

1. **First — Galois-qudit `GF(p^a)` (field).** Reuses the existing field machinery;
   qldpc/galois give construction + `k` for free, so only the `GF(p^a)`-valid
   distance path is new. (Roadmap **Phase 7a**.)
2. **Then — modular-qudit `Z_{p^a}` (ring).** The new Smith/Howell-normal-form ring
   backend that yields the *physical* clock-mod-`d` qudit and completes arbitrary
   `d` for non-square-free dimensions. (Roadmap **Phase 7b**.)

Net dimension coverage, in build order:
`prime fields` (MVP) → `square-free composite` via CRT (Phase 4.5, free) →
`GF(p^a)` field qudits (Phase 7a) → `Z_{p^a}` ring qudits (Phase 7b) ⇒ **every
integer dimension `d ≥ 2`**, with composite `d` always assembled by CRT from its
prime-power factors.
