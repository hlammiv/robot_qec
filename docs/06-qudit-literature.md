# 06 — Qudit QEC literature landscape (what's already known)

> Cited survey across **all** qudit dimensions, to ground novelty claims for the
> arbitrary-`d` discovery+certification pipeline. Produced by two web-grounded
> multi-agent surveys (qutrit-focused + cross-dimension), each with adversarial
> fact-checking, on 2026-06-14. Citations are arXiv / Error Correction Zoo (ECZoo)
> / codetables.de / journals. **Bottom line: qudit QEC *theory* is mature and
> dimension-generic; explicit, distance-*verified* qLDPC/BB *catalogs* are
> qubit-dominated and thin out fast with dimension — and are essentially empty for
> composite (modular `Z_d`) dimensions.**

## Foundations — mature, settled, dimension-generic

Nonbinary/Galois stabilizer codes are decades old and uniform in the field order:
- CSS construction (Calderbank–Shor; Steane, 1996) — generalizes to any field.
- **GF(4)↔qubit** correspondence (Calderbank–Rains–Shor–Sloane, 1998) — the binary
  template.
- **Qudit fault tolerance** for prime `d` (Gottesman, 1998).
- **Nonbinary stabilizer codes over GF(q)** for any prime power `q`
  (Ashikhmin–Knill, 2001) and the definitive "Galois theory" (Ketkar–Klappenecker–
  Kumar–Sarvepalli, 2006): trace-symplectic GF(q) ⇄ trace-Hermitian GF(q²)
  characterization, CSS, BCH/QR/Reed–Muller families.

Small qutrit codes are textbook: `[[5,1,3]]_q` perfect code (all `d`),
`[[9,1,3]]`/`[[7,1,3]]` over any field, the qutrit-native `[[11,1,5]]_3` ternary
Golay (used in magic-state distillation). ECZoo codifies the families cleanly:
`galois_stabilizer`, `galois_true_stabilizer`, `galois_css` (field) vs `qudit_css`,
`modular_qudit` (ring `Z_d`).

## Status by dimension class

| Class | Examples | Theory | Explicit/optimized codes | Best-known catalog | qLDPC instantiated? |
|---|---|---|---|---|---|
| **Qubit** `d=2` | 2 | complete | exhaustive (self-dual additive GF(4) ≤ len 12; uniqueness/nonexistence proofs) | codetables.de to n=256; Magma QECC to n=35 | **yes, mature** (all BB/GB/2BGA catalogs are binary) |
| **Prime qudit** `d=p` | 3,5,7,11,… | complete (`Z_p=GF(p)`) | good for small `p`; q=3 best | codetables.de advertises q=3,5,7 but **only q=3 populated** (q=5,7 "not available"); quantumcodes.info covers 3,5,7,11,13 | **yes — the active frontier** (Spencer; Halla; Liang–Chen, 2025–26) |
| **Prime-power Galois** `d=p^m, m≥2` | 4,8,9,16,25,… | complete as *construction* | GF(4) excellent; GF(8/9) moderate; GF(16)+ sparse | codetables.de: GF(4),GF(8) yes; **GF(9),GF(16) absent**; quantumcodes.info to 27 | **never searched natively** — only via extension-of-scalars |
| **Composite modular** `d` composite, ring `Z_d` | 6,10,12,15,… | **underdeveloped** (needs Smith/Howell normal form, torsion homology) | logical dim need not be a power of `d`; a few ring constructions (Gunderman 2025; Z₄ topological) | **none** — no `Z_d` qLDPC/BB catalog exists |

## Qutrit BB/qLDPC specifically — new (2025–26), partial, mostly uncertified

Field-generic BB/GB/2BGA machinery has existed for years (Lin–Pryadko **2BGA**,
arXiv:2306.16400; Panteleev–Kalachev GB), but **explicit qutrit instances appeared
only in the last ~8 months**, in three *restricted* families with *unproven*
distances:
- **Spencer et al.** (arXiv:2510.06495, Oct 2025) — generalizes BB/HGP/etc. to
  GF(q); selected q=3,5,7 codes; bespoke MIP decoder (notes the lack of open-source
  qudit decoders).
- **Halla** (arXiv:2602.04443, Feb 2026 — *cited in our `docs/05`*) — the most
  systematic GF(3) tables (~26–80 weight-6 twisted-torus BB codes, q=3,5,7), but
  distances are **QDistRnd upper bounds** (explicitly not exact).
- **Liang–Chen** (arXiv:2602.20158, Feb 2026) — generalized `Z_p` toric codes as a
  BB subclass, prime `p=3,5,7,11`.

All three are **restricted ansätze** (weight-6 twisted-torus / translation-invariant
generalized-toric), **not** a general GF(3) BB sweep over arbitrary monomial pairs +
the full GF(3)\* coefficient axis. `codetables.de` has **no GF(3) quantum table**
("currently not available for q=3").

## Databases / best-known-distance resources

| Resource | Coverage | Note |
|---|---|---|
| codetables.de (Grassl) | Galois only; verified q=2,3,4,8 (q=5,7 "not available", q=9 "wrong input") | gold standard but qubit-dominated; no qLDPC, no rings |
| quantumcodes.info (Aydin–Liu–Yoshino) | q∈{2,3,4,5,7,8,9,11,13,…} via Hermitian/CSS over GF(q²) | only multi-dimension qudit DB; small-n algebraic codes, not qLDPC |
| Magma QECC | binary only, to n=35 | built on Grassl + CRSS |
| Grassl–Roetteler MDS tables | Galois q≤32 (+a q=64 family) | MDS regime only |
| Error Correction Zoo | all families, Galois vs modular | encyclopedia of *families*, **not** an `[[n,k,d]]` lookup |

**qLDPC/BB codes are essentially absent from every best-distance lookup catalog, in
every dimension** — those catalogs are organized around small-n algebraic codes.

## Novelty implications for our pipeline (conservative, by dimension)

- **Qubit (`d=2`): not novel.** Ground-truth / regression only.
- **Prime qudit (`d=p`): low–moderate, bar rising.** q=3,5,7 now have three 2025–26
  efforts; qutrit is the best-covered higher dimension. Real openings: **(a) exact
  distance *certification*** where the literature has only QDistRnd *upper bounds*
  (our `[[72,6,8]]_3` is exactly this — converting a ceiling to a proof, valuable
  even at parameters they list); **(b)** primes `> 11`; **(c)** a general
  monomial-pair + coefficient sweep beyond the restricted weight-6/toric ansätze.
- **Prime-power Galois (`d=p^m, m≥2`): HIGH — the most actionable Galois gap.** *No
  qLDPC/BB code has ever been searched natively over GF(4)/GF(8)/GF(9)/GF(16)…* —
  they exist only via **extension-of-scalars**, which copies the base-prime
  `[[n,k,d]]` verbatim. A native GF(p^m) search (with **field** arithmetic, not
  integer mod-`q` — the correctness landmine in our `docs/05` Phase 7a) that **beats
  the extension-of-scalars baseline** would be clearly new; GF(9)/GF(16) aren't even
  in codetables.de. (Even-order `q=2^h` has an exact bijection to binary codes —
  Ball–Moreno–Simoens 2024 — so genuine novelty there requires escaping that map.)
- **Composite modular (`d` composite, ring `Z_d`): HIGHEST — essentially open.** No
  explicit qLDPC/BB/GB/toric code over any composite `Z_d` (6,10,12,15) exists; the
  family is "largely avoided due to ring-theoretic tools" (Gunderman 2025). This is
  our `docs/05` Phase 7b + CRT layer territory.

## Corrections folded back into our docs

1. **`distance = minᵢ dᵢ` (CRT, `docs/05`)** — the literature (ECZoo, attrib. Grassl
   priv. comm. 2024) states this as an **upper bound** in general. Our own
   derivation + kernel-factorization check give **equality for codes defined by a
   single `Z_d` check matrix reduced mod each factor** (our construction), but the
   honest general statement is `d ≤ minᵢ dᵢ`. `docs/05` updated to reflect this.
2. **CRT factors are prime *powers*, not prime *fields*** — confirmed: `Z_d` splits
   into `Z_{p^a}` (rings), reducing to prime *fields* only when `d` is **square-free**
   (already the basis of our Phase 4.5 vs Phase 7b split — we got this right).

## Key references

- Calderbank–Rains–Shor–Sloane, *QEC via codes over GF(4)* — arXiv:quant-ph/9608006 (1998)
- Gottesman, *FT QC with higher-dimensional systems* — arXiv:quant-ph/9802007 (1998)
- Ashikhmin–Knill, *Nonbinary quantum stabilizer codes* — arXiv:quant-ph/0005008 (2001)
- Ketkar–Klappenecker–Kumar–Sarvepalli, *Nonbinary stabilizer codes over finite fields* — arXiv:quant-ph/0508070 (2006)
- Lin–Pryadko, *Quantum two-block group algebra codes* — arXiv:2306.16400 (2023)
- Spencer et al., *Qudit low-density parity-check codes* — arXiv:2510.06495 (2025)
- Halla, *Qudit twisted-torus codes in the bivariate-bicycle framework* — arXiv:2602.04443 (2026)
- Liang–Chen, *Generalized `Z_p` toric codes as qudit LDPC codes* — arXiv:2602.20158 (2026)
- Ball–Moreno–Simoens, *Stabiliser codes over fields of even order* — arXiv:2401.06618 (2024)
- Aydin–Liu–Yoshino, *quantumcodes.info* database — arXiv:2108.03567
- Grassl–Roetteler, quantum MDS tables — arXiv:1502.05267
- Error Correction Zoo — https://errorcorrectionzoo.org (`galois_css`, `qudit_css`, `modular_qudit`)
- codetables.de (Grassl) — http://www.codetables.de/QECC.php

## Cautions

- **Galois (field GF(q)) vs modular (ring `Z_d`) is the load-bearing distinction**
  and is constantly conflated in casual usage. They are different code families.
- Existing qudit-BB distances (Halla, Liang–Chen) are **probabilistic upper bounds**,
  not certified — the same loose-incumbent failure mode we hit and corrected.
- "No composite-`Z_d` qLDPC catalog" is **absence-of-evidence** as of June 2026
  (checked against arXiv, ECZoo, codetables.de, quantumcodes.info), not a proof of
  absence.
- Scalable qudit decoders are immature (no open-source q-ary BP-OSD); FT/threshold
  analysis for qudit qLDPC is essentially absent.
