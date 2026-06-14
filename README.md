# robot_qec — LLM-guided discovery of **qudit** quantum error-correcting codes

> Goal: build an LLM workflow that discovers *new qudit* (GF(q)) error-correcting
> codes — and decides whether they are useful for fault-tolerant computation — by
> understanding and extending IBM Research's qubit code-discovery pipeline to
> higher-dimensional qudits.

This repository tracks the design **and implementation** of that effort: it starts
from one paper and its open-source code, comes to a precise understanding of how
that pipeline works, scopes the changes needed to search the **qudit** landscape,
and builds them — a trusted code-discovery pipeline, two magic-state-**distillation**
discovery arms, and a **universality** layer — all grounded against the literature
and adversarially verified (218 tests).

## Where this comes from

- **Paper** (ingested): *Evolutionary Discovery of Bivariate Bicycle Codes with
  LLM-Guided Search* — Cruz-Benito, Cross, Kremer, Faro (IBM Research),
  arXiv:2606.02418. See [`literature/`](literature/).
- **Reference pipeline** (studied): `qcode-discovery`
  (https://github.com/qiskit-community/qcode-discovery) — LLM (OpenEvolve /
  MAP-Elites) evolves Python generator ansätze for bivariate-bicycle (BB) and
  perturbed-BB (PBB) codes over `F₂[x,y]/(x^ℓ-1, y^m-1)`, validated by a staged
  cascade (GF(2) rank → BP-OSD → MILP exact distance) plus BLISS Tanner-graph
  dedup and Clifford-equivalence checks.

## The one-line thesis

The reference pipeline's code-construction library, **`qldpc`, is already
field-generic**: `BBCode(orders, A, B, field=q)` builds a *qudit* BB code over
GF(q) today. So the qudit extension is **not** a from-scratch rebuild — it is a
focused set of changes in (a) the evolutionary **genotype** (terms gain GF(q)
coefficients), (b) **distance/decoding** over GF(q) (a *trusted* GF(q) MILP
replaces the binary BP-OSD signal), and (c) threading `field=q` through the
cascade. On top of that "find codes by distance" core we add "is the code *useful*"
— magic-state distillation and gate-set universality.

## What's built — three capability pillars

### 1. Trusted code-discovery pipeline  (find codes by distance / FOM)
`qudit_qec/{genotype,construct,distance,distance_milp,distance_qudit,evaluator,structure,results,crt}.py`
- Field-aware genotype → `BBCode(field=q)` construction → **trusted distance**
  (prime-q mod-q MILP, with a **QDistRnd** second source and a weight-cut
  **lower-bound certificate**) → coefficient-aware Pareto catalog.
- **Arbitrary square-free `d`** via exact CRT factoring (`crt.py`): a `Z_d` code
  splits into prime-power reductions, distance = minᵢ dᵢ (trusted `Z_6` milestone).
- Results: **first MILP-certified qudit code `[[32,2,4]]₃`**; a GF(3) baseline sweep
  with adversarially re-verified certified baselines (e.g. `[[72,6,8]]₃`, FOM 5.33).
  Drivers: `scripts/discover.py`, `scripts/gf3_baseline_sweep.py`.

### 2. Distillation discovery arms  (find codes by magic-state yield / threshold)
`qudit_qec/{distillation,distill_discovery,distill_strange}.py` — see [`docs/08`](docs/08-distillation-arm-and-universality.md).
- **v0 suitability checker** (`distillation.py`): triorthogonality (cubic mod-p),
  transversal-gate level, weight enumerator + MacWilliams, yield `γ` / threshold.
- **d=2 triorthogonal T-gate arm** (general **odd** prime): family
  `T(p,m,k) → [[p²m−k,k,2]]_p`, reproduces arXiv:2403.06228 and generalizes to every
  odd prime; `γ → 1` asymptotically as `m` grows. Driver: `scripts/distill_discover.py`.
- **d>2 Reed-Muller triorthogonal** (`reed_muller_triortho`): reproduces the qubit
  **`[[15,1,3]]`** T-distiller (level-3 T, distance 3 via trusted MILP).
- **d>2 qutrit strange/QR sub-arm** (`distill_strange.py`): self-orthogonal cyclic
  F₃ codes; reproduces the **11-qutrit Golay `[[11,1,5]]₃`** (cubic `ν=3`, threshold
  ≈0.387). Driver: `scripts/distill_strange_search.py`.

### 3. Universality layer  (is `{Clifford + G}` universal at `d = p^m`?)
`qudit_qec/universality.py` — the decidable **p-adic sector-coverage checker** (does a
single-qudit diagonal gate couple all Pauli sectors? `O(d²)`). In-repo finding: CCZ's
single-qudit reductions are Pauli and **don't** couple sectors — necessary-condition
evidence that `{Clifford+CCZ}` is not native-`p^m` universal by those reductions.

> Every numeric claim above is reproduced by a test and was re-derived by independent
> adversarial review (see the `5a20373` verification round). Soundness is gated:
> distance carries a `trusted` flag, distillation is prime-only, and the strange route
> is qutrit-only — non-sound regimes refuse rather than return a number.

## Repository layout

```
robot_qec/
├── README.md                      ← you are here
├── qudit_qec/                     ← the implementation (14 modules + evolve/; see qudit_qec/README.md)
│   ├── genotype, construct, distance, distance_milp, distance_qudit,
│   │   evaluator, structure, results, crt, field_utils   ← trusted discovery pipeline
│   ├── distillation, distill_discovery, distill_strange  ← magic-state-distillation arms
│   ├── universality                                      ← gate-set universality layer
│   └── evolve/                                           ← OpenEvolve LLM-campaign scaffold
├── scripts/                       ← runnable drivers
│   ├── discover.py                  evaluate BB candidates over GF(q) (prime + Galois)
│   ├── gf3_baseline_sweep.py        the GF(3) certified baseline sweep
│   ├── distill_discover.py          d=2 triorthogonal distillation search + report
│   └── distill_strange_search.py    qutrit strange/QR distillation search + report
├── docs/                          ← understanding → scoping → build deliverables (01–08)
├── results/                       ← catalogs + reports (gf3_*, distill_*)
├── tests/                         ← 218 tests (run: `PYTHONPATH=. pytest`)
├── workflows/                     ← the multi-agent LLM workflow that produced docs/
└── literature/                    ← arXiv pointer (paper sources gitignored)
```

## How to run

```bash
# tests (the contract for every claim above)
PYTHONPATH=. pytest -q

# distillation discovery (light, local) — validates vs known codes, writes results/
PYTHONPATH=. python scripts/distill_discover.py          # d=2 triorthogonal arm
PYTHONPATH=. python scripts/distill_strange_search.py    # qutrit strange/QR (Golay)

# evaluate BB-code candidates over GF(q) (heavy distance work → run on a big box)
PYTHONPATH=. python scripts/discover.py --field 3 --in candidates.json
```

Heavy distance/enumeration work (large BB sweeps, distillation beyond the local enum
cap) is designed to run on a compute backend, not the local box; the local search
loops use cheap gates and refuse oversize work.

## Docs

| Doc | What it is |
|---|---|
| [`01`](docs/01-reference-architecture.md) | how `qcode-discovery` works, cluster by cluster |
| [`02`](docs/02-grounding-experiments.md) | live qldpc/galois grounding experiments (evidence) |
| [`03`](docs/03-qudit-extension-scope.md) | the math framing + per-module change inventory |
| [`04`](docs/04-implementation-roadmap.md) | phased plan to an MVP and beyond |
| [`05`](docs/05-arbitrary-dimension-crt.md) | spanning all dimensions via CRT factoring |
| [`06`](docs/06-qudit-literature.md) | cited qudit-QEC landscape, by dimension + novelty |
| [`07`](docs/07-magic-state-distillation-scope.md) | scoping a "can this code distill?" capability |
| [`08`](docs/08-distillation-arm-and-universality.md) | the distillation arms + universality layer |

## Honest findings & limitations

- The d=2 triorthogonal family has **no finite yield frontier** — `γ` decreases
  monotonically toward 1 as `m` grows; `[[13,5,2]]₃` (`γ=1.379`) is only the best in
  the default small sweep. Beating d=2 on *noise suppression* needs `d>2`.
- The small-`n` self-orthogonal **cyclic** F₃ strange family is **sparse** —
  essentially Golay-only locally; broader (non-cyclic, `n>13`) search is a
  compute-backend job.
- The small punctured-RM search finds **no `d>2` qutrit *triorthogonal* code** (only
  the d=2 family) — so the practical d>2 qutrit distiller is the Golay (strange/QR),
  not a triorthogonal T code.
- Deferred: non-CSS PBB codes, qudit-Clifford equivalence dedup, the prime-power
  `GF(p^a)` / modular `Z_{p^a}` distance backends (which complete arbitrary `d`), and
  running a live LLM evolution campaign (`evolve/` is built + tested, needs an LLM
  endpoint).

## License

Apache-2.0 (see `LICENSE`).
