# robot_qec — LLM-guided discovery of **qudit** quantum error-correcting codes

> Goal: build an LLM workflow that discovers *new qudit* (GF(q)) error-correcting
> codes, by understanding and extending IBM Research's qubit code-discovery
> pipeline to higher-dimensional qudits.

This repository tracks the design and implementation of that effort. It starts
from one paper and its open-source code, comes to a precise understanding of how
that pipeline works, and scopes — then builds — the changes needed to make it
search the **qudit** code landscape rather than only the qubit (F₂) one.

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
focused set of changes concentrated in (a) the evolutionary **genotype** (terms
gain GF(q) coefficients), (b) **distance/decoding** over GF(q) (the binary
BP-OSD path must yield to qldpc's qudit decoder + a GF(q) MILP), and (c)
threading a `field=q` parameter through the cascade and dedup. See the scoping
docs below.

## Repository layout

```
robot_qec/
├── README.md                     ← you are here
├── docs/                         ← understanding + scoping deliverables
│   ├── 01-reference-architecture.md   how qcode-discovery works, cluster by cluster
│   ├── 02-grounding-experiments.md    live qldpc/galois experiments (evidence)
│   ├── 03-qudit-extension-scope.md    the math framing + per-module change inventory
│   └── 04-implementation-roadmap.md   phased plan to an MVP and beyond
├── workflows/                    ← the multi-agent LLM workflow(s) used to produce the above
│   ├── README.md
│   └── qudit-qec-scope.js
├── qudit_qec/                    ← our code (scaffold; grows as we implement the roadmap)
│   ├── README.md
│   └── __init__.py
└── literature/                   ← arXiv pointer (paper source gitignored)
```

## The workflow (what "ingest → acquire → understand → scope" means here)

`workflows/qudit-qec-scope.js` is the reusable orchestration that produced the
`docs/`: it fans out parallel readers across the reference codebase to map every
GF(2)/qubit assumption, runs live grounding experiments against `qldpc`/`galois`,
scopes the extension by concern, adversarially verifies the algebra against the
literature, and synthesizes the result. See [`workflows/README.md`](workflows/README.md).

## Status

- [x] Ingest paper (arXiv:2606.02418) and acquire reference repo
- [x] Map the reference architecture + inventory every field assumption → [`docs/01`](docs/01-reference-architecture.md)
- [x] Ground qudit capabilities with live experiments → [`docs/02`](docs/02-grounding-experiments.md)
- [x] Scope the qudit extension (math + per-module changes) → [`docs/03`](docs/03-qudit-extension-scope.md)
- [x] Phased implementation roadmap → [`docs/04`](docs/04-implementation-roadmap.md)
- [x] Arbitrary-dimension strategy (CRT factoring) → [`docs/05`](docs/05-arbitrary-dimension-crt.md)
- [x] **Phase 0**: field substrate + genotype — `qudit_qec/{field_utils,genotype}.py` (round-trips GF(3)/GF(4) `[[72,·]]`, composite-`d` raises)
- [x] **Phase 1**: CSS construction over GF(q) — `qudit_qec/construct.py` (`build_bb_code(field=q)`; coeff flip collapses `k` 8→0)
- [x] **Phase 2** *(critical path)*: distance backend — `distance.py` (GUF pre-filter, field-gated), `distance_milp.py` (prime-q mod-q MILP, **trusted**), `distance_qudit.py` (trust-gated dispatcher). MILP verified vs brute-force (d=2,4) + qldpc exact; **first MILP-certified qudit code: `[[32,2,4]]₃`**. 83 tests.
- [x] **Phase 3**: evaluator cascade + dedup — `evaluator.py` (`evaluate_candidate(field=q)`), `structure.py` (decomposability), `results.py` (`CodeCatalog`, coeff-aware Pareto). End-to-end `[[32,2,4]]₃` `{n,k,d,fom}` trusted; coeff-differing codes kept distinct. 99 tests.
- [x] **Phase 4.5**: CRT layer — `qudit_qec/crt.py` (`evaluate_crt_candidate`). Arbitrary **square-free** `d` (6, 10, 15, 30…) via per-prime-factor reuse of the field pipeline; distance = minᵢ dᵢ. Milestone: trusted `Z_6` (qubit⊗qutrit) code. 124 tests.
- [~] **Phase 4**: OpenEvolve scaffolding **built & tested** — `qudit_qec/evolve/` (coeff-bearing seed, `adapter.evaluate`, `run_evolution --field`, GF(q) config + prompt). Ready to launch; needs `openevolve`+`litellm`+an LLM endpoint to run a live campaign.
- [ ] **Phase 7**: prime-power `GF(p^a)` (7a) then modular `Z_{p^a}` ring backend (7b) → completes arbitrary `d`
- [ ] **Phase 3–4**: cascade + qudit genotype/seeds/prompts; run a GF(3) campaign
- [ ] **MVP**: discover + MILP-verify a new `[[n,k,d]]₃` CSS code

### Key result of the scoping

The CSS qudit path is **nearly free** — `qldpc`'s `BBCode(field=q)` already builds
commuting GF(q) codes (auto-inserting the antipode+sign `H_Z=[Bᵀ,−Aᵀ]`) and gives
field-aware `k`. The real work is the **distance layer** (a GF(q) mod-q MILP is the
trusted signal and the ~2–3 day critical path) and a **coefficient-carrying
genotype**. MVP ≈ 1.5–2.5 weeks, CSS + prime `q`.

**Arbitrary dimension** (`docs/05`): a `Z_d` code factors *exactly* via CRT into
its prime-power reductions (verified), so **square-free `d` (6, 10, 15, …) reduces
to prime fields** and is delivered cheaply by a thin CRT layer on top of the
prime-`q` MVP (distance = min over factors). **Prime-power dimensions `p^a`** get
**both** backends, cheapest-first: Galois-qudit `GF(p^a)` (field; Phase 7a), then
modular-qudit `Z_{p^a}` (physical clock-mod-`d`, Smith-normal-form ring backend;
Phase 7b) — together with the CRT layer this **covers every integer dimension
`d ≥ 2`**. Build order: prime fields → square-free composite → `GF(p^a)` →
`Z_{p^a}`. Non-CSS PBB and qudit-Clifford equivalence remain deferred research
tracks. Full details in `docs/`.

## License

Apache-2.0 (see `LICENSE`). Chosen to match the surrounding QEC ecosystem;
adjust if you prefer otherwise.
