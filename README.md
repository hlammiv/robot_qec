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
- [x] Map the reference architecture + inventory every field assumption
- [x] Ground qudit capabilities with live experiments
- [x] Scope the qudit extension (math + per-module changes + roadmap)
- [ ] Implement Phase 0 (qudit construction + field-aware k + exact-distance MVP)
- [ ] Implement Phase 1 (GF(q) distance backend: qudit decoder + MILP)
- [ ] Implement Phase 2 (qudit genotype + seeds + prompt context; run a campaign)
- [ ] Discover and verify a new qudit code

## License

Apache-2.0 (see `LICENSE`). Chosen to match the surrounding QEC ecosystem;
adjust if you prefer otherwise.
