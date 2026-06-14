# Literature

This project is seeded by one paper. We do **not** vendor the paper source into
this (public) repository; we reference it by arXiv identifier and keep a local
copy under `.gitignore`.

## Source paper

**Evolutionary Discovery of Bivariate Bicycle Codes with LLM-Guided Search**
Juan Cruz-Benito, Andrew W. Cross, David Kremer, Ismael Faro (IBM Research).
arXiv: **2606.02418**.

LLM-guided evolutionary search (OpenEvolve / MAP-Elites, FunSearch-style) that
mutates Python *generator ansätze* producing bivariate-bicycle (BB) and
perturbed-bivariate-bicycle (PBB) quantum LDPC codes over `F_2[x,y]/(x^ℓ-1, y^m-1)`,
with a staged validation pipeline (GF(2) rank → BP-OSD → MILP exact distance,
BLISS Tanner-graph dedup, Clifford-equivalence).

## Magic-state-distillation papers (added for the MSD scoping, `docs/07`)

Qudit/qutrit magic-state-distillation cluster (source tarballs gitignored; referenced by arXiv ID):
- **2403.06228** — Low Overhead Qutrit Magic State Distillation (triorthogonal `[[9m−k,k,2]]₃`)
- **2408.00436** — A Search for High-Threshold Qutrit MSD Routines (weight-enumerator ⇒ performance)
- **2510.10852** — Sublogarithmic Distillation in all Prime Dimensions (punctured Reed–Muller)
- **2603.18560** — High-threshold MSD with quantum quadratic-residue codes
- **2605.30108** — Asymptotic magic state distillation with almost-linear rate

## Reference code repository (the repo the paper "talks about")

`qcode-discovery` — https://github.com/qiskit-community/qcode-discovery
(source + data for the paper). We clone it as a **sibling** directory for study;
it is never committed here.

Supporting libraries the pipeline builds on:
- `qldpc` — https://github.com/Infleqtion/qLDPC (field-generic; the key qudit enabler)
- `ldpc` — https://github.com/quantumgizmos/ldpc (binary BP-OSD)
- `OpenEvolve` — https://github.com/codelion/openevolve (evolutionary loop)

## To fetch the paper locally

```bash
# arXiv source tarball (placed in this folder during ingestion):
#   literature/arXiv-2606.02418v1.tar.gz   (gitignored)
mkdir -p literature/extracted
tar -xzf literature/arXiv-2606.02418v1.tar.gz -C literature/extracted
```
