# Workflows

Multi-agent LLM orchestrations used to produce the deliverables in `docs/`.

## `qudit-qec-scope.js`

The "ingest → acquire → understand → scope" pipeline for the qudit extension.
It is a [Claude Code `Workflow`](https://claude.com/claude-code) script (plain
JS) that runs many subagents in a deterministic structure:

| Phase | What it does | Fan-out |
|-------|--------------|---------|
| **Map** | Parallel readers over the reference codebase clusters (construction, evaluator cascade, BP-OSD, MILP, equivalence/dedup, CSS evolution, non-CSS evolution, prompts/entry/scripts). Each returns a structured module map **and an inventory of every GF(2)/qubit assumption** (the extension points), with `file:line` and a difficulty rating. | 8 readers |
| **Probe** | Parallel **live experiments** against `qldpc`/`galois` to ground the scoping in fact: qudit construction over GF(q), whether qldpc's GF(q) checks actually commute, the qudit decoder API/gap, field-aware `k` via `galois`, and BLISS/Clifford generalizability. | 4 probers |
| **Scope** | Parallel scoping by concern (algebra+genotype, distance+decoders, evolution loop, equivalence/post-campaign), each producing a file-level change plan with MVP-vs-full effort. | 4 scopers |
| **Verify** | Adversarial verification, distinct lenses: algebraic correctness vs. the **literature** (web search: 2BGA / generalized-bicycle / qudit BB codes), `qldpc` **implementation** reality (re-runs code), and **feasibility/effort** realism. Corrections override earlier claims. | 3 verifiers |
| **Synthesize** | One architect agent merges everything into the structured scoping document. | 1 |

Output (returned to the orchestrator) is written up into `docs/01..04`.

### Re-running / iterating

This is launched from within a Claude Code session via the `Workflow` tool.
To iterate on the script and re-run while reusing cached agent results:

```
Workflow({ scriptPath: "workflows/qudit-qec-scope.js", resumeFromRunId: "<run id>" })
```

Same script + same inputs ⇒ cached agents return instantly; only changed/new
agents re-run.

### Provenance

This file is a verbatim copy of the script that generated the current `docs/`.
The canonical run artifacts live under the session's workflow transcript dir.
