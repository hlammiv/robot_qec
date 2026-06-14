export const meta = {
  name: 'qudit-qec-scope',
  description: 'Deep-understand IBM qcode-discovery and scope its extension to qudit (GF(q)) code discovery',
  phases: [
    { title: 'Map' },
    { title: 'Probe' },
    { title: 'Scope' },
    { title: 'Verify' },
    { title: 'Synthesize' },
  ],
}

// ---------- shared context ----------
const REPO = '/home/hlamm/Desktop/QC/qcode-discovery'
const PAPER = '/home/hlamm/Desktop/QC/robot_qec/literature/extracted/paper.tex'

const FACTS = [
  'PROJECT GOAL: extend IBM qcode-discovery (qubit, F_2) into a pipeline that discovers QUDIT error-correcting codes over GF(q).',
  'Reference repo (read-only) at ' + REPO + '. Source paper at ' + PAPER + ' (arXiv 2606.02418).',
  'The paper: OpenEvolve/MAP-Elites evolves Python generator ansatze G(l,m) producing polynomial pairs (A,B) (CSS BB codes) or 4-tuples (A,B,C,D) (non-CSS PBB) over F_2[x,y]/(x^l-1,y^m-1). Cascade: Stage1 k via GF(2) rank; Stage2 BP-OSD distance; Stage3 MILP exact (HiGHS/scipy). Post-campaign: BLISS Tanner-graph dedup (python-igraph), decomposability, Clifford/LC equivalence. FOM=k*d^2/n.',
  'Installed and importable in base python (/home/hlamm/miniforge3/bin/python3): qldpc 0.2.1, ldpc 2.4.0, galois 0.4.6, sympy, scipy, igraph.',
  'ALREADY-GROUNDED FACTS (do not waste effort re-deriving; build on them):',
  ' (1) qldpc is FIELD-GENERIC. qldpc.codes.BBCode(orders, poly_a, poly_b, field=q) constructs a qudit BB code over GF(q) out of the box. Confirmed: BBCode({x:6,y:6}, x**3+y+y**2, y**3+x+x**2, field=3) -> num_qudits=72, dimension=8.',
  ' (2) qldpc.codes.CSSCode(code_x, code_z, field=...) and qldpc.codes.QuditCode(matrix, field=...) (non-CSS, from a symplectic stabilizer matrix) both accept field. code.dimension and code.num_qudits are field-aware.',
  ' (3) DECODER GAP: code.get_distance_bound_with_decoder(Pauli.X, T, bp_method="product_sum", osd_method=..., osd_order=...) WORKS on GF(2) (returns 12 for the gross code) but on GF(q>2) qldpc switches to a class GUFDecoder that raises TypeError on bp_method/osd kwargs. The ldpc library is BINARY ONLY. qcode-discovery evaluation/distance.py hardcodes those BP-OSD kwargs and will break over GF(q). qldpc.get_distance_exact() is field-generic but costly over GF(q).',
  'IMPORTANT for any agent running python: use OS-level timeout (e.g. shell: timeout 60 python3 script.py). Do NOT use in-process signal.alarm around qldpc C calls - it crashes the interpreter (observed exit 144).',
].join('\n')

// ---------- schemas ----------
const MODULE_SCHEMA = {
  type: 'object',
  required: ['cluster', 'summary', 'modules', 'field_assumptions'],
  properties: {
    cluster: { type: 'string' },
    summary: { type: 'string', description: '1-2 paragraphs: role of this cluster in the pipeline and how data flows through it' },
    modules: { type: 'array', items: { type: 'object', required: ['file', 'purpose', 'key_symbols', 'external_deps'], properties: {
      file: { type: 'string' }, purpose: { type: 'string' },
      key_symbols: { type: 'array', items: { type: 'string' }, description: 'key functions/classes with one-line role' },
      external_deps: { type: 'array', items: { type: 'string' }, description: 'qldpc/ldpc/scipy/igraph/sympy/galois calls relied upon' },
    } } },
    field_assumptions: { type: 'array', description: 'EVERY place that assumes GF(2)/qubits (the extension points)', items: { type: 'object', required: ['location', 'assumption', 'qudit_change', 'difficulty'], properties: {
      location: { type: 'string', description: 'file:line or file:function' },
      assumption: { type: 'string' },
      qudit_change: { type: 'string', description: 'concretely what must change for GF(q)' },
      difficulty: { type: 'string', enum: ['trivial', 'easy', 'moderate', 'hard', 'research'] },
    } } },
  },
}

const PROBE_SCHEMA = {
  type: 'object',
  required: ['question', 'commands_run', 'findings', 'verdict'],
  properties: {
    question: { type: 'string' },
    commands_run: { type: 'array', items: { type: 'string' } },
    findings: { type: 'array', items: { type: 'string' }, description: 'concrete observed outputs / signatures / error messages' },
    verdict: { type: 'string', description: 'grounded conclusion with evidence' },
  },
}

const SCOPE_SCHEMA = {
  type: 'object',
  required: ['concern', 'summary', 'changes', 'new_components', 'effort', 'risks', 'open_questions'],
  properties: {
    concern: { type: 'string' },
    summary: { type: 'string' },
    changes: { type: 'array', items: { type: 'object', required: ['target', 'change', 'difficulty'], properties: {
      target: { type: 'string', description: 'file/module to modify' },
      change: { type: 'string' },
      difficulty: { type: 'string', enum: ['trivial', 'easy', 'moderate', 'hard', 'research'] },
    } } },
    new_components: { type: 'array', items: { type: 'string' }, description: 'new files/modules to create' },
    effort: { type: 'string', description: 'rough effort + MVP vs full split' },
    risks: { type: 'array', items: { type: 'string' } },
    open_questions: { type: 'array', items: { type: 'string' } },
  },
}

const VERIFY_SCHEMA = {
  type: 'object',
  required: ['claim', 'verdict', 'confidence', 'evidence', 'corrections'],
  properties: {
    claim: { type: 'string' },
    verdict: { type: 'string', enum: ['confirmed', 'partially_correct', 'refuted', 'uncertain'] },
    confidence: { type: 'number' },
    evidence: { type: 'array', items: { type: 'string' } },
    corrections: { type: 'array', items: { type: 'string' }, description: 'specific fixes if not fully confirmed' },
  },
}

// ===================== PHASE 1 (Map) + PHASE 2 (Probe) — independent, run concurrently =====================
phase('Map')

const READERS = [
  { cluster: 'construction', files: ['evaluation/bb_code.py', 'evaluation/pbb_code.py', 'evaluation/mirror_code.py', 'evaluation/__init__.py'] },
  { cluster: 'evaluator-cascade', files: ['evaluation/evaluator.py', 'evaluation/results.py', 'evaluation/tracking.py'] },
  { cluster: 'distance-bposd', files: ['evaluation/distance.py', 'evaluation/distance_bposd_noncss.py'] },
  { cluster: 'distance-milp', files: ['evaluation/distance_milp.py'] },
  { cluster: 'equivalence-dedup', files: ['evaluation/tanner_equivalence.py', 'evaluation/clifford_equivalence.py'] },
  { cluster: 'evolve-css', files: ['evolve/seed_solution.py', 'evolve/seed_solution_ansatz.py', 'evolve/openevolve_evaluator.py', 'evolve/run_evolution.py'] },
  { cluster: 'evolve-noncss-config', files: ['evolve/seed_solution_noncss.py', 'evolve/openevolve_evaluator_noncss.py', 'evolve/_noncss_distance_worker.py', 'evolve/config.yaml', 'evolve/config_noncss.yaml', 'evolve/config_ansatz.yaml'] },
  { cluster: 'prompts-entry-scripts', files: ['evolve/prompt_context.md', 'evolve/prompt_context_ansatz.md', 'evolve/prompt_context_noncss.md', 'main.py', 'scripts/verify_publication.py', 'scripts/verify_from_scratch.py'] },
]

const readerThunks = READERS.map((r) => () => agent(
  'You are mapping one cluster of the IBM qcode-discovery codebase to plan a qudit (GF(q)) extension.\n\n' + FACTS +
  '\n\nYOUR CLUSTER: "' + r.cluster + '". Read these files (under ' + REPO + '): ' + r.files.join(', ') +
  '\n\nFor each file: capture purpose, key functions/classes, and which qldpc/ldpc/scipy/igraph/galois/sympy facilities it relies on. ' +
  'THEN, most importantly, hunt for EVERY place that assumes binary/GF(2)/qubits: e.g. mod-2 arithmetic, np.int8/uint8 over GF(2), GF(2) rank/nullspace, ldpc BP-OSD with bp_method/osd kwargs, Pauli.X/Z only (no qudit clock/shift powers), coefficient-free genotype (exponent tuples with implicit coeff 1), n=2*l*m, k=2lm-2rank, symplectic-over-F2, BLISS coloring that ignores edge coefficients, hardcoded {I,S,H} qubit Cliffords, etc. Give file:line for each. For each, state concretely what must change for GF(q) and rate difficulty. Use Read and Grep. Be exhaustive and precise; cite line numbers.',
  { label: 'map:' + r.cluster, phase: 'Map', schema: MODULE_SCHEMA }
))

const PROBES = [
  { id: 'construct', q: 'Exactly how does qldpc build qudit BB/CSS/non-CSS codes over GF(q), and are the stabilizer checks genuinely commuting?',
    body: 'Run python (via: timeout 90 python3 file.py). (a) Build qldpc.codes.BBCode over GF(q) for q in {2,3,4,5} at small lattices (e.g. (6,6),(3,6)); print num_qudits, dimension, field.order. (b) CRITICAL: verify the stabilizer checks COMMUTE over GF(q) - obtain the symplectic stabilizer matrix (try code.matrix, code.matrix_x/matrix_z, or code.get_stabilizer / code.canonicalized) and check the symplectic inner product of all stabilizer pairs is 0 mod q. This decides whether qldpc auto-handles the CSS commutativity (antipode/sign) for odd q. (c) Inspect what H_X and H_Z look like vs A,B for q=3: does qldpc use H_Z=(B*,-A*) (antipode + sign)? Compare matrix_z to the transpose/conjugate of A,B blocks. (d) Try QuditCode(matrix, field=q) from a hand-built symplectic matrix to confirm non-CSS qudit construction. Report concrete numbers and any API used.' },
  { id: 'distance', q: 'What distance facilities does qldpc expose over GF(q), what is the exact qudit decoder API, and how costly is exact distance?',
    body: 'Run python with OS timeouts. (a) For a small GF(3) BB code, call get_distance_exact() under timeout 60 and report whether it returns and how long. (b) Determine the GF(q) decoder API: inspect qldpc source for GUFDecoder and for get_distance_bound / get_distance_bound_with_decoder; find the call signature that DOES work over GF(q) (the binary one with bp_method fails). Actually invoke a working qudit distance-bound call and report the value+signature. (c) Test get_distance() and get_distance_bound() on the GF(3) code. (d) Note how cost scales with q (alphabet grows q-fold). Locate qldpc install path via python -c "import qldpc,os;print(os.path.dirname(qldpc.__file__))" and grep its source for the decoder classes. Report exact signatures and evidence.' },
  { id: 'linalg', q: 'How is field-aware k computed, and can galois do the GF(q) linear algebra the pipeline needs?',
    body: 'Run python. (a) Show that code.dimension equals n - rank(H_X) - rank(H_Z) over GF(q): extract H_X,H_Z (or the CSS matrices) from a GF(3) BBCode and compute galois ranks; compare to code.dimension. (b) Demonstrate galois GF(q) operations the pipeline will need: GF=galois.GF(3); matrix rank, null_space/left_null_space, row_reduce. (c) Briefly assess how a GF(q) minimum-weight-logical MILP differs from GF(2): mod-q linear constraints need integer-multiple slack variables; symplectic weight = (X-support OR Z-support) per qudit. State the formulation delta concretely (you do NOT need to implement it). Report code snippets that worked and their output.' },
  { id: 'dedup-clifford', q: 'Can the BLISS Tanner-graph dedup and the Clifford/LC-equivalence machinery generalize to GF(q)?',
    body: 'Read ' + REPO + '/evaluation/tanner_equivalence.py and skim ' + REPO + '/evaluation/clifford_equivalence.py. (a) Describe exactly how the colored Tanner graph is built with python-igraph (vertex colors, edges) and whether edge COEFFICIENTS (GF(q) values, not just 0/1) are representable - igraph BLISS supports vertex colors; can coefficients be encoded (e.g. as extra colored vertices per (check,qubit,coeff) or edge subdivision)? Propose the concrete encoding. (b) For Clifford/LC equivalence: identify everything that is qubit-Pauli-specific ({I,S,H,HS,...}, binary symplectic) and state what the qudit-Clifford generalization requires (GF(q) symplectic group, qudit phase/multiplication gates). (c) Check qldpc for qudit logical-op / symplectic support (code.get_logical_ops over GF(q), qldpc.objects). Report findings with file:line where relevant.' },
]

const proberThunks = PROBES.map((p) => () => agent(
  'You are grounding a qudit-QEC scoping effort with LIVE experiments and source reading.\n\n' + FACTS +
  '\n\nYOUR QUESTION: ' + p.q + '\n\nDO THIS: ' + p.body +
  '\n\nWrite throwaway python to /tmp and run via "timeout N python3 /tmp/...py". Report exactly what you observed (numbers, signatures, errors). Prefer evidence over speculation; if something fails, report the precise error.',
  { label: 'probe:' + p.id, phase: 'Probe', schema: PROBE_SCHEMA }
))

const [maps, probes] = await Promise.all([
  parallel(readerThunks),
  parallel(proberThunks),
])

const mapsOk = maps.filter(Boolean)
const probesOk = probes.filter(Boolean)
log('Map: ' + mapsOk.length + '/' + READERS.length + ' clusters; Probe: ' + probesOk.length + '/' + PROBES.length + ' experiments')

// Consolidated context for scoping
const MAP_CTX = mapsOk.map((m) => '### Cluster: ' + m.cluster + '\n' + m.summary + '\nField assumptions:\n' +
  (m.field_assumptions || []).map((a) => '- [' + a.difficulty + '] ' + a.location + ': ' + a.assumption + ' -> ' + a.qudit_change).join('\n')).join('\n\n')
const PROBE_CTX = probesOk.map((p) => '### Probe: ' + p.question + '\nVERDICT: ' + p.verdict + '\nFindings: ' + (p.findings || []).join(' | ')).join('\n\n')
const GROUND = '=== ARCHITECTURE MAP (field-assumption inventory) ===\n' + MAP_CTX + '\n\n=== GROUNDED EXPERIMENTS ===\n' + PROBE_CTX

// ===================== PHASE 3 (Scope) — by concern, concurrent =====================
phase('Scope')

const CONCERNS = [
  { id: 'algebra-genotype', focus: 'Algebra, code construction, and the evolutionary GENOTYPE. Cover: threading field=q through bb_code.py/pbb_code.py construction; the 2-block group-algebra (2BGA) framing of qudit BB codes and whether qldpc already handles CSS commutativity (antipode/sign) for odd q; the genotype change from coefficient-free exponent tuples (x_exp,y_exp) to coefficient-carrying terms (x_exp,y_exp,coeff in GF(q)*); terms_to_poly / validate_terms updates; whether to restrict to prime q (field) vs prime-power GF(q) vs composite Z_d (ring, no field - Smith/Howell normal form); seed programs.' },
  { id: 'distance-decoders', focus: 'Distance estimation and decoders over GF(q). Cover: evaluation/distance.py (binary BP-OSD kwargs break on GF(q) -> must branch to qldpc GUFDecoder / qudit path); exact-distance Tier-1 enumeration blow-up over GF(q); distance_milp.py reformulation for GF(q) (mod-q constraints, symplectic weight); distance_bposd_noncss.py achievable-syndrome sampling generalization; the trust filter. Define an MVP distance backend (likely exact-enum for small n + qldpc qudit decoder-bound + MILP-GF(q)).' },
  { id: 'evolution-loop', focus: 'The OpenEvolve evolutionary loop. Cover: openevolve_evaluator.py and openevolve_evaluator_noncss.py changes (field param, fitness, MAP-Elites behavioral features, trust filter); run_evolution.py CLI (add --field / --qudit-dim); the config_*.yaml; prompt_context*.md rewrites teaching the LLM qudit BB algebra (clock/shift Paulis, GF(q) coefficients, antipode, qudit-specific code families); ensemble/model config. What new domain knowledge must the prompt convey so mutations explore the qudit search space well?' },
  { id: 'equivalence-postcampaign', focus: 'Dedup and post-campaign verification over GF(q). Cover: tanner_equivalence.py BLISS coloring to encode GF(q) edge coefficients; decomposability analysis; clifford_equivalence.py for qudit Cliffords (GF(q) symplectic, qudit phase gates) - likely the hardest, possibly deferred for MVP; novelty assessment vs known qudit code catalogs. Mark what is MVP vs research.' },
]

const scopes = await parallel(CONCERNS.map((c) => () => agent(
  'You are scoping ONE concern of extending IBM qcode-discovery (qubit, F_2) to QUDIT (GF(q)) code discovery. Produce an actionable, file-level change plan.\n\n' +
  FACTS + '\n\nYou may Read any file under ' + REPO + ' and run grounding python (timeout N python3 ...). Build on this consolidated team context:\n\n' + GROUND +
  '\n\nYOUR CONCERN: ' + c.focus +
  '\n\nDeliver concrete, file-level changes (each with a difficulty), the new components to create, an effort estimate split into MVP vs full, risks, and open questions. Be specific and realistic; prefer reusing qldpc field-generic facilities over reimplementing. Distinguish "prime q (field, MVP)" from "composite-d qudits (ring, research)".',
  { label: 'scope:' + c.id, phase: 'Scope', schema: SCOPE_SCHEMA }
)))

const scopesOk = scopes.filter(Boolean)
const SCOPE_CTX = scopesOk.map((s) => '### Concern: ' + s.concern + '\n' + s.summary + '\nChanges:\n' +
  (s.changes || []).map((ch) => '- [' + ch.difficulty + '] ' + ch.target + ': ' + ch.change).join('\n') +
  '\nNew: ' + (s.new_components || []).join('; ') + '\nEffort: ' + s.effort + '\nRisks: ' + (s.risks || []).join('; ') + '\nOpen: ' + (s.open_questions || []).join('; ')).join('\n\n')

// ===================== PHASE 4 (Verify) — adversarial, concurrent =====================
phase('Verify')

const VERIFY = [
  { id: 'math-literature', lens: 'algebraic correctness against the published literature', task:
    'CLAIM TO TEST: "Qudit BB codes are the GF(q) case of two-block group-algebra (2BGA) / generalized-bicycle (GB) codes over GF(q)[Z_l x Z_m]; with H_X=(A|B), the CSS condition needs H_Z=(B*,-A*) where * is the antipode (x->x^-1,y->y^-1) transpose and the minus sign matters for odd q (over F_2 it reduces to the paper s H_Z=(B^T,A^T)); for ABELIAN groups commutativity is then automatic; the genotype must carry GF(q) coefficients; prime-power q gives a field (GF(q)) while composite d gives a ring Z_d where rank/distance need Smith/Howell normal form." Use WebSearch/WebFetch to check against the literature (search: 2-block group algebra codes Lin Pryadko; generalized bicycle codes Panteleev Kalachev; bivariate bicycle qudit / GF(q); qudit CSS commutativity antipode). Confirm or correct the signs, the field-vs-ring distinction, which q are valid, and the coefficient claim. Give citations.' },
  { id: 'implementation', lens: 'what qldpc actually does (run code)', task:
    'CLAIM TO TEST: "qldpc already constructs COMMUTING qudit BB/CSS/non-CSS codes over GF(q) (so the hard commutativity algebra is solved inside qldpc), gives field-aware k, exposes a working but different qudit decoder (GUFDecoder) and a field-generic get_distance_exact; galois provides the GF(q) linear algebra." Re-verify by running python (timeout N python3): build a GF(3) BBCode, prove its stabilizers commute over GF(3), confirm dimension==n-rank_X-rank_Z via galois, confirm the qudit decoder-bound call signature works and the binary one fails, and confirm get_distance_exact returns on a tiny GF(3) code. Report any claim that does NOT hold.' },
  { id: 'feasibility', lens: 'feasibility, critical path, and effort realism', task:
    'Critically assess the consolidated scoping below for an MVP that can actually evolve and verify NEW qudit codes. Identify the true critical path, the single biggest risk, any effort estimate that looks optimistic, and anything mis-scoped as MVP that is really research (e.g. qudit Clifford/LC equivalence, composite-d distance, GF(q) MILP performance, non-CSS qudit distance). State the minimal believable MVP. Here is the scoping:\n\n' + SCOPE_CTX },
]

const verifies = await parallel(VERIFY.map((v) => () => agent(
  'You are an adversarial verifier. Lens: ' + v.lens + '. Default to skepticism; only mark "confirmed" with concrete evidence.\n\n' + FACTS +
  '\n\nTASK: ' + v.task +
  '\n\nReturn a verdict (confirmed/partially_correct/refuted/uncertain), confidence 0-1, evidence, and SPECIFIC corrections where the claim is wrong or incomplete.',
  { label: 'verify:' + v.id, phase: 'Verify', schema: VERIFY_SCHEMA }
)))

const verifiesOk = verifies.filter(Boolean)
const VERIFY_CTX = verifiesOk.map((v) => '### ' + v.claim.slice(0, 80) + '\nVerdict: ' + v.verdict + ' (conf ' + v.confidence + ')\nEvidence: ' + (v.evidence || []).join(' | ') + '\nCorrections: ' + (v.corrections || []).join(' | ')).join('\n\n')

// ===================== PHASE 5 (Synthesize) =====================
phase('Synthesize')

const SYNTH_SCHEMA = {
  type: 'object',
  required: ['executive_summary', 'math_framing_md', 'architecture_map_md', 'change_inventory_md', 'decoder_strategy_md', 'phased_roadmap_md', 'mvp_definition', 'risks_md', 'open_questions'],
  properties: {
    executive_summary: { type: 'string', description: '4-8 sentence summary of the whole scoping conclusion' },
    math_framing_md: { type: 'string', description: 'Markdown: the qudit BB math (2BGA/GF(q)), corrected per the verifiers' },
    architecture_map_md: { type: 'string', description: 'Markdown: how qcode-discovery works, cluster by cluster' },
    change_inventory_md: { type: 'string', description: 'Markdown table/list: per-module changes with difficulty, MVP-flag' },
    decoder_strategy_md: { type: 'string', description: 'Markdown: the distance/decoder strategy over GF(q)' },
    phased_roadmap_md: { type: 'string', description: 'Markdown: phased implementation roadmap (Phase 0..N) with concrete milestones' },
    mvp_definition: { type: 'string', description: 'The minimal believable MVP that evolves+verifies a new qudit code' },
    risks_md: { type: 'string', description: 'Markdown: top risks and mitigations' },
    open_questions: { type: 'array', items: { type: 'string' } },
  },
}

const synthesis = await agent(
  'You are the lead architect. Synthesize a single coherent scoping document for extending IBM qcode-discovery (qubit, F_2) to QUDIT (GF(q)) error-correcting code discovery. Be technically precise, file-specific, and honest about what is MVP vs research. Incorporate the verifiers corrections (they override earlier claims on conflict).\n\n' +
  FACTS + '\n\n=== ARCHITECTURE + ASSUMPTION MAP ===\n' + MAP_CTX + '\n\n=== GROUNDED EXPERIMENTS ===\n' + PROBE_CTX +
  '\n\n=== PER-CONCERN SCOPING ===\n' + SCOPE_CTX + '\n\n=== ADVERSARIAL VERIFICATION (authoritative) ===\n' + VERIFY_CTX +
  '\n\nProduce the structured document. The phased roadmap must let us reach an MVP that actually evolves and MILP/exact-verifies at least one NEW qudit code, then expands toward the full pipeline.',
  { label: 'synthesize', phase: 'Synthesize', schema: SYNTH_SCHEMA }
)

return {
  counts: { maps: mapsOk.length, probes: probesOk.length, scopes: scopesOk.length, verifies: verifiesOk.length },
  maps: mapsOk,
  probes: probesOk,
  scopes: scopesOk,
  verifies: verifiesOk,
  synthesis,
}
