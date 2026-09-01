# UI expansion — implementation plan (2026-08-31)

Scope: presentation/UI layer only. No model fitting, no inference recompute, no
test-split access. `case_study.json` and `random_projection_results.json` are
read-only; their hashes must not change. `phase2_report.html` hash changes
intentionally (generator redesign).

Environment note: Node.js v22.23.2 (LTS) + npm 10.9.8 were installed as a
portable extract at `C:\Users\Leo He\tools\node` (no Node was present).

## Stage 1 — `feat: add modular capstone research interface`  (`ui/**`)

Toolchain: Vite 7 + React 19 + TypeScript 5 + Vitest 3 + @testing-library/react
+ ESLint 9 (flat config) + jsdom. Pinned via `ui/package-lock.json`.

Data flow (single source of truth):
- `data/processed/case_study.json` is the ONLY scientific source.
- A Node sync script (`ui/scripts/sync-case-study.mjs`) copies it verbatim to
  `ui/src/data/case_study.generated.json` and writes its SHA-256 to
  `ui/src/data/case_study.sha256.txt`. Runs on `predev`/`prebuild`/`pretest`.
- `ui/src/data/schemas/caseStudy.ts` — zod-free hand-written type guard + parse.
- A test asserts the synced copy's SHA-256 == the value committed in
  `capstone/data-integrity-hashes.md` (kept in `ui/src/data/case_study.sha256.txt`).

Architecture:
- `src/data/CapstoneDataSource.ts` — interface: `getProjectSummary`, `getSamples`,
  `getModelRanking(sample, model)`, `getGeneEvidence(entrez)`, `getSampleMetadata(sample)`,
  `getStructureMetadata(entrez)`.
- `src/data/adapters/StaticCaseStudyDataSource.ts` — current impl (reads the
  generated JSON).
- `src/data/adapters/ApiDataSource.contract.md` — documented future contract; NOT
  implemented (no fake endpoints).
- `src/data/providers/` — `uniprot.ts`, `rcsb.ts`, `alphafold.ts` (live structure
  metadata; identifiers only; AbortSignal support).
- `src/lib/format.ts` — scientific formatting (GeneEffect sign, ranks, deltas).
- `src/lib/export.ts` — CSV + JSON export of the displayed view.
- `src/lib/a11y.ts` — focus-trap, reduced-motion, id helpers.
- `src/app/` — providers, routing, state (URL-synced selection).
- `src/pages/` — Overview, DependencyExplorer, ModelComparison, ProteinStructure,
  Methods.
- `src/components/<area>/` — navigation, overview, samples, models, rankings,
  genes, evidence, structures, safety, export, common.
- `src/styles/tokens.css` — green/white/black design system.

Nav: Overview · Dependency Explorer · Model Comparison · Protein Structure ·
Methods & Limitations. Desktop bar + mobile disclosure menu. Selection
(sample/model/gene/search/evidence filter) kept in the URL query string.

Functions 1–14 map to components; rank order 1–25 is never re-sorted (search/
filter only toggle visibility). Tests cover schema load, sample/model select,
rank-order invariance, search, evidence filter, reset, observed vs unavailable,
independent side-by-side rankings, gene detail open, CSV/JSON export.

## Stage 2 — `feat: add protein structure explorer`  (`ui/**`)

`molstar` pinned npm dep, lazy-loaded, error boundary, WebGL guard.
Provider chain: Entrez → UniProt (reviewed, taxonomy 9606) → RCSB PDB candidates
→ AlphaFold fallback. User picks an experimental candidate or falls back to
AlphaFold (clearly labelled "predicted structure"). Controls: rotate/zoom/pan/
reset/fullscreen, representation (cartoon/surface/ball-and-stick), colour by
chain, ligands, water, background toggle, screenshot.
Exploded view: implement the **component-isolation panel** (focus/hide chains &
ligands) — true translated exploded mode is not reliably doable in scope; report
honestly. All error/empty/offline/no-mapping/no-structure/WebGL states explicit.
Tests mock all network. One documented live smoke test with a real ranked gene.

## Stage 3 — `style: redesign offline Phase 2 report`  (`report.py`, `phase2_report.html`, `config.py`, docs)

Rewrite `_CSS` / `_JS` / `render_html` in `report.py` with the same design
system: hero + plain-language intro, sticky section nav, result cards, sample/
model segmented controls, search + evidence filter + reset, ordered dependency
table, gene/evidence expansion, side-by-side model comparison section, sample
information, methods/limitations, print/PDF, responsive, visible warnings, empty
states. A short note that interactive structures need the `ui/` app (no viewer,
no fetch, no CDN).
Embedded `case_study.json` stays byte-equal to the committed file.
Update `validate()` structural checks + `_SMOKE_HARNESS` deliberately for any
markup/id changes; keep all scientific-claim and offline-guarantee checks.
Rebuild twice (byte-identical) → new SHA-256 → update ONLY
`config.REPORT_HTML_SHA256`, plus the hash references + history in `CLAUDE.md`
and `capstone/data-integrity-hashes.md` (previous hash kept as provenance).

## Cross-cutting

- `capstone/scope-decisions.md`: dated 2026-08-31 entry.
- `README.md`: run instructions for both interfaces.
- `.gitignore`: `ui/node_modules/`, `ui/dist/`, `ui/coverage/`, Vite caches,
  temp screenshots.
- Python validations after Stage 3: case_study 43/43, evidence 34/34,
  E1 self-test 10/10, E1 artifact 24/24, checks.py 55/55, report.py --validate
  all structural + browser smoke.
