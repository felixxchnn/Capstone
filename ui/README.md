# Capstone Research Interface (`ui/`)

A connected, modular presentation layer for the DepMap perturbation-response
capstone. It shows the project's **committed model predictions** and provenance,
and it **performs no model inference**. It is the interactive counterpart of the
fully static offline `../phase2_report.html`.

---

## What this is (and is not)

| It is | It is not |
| --- | --- |
| A read-only view over the frozen `data/processed/case_study.json` | A place where models run |
| Two **independent** model rankings, side by side | A consensus / merged ranking |
| Drug–gene interaction **evidence retrieval**, grouped after ranking | Treatment-efficacy or clinical guidance |
| A live protein-structure explorer (UniProt → RCSB → AlphaFold → Mol\*) | A claim that structure proves function or drug response |
| Deterministic: the scientific data is a byte-verified copy of the committed artifact | A second, hand-maintained copy of the science |

No test-split features, outcomes, predictions, or performance numbers are read or
shown anywhere.

---

## Architecture

```
src/
  app/            Providers (the one data-source seam), routes, shell (App.tsx)
  pages/          Overview · DependencyExplorer · ModelComparison · ProteinStructure · Methods
  components/
    navigation/   top nav + mobile menu + DNA-helix mark
    overview/     hero, helix SVG, result cards
    samples/      sample selector, role badge, sample-info panel
    rankings/     ordered table, mobile cards, search + filters, results summary
    genes/        gene detail drawer, predicted-vs-observed status
    evidence/     evidence group, evidence record, status pill
    comparison/   side-by-side independent-ranking comparison
    methods/      (rendered inline on the Methods page)
    structures/   provider status, candidate list, Mol* viewer (lazy), text summary, WebGL guard
    export/       CSV / JSON / print buttons
    safety/       non-clinical banner
    common/       Card, Callout, Badge, Disclosure, Drawer, SegmentedControl, ErrorBoundary, EmptyState, Spinner
  data/
    CapstoneDataSource.ts          the interface every component consumes
    adapters/
      StaticCaseStudyDataSource.ts  the ONLY implementation today (reads the committed JSON)
      ApiDataSource.contract.md     documented future backend contract (not implemented)
    providers/                      uniprot.ts · rcsb.ts · alphafold.ts · structureProvider.ts · http.ts
    schemas/caseStudy.ts            runtime validation of the committed artifact
    case_study.generated.json       byte-for-byte copy of ../data/processed/case_study.json (build artifact)
    case_study.expected-sha256.txt  pinned hash (mirrors capstone/data-integrity-hashes.md)
  lib/            format.ts (scientific wording) · export.ts · selection.ts · a11y.ts
  hooks/          useAsync (abortable) · useSelection (URL-synced) · useMediaQuery · useCaseStudyViews
  styles/         tokens.css (green/white/black design system) · base.css · components.css
  tests/          Vitest + React Testing Library (network mocked)
scripts/
  sync-case-study.mjs   copies + hashes the committed case study into src/data/
  browser-smoke.mjs     headless-Chrome screenshots of every route (dev aid; not a unit test)
  structure-smoke.mjs   live UniProt/RCSB/AlphaFold integration smoke (records IDs + outcome)
```

### How committed model results reach the UI

1. `data/processed/case_study.json` is the single scientific source of truth. It
   is produced by `case_study.py` from hash-pinned Phase 1 artifacts and the
   reconstructed fitted state; nothing in this app changes it.
2. `npm run sync:data` (auto-run on `predev` / `prebuild` / `pretest`) copies it
   **verbatim** to `src/data/case_study.generated.json` and writes its SHA-256 to
   `src/data/case_study.sha256.txt`.
3. `src/tests/caseStudyLoad.test.ts` asserts the copy is byte-identical to the
   committed file **and** matches the hash pinned in
   `src/data/case_study.expected-sha256.txt`
   (`a962c01a…`, the value recorded in `capstone/data-integrity-hashes.md`). If
   the committed JSON ever changes, this test tells you to review and bump the
   pin deliberately.
4. `StaticCaseStudyDataSource` parses the copy through the schema guard and
   exposes it via `CapstoneDataSource`. Every page/component depends only on that
   interface.

### Why no model inference happens in the frontend

The Phase 1 comparison is frozen and already quantified. The reconstructed
fitted state reproduces the committed validation statistics exactly, and the
`case_study.json` rankings are the deterministic output of running those models
once. Re-running a model in the browser could only reproduce the same numbers at
best, and at worst drift from them. The UI's job is to make the *committed*
result legible, not to recompute it.

### Future Python backend

`src/data/adapters/ApiDataSource.contract.md` specifies the endpoints and
invariants a backend must satisfy. To switch, add
`src/data/adapters/ApiDataSource.ts` implementing `CapstoneDataSource` and pass it
to `<DataSourceProvider value={…}>` in `src/main.tsx`. **No component changes.**
The backend must still serve the same committed predictions + provenance, never
merge the two models, never expose test-split data, and stamp every response
with the `case_study.json` SHA-256 it derived from.

---

## Protein structure providers

The structure explorer is **independent of the data source** and is contacted
only with identifiers (Entrez Gene ID, UniProt accession, PDB id) and human
taxonomy `9606` — never expression values or predictions.

| Step | Service | Endpoint (verified 2026‑08‑31) |
| --- | --- | --- |
| 1 | UniProt ID mapping | `https://rest.uniprot.org/uniprotkb/search?query=xref:GeneID-<entrez> AND organism_id:9606 AND reviewed:true` |
| 2 | RCSB Search API | `https://search.rcsb.org/rcsbsearch/v2/query` (by `reference_sequence_identifiers.database_accession`, `results_content_type: experimental`, sorted by resolution) |
| 2 | RCSB Data API | `https://data.rcsb.org/rest/v1/core/entry/<id>` (method, resolution, release date, primary citation) |
| 2 | RCSB model file | `https://files.rcsb.org/download/<ID>.cif` (loaded by Mol\* via `viewer.loadPdb`) |
| 3 | AlphaFold DB API | `https://alphafold.ebi.ac.uk/api/prediction/<accession>` — the current `cifUrl` / `pdbUrl` and `globalMetricValue` (mean pLDDT) are **read from the response**, never guessed |
| 4 | Mol\* viewer | `molstar` npm package, `Viewer` app from `molstar/lib/apps/viewer/app`, lazy-loaded, WebGL-guarded, behind an error boundary |

Experimental structures are preferred and clearly labelled; the AlphaFold model
is the labelled **predicted** fallback. Candidates are listed (with method +
resolution) and the user chooses — the viewer never silently picks an
unexplained "best".

**Exploded / component view:** true translated exploded assembly is not attempted
(it cannot be done reliably and tested within scope). Instead the viewer exposes
Mol\*'s own maintained **Components** panel for isolating individual chains,
ligands and water. Nothing fakes a separated structure or edits coordinates.

### Network requirements

- **Static scientific data:** none — bundled at build time.
- **Protein structure:** live HTTPS to `rest.uniprot.org`, `search.rcsb.org`,
  `data.rcsb.org`, `files.rcsb.org`, `alphafold.ebi.ac.uk`. Offline / blocked /
  API-error / no-mapping / no-structure / WebGL-unavailable all have explicit UI
  states; the rest of the app is unaffected.
- No analytics, trackers, remote fonts, or CDN scripts anywhere.

---

## Run · test · build

```bash
cd ui
npm install            # uses package-lock.json (exact versions)
npm run dev            # Vite dev server (syncs the case study first)
npm run lint           # ESLint (flat config)
npm run typecheck      # tsc --noEmit
npm run test -- --run  # Vitest + RTL, network mocked
npm run build          # sync + typecheck + production build to dist/
npm run preview        # serve the production build

# optional dev aids (need a local Chrome/Edge; not part of CI)
npm run smoke:structure       # live UniProt/RCSB/AlphaFold integration check
node scripts/browser-smoke.mjs [--mobile]   # route screenshots -> ui/.smoke/
```

Node ≥ 20.19 (`< 23`). This repo's Node is a portable install at
`C:\Users\Leo He\tools\node` (see `capstone` docs); add it to `PATH` first.

### Pinned toolchain (from `package-lock.json`)

| Package | Version | | Package | Version |
| --- | --- | --- | --- | --- |
| react / react-dom | 19.2.0 | | vite | 8.2.2 |
| react-router-dom | 7.9.4 | | @vitejs/plugin-react | 6.1.1 |
| molstar | 5.11.0 | | vitest | 4.1.11 |
| typescript | 5.9.3 | | @vitest/coverage-v8 | 4.1.11 |
| eslint | 9.39.5 | | jsdom | 30.0.1 |
| typescript-eslint | 8.69.0 | | @testing-library/react | 16.3.3 |
| eslint-plugin-react-hooks | 7.1.1 | | @testing-library/user-event | 14.6.6 |
| fp-ts | 2.16.11 (molstar peer) | | @types/node | 22.19.0 |

---

## Scientific limitations carried in the UI

- Both models are **reconstructed** fitted state, not the unavailable historical
  Phase 1 objects (they reproduce the committed statistics exactly).
- `ACH-000364` is one held-out validation cell line — a pipeline check, not a
  performance estimate. Its observed CRISPR values were attached **after**
  ranking.
- `BG003082` is a real primary tumour with **no CRISPR outcome**; bulk tumour
  tissue is a domain shift from the cultured cell lines the models were trained
  and validated on. Every number for it is exploratory.
- Drug–gene evidence is **retrieval** from a licence-filtered 2026‑06b DGIdb
  snapshot (interaction data vintage Dec‑2023). It never establishes efficacy.
- A predicted (AlphaFold) structure is **not** an experimental measurement and
  proves nothing about function, drug response, or therapeutic relevance.
- The frozen Phase 1 headline (ridge_pca 0.2356 vs ridge_head 0.2047, Δ −0.0308)
  is the authoritative result. The osteosarcoma-cohort aggregate (n = 5) is
  descriptive and unstable and is not a substitute.
