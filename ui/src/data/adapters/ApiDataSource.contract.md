# `ApiDataSource` — future backend contract (NOT implemented)

This file documents the contract a future Python-backed `ApiDataSource` will
satisfy so it can replace `StaticCaseStudyDataSource` with **no change to any UI
component**. No endpoints are implemented today; the UI ships with the static
adapter only. There are no fake API results and no non-functional "predict"
buttons anywhere in the app.

## Rules the backend must keep

1. **No new inference in response to the UI.** The frozen Phase 1 models and the
   reconstructed fitted state are the source of every prediction. A backend
   serves the *same committed predictions and provenance*; it does not fit,
   fine-tune, or re-run a model because a user clicked something.
2. **No test-split anything.** The backend must never expose test-split
   expression features, test outcomes, test predictions, or test performance
   numbers. Split-assignment labels may be read only for integrity / role
   assertions.
3. **`ridge_pca` and `ridge_head` stay independent.** No endpoint returns a
   merged or consensus ranking.
4. **Evidence is retrieval, after ranking.** Drug–gene evidence is returned from
   the licence-filtered DGIdb snapshot and is never framed as efficacy.
5. **Determinism / provenance.** Every response carries the `case_study.json`
   SHA-256 (or an equivalent version id) it was derived from.

## Method → transport mapping

| `CapstoneDataSource` method | Method + path (proposed) | Notes |
| --- | --- | --- |
| `getProjectSummary()` | `GET /api/v1/summary` | frozen Phase 1 headline, osteo aggregate, coverage, env, hashes |
| `getSamples()` | `GET /api/v1/samples` | two samples, explicit roles |
| `getModelRanking(sample, model)` | `GET /api/v1/rankings/{sample}/{model}` | frozen top-25; `model ∈ {ridge_pca, ridge_head}` only |
| `getGeneEvidence(entrez)` | `GET /api/v1/evidence/{entrez}` | one bucket or `404`; retrieval-only |
| `getSampleMetadata(sample)` | `GET /api/v1/samples/{sample}/metadata` | reconciliation, imputation, domain-shift |
| `getStructureMetadata(entrez)` | `GET /api/v1/genes/{entrez}/structure-hint` | Entrez + `taxonomy=9606`; **no** expression/prediction data |
| `getRawCaseStudy()` | `GET /api/v1/case-study` | whole immutable artifact, read-only |

## Where the seam is in the code

`src/app/Providers.tsx` constructs exactly one `CapstoneDataSource` and puts it in
context via `DataSourceProvider`. To switch to a live backend, add
`src/data/adapters/ApiDataSource.ts` implementing `CapstoneDataSource`, then swap
the single `new StaticCaseStudyDataSource()` call. Components use
`useDataSource()` and never import an adapter directly.

## Structure services are separate

`src/data/providers/{uniprot,rcsb,alphafold}.ts` talk to public structural
databases directly from the browser and are **independent of the data source**.
They:

- are contacted only with gene / protein identifiers (Entrez, UniProt accession,
  PDB id) and human taxonomy `9606`;
- never receive sample expression values or model predictions;
- do not influence rankings, metrics, or any scientific conclusion;
- return *experimental or predicted structural evidence*, never drug-response
  evidence.
