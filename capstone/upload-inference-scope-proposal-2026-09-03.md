# Upload-enabled exploratory inference mode — approved scope

> **APPROVED APPLICATION-LAYER EXTENSION — NOT YET IMPLEMENTED**

Author: Claude Code session, 2026-09-03. Repository HEAD at drafting:
`5fbf342f1536a125d80f5a3b0ed3c8f95dc58ba7` (`main`, level with `origin/main`).

**Approval:** Felix explicitly approved this narrowly defined application-layer scope
extension on **2026-09-03**. The decision is recorded in
[`capstone/scope-decisions.md`](scope-decisions.md) under the CLAUDE.md §13
material-change checklist. This document is the detailed contract that approval refers to.

**Nothing here is implemented.** No upload backend, no upload UI, no `ApiDataSource`, no
`runInference` method, and no `/api/**` endpoint exists. Upload inference is **not
working** and must not be described as working anywhere. As of this document's date the
only change on disk is this file plus the scope-decisions entry; no model code, UI code,
frozen data, result artifact, manifest, or hash was modified.

**Approved objective (verbatim intent):** allow a user to upload **one** supported
research gene-expression profile, validate and align it to the frozen feature space, run
the frozen `ridge_pca` dependency model **without retraining**, and return **traceable
exploratory dependency rankings with evidence and limitations**. The output is **research
decision support** — it is **not** patient treatment-response prediction, treatment
recommendation, clinical validation, or medical advice.

**Runtime figures in this document are unbenchmarked estimates.** Any "seconds" / "< 1 s" /
"~1–2 s" is a design expectation, not a measurement, and must be labelled as an estimate
until a committed benchmark exists (§F.1).

---

## 0. Authoritative current scope (unchanged, restated)

**Phase 1 research question (validated, authoritative):** do frozen Geneformer embeddings
beat a ridge-on-PCA-expression baseline at predicting CRISPR gene dependency across
held-out cancer cell lines, measured by per-target Spearman correlation?

**Phase 1 result (authoritative):** ridge on PCA-200 expression = 0.2356 mean per-target
Spearman; ridge on frozen 768-dim Geneformer embeddings = 0.2047; head − baseline
= −0.0308, paired 95% bootstrap CI [−0.0365, −0.0255]. A rigorous negative finding.

**Phase 2 application layer (implemented, validated):** two named samples (ACH-000364
val-split anchor; BG003082 exploratory external osteosarcoma tumour), two **independent**
frozen models (`ridge_pca`, `ridge_head`), drug–gene interaction **evidence retrieval**
(DGIdb, after rankings freeze), a self-contained offline report, and a modular
React/TypeScript interface with a real protein-structure viewer. It predicts **no**
treatment efficacy, patient response, adverse reactions, or clinical recommendations. The
interfaces **display** committed predictions from `data/processed/case_study.json`; they
run **no** inference. `ApiDataSource.contract.md` describes a future backend seam; no
`ApiDataSource`, upload endpoint, or `runInference` method exists.

**This extension does NOT change** the Phase 1 research question, dataset, prediction
target, evaluation metric, frozen train/val/test split, the validated negative result, or
its scientific interpretation. It changes **the application's accepted input** (one
user-supplied expression profile) and **its runtime behaviour** (from static display to
on-request local inference on that file). Under CLAUDE.md §13 that is a material change to
*input data*, *intended user/use case*, and *runtime behavior* — hence the scope-decisions
entry.

---

## A. Current vs. approved scope (full mapping)

| Dimension | Current (Phase 1 + Phase 2 as shipped) | Approved upload-inference mode |
|---|---|---|
| **Input** | DepMap Public 26Q1 (frozen matrices) + two committed named samples (BG003082 GCT, ACH-000364 DepMap row). No user-supplied data anywhere. | Adds **exactly one** user-supplied research gene-expression profile per request: **GCT v1.2** in `.gct` or `.gct.gz` form, **Ensembl** gene identifiers, **linear TPM** values that the user **explicitly declares** (never inferred). Nothing else in the approved first slice. |
| **Unit of analysis** | Cancer cell line (Phase 1); two demonstration samples (Phase 2). | One uploaded research expression profile ("sample"), handled exactly as BG003082 is: exploratory, external, no ground truth. Not a patient, not a cohort. |
| **Prediction target** | Per-target CRISPR gene effect (Chronos), the **4,297** selective targets in their frozen order. | **Unchanged.** Same 4,297 targets, same frozen order, same meaning. No new target, no re-selection. |
| **Model inputs** | `ridge_pca`: 18,460-dim `log2(TPM+1)` expression vector in canonical order. `ridge_head`: 768-dim Geneformer CLS embedding. | **Unchanged contracts.** The upload is aligned to the *same* 18,460-dim canonical `log2(TPM+1)` vector and run through **frozen `ridge_pca` only**. `ridge_head` is **not** part of the approved first slice (§C.2, §B.3). |
| **Output** | Committed `case_study.json`: top-25 ranked predicted dependencies per model per sample + evidence + osteo aggregate. | A **session-scoped**, non-committed result: top-25 ranked predicted dependencies from `ridge_pca` for the uploaded sample, + full provenance + evidence-retrieved-after-ranking + limitations. Never written to `data/processed/`, never merged into `case_study.json`. |
| **Intended user** | The capstone author; university reviewers cloning the repo; presentation audience. | The same audience **plus** a researcher who wants the frozen model's exploratory output on their **own research** expression profile. **Not** clinicians, **not** patients, **not** for clinical decision-making. |
| **Intended use** | Demonstrate and document a validated negative finding and a possible future workflow. | Same, extended: let a viewer run that possible-future-workflow demonstration on a profile of their choosing. It is **research decision support** and a demonstration of a workflow shape — never a service, never advice. |
| **Evaluation** | Per-target Spearman vs. held-out CRISPR (Phase 1); descriptive n=5 osteo aggregate (Phase 2). No evaluation of BG003082 (no ground truth). | **No evaluation.** An uploaded sample has no CRISPR ground truth; nothing is scored, no metric is computed or shown for it. `outcome_status = unavailable`, exactly as BG003082. |
| **Supported claims** | "Frozen embeddings do not beat the PCA-expression baseline (quantified)." "Here is what the frozen models rank for these two samples." "Here is cited interaction evidence for those genes." | Adds only: *"Here is what the frozen `ridge_pca` model ranks as predicted dependencies for your research expression profile, under a documented domain shift, with no ground truth and no validation for your sample."* |
| **Unsupported claims** (explicitly out) | Treatment efficacy; patient response; clinical recommendation; drug safety / pharmacogenomics; that embeddings ≥ baseline. | **All of the above, unchanged**, plus: that the uploaded ranking is validated, accurate, or generalises; that the sample is comparable to DepMap training data; that a ranked gene is a therapeutic target; that mapping expression into canonical space produces a Geneformer embedding; any clinical-grade interpretation. |
| **Domain-shift limitations** | BG003082 carries a bulk-tumour-vs-cultured-cell-line caveat and `exploratory_external_prediction`. | Every uploaded result carries the **same class of caveat, made explicit**: the model was trained and validated only on cultured DepMap cell lines; an arbitrary research profile (different tissue, protocol, quantifier, filtering) is a domain shift; missing genes are training-mean-imputed; no ground truth exists for the uploaded sample; the ranking is exploratory and uncalibrated. |
| **Privacy boundary** | None — all data committed & public. | The upload may contain human research data. It is processed **locally only** (in memory / an ephemeral temp file), **transmitted to no third party** (the DGIdb snapshot is local/offline; the protein viewer receives only gene identifiers + taxonomy 9606, as today), and a pre-upload warning **prohibits protected health information and identifying patient information**. No expression values are logged. |
| **Persistence policy** | All data committed. | **No persistence by default.** The uploaded bytes and any derived vector live only for the request; the temp file is deleted in a `finally`; nothing is written under `data/processed/`; no database, no cache of user data, no user accounts, no cloud storage. |
| **Current implementation status** | Phase 1 + Phase 2 complete and validated. | **NOT IMPLEMENTED.** Approved for implementation; design only. |
| **Runtime behavior** | Static: both interfaces read a byte-verified copy of a committed JSON. Zero server-side computation. Structure viewer makes live HTTPS calls with identifiers only. | Adds a **local** inference path: a small Python backend validates the upload, aligns it via the **same** `sample_profile` code the committed pipeline uses, runs `fitted_artifacts` **closed-form `predict()`** (no `fit()`, no `fit_transform()`, no hyperparameter selection, no test-split access), ranks, retrieves evidence after the ranking freezes, and returns a session result. The static two-sample demo is unchanged and is the **default / fallback** mode. |

**Deferred capabilities (documented, NOT in the approved first slice):**

| Capability | Status | Condition to revisit |
|---|---|---|
| Precomputed 768-dim Geneformer embedding upload → `ridge_head` inference | **Optional second slice, subordinate** | Only after the `ridge_pca` expression slice is implemented **and** validated (incl. §F gate G1). Documented in §B.3 / §C.2. Not required for the first implementation. |
| CSV / TSV expression matrix (`SYMBOL (ENTREZ)` or unambiguous Entrez) | **Deferred** | A later additive slice; even then it is a *constrained* schema. **Arbitrary CSV / arbitrary datasets are permanently out of scope.** (§B.4) |
| Multi-sample / batch dataset upload | **Future work only** | An additive `:batch` endpoint after the single-sample slice ships. **Not an October capability.** (§B.5) |
| Raw expression → Geneformer embedding generation (local) | **Out of scope** | Needs `geneformer` + GPU; a separate async/GPU extension with its own §13 entry. (§C.2) |

**Permanently out of scope (this extension reopens none of them):** model fitting /
fine-tuning; hyperparameter / alpha selection; changing the target set or its order; E2;
E3 beyond its committed state; **any** test-split evaluation or test-set-driven
development; treatment-response prediction; drug-efficacy scoring; treatment
recommendations; pharmacogenomic safety prediction; patient-specific clinical guidance;
clinical-record ingestion; user accounts / cloud storage; a free-form medical chatbot;
arbitrary dataset support; silently merged model rankings; automatic interpretation beyond
the existing evidence-retrieval layer.

**Plain statement (required):** this extension changes *what the application accepts as
input* and *whether it computes at request time*. It does **not** change the Phase 1
scientific question, the dataset it is answered on, the prediction target, the evaluation
metric, the frozen split, or the validated negative result. The frozen models, their
frozen alphas (`ridge_pca` α=100000.0, `ridge_head` α=3162.0), their frozen preprocessing,
and every committed artifact and hash remain untouched.

---

## B. Supported input contracts

### B.0 Parsers that already exist (inspected)

| Parser | Location | What it does | Reuse verdict |
|---|---|---|---|
| `sample_profile.parse_gct` | [sample_profile.py:155](../sample_profile.py#L155) | Full GCT v1.2 schema validation (`#1.2` tag, `<rows>\t<cols>` dim line, header width, `Name`/`Description` cols, per-row field count, declared-vs-actual row count). Reads gzip **or plain** by suffix. Raises `GCTFormatError` on any violation. | **Reuse verbatim.** Reference GCT reader for both `.gct` and `.gct.gz`. |
| `sample_profile.load_external_sample` | [sample_profile.py:257](../sample_profile.py#L257) | GCT → linear-TPM column → strip Ensembl version → collapse duplicate Ensembl (sum linear TPM) → join to canonical Entrez via `ensembl_map.csv` **only** → restrict to canonical space → detect Entrez collisions (sum) → reindex to `gene_columns.json` order → `log2(TPM+1)` → **missing left as explicit `NaN`** → `(Series, provenance)`. No symbol fallback (deliberate; see its docstring). Rejects non-finite / negative TPM; rejects non-unique `ensembl_id` map. | **Reuse verbatim** for the approved contract (B.1). It is the **exact** path `case_study.py` uses for BG003082 — this is what makes the §F gate G1 (byte-exact parity) achievable and mandatory. |
| `geneformer_sample_input.build_geneformer_frames` | [geneformer_sample_input.py:200](../geneformer_sample_input.py#L200) | Canonical linear-TPM vector → Geneformer input frames (`X`, `var`, `obs`), Entrez→Ensembl re-key, "34th drop" guard, no symbol fallback. **Does not tokenise or embed.** | Not in scope for the approved first slice. Relevant only if a future GPU embedding path is approved. |
| `gene_ids.parse_gene_label` | [gene_ids.py:43](../gene_ids.py#L43) | Parse one `SYMBOL (ENTREZ)` label; returns `None` on non-match. | Reuse **only if** the deferred CSV contract (B.4) is later approved. Not first slice. |
| `gene_ids.map_external_matrix` | [gene_ids.py:228](../gene_ids.py#L228) | Reindex an external symbol/Entrez/`SYMBOL (ENTREZ)` matrix into canonical columns. **`fill_value` defaults to `0.0`** (missing vs measured-zero collide) and its symbol pass uses `lookup.setdefault` (silently keeps first of duplicate symbols). CLAUDE.md §7 records why `sample_profile.py` does **not** call it. | **Do NOT reuse.** Its missing/zero conflation and silent duplicate handling are exactly what the approved contract forbids. |

### B.1 Approved first-slice contract — GCT v1.2 (`.gct` / `.gct.gz`), linear TPM, Ensembl IDs

| Property | Requirement |
|---|---|
| File form | GCT v1.2, either plain `.gct` or gzip `.gct.gz`. Both handled by `parse_gct` (suffix-dispatched). No other container. |
| Orientation | GCT native: rows = genes; columns = `Name`, `Description`, then one column **per sample**. Fixed by the format — not a user choice. |
| Identifier column | `Name` column = Ensembl gene IDs, versioned (`ENSG00000123456.7`) or bare. Non-`ENSG` values are counted as `unexpected_identifiers` and excluded from the join (not fatal on their own; see rejection row for the zero-mapped case). |
| Sample identifier | The uploader **names** the sample column to read. If the GCT has exactly one sample column it is used and echoed back for confirmation. Extra sample columns are ignored with an explicit notice (one sample per request — §B.5). |
| Expression scale | **Linear TPM only, explicitly declared by the uploader.** Never inferred from value ranges (§B.6). `load_external_sample` applies `log2(TPM+1)` to match `data/processed/expression.npz`. |
| Alignment | Reindex to the **frozen 18,460-feature order** from `gene_columns.json` (`feature_order.sha256 = b9ca520c…`), independent of GCT row order. |
| Duplicate genes | Two GCT rows with equal stripped Ensembl ID → **summed as linear TPM** (transcript-additive), reported under `duplicate_external_ids`. Two Ensembl IDs → same canonical Entrez → **summed**, reported under `canonical_id_collisions`. Never merged on shared symbol. |
| Missing / unmapped genes | A canonical gene the sample never resolves (absent, or Ensembl ID not in `ensembl_map.csv`) → **explicit `NaN`**, counted (`canonical_genes_missing`). At inference the reconstructed baseline applies its **stored training-mean impute vector** to those positions (the identical impute step every Phase 1 row uses). The imputed-feature count is reported separately. |
| Measured zero | A canonical gene measured at TPM 0 → `log2(0+1) = 0.0`, counted separately (`canonical_genes_measured_zero`). **Never conflated with missing/unmapped/imputed.** These three are distinct, separately reported counts. |
| Non-finite / negative | Any non-finite or negative value in the chosen sample column → **hard reject** (`ExternalSampleError`), offending rows listed, no guessing. |
| Mapping provenance | Ensembl→Entrez via `data/processed/ensembl_map.csv` **only** (NCBI `gene2ensembl`, single-provenance). **No symbol fallback.** A non-unique `ensembl_id` column → hard reject (ambiguous join). |
| Model | Frozen `ridge_pca` **only** (α 100000.0, read not selected). Closed-form `fitted_artifacts` predict. |
| Targets | The frozen 4,297 selective targets in frozen order (`target_order.sha256 = cc49219e…`). |
| Rejection conditions | Wrong version tag; malformed / mismatched dimension line; header not `Name`,`Description`,…; ragged row; named sample column absent; non-finite / negative value; zero canonical genes mapped; ambiguous `ensembl_map`; file over `MAX_UPLOAD_BYTES`; decompresses over `MAX_DECOMPRESSED_BYTES`; unsupported compression; hostile filename. Each returns a specific structured error. |

### B.2 (reserved — was "CSV contract"; see B.4, deferred)

### B.3 Precomputed 768-dim Geneformer embedding + strict sidecar manifest — **OPTIONAL SECOND SLICE, subordinate**

Documented for completeness. **Not part of the approved first implementation.** May be
built **only after** the B.1 `ridge_pca` slice is implemented and passes every §F gate.

| Property | Requirement |
|---|---|
| Orientation | CSV: one data row per sample, first column = sample ID, then columns `0`…`767` in **exact** order (matches `head_ridge_head/feature_names.json`). |
| Manifest (required alongside) | JSON declaring: `sample_id`; `n_dims == 768`; `column_order == "0..767"`; SHA-256 of the embedding CSV; the Geneformer config that produced it (`model_repo`/`model_subdir` = `Geneformer-V2-104M_CLcancer`, `context_length` 4096, `special_token` true, `pooling` CLS layer −1, `geneformer_revision`); the input provenance (source expression file, Entrez→Ensembl map, pseudo-count basis); `all_finite: true`. Mirrors `case_study._bg003082_embedding_block` ([case_study.py:558](../case_study.py#L558)). |
| Framing | Result labelled `exploratory_user_upload` **and** carries the same "commensurability not proven" disclosure BG003082 carries (bulk vs cultured; NCBI-not-`mygene` map; revision pin not proven equal to Phase 1's). The app never claims the embedding is comparable to the Phase 1 training embeddings. |
| Rejection conditions | Not 768 columns; column order ≠ `0..767`; CSV SHA-256 ≠ manifest; any required manifest field missing; model config ≠ the frozen Geneformer config; non-finite value; oversize. |

### B.4 CSV / TSV expression matrix — **DEFERRED (not approved for October)**

A later additive slice **may** support a *constrained* CSV/TSV: identifiers restricted to
`SYMBOL (ENTREZ)` or unambiguous bare Entrez, orientation declared by the user, scale
declared, duplicate Entrez a hard reject, missing → explicit `NaN`, written to
`sample_profile`'s contract (a new `expression_matrix_profile.py`, **never**
`map_external_matrix`). **Arbitrary CSV and arbitrary datasets are permanently out of
scope** and must be rejected at the schema gate. This slice is future work, not an October
promise.

### B.5 Multi-sample / batch — **FUTURE WORK**

The approved slice is **one sample per request**. A batch path, if ever approved, is an
**additive** endpoint (`/inference/expression:batch`, 2 ≤ N ≤ a config cap) that does not
change the single-sample contract. It is **not** an October capability and must not be
presented as one.

### B.6 Scale is declared, never inferred

The uploader **must** select the expression scale from an explicit control. In the
approved first slice the only accepted value is **`linear TPM`**; anything else is
rejected. The backend records the declared scale in the result and in the input-hash
context. It does **not** sniff value ranges to guess. Rationale: the training matrix is
`log2(TPM+1)`; feeding linear TPM as if it were log space (or vice-versa) silently
produces confident nonsense, and range heuristics fail exactly on edge cases (low-depth
samples, heavily filtered gene sets).

### B.7 Hard exclusions (all slices)

Raw FASTQ / BAM / CRAM; VCF / gVCF; WGS / WES / panel variant calls; methylation / ATAC /
single-cell matrices; imaging; clinical records, treatment histories, outcomes, labs;
arbitrary tabular data; bare-symbol-only identifier columns (ambiguous). Anything
requiring alignment, quantification, or germline/somatic calling. Rejected at the schema
gate with a specific message.

---

## C. Model feasibility (the two models treated independently)

### C.1 `ridge_pca` (expression baseline) — **approved, feasible locally**

- **Path exists and is proven.** `fitted_artifacts.load_baseline_ridge_pca()`
  ([fitted_artifacts.py:260](../fitted_artifacts.py#L260)) returns a
  `ReconstructedBaselineRidgePCA` whose `.predict(X)`
  ([fitted_artifacts.py:231](../fitted_artifacts.py#L231)) is **closed-form array
  arithmetic**: `impute(non-finite → stored train mean) → (X − scaler_mean)/scaler_scale
  → (Xs − pca_mean) @ pca_components.T → Z @ ridge_coef.T + ridge_intercept`. **No
  scikit-learn import, no `fit()`, no `fit_transform()`, no hyperparameter selection.**
  Every `.npy` is SHA-256-verified against `manifest.json` on load.
- **Input it needs:** a length-18,460 vector in `gene_columns.json` order, `NaN` allowed
  for missing (imputed by the stored mean). Exactly what
  `sample_profile.load_external_sample()` produces. `.assert_feature_order(names)` /
  `.assert_target_order(names)` enforce alignment before inference; `transform()` also
  hard-checks `X.shape[1] == 18460`.
- **Output:** `(1, 4297)` predicted GeneEffect. Ranking reuses `case_study._full_ranking`
  semantics ([case_study.py:267](../case_study.py#L267)): sort ascending by **raw float64**
  predicted value (most-negative = strongest predicted dependency), numeric Entrez breaks
  exact ties only, round for display **after** top-25 is frozen.
- **Runtime (ESTIMATE, not measured):** a `(1×18460)` scale + `(18460×200)` matmul +
  `(200×4297)` matmul — order a few million FLOPs, expected sub-second on CPU. Artifact
  load (hashing ~37 MiB of `.npy`, incl. a 29 MiB `pca_components.npy`) is a one-time
  cost per process, expected ~1–2 s; cache the loaded model. **These figures are
  unbenchmarked (§F.1).**
- **`fit()` / selection:** none. α (100000.0) is read from the manifest, which read it
  from `baseline_results.json`. Verified this session: `reconstruct_fitted.py --validate`
  reproduces baseline `spearman_mean` `0.2356` exactly; `fitted_artifacts` vs the
  gitignored prediction matrix max-abs-diff `4.44e-16`.
- **Verdict:** an aligned upload can return an exploratory top-25 `ridge_pca` ranking
  locally, prediction-only. **This is the approved first slice.**

### C.2 `ridge_head` (Geneformer embedding) — **not in the approved first slice**

- **What it consumes:** a **768-dim Geneformer CLS embedding**, not an expression vector.
  `ReconstructedHeadRidge.predict(X)` ([fitted_artifacts.py:242](../fitted_artifacts.py#L242))
  is closed-form (`impute → scale → Xs @ ridge_coef.T + intercept`), `X.shape[1] == 768`
  enforced — that half is as cheap and safe as `ridge_pca`.
- **The gap:** producing the embedding needs the `geneformer` package **and a GPU**.
  `geneformer_sample_input.py` builds tokeniser input frames but **does not tokenise or
  embed**; `anndata` is absent locally; the BG003082 embedding was made on **Kaggle**
  (`capstone/kaggle_bg003082_embedding.py`, Tesla T4). **No practical local embedding
  runtime exists.**
- **Aligning expression into the canonical 18,460-gene space does NOT create a Geneformer
  embedding.** Different representations; the app must never imply otherwise.
- **Approved handling for the first slice:** run `ridge_pca` only; the UI shows, for the
  second model, `ridge_head unavailable — no compatible embedding supplied` (a real
  state, not a disabled button implying a hidden capability). **Never fabricate a second
  ranking; never merge the two models.**
- **Later, subordinate:** a precomputed-embedding upload (B.3) → `ridge_head` closed-form,
  with the BG003082-style commensurability caveat. Only after the `ridge_pca` slice ships
  and passes §F. Automatic GPU embedding generation is a separate §13-gated extension.

### C.3 Feasibility summary

| Capability | In approved first slice? | Mechanism | Runtime (ESTIMATE) | Remaining work |
|---|---|---|---|---|
| `ridge_pca` on an uploaded GCT expression profile | **Yes** | `sample_profile.load_external_sample` → `fitted_artifacts.load_baseline_ridge_pca().predict` → `_full_ranking` semantics → evidence after ranking | sub-second inference; ~1–2 s one-time model load *(unbenchmarked)* | backend (validation, alignment reuse, ranking reuse, serializer, `/inference/expression`), `ApiDataSource` + `UploadInferenceClient`, UI "Upload Research Profile" page, the §F test matrix + gate G1 |
| `ridge_head` on a validated precomputed embedding | No — optional second slice | `fitted_artifacts.load_head_ridge().predict` on `(n,768)` | sub-second *(unbenchmarked)* | embedding+manifest validator, `/inference/embedding`, commensurability-caveat plumbing, tests — **only after the first slice passes §F** |
| `ridge_head` from raw uploaded expression (local) | No — out of scope | needs `geneformer` + GPU | n/a | separate async/GPU extension, own §13 entry |
| Single-sample upload | **Yes** (`ridge_pca`) | as row 1 | *(unbenchmarked)* | as row 1 |
| Batch / multi-sample upload | No — future work | additive `:batch` endpoint | n/a | per-sample schema, progress, abort, config cap |

---

## D. Approved October MVP boundary

1. **Keep the current static two-sample demonstration unchanged and available** as the
   **default / fallback** mode. `phase2_report.html` is untouched; `StaticCaseStudyDataSource`
   stays the default adapter; ACH-000364 and BG003082 render exactly as today.
2. **Add a separate "Upload Research Profile" workflow** to the connected React app only
   (not the offline report, which stays a single self-contained file).
3. **One sample per request.**
4. **Input:** GCT v1.2 (`.gct` / `.gct.gz`), Ensembl identifiers, **linear TPM explicitly
   declared** by the user, never inferred.
5. **Alignment** to the frozen 18,460-feature baseline order via the **same**
   `sample_profile.load_external_sample` code the committed pipeline uses (not a
   re-implementation).
6. **Explicit three-way distinction** among measured zero, unmapped/missing, and
   training-mean-imputed values — separately counted and displayed.
7. **Frozen `ridge_pca` inference only** — `fitted_artifacts` closed-form `predict`, no
   `fit()`, no `fit_transform()`, no alpha selection, no target change, **no test-split
   access**.
8. **The frozen 4,297 CRISPR targets in frozen order.**
9. **Session-scoped, local, no persistence by default.** Ephemeral temp file deleted in a
   `finally`; nothing written under `data/processed/`; no accounts, no cloud storage, no
   caching of user data; no expression values in logs.
10. **Pre-upload warning** prohibiting protected health information and identifying patient
    information, alongside the existing non-clinical banner.
11. **`prediction_status = exploratory_user_upload`**; **`outcome_status = unavailable`** —
    always.
12. **DGIdb evidence retrieved only after** the dependency ranking is produced, via
    `evidence.get_evidence_for_gene` against the committed offline snapshot — identical to
    `case_study.collect_evidence` ([case_study.py:390](../case_study.py#L390)). Retrieval
    only, never efficacy.
13. **Reuse the existing protein-structure viewer** for a selected ranked gene (Entrez ID
    + taxonomy 9606 only; unchanged).
14. **Display, per result:** mapping coverage (mapped / 18,460), missing count, imputed
    count, measured-zero count, model ID + frozen alpha ("read, not re-selected"),
    model-artifact hashes, `feature_order_sha256`, `target_order_sha256`, code commit,
    input file SHA-256, declared schema + scale, domain-shift warnings, fixed non-clinical
    disclaimer.
15. **`case_study.json` and every frozen artifact are read-only** to the backend;
    `checks.py` stays green; no protected hash changes (§F row 21).

**The October MVP does NOT promise:** arbitrary CSV / arbitrary dataset compatibility;
multi-sample or batch inference; raw-expression-to-Geneformer inference; patient
treatment-response prediction; drug-efficacy prediction; treatment recommendations;
clinical-grade interpretation; upload persistence or user accounts; retraining or model
selection; any test-set-driven development.

**Why this boundary:** it delivers a genuine end-to-end "bring your own profile" slice
using only code paths already written and validated (`sample_profile`, `fitted_artifacts`,
`evidence`), adds one small backend and one UI workflow, and introduces **zero** new
scientific claims. Everything riskier (embedding upload, CSV, batch, local Geneformer,
data-tuned thresholds) is an explicit fast-follow or deferred, not a launch blocker.

---

## E. Backend and frontend architecture (design — not implemented)

### E.1 Modular backend

A small local Python service (FastAPI or stdlib `http.server` + a thin router; avoid a
new heavyweight dependency if practical), structured so upload logic never leaks into
React:

```
server/
  app.py                 # routing, localhost CORS, request-id middleware, error mapper
  config_upload.py       # ALL upload limits & thresholds — one place (§F.1)
  validation/
    gct.py               # wraps sample_profile.parse_gct + schema gates + declared-scale gate
    embedding.py         # B.3 contract (mirrors case_study._bg003082_embedding_block) — 2nd slice
    filenames.py         # path-traversal / extension / compression allow-list
    limits.py            # size + streaming-decompression caps in front of parse_gct
  alignment/
    to_canonical.py      # calls sample_profile.load_external_sample UNMODIFIED; surfaces the provenance/coverage report
  inference/
    run.py               # fitted_artifacts load (cached, hash-verified) + closed-form predict; NO fit()
    ranking.py           # reuse case_study._full_ranking semantics verbatim
  evidence/
    lookup.py            # reuse evidence.get_evidence_for_gene against the committed snapshot
  serialize/
    result.py            # the inference-response schema (E.4); deterministic field order
  tests/                 # the §F matrix, incl. gate G1
```

- **`fitted_artifacts` is imported; `reconstruct_fitted` is not** — the same rule
  `case_study.py` enforces by AST check; add the equivalent guard here.
- **`sample_profile.load_external_sample` is called unmodified.** The backend does not
  re-implement parsing, mapping, collapsing, log-transform, or NaN handling. This is what
  makes gate G1 (§F) meaningful.
- Model artifacts loaded **once** per process, reused, hash-verified on load.
- The service opens frozen files **read-only**; it never opens `data/processed/` for
  write. A smoke test asserts frozen hashes are unchanged after a request.

### E.2 New `ApiDataSource`, static adapter retained

- Add `ui/src/data/adapters/ApiDataSource.ts implements CapstoneDataSource` per
  `ApiDataSource.contract.md`: it serves the **same committed** summary / samples /
  rankings / evidence / metadata for the two existing samples (contract rules 1–5),
  unchanged.
- `StaticCaseStudyDataSource` stays the **default and fallback**. If the backend is
  unreachable, the app still fully works for the committed demo.
- **Provider seam:** `DataSourceProvider` in
  [ui/src/app/Providers.tsx](../ui/src/app/Providers.tsx) already accepts an injected
  `value`. Add a `static` (default) vs `api` mode there. No component changes for existing
  pages.
- **Upload is a new capability, not a `CapstoneDataSource` method.** Define a separate
  `UploadInferenceClient` interface (`validateUpload`, `runExpressionInference`; later
  `runEmbeddingInference`) consumed only by the new page. This keeps `CapstoneDataSource`
  meaning "read the frozen case study" and preserves the "no `runInference` on the data
  source" rule the contract states.

### E.3 Endpoints (versioned; proposed, not implemented)

| Method + path | Slice | Purpose |
|---|---|---|
| `POST /api/v1/uploads/validate` | first | Schema + declared-scale + size + identifier checks only; returns a mapping/coverage preview + a normalized upload token. **No inference.** |
| `POST /api/v1/inference/expression` | first | Validate → align (`sample_profile`) → `ridge_pca.predict` → rank → evidence-after-ranking → session result. The second model is reported `unavailable`. |
| `POST /api/v1/inference/embedding` | second (subordinate) | Validate 768-dim embedding + sidecar manifest → `ridge_head.predict` → rank → evidence → result, with commensurability caveat. Independent; never synthesizes the other model. |
| `GET /api/v1/summary`, `/samples`, `/rankings/{sample}/{model}`, `/evidence/{entrez}`, `/samples/{sample}/metadata`, `/genes/{entrez}/structure-hint`, `/case-study` | first | The read-only routes from `ApiDataSource.contract.md`, serving the committed two-sample data unchanged; every response stamped with the `case_study.json` SHA-256. |

- **Config-driven limits** in `config_upload.py`: `MAX_UPLOAD_BYTES`,
  `MAX_DECOMPRESSED_BYTES`, `MAX_SAMPLE_COLUMNS`, `MAX_GENE_ROWS`, `REQUEST_TIMEOUT_S`,
  `MIN_MAPPING_COVERAGE`, concurrency cap — none guessed (§F.1).
- **Structured errors:** `{error_code, message, field, limit?}` with a stable
  `error_code` enum; the UI maps codes to friendly copy.
- **Abortable requests:** the existing `useAsync` hook is abortable; wire an
  `AbortController` through the upload client; the backend honours `REQUEST_TIMEOUT_S`
  and returns a structured timeout.
- **Progress / failure states:** upload → validating → aligning → predicting → retrieving
  evidence → done / failed(code). One labelled-stage spinner; explicit empty/failed
  panels (match the report's `#no-results` pattern).

### E.4 Inference-response schema (design target — every field required)

```
request_id                 # uuid, echoed in logs (no PII)
input_sha256               # SHA-256 of the exact uploaded bytes
declared_schema            # "gct-v1.2/ensembl"  (first slice)
declared_scale             # "linear_tpm"        (first slice; only accepted value)
sample_id                  # from the file
mapping_report             # canonical_genes, mapped, missing, measured_zero,
                           #   measured_nonzero, duplicate_external_ids,
                           #   canonical_id_collisions, resolved_outside_canonical
missing_count / imputed_count / measured_zero_count   # the three-way distinction, separate
model_id                   # "ridge_pca"   (first slice)
model_provenance_status    # fitted_artifacts.RECONSTRUCTION_STATUS string, verbatim
model_artifact_hashes      # manifest sha + each .npy sha (from the loaded artifact)
feature_order_sha256       # from the loaded artifact  (b9ca520c…)
target_order_sha256        # from the loaded artifact  (cc49219e…)
frozen_alpha               # 100000.0, annotated "read from baseline_results.json, not re-selected"
code_commit                # git HEAD at serve time
ranked_outputs             # top-25: rank, entrez_id, symbol, predicted_geneeffect(raw + display)
                           #   ranking_rule string (raw float64, ascending, numeric Entrez tie-break)
evidence                   # by_entrez, retrieved AFTER ranking; retrieval-only framing;
                           #   coverage counts; DGIdb snapshot + manifest sha
evidence_retrieval_status  # "ok" | "snapshot_unavailable"
prediction_status          # "exploratory_user_upload"
outcome_status             # "unavailable"
domain_shift_warnings      # list (trained/validated on cultured cell lines only; your sample
                           #   is a domain shift; N genes training-mean-imputed; no ground
                           #   truth; ranking is exploratory and uncalibrated)
non_clinical_disclaimer    # config.DGIDB_EVIDENCE_DISCLAIMER + the fixed report disclaimer
```

---

## F. Safety and integrity requirements (test matrix — required before release)

"Reject" = structured error, no inference, no file retained.

| # | Scenario | Expected behaviour |
|---|---|---|
| 1 | Malformed GCT dimensions (declared rows/cols ≠ actual; header width wrong) | Reject via `GCTFormatError`; message names the mismatch. |
| 2 | Malformed schema (header not `Name`,`Description`,…; ragged row) | Reject; message names the offending line. |
| 3 | Unsupported / undeclared expression scale (anything but declared `linear TPM`) | Reject; state the only accepted value; no range-sniffing fallback. |
| 4 | Invalid identifiers (non-`ENSG` in `Name`) | Count + exclude non-joinable rows; if **zero** canonical genes map → reject. |
| 5 | Duplicate identifiers — same stripped Ensembl ID twice | Sum linear TPM; report `duplicate_external_ids`; proceed. |
| 6 | Ambiguous mapping (`ensembl_map.csv` non-unique `ensembl_id`; two Ensembl → one Entrez) | Non-unique map → hard reject. Two→one Entrez in the upload → sum, report `canonical_id_collisions`. |
| 7 | Negative expression values in the chosen column | Reject (`ExternalSampleError`), row indices listed. |
| 8 | `NaN` / `Infinity` in a **provided** cell | Reject; "refusing to guess a replacement". |
| 9 | Measured zero vs unmapped/missing vs imputed | `0.0` provided → `log2(1)=0.0`, counted `measured_zero`. Gene absent / Ensembl unmapped → `NaN` → training-mean-imputed at inference; counted `missing` and `imputed`. **All three counts are distinct and never combined.** |
| 10 | Low mapping coverage (below `MIN_MAPPING_COVERAGE`) | Reject **or** proceed-with-prominent-warning per the config threshold; the threshold is **benchmarked, not invented** (§F.1). |
| 11 | Incorrect feature order into the model | `fitted_artifacts.assert_feature_order` raises `ReconstructedArtifactError`; request fails closed. |
| 12 | (2nd slice) Incorrect embedding dimensions (≠ 768; columns not `0..767`) | Reject before `predict`. |
| 13 | Hash mismatch (a model `.npy` ≠ its manifest; 2nd slice: embedding CSV ≠ sidecar sha) | Fail closed; never predict on unverified arrays. |
| 14 | Oversized input (bytes > `MAX_UPLOAD_BYTES`) | Reject at the transport boundary before reading into memory. |
| 15 | Excessive decompressed size (gzip bomb > `MAX_DECOMPRESSED_BYTES`) | Streaming gunzip limiter aborts + rejects. **Note:** `parse_gct` currently `read()`s the whole file — the cap must sit in front of it (`validation/limits.py`). |
| 16 | Path-traversal / hostile filename (`../`, absolute, NUL, control chars, wrong extension) | Reject; never use the client filename for any filesystem path — assign a server UUID. |
| 17 | Unsupported compression (`.zip`, `.bz2`, `.xz`, double-gzip) | Reject; allow-list is `{plain .gct, .gct.gz}` only. |
| 18 | Multiple simultaneous requests | Each request isolated (own temp; a guarded shared read-only model or a per-request reference); no shared mutable state; a concurrency cap returns a structured "busy" error. |
| 19 | Cancellation / timeout | UI `AbortController` aborts the fetch; backend enforces `REQUEST_TIMEOUT_S`, returns a structured timeout; partial temp files removed. |
| 20 | No upload persistence | After the response (success or failure) the temp file is gone; assert the upload dir is empty; no expression values in logs. |
| 21 | No mutation of frozen artifacts | SHA-256 of every frozen input (`expression.npz`, `case_study.json`, both reconstruction manifests + their `.npy`, DGIdb snapshot/manifest, `phase2_report.html`, `baseline_results.json`) identical before/after a full request cycle; `py checks.py` still 55/55. |
| 22 | No test-split access | The backend never reads test-split rows for features or outcomes; it does not load `splits.json` at all in the first slice. A grep/AST test asserts no `"test"` split selection anywhere in `server/`. |
| 23 | **Gate G1 — BG003082 round-trip parity (MANDATORY, blocks release)** | Feeding the **committed** `data/external/sid_osteosarc/BG003082.gene_tpm.gct.gz` through `POST /api/v1/inference/expression` (declared `linear TPM`, sample `BG003082`) reproduces `case_study.json`'s `rankings["BG003082"]["ridge_pca"]` **exactly**: all 25 rows identical in `rank`, `entrez_id`, `symbol`, and `predicted_geneeffect` at the committed display precision, and the same `ranking_rule` semantics. **Any difference fails validation and blocks release.** |
| 24 | Gate G2 — same-code alignment | The canonical 18,460-vector the backend feeds to `ridge_pca` for a given GCT is byte-identical to `sample_profile.load_external_sample()` on the same file. The backend must **call** that function, not re-implement it. Proving matched shape / feature order is **not** sufficient — exact preprocessing compatibility must be proven. |
| 25 | Determinism | Same upload → identical `ranked_outputs` and `mapping_report` across repeated requests and process restarts. |
| 26 | Evidence-after-ranking | The evidence lookup receives the already-frozen gene list and cannot change ranks; the result asserts `evidence_availability_did_not_affect_ranking = true`. |
| 27 | Model independence | The response never contains a merged/consensus list; when the second model is `unavailable` there is no cross-fill or fabricated ranking. |

### F.1 Thresholds must be benchmarked, not guessed

The repository currently enforces **no** minimum mapping coverage and **no** upload-size
limit (BG003082 maps 18,427 / 18,460 ≈ 99.8% and is ~778 KB compressed). New knobs and
how to set them:

| Knob | Home | How to set it (evidence-based) |
|---|---|---|
| `MIN_MAPPING_COVERAGE` | `server/config_upload.py` | Benchmark on the 170 **val** expression rows (allowed — not test): randomly mask increasing fractions of genes to `NaN`, run `ridge_pca.predict`, measure the per-target Spearman drop vs. the unmasked val result. Set the threshold where degradation becomes material (e.g. mean Spearman falls to the `lineage_mean` control, 0.15); default to **warn-not-reject** above it, reject below. Commit the curve as a short `capstone/` note. |
| `MAX_UPLOAD_BYTES` | `server/config_upload.py` | From the largest legitimate single-sample GCT (~3 MB uncompressed for ~75k rows); ~10× headroom (e.g. 32 MB); documented. |
| `MAX_DECOMPRESSED_BYTES` | `server/config_upload.py` | ~5–10× `MAX_UPLOAD_BYTES`; enforced by a streaming gunzip limiter before `parse_gct`. |
| `MAX_SAMPLE_COLUMNS` / `MAX_GENE_ROWS` | `server/config_upload.py` | Rows: a little above the human transcriptome (~80k). Sample columns: 1 in the first slice. |
| `REQUEST_TIMEOUT_S` | `server/config_upload.py` | Measure p99 of a full request on the target host; set ~5×. Parsing/hashing a large GCT dominates; inference is expected sub-second. **Both are estimates until measured.** |
| Concurrency cap | `server/config_upload.py` | 1–2 for a laptop demo; a structured "busy" response above it. |

Any runtime number quoted to a reviewer or in the UI must be sourced from a **committed
benchmark**, not from this document's estimates.

---

## G. Implementation gates (all must pass before the upload mode is enabled by default)

1. **G1** — BG003082 round-trip parity (§F row 23): exact reproduction of the committed
   `ridge_pca` top-25 for BG003082 through the upload path. A blocker.
2. **G2** — same-code alignment (§F row 24): the backend calls
   `sample_profile.load_external_sample` unmodified; exact preprocessing compatibility is
   proven, not inferred from dimensionality.
3. **G3** — `py checks.py` → 55/55, and every protected-artifact hash unchanged after a
   full request cycle (§F row 21).
4. **G4** — no test-split code path anywhere in `server/` (§F row 22), asserted by test.
5. **G5** — the static two-sample demo still renders and still passes its existing tests
   with the backend both present and absent.
6. **G6** — every §F row has a passing automated test.

Until G1–G6 pass, the upload workflow ships (if at all) behind an explicit
non-default toggle labelled experimental.

---

## H. Approval record

| Item | Decision (2026-09-03, Felix) |
|---|---|
| Scope entry in `capstone/scope-decisions.md` with the §13 checklist | **Approved** — application input + runtime change recorded; Phase 1 science explicitly unchanged. |
| First slice = GCT / Ensembl / declared-linear-TPM → frozen `ridge_pca` only, one sample per request, session-scoped, no persistence | **Approved.** |
| Precomputed 768-dim embedding upload → `ridge_head` | **Approved as an optional, subordinate second slice only**, after the first slice is implemented and validated. Not required for the first implementation. |
| CSV / arbitrary dataset support | **Not approved** — arbitrary CSV / datasets are permanently out; a constrained CSV schema is possible future work only. |
| Batch / multi-sample upload | **Deferred to future work** — not an October capability. |
| Raw-expression → Geneformer inference (local) | **Out of scope** — separate GPU extension, own §13 entry. |
| Clinical / treatment-response / efficacy / recommendation claims | **None introduced.** Output is research decision support only. |
| Implementation in this session | **Not authorized** — scope recording, E3 decision, and the model freeze only. |

This document remains **APPROVED — NOT YET IMPLEMENTED** until an implementation session is
separately authorized and gates G1–G6 pass.
