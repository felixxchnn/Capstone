# Scope decisions

Dated record of material scope changes, per CLAUDE.md §13. Append only; do not edit or
remove a past entry if a later decision supersedes it — add a new entry instead.

---

## 2026-08-25 — Phase 2: precision-oncology proof-of-concept demo

**Approved by:** Felix, in review of the draft implementation plan (two `AskUserQuestion`
prompts, both answered "Approve as specified").

**Current approved scope before this change:** the capstone answers one question — do frozen
Geneformer embeddings beat a ridge-on-PCA-expression baseline at predicting CRISPR gene
dependency across held-out DepMap cancer cell lines, measured by per-target Spearman
correlation. That comparison (Phase 1) is complete for its pre-specified headline (ridge_pca
vs ridge_head) and is not reopened, discarded, or rewritten by anything below. E1 (random-
projection control), E3's confirmatory wider-grid run, the model-set freeze, and the one-time
test-split evaluation remain open, independent of this decision.

**Change:** the October deliverable also includes a proof-of-concept application layer —
using Phase 1's frozen models to rank predicted CRISPR dependencies for two samples, then
connecting the top-ranked genes to cited drug–gene interaction evidence (not a trained
model), shown in a self-contained offline HTML report.

**Why:** Felix's stated goal is a substantial, honest, end-to-end prototype the team can
continue developing after the presentation, demonstrating how genomic dependency prediction
could eventually support individualized cancer-treatment research — not a claim that it does
so today.

**What changes, precisely (material-change checklist, CLAUDE.md §13):**

- *Population/unit of analysis:* Phase 1 is 1,140 DepMap cell lines. Phase 2 adds two
  individually-named samples used for demonstration, not for any statistic Phase 1 reports:
  - **ACH-000364** (DepMap ID, cell line "U-2 OS", val split) — internal verification anchor.
    Selected for osteosarcoma identity and having every required artifact already committed
    (expression, CRISPR ground truth, Geneformer embedding, all four saved val-split
    predictions) — explicitly **not** selected for a favorable predicted score. Shown
    alongside the aggregate result across all 5 val-split osteosarcoma lines so this one line
    cannot be mistaken for cherry-picked validation.
  - **BG003082** — a real RNA-seq sample from Sid Sijbrandij's self-released osteosarcoma
    dataset (`osteosarc.com`, CC0 1.0), primary tumor, resected 2022-12-16. This is **bulk
    tumor tissue**, not a cultured cell line — a real domain shift from everything the model
    was trained and validated on. No CRISPR screen exists for this tissue.
- *Input data:* adds one committed external RNA-seq profile (BG003082's gene-level TPM,
  778 KB, CC0) and one committed drug–gene interaction evidence snapshot (DGIdb, license
  verified and filtered before committing, per-source). Also required: a static Entrez↔
  Ensembl cross-reference to fill a previously-undiscovered gap — `ensembl_map.csv`, which
  `prepare_geneformer_input.py` already expects to cache but which has never existed anywhere
  in this repository (the original 18,460/18,460 zero-attrition Ensembl mapping was produced
  and consumed entirely on Kaggle, via `mygene`, and never carried back here).
- *Prediction target/label:* unchanged (per-target CRISPR gene effect). No new prediction
  target is introduced; the evidence layer retrieves existing, cited gene–drug interaction
  records, it does not predict anything.
- *Evaluation metric:* unchanged for Phase 1. The two Phase 2 samples are explicitly **not**
  scored against Phase 1's metric — BG003082 has no ground truth to score against at all, and
  ACH-000364's role is demonstration, not a new evaluation.
- *Intended user/use case:* unchanged core research question; adds an illustrative
  application-layer demonstration of a possible future workflow, labeled as such throughout.
- *Scientific claim:* no new claim about model performance. Every Phase 2 output carries an
  explicit, machine-readable status:

  | Field | ACH-000364 | BG003082 |
  |---|---|---|
  | `prediction_status` | `held_out_prediction` | `exploratory_external_prediction` |
  | `outcome_status` | `measured_crispr` | `unavailable` |

  BG003082's prediction is never described as validated, never compared to a measurement
  (none exists), and always carries the cell-line-vs-bulk-tumor domain-shift caveat.
- *Distinction between experimental, patient-level, and clinical evidence:* BG003082 is a
  real, named, identifiable person's real tumor data, self-released specifically for reuse —
  it is accurate to call it a real patient sample, and the report does so, but never implies
  a treatment recommendation, an efficacy claim, or clinical validation. The drug–gene
  evidence section is titled "Drug–gene interaction evidence," never "candidate treatments"
  or "recommended treatments," and every retained record carries its source, interaction
  direction, evidence tier, licence, and a fixed efficacy disclaimer. This does not reopen
  the CPIC/PharmGKB pharmacogenomic *safety* layer (germline adverse-reaction prediction),
  which remains cut — the evidence-retrieval layer approved here is a narrower, different
  thing: cited interaction evidence, no model, no safety claims.

**Alternatives considered and rejected:**
- DepMap val-split osteosarcoma line alone (zero new engineering risk, no real-patient
  narrative value) — rejected in favor of pairing it with BG003082 for higher narrative
  value, accepting the added engineering risk.
- Installing an interactive web framework (e.g. Streamlit) for the UI — rejected in favor of
  a dependency-free static HTML report with embedded vanilla JavaScript, to avoid a new
  failure surface six weeks before the presentation.

**Record:** approved 2026-08-25. E2 (concatenation control) is deprioritized past October
given the added Phase 2 workload, or reduced to its minimum-honest version if time remains —
this is a scheduling call, not a scope change, and does not require a separate entry here.

See `C:\Users\Leo He\.claude\plans\moonlit-dazzling-dream.md` for the full implementation
plan (gene-ID reconciliation order, log-scale handling, no-leakage guarantees, the evidence
schema, and the work sequence).

---

## 2026-08-29 — Phase 2: reconstructed fitted artifacts (no case_study.py yet)

**Approved by:** Felix, this session, "Choose Option 2. First generate and commit
deterministic reconstructed fitted artifacts. Do not implement case_study.py or report.py
in this pass."

**Context / why.** `case_study.py` needs to run the two frozen Phase 1 linear models
(`ridge_pca`, `ridge_head`) on the two demo samples. The original Phase 1 fitted
`StandardScaler` / `PCA` / `Ridge` objects were **never serialised** and cannot be
recovered; `baseline_results.json` / `head_results.json` hold only hyper-parameters and
metrics. A prior pass established that inference therefore could not proceed without either
refitting or first committing fitted state. Option 2 (commit reconstructed fitted state,
then have `case_study.py` load it with no `fit()`) was chosen over refitting inside
`case_study.py`.

**What was built (this decision, implemented same day):**
`reconstruct_fitted.py` fits `impute(train-mean) -> StandardScaler -> PCA(200) -> Ridge`
(baseline) and `impute -> StandardScaler -> Ridge` (head) on **exactly the committed
`train` split**, at the alpha **read from** `baseline_results.json` / `head_results.json`
(100000.0 / 3162.0 — no hyper-parameter selection re-run), and serialises plain `.npy`
arrays + a `manifest.json` under `data/processed/reconstructed_fitted/`.
`fitted_artifacts.py` loads them and predicts with array arithmetic only — no sklearn
import, no `fit()`. Verified to reproduce every committed Phase 1 validation statistic
exactly at the recorded 4-dp precision; byte-identical on rebuild.

**These are described everywhere, without exception, as:**

> "reconstructed fitted state at the frozen Phase 1 alpha from the unchanged frozen
> training data"

They are **not** the unavailable original fitted Phase 1 objects. This is a
reproducibility convenience, not a new result: no Phase 1 metric, split, target ordering,
or committed results file changes.

### Osteosarcoma five-line validation aggregate — definition locked (NOT computed here)

`case_study.py` (a later pass) will report one descriptive aggregate. Its full definition,
fixed now so it cannot be improvised later:

- **cohort:** the validation-split cell lines selected by the committed osteosarcoma
  predicate `config.osteosarcoma_mask(model_metadata)` (keyword match on
  `OncotreePrimaryDisease` / `OncotreeSubtype`) intersected with `splits.json` assignment
  `== "val"`. That set has **five** cell lines (`checks.py` §8: val osteosarcoma = 5).
  The five DepMap IDs are whatever that predicate + split selects — not hand-listed.
- **split:** `val` only. The `test` split is never touched.
- **models:** `ridge_pca` and `ridge_head` **only** (the two reconstructed linear models).
  Not `lineage_mean`, not `global_mean`, not `mlp_head`.
- **targets:** all **4,297** selective targets (`selective_genes.json` filtered to
  `crispr_effect` columns), same ordering the reconstruction uses.
- **per-target statistic:** for each target, the Spearman correlation between predicted and
  observed CRISPR GeneEffect **across the five cohort cell lines** (via
  `baseline.per_target_spearman`, its own `min_samples`/constant-column rules unchanged).
- **common target set:** the intersection of targets whose per-target Spearman is **finite
  for both models**. Report the included count and the excluded count (and why excluded:
  too few finite pairs, or a constant column, at n = 5).
- **aggregate:** the plain arithmetic mean of the per-target Spearman values over that
  identical common target set, computed **separately for each model**.
- **delta:** `mean(ridge_head) - mean(ridge_pca)` over the common set.
- **direction:** more negative predicted GeneEffect = stronger predicted dependency
  (`config.DEPENDENCY_THRESHOLD = -0.5`, "dependent below this score"). The aggregate is a
  correlation, so direction only matters for the per-sample ranking, not for this number.
- **status:** explicitly **descriptive and unstable because n = 5**. No confidence
  interval, no significance test. It does **not** replace or restate the frozen Phase 1
  primary result (val ρ 0.2356 vs 0.2047 over 170 held-out lines); it is shown only so
  ACH-000364's single-line result does not read as cherry-picked.

**Not computed in this pass.** Only the definition is recorded here.
