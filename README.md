# DepMap Perturbation-Response Engine — Data Spine

**Research question (Phase 1, validated):** do frozen Geneformer embeddings beat a
ridge-on-PCA-expression baseline at predicting CRISPR gene dependency across held-out cancer
cell lines, measured by per-target Spearman correlation?

**Answer: no, and the deficit is quantified.** A ridge on 200 PCA components of expression
scores 0.2356 mean per-target Spearman (4,297 held-out targets, 170 val-split cell lines); a
ridge head on frozen 768-dim Geneformer embeddings scores 0.2047 — a paired deficit of
−0.0308 (95% CI [−0.0365, −0.0255], 1,000-resample cell-line bootstrap). This is a clean,
quantified negative result, not a failed experiment. See "Reading the results" below and
`data/processed/{baseline_results,head_results,analysis_results}.json` for every number and
its provenance.

**Phase 2 (approved scope, not yet implemented):** a proof-of-concept application layer that
uses Phase 1's frozen models to rank predicted CRISPR dependencies for two samples and
connects the top-ranked genes to cited drug–gene interaction evidence. See "Phase 2" below
and `capstone/scope-decisions.md` for the full approved scope.

## Data and licences

| Source | Version | Licence |
|---|---|---|
| [DepMap](https://depmap.org/) — expression, CRISPR gene effect, model metadata | DepMap Public **26Q1**, downloaded 2026-07-22 | CC BY 4.0 |
| [Geneformer](https://huggingface.co/ctheodoris/Geneformer) — `Geneformer-V2-104M_CLcancer` | V2, Dec 2024 | Apache 2.0 |
| [Osteosarcoma dataset](https://osteosarc.com) (Sid Sijbrandij) — Phase 2, planned, not yet committed to this repo | — | CC0 1.0 |

**Data citation.** DepMap, Broad (2026). DepMap Public 26Q1. Dataset. https://depmap.org

Files in `data/processed/` are derived from DepMap Public 26Q1 and remain under CC BY 4.0.
Code in this repository is MIT-licensed — see `LICENSE`.

---

## Quick start

```bash
pip install pandas numpy scipy scikit-learn
pip install pyarrow          # optional; smaller/faster files if present
```

Put the three DepMap CSVs in `data/raw/` (or set `DEPMAP_DATA_DIR`), then:

```bash
python build_dataset.py     # join, filter, freeze the gene space
python splits.py            # leakage-safe train/val/test
python checks.py            # 30+ integrity assertions
python baseline.py          # the number to beat
```

Each script prints what to run next. Run them in that order; each depends on
the previous one's output.

### Required files

| Logical name | Expected filename |
|---|---|
| expression | `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` |
| crispr | `CRISPRGeneEffect.csv` |
| model | `Model.csv` |

PRISM is optional — the pipeline runs without it and picks it up automatically
if you drop it in later. Filenames vary between releases; if one isn't found
you get a message listing what *is* in the directory and how to fix it. Add
new names to `FILE_ALIASES` in `config.py`.

### Try it without any data

```bash
python make_fixture.py --out data/fixture
DEPMAP_DATA_DIR=data/fixture DEPMAP_PROCESSED_DIR=data/fixture_out python build_dataset.py
```

The fixture reproduces the real file structure, including both traps described
below, and plants a known signal so you can confirm the pipeline recovers it.

### On Kaggle or Colab

```python
import os
os.environ["DEPMAP_DATA_DIR"] = "/kaggle/input/depmap"
os.environ["DEPMAP_PROCESSED_DIR"] = "/kaggle/working/processed"
```

Set `DEBUG_NROWS = 200` in `config.py` while iterating; `None` for real runs.

---

## The two traps this defuses

**1. The expression matrix isn't indexed by ModelID.** Its first five columns
are sequencing metadata. `read_csv(..., index_col=0)` silently grabs a row
counter, and naive column slicing drags string columns into your features.

**2. There are multiple sequencing profiles per cell line.** Without filtering
on `IsDefaultEntryForModel`, ModelID isn't unique. Merging on a duplicated key
multiplies rows — the dataset looks bigger, near-identical rows land on both
sides of the train/test boundary, and the score comes out inflated. **No error
is raised at any point.** This is the most likely way this project produces a
confident wrong answer, which is why `checks.py` asserts uniqueness
independently of the code that built the file.

---

## Files

| File | Role |
|---|---|
| `config.py` | Paths, filename resolution, every threshold. Change things here. |
| `io_utils.py` | Parquet with an automatic `.npz` fallback when pyarrow is absent. |
| `gene_ids.py` | Parses `SYMBOL (ENTREZ)`, reconciles gene spaces, maps external data in. |
| `build_dataset.py` | The join. Writes matrices + a full audit report. |
| `splits.py` | Patient-grouped, lineage-stratified split. Frozen to disk. |
| `checks.py` | Independent integrity assertions. Run after any change. |
| `baseline.py` | Global mean → lineage mean → ridge on PCA of expression. |
| `train_head.py` | Ridge / MLP heads on frozen Geneformer embeddings. |
| `analysis.py` | Bootstrap CI, paired per-target comparison, effective degrees of freedom. |
| `make_fixture.py` | Synthetic data mirroring the real formats. |

`prepare_geneformer_input.py` and `run_geneformer_embeddings.py` produced
`data/processed/geneformer_embeddings.csv`, but ran on Kaggle and are not
re-runnable from a plain clone of this repository — `train_head.py` can still
refit and rescore the head exactly, since the embeddings ship as a committed
data artifact, but the embeddings themselves cannot currently be regenerated
locally. This asymmetry is intentional and disclosed, not an oversight.

---

## Outputs

Written to `data/processed/`:

- `expression`, `crispr_effect`, `prism_response` — aligned matrices
- `model_metadata` — curated cell line metadata
- **`gene_columns.json`** — the canonical gene order
- `gene_id_map` — entrez / symbol / label
- `selective_genes.json` — CRISPR targets worth predicting
- `splits.json` — frozen train/val/test assignment
- `join_report.txt` / `.json` — every filtering step, with counts
- `geneformer_embeddings.csv` — frozen 768-dim embeddings, one row per cell line
- `baseline_results.json` — expression-baseline metrics (`baseline.py`)
- `head_results.json` — Geneformer-head metrics (`train_head.py`)
- `analysis_results.json` — bootstrap CI, paired comparison, effective df (`analysis.py`)

`join_report.txt` is methods-section material. "We intersected three DepMap
tables, reconciled 18,463 genes across two matrices by Entrez ID, and retained
4,297 selective dependencies" is real provenance, and it's already written for you.

---

## Three design decisions worth defending

**Selective targets only.** Chronos scores 0 as no effect, −1 as the median
pan-essential gene. Pan-essential genes are predictable by learning a constant;
never-essential genes are noise. Neither says anything about context-specific
biology, but both inflate aggregate metrics. Only genes essential in *some*
lines and not others are kept. Thresholds are in `config.py`.

**Grouped by patient, not just cell line.** Several DepMap lines can derive
from one donor (`PatientID`). Sibling lines are far more similar to each other
than to the rest of the panel, so letting them straddle the split is quiet,
real leakage. This is the guard most projects miss. `checks.py` fails the run
if any group straddles — verified by deliberately corrupting a split and
confirming it's caught.

**Per-target Spearman, not pooled.** Correlation over all (line, gene) pairs
is dominated by between-gene differences in mean effect and looks impressive
for a model that has learned nothing about individual cell lines. Correlation
is computed per target across held-out lines, then summarised.

---

## Reading the results

Five models, increasing information, so any gain can be attributed (all scored on the same
170 val-split cell lines and 4,297 selective targets):

1. **Global mean** — the null model. Its Spearman is *undefined*, not zero: a
   constant prediction has no variance to rank. That's the point, and it's why
   R² is reported alongside. R² mean −0.0066.
2. **Lineage mean** — the control that matters. Much of what looks like learned
   biology is "this is a bone tumour, and bone tumours behave like this." Any
   model must beat this to have learned anything beyond tissue identity. Spearman mean 0.1500.
3. **Ridge on PCA of expression (the baseline)** — Spearman mean **0.2356**, R² mean 0.0690.
4. **Ridge on frozen Geneformer embeddings (the head)** — Spearman mean **0.2047**, R² mean 0.0523.
5. **MLP on the same embeddings** — Spearman mean 0.1274, R² mean −0.1217. A confirmatory
   control, not part of the pre-specified headline comparison.

**The pre-specified headline comparison is (3) vs (4), and the embeddings lose.** Paired over
the same 4,297 targets: head − baseline = −0.0308 (95% CI [−0.0365, −0.0255], SE 0.0028,
1,000 bootstrap resamples over held-out cell lines). The two models succeed and fail on
largely the same targets (per-target correlation of the two rho vectors: 0.9283 Pearson,
0.9189 Spearman), consistent with the embeddings carrying a compressed, not orthogonal,
version of the same expression signal. Five models were scored on val; only the (3)-vs-(4)
comparison was fixed before any of them were run, and it is reported here regardless of how
it came out — a rigorous negative result is a stronger capstone than a vague positive one.
Full bootstrap, Wilcoxon, and effective-degrees-of-freedom analysis: `analysis.py` and
`data/processed/analysis_results.json`.

### The test set

`splits.json` marks it; nothing reads it except `baseline.py --split test`.
Use `val` for every decision — features, hyperparameters, whether the
transformer won. Touch `test` once, at the end, and report what it says. Every
extra look quietly converts it into a validation set.

The ridge penalty is never tuned on the split being reported: on `val` it's
selected by grouped 5-fold CV inside the training set; on `test` it uses `val`.

---

## What comes next

`gene_columns.json` is the socket. Any external expression matrix — including a
patient tumour profile quantified by a completely different pipeline — must be
reindexed into exactly those columns, in exactly that order, or the trained
model receives scrambled features and returns confident nonsense.

`gene_ids.map_external_matrix()` does this, handling bare symbols, bare Entrez
IDs, or `SYMBOL (ENTREZ)`, and reports how many genes matched. Built now
because retrofitting it in September is genuinely painful.

---

## Phase 2 — proof-of-concept demo (approved scope, not yet implemented)

The October deliverable also includes a proof-of-concept application layer: using Phase 1's
frozen models to rank predicted CRISPR dependencies for two samples, then connecting the
top-ranked genes to cited drug–gene interaction evidence (retrieved from existing sources,
not predicted by any model), shown in a self-contained offline HTML report. Full approved
scope, rationale, and material-change record: `capstone/scope-decisions.md`.

Two samples, two different, explicitly labeled roles:

- **ACH-000364** (DepMap ID, cell line "U-2 OS", val split) — internal verification anchor.
  Selected for osteosarcoma identity and complete committed artifacts, not for a favorable
  predicted score; shown alongside the aggregate result across all 5 val-split osteosarcoma
  lines so it cannot be mistaken for cherry-picked validation.
  `prediction_status = held_out_prediction`, `outcome_status = measured_crispr`.
- **BG003082** — a real primary-tumor RNA-seq sample from Sid Sijbrandij's self-released
  osteosarcoma dataset (CC0 1.0), the primary demo sample, exploratory only.
  `prediction_status = exploratory_external_prediction`, `outcome_status = unavailable` — no
  CRISPR screen exists for this tissue, and bulk tumor tissue is a real domain shift from the
  cultured cell lines the model was trained and validated on.

**As of this commit, no Phase 2 code or external data lives in this repository** — no
`data/external/`, no evidence-retrieval or report-generation modules. This section describes
approved, not completed, work.

**Evidence boundaries.** Ranking a predicted CRISPR dependency and citing drug–gene
interaction evidence for it is not a treatment-efficacy estimate, not a patient-response
prediction, and not a clinical recommendation — none of those are modeled or claimed anywhere
in this project. The evidence layer itself is not a drug-response prediction either: it
retrieves existing, cited interaction records and does not predict anything. This is distinct
from `prism_response` — real drug-response *prediction* infrastructure that does exist
elsewhere in this codebase (`baseline.py`/`train_head.py`'s `task="prism"` path, `checks.py`
§7), but has never been run because no raw PRISM file has ever been downloaded; it is no more
part of Phase 2's evidence layer than it is part of Phase 1's headline result. BG003082's
prediction is never described as validated or compared to a measurement, because none exists
for that tissue.

---

## Data provenance

Record release version and access date for every file in a `LICENSES.md`.
DepMap is CC BY 4.0; Sanger-derived data carries separate terms. DepMap data is
generated for research purposes and is not intended for clinical use.
