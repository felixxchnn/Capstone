# DepMap Perturbation-Response Engine — Data Spine

Week 1–2 of the capstone: a verified, leakage-checked dataset and a baseline
number. Everything later is reported as a delta against that number.

## Data and licences

| Source | Version | Licence |
|---|---|---|
| [DepMap](https://depmap.org/) — expression, CRISPR gene effect, model metadata | DepMap Public **26Q1**, downloaded 2026-07-22 | CC BY 4.0 |
| [Geneformer](https://huggingface.co/ctheodoris/Geneformer) — `Geneformer-V2-104M_CLcancer` | V2, Dec 2024 | Apache 2.0 |
| [Osteosarcoma dataset](https://osteosarc.com) (Sid Sijbrandij) | — | CC0 1.0 |

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
| `baseline.py` | Global mean → lineage mean → ridge. Produces the number. |
| `make_fixture.py` | Synthetic data mirroring the real formats. |

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
- `baseline_results.json` — all metrics

`join_report.txt` is methods-section material. "We intersected three DepMap
tables, reconciled 19,221 genes across two matrices by Entrez ID, and retained
N selective dependencies" is real provenance, and it's already written for you.

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

## Reading the baseline

Three models of increasing information, so any gain can be attributed:

1. **Global mean** — the null model. Its Spearman is *undefined*, not zero: a
   constant prediction has no variance to rank. That's the point, and it's why
   R² is reported alongside.
2. **Lineage mean** — the control that matters. Much of what looks like learned
   biology is "this is a bone tumour, and bone tumours behave like this." Your
   model must beat this to have learned anything beyond tissue identity.
3. **Ridge on PCA of expression** — a real, simple learned model.

**If ridge barely beats the lineage mean, that is a finding, not a failure.**
Report it plainly. There's published work arguing single-cell foundation models
sometimes fail to beat simple baselines on downstream tasks — a rigorous
negative result is a stronger capstone than a vague positive one, and it's
honest.

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

## Data provenance

Record release version and access date for every file in a `LICENSES.md`.
DepMap is CC BY 4.0; Sanger-derived data carries separate terms. DepMap data is
generated for research purposes and is not intended for clinical use.
