# CLAUDE.md

Operating instructions for Claude Code in this repository. Read this before touching
anything. Last updated 2026-08-24 at commit `9a690d0`.

---

## 1. What this project is

A high-school capstone testing one question: **do frozen Geneformer embeddings beat a
ridge-on-PCA-expression baseline at predicting CRISPR gene dependency across held-out
cancer cell lines, measured by per-target Spearman correlation?**

The answer is no, and the result is now quantified. This is a defensible negative finding
with identified mechanisms, not a failure. The presentation is early October 2026.

The public repository is itself a deliverable — it supports university applications.
Reviewers may clone it. Code quality, honesty of documentation, and reproducibility matter
as much as the result.

---

## 2. Environment — non-negotiable

| Constraint | Detail |
|---|---|
| Python | Use **`py`**, never `python`. A Microsoft Store alias shim intercepts `python`. Python 3.14.6. |
| Git | Installed at `C:\Program Files\Git\cmd\git.exe`, **not on PATH**. Prefix each session with `$env:Path = "C:\Program Files\Git\cmd;$env:Path"` |
| Shell | Windows PowerShell. Multi-line `if / else` breaks when pasted interactively — keep `} else {` on one line. |
| **`pyarrow`** | **Absent. Do not install it.** `io_utils` falls back to `.npz` + `.labels.json` pairs. Installing pyarrow would change the on-disk format mid-project and invalidate every recorded hash. |
| Absent packages | `mygene`, `anndata`, `scanpy`. The Geneformer prep ran on Kaggle for this reason. |
| Paths | Work only in `C:\Dev\Capstone`. Deliberately outside OneDrive — `.git` corruption under sync is a documented hazard. Other copies exist under `OneDrive\` and are superseded; never edit them. |
| Line endings | Working tree is CRLF; `.gitattributes` sets `data/processed/** -text`. `git diff` on those files needs `--text` to show content. |

---

## 3. Invariants — do not break these

These are the things that make the project defensible. Breaking one silently is worse than
not making the change at all.

1. **`baseline_results.json`'s model-behaviour values must stay byte-identical.**
   Every existing key — selected alpha, the alpha sweep, all spearman and r2
   statistics — must reproduce exactly after any change to `baseline.py`. A fresh
   clone reproducing those values is the strongest reproducibility claim the
   project makes.

   Adding a *new* key that records something previously unrecorded is permitted,
   but only when every existing value is unchanged, and it must be called out in
   the commit message. Verify with:

       py baseline.py --split val
       git diff --text data/processed/baseline_results.json

   The `--text` is required; `.gitattributes` marks the directory binary. Read the
   diff. Additions only, no modifications to existing lines, or stop and diagnose.

   Last permitted schema change: `alpha_at_grid_boundary` on `ridge_pca`, commit
   after `9a690d0`.

2. **The test split has been touched zero times. Keep it that way.** Do not run
   `--split test` on anything. It runs exactly once, with all models together, only after
   the model set is frozen in writing. This is the single irreversible step in the project.

3. **The metric is imported from `baseline.py`, never reimplemented.** `per_target_spearman`
   and `evaluate` are the single source of truth for what "performance" means. If you write
   a faster path, it must ship with an assertion that it matches `per_target_spearman`
   exactly on the unresampled data before it is allowed to run.

4. **All constants live in `config.py`.** Its own docstring: "No other module should
   hard-code a filename, a path, or a threshold." `inspect_data.py` already violates this
   and is the known exception. Do not add more.

5. **`.npz` and `.labels.json` files are pairs.** Never move, copy or delete one without
   the other.

6. **Never quote a number in a comment or document that no code produced.** This has failed
   twice — the truncated alpha grid, and a hand-computed effective-df table that was wrong
   by up to 33%. Every reported figure must be traceable to a committed artifact.

7. **New feature spaces need their own alpha grid, and the selected alpha must be interior.**
   The baseline regularises 200 PCA components (variances = eigenvalues, up to thousands);
   the head regularises 768 standardised dimensions (variance 1.0 each). They independently
   selected alphas 32× apart. Any new representation is a third spectrum.

---

## 4. Verification protocol

Run this after any change to `baseline.py`, `train_head.py`, or `analysis.py`.

```powershell
cd C:\Dev\Capstone
$env:Path = "C:\Program Files\Git\cmd;$env:Path"

py baseline.py --split val          # must print alpha 100000, spearman 0.2356
py train_head.py --split val        # must print alpha 3162, ridge 0.2047, mlp 0.1274
git status --short                  # neither results JSON may appear
```

If either results file shows as modified, the change altered pipeline behaviour. Stop and
diagnose with `git diff --text data/processed/<file>.json`.

Before running anything long, verify the module imports and the constant you edited is what
you think it is:

```powershell
py -c "import train_head; print(len(train_head.HEAD_RIDGE_ALPHAS))"   # 13
py -c "import config; print(config.PCA_COMPONENTS, config.RANDOM_SEED)"  # 200 20260722
```

`make_fixture.py` exists in the repo for building small test data. I have not read it —
inspect before relying on it.

---

## 5. Known footguns

**`py train_head.py --no-mlp` silently destroys the MLP results.** It rewrites the whole of
`head_results.json` with `"mlp": null` and no `mlp_head` block. Recoverable only because the
file is committed: `git restore data/processed/head_results.json`. This has already happened
twice. Prefer running without `--no-mlp` unless you specifically need the shorter run, and
check `git status` afterwards.

**`run_ridge_pca` has no grid-boundary guard.** `run_ridge_head` sets
`alpha_at_grid_boundary` and warns; `run_ridge_pca` does not. The function that actually
suffered the truncated-grid bug is the one without the check. Adding the same guard to
`run_ridge_pca` is a small, safe improvement worth making before E1 and E2.

**`evaluate()` discards the per-target correlation vectors.** It returns summary statistics
only. Anything needing the vectors must recompute them via `per_target_spearman`, or read
the saved prediction matrices.

**`run_task` in both modules discards predictions** unless handed a `predictions_out` dict.
The model runner functions (`run_ridge_pca`, `run_ridge_head`) do return them.

**A latent no-op in both `run_task` functions:** `train_groups` is built with
`.astype(str).fillna(...)`, but `.astype(str)` converts NaN to the string `"nan"` first, so
the `fillna` never fires. Lines with no PatientID are pooled into one group. This is
conservative, not leaky, and it is how every published number was produced. **Do not "fix"
it** — that would change the results.

**Downloads to `C:\Users\Leo He\Downloads` frequently fail to arrive.** Prefer writing files
directly, or editing in place.

---

## 6. Numbers that must reproduce

Evaluated on **val**: 170 held-out cell lines, 4,297 targets.

| Model | Selected α | Spearman mean | R² mean |
|---|---|---|---|
| `ridge_pca` | 100,000 | **0.2356** | 0.0690 |
| `ridge_head` | 3,162 | **0.2047** | 0.0523 |
| `lineage_mean` | — | 0.1500 | 0.0124 |
| `mlp_head` | 0.1 | 0.1274 | −0.1217 |
| `global_mean` | — | undefined | −0.0066 |

`global_mean`'s Spearman is undefined for all targets because a constant prediction has no
rank variance. This is correct behaviour, not a bug.

**Analysis results** (`data/processed/analysis_results.json`, 1,000 bootstrap resamples):

```
Delta (head - baseline)     -0.0308   95% CI [-0.0365, -0.0255]   SE 0.0028
Share of expression gain     0.6373   95% CI [ 0.5561,  0.7080]
Per-target rho correlation   0.9283 Pearson, 0.9189 Spearman
Targets where head wins      24.65%
Wilcoxon p                   1.25e-288  (optimistic: targets are correlated)
Effective df, baseline       49.71  at alpha=100000, design 800 x 200
Effective df, head           51.78  at alpha=3162,   design 800 x 768
```

**Dataset constants:** 1,140 cell lines (800/170/170), 18,460 protein-coding genes, 4,297
selective targets, 768 embedding dims, PCA 200 at 75.56% variance, `RANDOM_SEED = 20260722`.
Splits are patient-grouped and lineage-stratified.

**Geneformer:** checkpoint `Geneformer-V2-104M_CLcancer`, 4,096 context, `special_token=True`,
CLS pooling from layer −1. `count_source = reconstructed_from_log_tpm`. Ensembl-mapping
attrition zero.

---

## 7. Module map

Thirteen modules. I have read four in full: `config.py`, `io_utils.py`, `baseline.py`,
`train_head.py`. The rest I know only from project documentation — **inspect them before
relying on any claim about their internals.**

| Module | Role | Read? |
|---|---|---|
| `config.py` | All constants, paths, file alias resolution, thresholds | yes |
| `io_utils.py` | `save_matrix` / `load_matrix` / `save_table` / `save_json`, parquet-or-npz fallback | yes |
| `baseline.py` | `prepare_task`, `impute_with_train_mean`, `per_target_spearman`, `evaluate`, `run_global_mean`, `run_lineage_mean`, `run_ridge_pca`, `_select_alpha_inner_cv`, `save_prediction_bundle`, `verify_prediction_bundle` | yes |
| `train_head.py` | `load_embeddings`, `prepare_task`, `run_ridge_head`, `run_mlp_head`, `HEAD_RIDGE_ALPHAS` | yes |
| `analysis.py` | A1–A4: bootstrap CI, Wilcoxon, per-target correlation, effective df | yes (authored) |
| `gene_ids.py` | Entrez/Ensembl handling | no |
| `build_dataset.py` | Joins DepMap CSVs into the processed matrices | no |
| `splits.py` | Patient-grouped, lineage-stratified split generation | no |
| `checks.py` | Integrity assertions; fails the run on group straddling | no |
| `make_fixture.py` | Small test data generation | no |
| `prepare_geneformer_input.py` | Tokeniser input prep; ran on Kaggle | no |
| `run_geneformer_embeddings.py` | Embedding extraction; ran on Kaggle | no |
| `inspect_data.py` | Ad-hoc inspection. Hard-codes filenames, violating `config.py`'s convention. Unrunnable from a clone since raw CSVs are gitignored. Ship-or-cut undecided. | no |

**`data/processed`** holds 15 tracked files. `data/processed/predictions/` is gitignored —
12 files, ~12 MB, regenerable in ~3 minutes from `--save-predictions` on both modules.

---

## 8. Code conventions observed in this repo

- Module docstrings are long and explain *why*, not just what. Match that register.
- Functions return `(predictions, info_dict)` where `info_dict` carries hyperparameters and
  selection provenance, which is then merged into the results JSON.
- New behaviour is added opt-in via a flag, defaulting off, so the existing command remains
  the one that produced the published artifact. `--save-predictions` is the model for this.
- Anything written to disk should be verifiable. `save_prediction_bundle` is paired with
  `verify_prediction_bundle`, which reloads and rescores from disk.
- Print output is structured with `=` rules and indented blocks. Match it.
- Warnings about methodological hazards go in comments next to the code that could trip
  them, not in a separate document.

---

## 9. Work queue

### 9.1 Documentation pass — currently the bottleneck

`capstone/handoff-review-2026-08-20.md`, `capstone/status-report-2026-08-20.md`,
`capstone/RESUME-HERE.md`, `capstone/geneformer-provenance-findings.md`, and
`capstone/data-integrity-hashes.md` all contain superseded figures. Slides drafted from them
would inherit the errors. Corrections needed:

| Stale claim | Correct value |
|---|---|
| Head effective df 64.2 / 64.3 | **51.78** |
| Baseline effective df "unknown" | **49.71** |
| "No uncertainty estimate anywhere" | Δ = −0.0308, CI [−0.0365, −0.0255] |
| Per-target ρ correlation "not computed" | 0.9283 |
| MLP head ρ "not run" | **0.1274, already committed** |
| "n=170 means a wide interval" | Falsified — pairing gives SE 0.0028 |
| `data/processed` holds 14 files | 15 |
| Ensembl attrition unknown | Zero, 18,460/18,460 |
| `count_source` unknown | `reconstructed_from_log_tpm` |
| Repo at `eef71e8` | `9a690d0` |
| 11 or 12 modules | 13 |

Also: `data-integrity-hashes.md` needs `analysis_results.json` added as the 15th file with
its hash, and the determinism claim pinned to commit `9a690d0`.

The README rewrite is larger: it still describes the project as "Week 1–2," lists 8 of 13
modules, and predates `train_head.py`, `analysis.py` and the entire Geneformer arm.

### 9.2 Vectorised bootstrap — do before E1

`analysis.py`'s bootstrap loops over 4,297 targets calling `spearmanr` on 170 values each:
~12.9 million calls per 1,000-resample run, single-threaded, **30–45 minutes**. E1, E2 and
the test run mean another 2–3 hours of waiting.

Available speedup: ridge predictions contain no NaNs, so per-column missingness is determined
entirely by `y_true` and is **fixed across resamples**. Fully-observed target columns can be
ranked with `rankdata(axis=0)` and correlated as one vectorised block; partially-missing
columns fall back to the loop. Realistically 8–15× faster.

**Requirement:** must assert exact agreement with `per_target_spearman` on the unresampled
data before running. See invariant 3.

### 9.3 E1 — random projection control

`GaussianRandomProjection(n_components=768, random_state=config.RANDOM_SEED)` on standardised
expression, then the identical ridge head and its own 13-point alpha grid. Asks whether
pretraining bought anything over a random bottleneck of the same width.

Scale note: sklearn draws components from N(0, 1/n_components), so projected columns land at
variance ≈ 18,460/768 ≈ 24 — between the head's 1.0 and the PCA block's eigenvalues. Confirm
the selected alpha is interior.

**Expect the outcome older documents treat as unlikely.** Under heavy shrinkage a random
768-dim linear sketch may match or beat PCA-200, and therefore beat Geneformer. That is a
stronger result, not a broken one: "a random linear compression of expression outperforms a
pretrained cancer transformer at matched width."

### 9.4 E2 — concatenation, with a required design fix

`[expression PCs ‖ embeddings]` through the same ridge head.

**⚠ As naively specified this will return an uninterpretable result.** In `run_ridge_pca`
the PCA scores are **not** re-standardised, so those 200 columns carry their eigenvalues as
variances (up to thousands). In `run_ridge_head` the 768 embedding columns are standardised
to unit variance. Under a single scalar alpha the high-variance PC directions are shrunk
proportionally far less, so the likely outcome is ρ_concat ≈ ρ_baseline with no way to
distinguish "the embeddings added nothing" from "the embeddings were shrunk out of
existence." That the two models independently selected alphas 32× apart is direct evidence
the spectra are incompatible under one alpha.

Fixes, ascending cost:
- (a) standardise both blocks to comparable scale, and say so in the write-up
- (b) report the fraction of fitted coefficient norm attributable to each block, as a diagnostic
- (c) block-wise alpha — the correct experiment, more work

Minimum honest version is (a) + (b). Note that (a) changes the baseline's inductive bias, so
it wants a whitened-PCs-only control, which adds a model to val.

### 9.5 E3 — already complete

MLP head scored 0.1274, below the tissue-identity control, 83% of targets at negative R².
Committed in `head_results.json`. **Do not re-run it as new work.**

One open caveat: its record carries `"alpha_at_grid_boundary": true` — it selected 0.1, the
ceiling of `MLP_ALPHAS = [1e-3, 1e-2, 1e-1]`. The sweep was flat (0.1134 / 0.1123 / 0.1135,
range 0.0012), unlike the ridge case which was monotonically accelerating at its ceiling, so
extending the grid almost certainly will not rescue it. A single confirmatory run with a
wider grid would let the slide say "we checked" rather than "we think." ~20 minutes.

### 9.6 Then

Freeze the model set in writing, dated. Then F1 (test split, once, all models together),
F2 (osteosarcoma — `analysis.py`'s saved-prediction machinery gives this nearly free by
resampling osteosarcoma rows instead of all 170; expect ~4 lines in val, state n explicitly),
F3 (narrative and slides — the actual deliverable, draft complete two weeks before).

---

## 10. Do not

- Do not run anything on the test split.
- Do not install `pyarrow`.
- Do not add scope. Fine-tuning, alternative pooling, and a CPIC/PharmGKB pharmacogenomic
  safety layer are all cut. Six weeks does not accommodate them alongside a proper test run
  and a rehearsed presentation.
- Do not "fix" the `train_groups` fillna no-op.
- Do not edit copies under `OneDrive\`.
- Do not rename the Windows profile folder.
- Do not present an inference as a measurement. If a number was not computed, say so.

---

## 11. Standing methodological cautions

**Val-set optimism.** Five models have been scored on val. E1 and E2 would make seven.
Selecting among seven on 170 lines inflates the reported numbers. The defence: the headline
comparison (baseline vs head) was fixed before E1 and E2 existed. Label it the pre-specified
comparison and E1/E2 as exploratory controls.

**Reproducibility is asymmetric and must be stated, not discovered.** The baseline arm is
reproducible end to end from the public repository. The Geneformer arm is not — the
embeddings came from an interactive Kaggle GPU session and ship as a data artifact, so the
head can be refitted and rescored exactly but the embeddings cannot currently be regenerated
from the repo. Checkpoint, context length, pooling and input provenance are all recorded.

**Disclosures that must travel with their numbers.** Bootstrap resampling with replacement
duplicates cell lines, introducing ties into the rank correlation; `spearmanr` averages tied
ranks, so the computation is correct, but the attenuation is real (measured at ~0.001 per
model) and cancels in the paired difference. The Wilcoxon p-value is optimistic because gene
dependencies are correlated. The effective-df convention is the hat-matrix trace using raw
singular values, **not** `PCA.explained_variance_`, which would rescale alpha by n−1.

---

## 12. Working style

Felix executes commands as given and pastes raw terminal output back. Give explicit,
sequentially numbered steps rather than conditional branches. **State which program each
block goes in** — PowerShell, VS Code editor, browser — before the block, and make clear
when something is text to paste into a file rather than a command to run. Write code in full
with no placeholders. Give critical, objective feedback; do not agree by default. Preserve
everything in existing files except what is being changed.
