# Capstone — Resume Here

**Updated 2026-08-23. Provenance investigation is CLOSED.**
Presentation early October 2026 — about 6 weeks out.
Supersedes `RESUME-HERE.md` (2026-08-20).

---

## Where to work

**`C:\Dev\Capstone`. Always. Nowhere else.**

It is a clean clone, outside OneDrive, with a healthy `.git`. Every module, all 14 data
files, and the recovered Kaggle notebook live here.

**Do not work in `OneDrive\Desktop\Capstone`.** It is superseded. It contains a dead
617 MB `.git` with no remote, and it sits inside OneDrive sync — the documented `.git`
corruption hazard that caused the 20 Aug migration in the first place.

**One thing that folder still holds, and must keep holding:** the three raw DepMap CSVs.

| File | Size |
|---|---|
| `CRISPRGeneEffect.csv` | 440,646,050 B |
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | 305,007,605 B |
| `Model.csv` | 697,455 B |

They are gitignored and are **not** on GitHub. They are only needed to re-run
`build_dataset.py`, which you have no reason to do — the processed outputs are committed.
If you ever do, `RAW_DIR` will not find them from `C:\Dev\Capstone`; set the environment
variable first:

```powershell
$env:DEPMAP_DATA_DIR = "C:\Users\Leo He\OneDrive\Desktop\Capstone"
```

**Do not edit files on the GitHub website.** Twice this session a web edit created a
commit the local clone didn't have, forcing a rebase. Edit in VS Code, push from
PowerShell. GitHub only receives.

---

## Current state

**Repository:** `github.com/felixxchnn/Capstone`, `main` @ `759583f`

Commits this session:

| SHA | Message |
|---|---|
| `759583f` | Add MIT licence, DepMap 26Q1 attribution, and licence documentation ← HEAD |
| `d7cca57` | Delete LISCENCES.md *(web UI)* |
| `75c7b7e` | Add Licences and attribution section *(web UI)* |
| `3de9590` | Add Kaggle notebook v2 snapshot: provenance record |

**New files:** `LICENSE` (MIT), `LICENSES.md`, `capstone/kaggle_notebook_v2_2026-08-06.ipynb`

---

## What was resolved this session

### The provenance question — answered directly, not inferred

**`count_source = "reconstructed_from_log_tpm"`.** Confirmed from the execution log
preserved in Kaggle notebook `notebook04de105da7`, Version 2 (Quick Save, 2026-08-06
3:25 PM), cell `In [10]`:

```
[1/4] Expression counts
No expected-counts file found -- reconstructing pseudo-counts from log-TPM.
      source     : reconstructed_from_log_tpm
      cell lines : 1140
...
  count source : reconstructed_from_log_tpm
  NOTE: pseudo-counts were used.
```

DepMap's gene-level expected-counts file was never downloaded, never uploaded to Kaggle,
and never present in any attached dataset. Only `dep-map processed` (your `processed.zip`)
and `dep-map scripts` (five `.py` files) were mounted.

### Mechanism 3c is eliminated

Same cell output:

```
[2/4] Ensembl mapping    mapped: 18460 / 18460    unmapped: 0
[3/4] Building AnnData   genes kept: 18460        genes dropped: 0
                         median n_counts: 782190
```

**Zero attrition.** Ensembl-mapping loss contributed nothing to the −0.0309. Remove it
from the mechanism list. This is a gain: the 18,460 → ≤4,096 truncation is now cleanly
attributable to Geneformer's context limit alone, with no confound.

### The `KeyError` mystery — solved

`prepare_geneformer_input.py` line 93 calls `config.resolve_file("expression_counts",
required=False)`, but `expression_counts` is not a key in `FILE_ALIASES` in any committed
version of `config.py`. `resolve_file` raises `KeyError` *before* checking `required`.

Kaggle **cell 7** patched the alias in at runtime, splicing this into the working copy:

```python
"expression_counts": (
    [
        "OmicsExpressionGenesExpectedCountProfile.csv",
        "OmicsExpressionGenesExpectedCount.csv",
    ],
    "OmicsExpression*ExpectedCount*.csv",
),
```

With the alias registered, `resolve_file` returned `None` (no matching file) → pseudo-count
branch. **`prepare_geneformer_input.py` itself ran unmodified.** The only divergence between
published and executed code was those seven lines in `config.py`.

*Caveat:* those two filenames were guessed at the time and have never been checked against
a real DepMap release. Verify before any future re-run.

### Reconstructed timeline

| When | What |
|---|---|
| 2026-07-22, 7:40–9:03 AM | DepMap 26Q1 CSVs downloaded (`RANDOM_SEED = 20260722` encodes this date) |
| 2026-08-06, 11:37:29–30 AM | Both Geneformer scripts written locally, 1 s apart |
| 2026-08-06, ~11:37 AM – 3:19 PM | Kaggle interactive session, ~3h 42m |
| — | Version 1 (Save & Run All) **failed** at 1m 44s |
| 2026-08-06, 3:19:39 PM | `geneformer_embeddings.csv` downloaded — **the only file retrieved** |
| 2026-08-06, 3:25 PM | Version 2 Quick Save — the snapshot that preserved everything |
| 2026-08-19, 4:40 PM | All 15 files browser-uploaded in one commit (`634b2c2`) |
| 2026-08-20 | Consolidation to `C:\Dev\Capstone`, data published |
| 2026-08-23 | This session |

### Other findings

- **Version 1 failed on an environment bug**, not a code bug: `AttributeError: module
  'numpy._core._multiarray_umath' has no attribute '_blas_supports_fpe'` → numpy/scipy ABI
  mismatch, at `import geneformer`. **The notebook is not re-runnable as saved.** A clean
  Save & Run All fails immediately. This is a reproducibility caveat for the methods section.
- **The embeddings came from an uncommitted interactive session.** Nothing landed in the
  Kaggle Output tab; only the CSV was downloaded by hand.
- **`notebook.ipynb` in git is 0 bytes** — uploaded empty. `RESUME-HERE.md` §4's claim that
  it is "recoverable" is technically true and practically worthless.
- **`RAW_DIR` silently falls back to `PROJECT_ROOT`** when `data/raw/` is absent — which is
  every clean clone. It doesn't warn; it just scans the wrong folder. Third silent-fallback
  bug in the project, after the alpha grid and `_select_alpha_inner_cv`.
- **1,089 vs 1,140 reconciled.** The portal's 1,089 counts Broad-only screens;
  `CRISPRGeneEffect.csv` has 1,208 cell lines (Chronos-integrated, includes Sanger). The
  1,140 intersection is sound.
- **`evaluate()` discards the per-target ρ vectors** — only summaries reach the JSONs, and
  `io_utils` has no prediction-saving path. A1/A2/A3 therefore require refitting both models.

### DepMap citation — the format changed

DepMap no longer requests a figshare DOI. The portal's Data Usage panel now specifies:

> DepMap, Broad (2026). DepMap Public 26Q1. Dataset. https://depmap.org

Release confirmed as **DepMap Public 26Q1** (April 2026). Stop hunting for a DOI — it is
not part of the requested citation.

---

## Next actions, in priority order

### 1. Two loose ends (~20 min)

- **README TODO line.** The stale instruction `Record release version and access date for
  every file in a LICENSES.md.` is still in `README.md` — the edit was never saved
  (`git commit` reported "working tree clean"). Delete it in VS Code, keep the two
  sentences after it, commit, push.
- **Delete the dead 617 MB `.git`.** Unblocked now that provenance is answered.

```powershell
$env:Path = "C:\Program Files\Git\cmd;$env:Path"
git --git-dir "C:\Users\Leo He\OneDrive\Desktop\Capstone\.git" ls-tree -r --name-only HEAD
Rename-Item -LiteralPath "C:\Users\Leo He\OneDrive\Desktop\Capstone\.git" -NewName "_git_disabled"
```

Leave it a day, then `Remove-Item ... -Recurse -Force`. **Keep the three raw CSVs.** Empty
the OneDrive *web* recycle bin afterwards or the quota won't drop.

### 2. A1–A3 — bootstrap CI, Wilcoxon, per-target ρ correlation

**The largest remaining scientific gap.** The headline −0.0309 has no standard error, no
confidence interval, and no significance test anywhere in the project.

- **A1** — resample the 170 val cell lines with replacement, recompute per-target ρ for
  both models, take the difference of means, ×1,000. Report `−0.031 [95% CI …]`.
- **A2** — Wilcoxon signed-rank on the 4,297 paired per-target ρ differences. State the
  caveat: gene dependencies are correlated (co-essential modules), so effective n is well
  below 4,297 and the p-value is optimistic.
- **A3** — Pearson/Spearman between the two per-target ρ vectors. Settles lossy-compression
  vs different-signal.

**Blocked on:** uploading `baseline.py` and `train_head.py` so the script can reuse
`prepare_task`, `run_ridge_pca`, `run_ridge_head` and `run_task` exactly — the refit must
reproduce ρ = 0.2356 and 0.2047, not something subtly different.

Bootstrap over **cell lines**, not targets. With n = 170 the interval will be wide, and
that width *is* the finding.

### 3. A4 — baseline effective df

`Σ λᵢ/(λᵢ+α)` over the 200 retained PCs at α = 100,000. The head's is 64.2 of 768. Without
both, "the two models had comparable capacity" is unsupported. ~1 hour.

### 4. E1 → E2 → E3, then freeze the model set

1. **E1 — random projection control** at 768 dims. `GaussianRandomProjection(n_components=768,
   random_state=RANDOM_SEED)` on standardised expression, identical ridge head, identical
   13-point α grid. Asks *was the representation ever informative?* — more fundamental than
   the MLP, same cost.
2. **E2 — concatenation** `[expression PCs ‖ embeddings]` through the same ridge head.
3. **E3 — MLP head.** Possibly already implemented; verify before writing anything.

**Freeze the model set after E3.** Every additional model evaluated on val makes the
reported val numbers more optimistic.

### 5. D2 — README rewrite

Still says "Week 1–2," lists 8 of 12 modules, predates `train_head.py` and the entire
Geneformer arm. Now unblocked (D1 closed). Must add: the Geneformer checkpoint and settings,
`RANDOM_SEED`, the patient-grouped CV design, the results, and the scope statement below.

### 6. D4 — corrections queue

| # | Correction | Where |
|---|---|---|
| 1 | `count_source = reconstructed_from_log_tpm` — confirmed, with source | Status report §7.1, §9 |
| 2 | Mechanism 3c **eliminated** — 18,460/18,460, zero dropped | Status report §3; provenance findings §3c |
| 3 | `median n_counts` = 782,190 (78.2% of nominal 1e6) | New, methods |
| 4 | Notebook not re-runnable — Save & Run All fails on numpy/scipy ABI | New, reproducibility caveat |
| 5 | Embeddings came from an uncommitted interactive session; only the CSV was retrieved | Methods |
| 6 | `notebook.ipynb` in git is 0 bytes | `RESUME-HERE.md` §4 |
| 7 | Reproducibility claim covers `baseline.py` only, not the Geneformer arm | Status report §1, §4.2 item 8 |
| 8 | `RAW_DIR` silently falls back to `PROJECT_ROOT` in a clean clone | Code / known defects |
| 9 | `expression_counts` alias absent; docstring line 22 cites a key that doesn't exist | Code |
| 10 | `_self_test()` doesn't reach `load_expression_counts()` — the gap that hid #9 | Methods |
| 11 | 512→768 docstring fix; reconcile 64.2 vs 64.3; git email in `RESUME-HERE.md` §6 is `felixxchnn@`, actual is `felixxchn@` | Code / docs |

### 7. F1–F4 — final phase

- **F1** — test set, run **once**, **all models together**. Only after E1–E3 and the freeze.
- **F2** — osteosarcoma val-restricted results (~4 of 26 lines; state the n explicitly),
  then the n=1 case study as qualitative application only.
- **F3** — narrative and slides. Draft complete **two weeks before**, not two days.
- **F4** — cleanup (item 1 above).

---

## Scope statement — current, ready to use

> Frozen Geneformer CLS embeddings (`Geneformer-V2-104M_CLcancer`, 4,096-token context,
> CLS pooling from layer −1) with a linear ridge head do not outperform ridge-on-PCA-expression
> for per-target Spearman prediction of CRISPR dependency across held-out DepMap cell lines
> (ρ = 0.2047 vs 0.2356).
>
> The comparison is **not input-matched**. The baseline sees all 18,460 protein-coding genes;
> Geneformer's rank-value encoding sees at most 4,096 per cell line.
>
> Geneformer's input was **reconstructed pseudo-counts** (`TPM = 2^log_TPM − 1`,
> `count_source = reconstructed_from_log_tpm`, verified from the execution log), which
> preserve within-sample relative abundance but carry TPM's gene-length normalisation — a
> systematic distortion of the ordering the model was pretrained to read. Median `n_counts`
> was 782,190 and near-constant across lines by construction, where real single-cell input
> varies with sequencing depth.
>
> Gene-mapping attrition was **not** a factor: all 18,460 genes mapped to Ensembl with none
> dropped before tokenisation.
>
> Untested: MLP head, concatenation, fine-tuning, alternative pooling, and re-running with
> DepMap expected counts.

**On re-running with real counts:** it is now a well-defined experiment and the obvious
"does mechanism 2 matter?" test. **Do not attempt it before October.** It needs a multi-GB
download, a Kaggle dataset upload, a fix for the numpy/scipy breakage, a fresh GPU
tokenise-and-extract, and a head re-run — several days with a live chance of environment
hell. Make it **future work item #1**, stated explicitly. That is a strong closing slide.

---

## Verified constants

```
RANDOM_SEED                20260722        (= the DepMap download date)
Cell lines                 1,140           (800 train / 170 val / 170 test)
Expression features        18,460 protein-coding genes
Selective CRISPR targets   4,297
PCA components             200  (75.56% of expression variance)
Embedding dimensions       768
Split design               Patient-grouped (PatientID), lineage-stratified, frozen

DepMap release             DepMap Public 26Q1 (April 2026), downloaded 2026-07-22
CRISPRGeneEffect.csv       1,208 cell lines (1,209 lines incl. header)

Geneformer checkpoint      Geneformer-V2-104M_CLcancer
Geneformer context         4,096 tokens, special_token=True, CLS pooling, layer -1
Geneformer count_source    reconstructed_from_log_tpm   ← CONFIRMED 2026-08-23
Ensembl mapping            18,460 / 18,460 mapped, 0 dropped
median n_counts            782,190

Alpha grid (both models, 13 points):
  1.0, 3.16, 10.0, 31.6, 100.0, 316.0, 1000.0, 3162.0,
  10000.0, 31623.0, 100000.0, 316228.0, 1000000.0

Baseline  ridge_pca    α = 100,000  ρ = 0.2356   median 0.2302   q25/q75 0.1439/0.3230
Head      ridge_head   α = 3,162    ρ = 0.2047   median 0.1980   q25/q75 0.1136/0.2884
                                    effective df 64.2 of 768
Control   lineage_mean              ρ = 0.1500
Control   global_mean               undefined (constant prediction, no rank variance)

Δ head − baseline    −0.0309      ← still no CI. This is A1.
Δ head − lineage     +0.0547
Δ baseline − lineage +0.0856

Kaggle notebook   kaggle.com/code/felixxchn/notebook04de105da7
                  Version 2, Quick Save, 2026-08-06 3:25 PM  (Version 1 failed)
Repository        github.com/felixxchnn/Capstone   main @ 759583f
Working copy      C:\Dev\Capstone
```

---

## Environment notes

- **`py`, not `python`** — a Microsoft Store alias shim intercepts `python`. Python 3.14.6.
- **Git is not on PATH.** Either call it in full — `& "C:\Program Files\Git\cmd\git.exe"` —
  or prepend for the session:
  ```powershell
  $env:Path = "C:\Program Files\Git\cmd;$env:Path"
  ```
- **`git --no-pager`** on every `log` / `show`, or you land in the pager (`q` to exit).
- **`git add` aborts entirely if any pathspec fails** — one typo and *nothing* stages.
  Run `git --no-pager status --short` after `add` and before `commit`.
- **Quote paths containing spaces** — `C:\Users\Leo He\...`.
- **Multi-line `if / else` breaks when pasted** into an interactive PowerShell console.
  Keep `} else {` on one line, or use single-line statements.
- **PowerShell 5.1 mangles UTF-8 on read.** `Get-Content` shows em-dashes as `â€"`. Use
  `-Encoding UTF8`. The files themselves are fine.
- **Red text from `git clone` / `git push` is usually not an error** — git writes progress
  to stderr. Read the result, not the colour.
- **`$env:USERPROFILE\Desktop` is wrong on this machine** — Desktop and Documents are
  redirected into OneDrive. Use `[Environment]::GetFolderPath('Desktop')`.
- **`pyarrow` is absent**; `io_utils` falls back to `.npz` + `.labels.json`. Both files are
  needed together. **Do not install pyarrow** — it would change the format mid-project.
- **The Windows profile folder is `C:\Users\Leo He`.** Do not rename it. Unsupported by
  Microsoft, breaks registry paths, purely cosmetic. Defer past October.
- **Git identity is correct** as `felixxchnn <felixxchn@gmail.com>`. The GitHub *username*
  has two n's; the *email* has one. Both right. `RESUME-HERE.md` §6 has the typo, not the config.

---

## Other locations — do not work in these

| Path | Status |
|---|---|
| `OneDrive\Desktop\Capstone` | Superseded. **Holds the three raw CSVs — keep those.** Dead 617 MB `.git` pending deletion |
| `OneDrive\Documents\GitHub\Capstone` | Superseded — inside OneDrive |
| `C:\Dev\_capstone_backup_2026-08-20` | Safety copy, 4 files, same physical disk |
| `C:\Dev\baseline_results.PRERUN.json`, `head_results.PRERUN.json` | Pre-re-run snapshots |
| `C:\Dev\Capstone_hashes.txt` | The 14 SHA-256 hashes |
| `C:\Users\Leo He\.copilot\repos\Capstone` | Stale June clone. Ignore |

---

## Calendar

Six weeks to early October, and school starts in roughly two. **August has materially more
available hours than September will.**

F3 (narrative and slides) is the actual deliverable — everything else is input to it.
Working backward from "draft complete two weeks before": slides underway by mid-September →
F1 done by then → E1–E3 done by early September → A1–A4 **this week**.

The hard bottleneck is calendar, not compute. Every remaining analysis item runs in minutes
on data that already exists. What takes time is writing, interpreting, and rehearsing.

**Do not add** the CPIC/PharmGKB safety layer, Geneformer fine-tuning, or alternative
pooling. Depth on a defensible null beats breadth on unfinished threads.

---

## Reference docs

- `capstone/kaggle_notebook_v2_2026-08-06.ipynb` — **the executed code.** Only surviving record
- `capstone/handoff-review-2026-08-20.md` — critical review of the original handoff
- `capstone/geneformer-provenance-findings.md` — checkpoint, truncation, pseudo-counts
  *(§3c is now superseded — attrition was zero)*
- `capstone/data-integrity-hashes.md` — SHA-256 for all 14 data files
- `capstone/status-report-2026-08-20.md` — full project status
