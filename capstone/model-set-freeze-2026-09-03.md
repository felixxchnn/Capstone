# Model-set freeze — 2026-09-03

**This is the written, dated freeze of the capstone's model set, required by CLAUDE.md
§9.6 before the one-time test-split evaluation (F1).**

After this freeze, **no further model selection, no new model, no re-tuning, and no
change to any preprocessing pipeline, feature set, target set, or alpha may occur.** The
only remaining modelling step is F1: run the frozen set on the **test** split exactly
once, all eligible models together, and report.

- **Freeze date:** 2026-09-03
- **Frozen as of source commit:** `5fbf342f1536a125d80f5a3b0ed3c8f95dc58ba7` (`main`).
  All frozen result artifacts named below already exist at this commit; this document is
  added in a child commit and changes none of them.
- **Environment that produced (and reproduces) the frozen numbers:** Python **3.14.6**;
  numpy **2.5.0**; scipy **1.18.0**; pandas **3.0.3**; scikit-learn **1.9.0**. (Confirmed
  live on the freeze machine 2026-09-03; identical to the reconstructed-fitted and E1
  pinned environments. `requirements-e1.txt` pins the E1 stack.)

---

## 1. Pre-specified headline comparison

**`ridge_pca` vs `ridge_head`**, per-target Spearman across held-out cell lines, mean over
the 4,297 selective targets. This comparison was fixed in writing **before** E1 and E3
existed and is reported regardless of outcome.

| | Selected α | α interior? | val ρ (mean per-target Spearman) | val R² mean |
|---|---|---|---|---|
| `ridge_pca` (expression / PCA-200 baseline) | **100000.0** | yes (`alpha_at_grid_boundary: false`) | **0.2356** | 0.0690 |
| `ridge_head` (frozen 768-d Geneformer embeddings) | **3162.0** | yes (`alpha_at_grid_boundary: false`) | **0.2047** | 0.0523 |

**Result (authoritative, unchanged):** head − baseline = **−0.0308**, paired 95 % bootstrap
CI **[−0.0365, −0.0255]**, SE 0.0028, 1,000 cell-line resamples
(`analysis_results.json::A1_bootstrap.delta_head_minus_baseline`). The two ρ vectors
correlate 0.9283 Pearson / 0.9189 Spearman. **A rigorous, quantified negative finding:
frozen Geneformer embeddings do not beat the ridge-on-PCA-expression baseline.**

---

## 2. The frozen model set

### 2.1 Every model and control

| Model | Role | Split scored | Selected α | α at grid boundary? | val ρ mean | val R² mean | Source artifact |
|---|---|---|---|---|---|---|---|
| `ridge_pca` | **headline A** (pre-specified) | val | 100000.0 | no | 0.2356 | 0.0690 | `baseline_results.json` |
| `ridge_head` | **headline B** (pre-specified) | val | 3162.0 | no | 0.2047 | 0.0523 | `head_results.json` |
| `lineage_mean` | tissue-identity control | val | — | — | 0.1500 | 0.0124 | `baseline_results.json` |
| `global_mean` | null model | val | — | — | undefined¹ | −0.0066 | `baseline_results.json` |
| `mlp_head` | **exploratory** confirmatory control (E3) | val | 0.1 | **yes** (`MLP_ALPHAS = [1e-3, 1e-2, 1e-1]`, selected the ceiling) | 0.1274 | −0.1217 | `head_results.json` |
| E1 random-projection ridge head | **exploratory** control (§9.3) | val **only** | 3162.0 | no | 0.2104 | 0.0477 | `random_projection_results.json` |

¹ `global_mean`'s per-target Spearman is *undefined* for every target (a constant
prediction has no rank variance). This is correct behaviour; R² is reported alongside for
this reason.

### 2.2 Eligible for the one-time test evaluation (F1)

**F1-eligible (run together, once, on `test`):** `ridge_pca`, `ridge_head`,
`lineage_mean`, `global_mean`, `mlp_head`.

**Validation-only, NOT run on test:** the **E1 random-projection control**.
`random_projection.py` has **no `--split test` code path** and, per CLAUDE.md §9.3 / §11
and the E1 artifact's own `no_test_evaluation` statement, **one must not be added**. E1
stays a validation-set exploratory control. (This resolves the wording "all models
together" in CLAUDE.md §9.6: it means all models that have a test path; E1 does not and
remains val-only.)

### 2.3 Exploratory labels (val-set-optimism, CLAUDE.md §11)

Six models have now been scored on the same 170 val lines. Selecting among six on 170
lines inflates reported numbers. Mitigation, fixed here:

- The **headline comparison `ridge_pca` vs `ridge_head` is pre-specified** — fixed before
  E1 and E3 existed.
- **E1** (random-projection control) and **E3** (`mlp_head`) are **exploratory** and are
  reported regardless of outcome. They do not alter, qualify, or override the headline
  result.
- `lineage_mean` and `global_mean` are controls, not competitors.

---

## 3. E3 resolution (recorded before the freeze)

**E3 = the `mlp_head` experiment. It is complete and committed** in
`head_results.json::tasks.crispr.models.mlp_head`.

**Was an exact wider MLP alpha grid predeclared anywhere in the authoritative documents?
No.** Checked 2026-09-03:

- `train_head.py:102` defines only `MLP_ALPHAS = [1e-3, 1e-2, 1e-1]`.
- CLAUDE.md §9.5 and `capstone/CONTEXT-2026-08-24.md` §7.2 mention "a wider grid" **only in
  the abstract** — no specific alpha values, no command, no predeclared grid.
- `capstone/scope-decisions.md` (2026-08-25) lists "E3's confirmatory wider-grid run" as
  *open but independent of that decision* — never as a freeze prerequisite.
- `capstone/RESUME-HERE-2026-08-23.md` §4 says "Freeze the model set after E3" and treats
  E3 as an already-implemented experiment.

**Per the instruction "if no exact wider grid was predeclared, do not invent one now":
`mlp_head` is frozen as-is.** No wider-grid run was performed in this session. The
`RESUME-HERE-2026-08-23.md` to-do item "reconcile 64.2 vs 64.3" and the abstract
"wider grid would let the slide say 'we checked'" note are **not** blocking prerequisites
and are explicitly closed here as optional, not required.

**`mlp_head` frozen record:**

| Field | Value |
|---|---|
| Selected alpha | **0.1** |
| `alpha_at_grid_boundary` | **true** (selected the ceiling of `MLP_ALPHAS = [1e-3, 1e-2, 1e-1]`) |
| Observed inner-CV alpha sweep | 0.001 → 0.1134, 0.01 → 0.1123, 0.1 → 0.1135 — **flat**, total range 0.0012 |
| val Spearman mean | **0.1274** |
| val R² mean | −0.1217 (83 % of targets at negative R²) |
| Status | **exploratory confirmatory control**, not part of the pre-specified headline comparison |

**Why the boundary selection cannot change the headline conclusion:**

1. `mlp_head` is **not** in the pre-specified comparison (`ridge_pca` vs `ridge_head`).
2. It scores **0.1274**, *below* the `lineage_mean` tissue-identity control (0.1500) — a
   nonlinear head on the frozen embeddings is dramatically worse than the linear ridge
   head, which **strengthens** the linear ridge as the fair test of the representation.
3. Unlike the truncated ridge grid (which was **monotonically accelerating** at its
   ceiling, a genuine truncation artifact that the project's methodological narrative is
   built on), the MLP sweep is **flat** (range 0.0012). Extending the grid changes nothing
   material; the boundary hit is a plateau, not a cut-off.
4. F1 will still run `mlp_head` on `test` alongside the others, labelled exploratory, so
   the record shows it was carried through, not quietly dropped.

---

## 4. Frozen preprocessing pipelines

Both pipelines are fit on the **800 training rows only**. The original Phase 1
`StandardScaler` / `PCA` / `Ridge` objects were never serialised; the committed
`baseline_results.json` / `head_results.json` hold the model definitions, selected alphas,
alpha sweeps, and every val statistic, and are byte-stable (CLAUDE.md invariant 1). The
`data/processed/reconstructed_fitted/` subtree reproduces every committed val statistic
exactly at the recorded 4-dp precision (`reconstruct_fitted.py --validate`, re-run
2026-09-03: `ridge_pca` 0.2356 == 0.2356, `ridge_head` 0.2047 == 0.2047).

### 4.1 `ridge_pca`

```
impute(train-mean)                     # per-feature training mean; expression.npz has no NaN,
                                       #   so this vector equals the StandardScaler mean
StandardScaler(with_mean=True, with_std=True)          # fit on 800 train rows
PCA(n_components=200, whiten=False,
    svd_solver="auto" -> resolved "randomized",
    random_state=20260722)             # 200 comps, 75.56% cumulative explained variance
Ridge(alpha=100000.0, fit_intercept=True,
      solver="auto" -> resolved "cholesky", tol=1e-4)   # multi-output, 4,297 targets
```

- Feature space: **18,460** protein-coding genes, canonical `SYMBOL (ENTREZ)` order from
  `gene_columns.json` / `expression.labels.json`.
- Ordering fingerprint (SHA-256 of the newline-joined feature label list, per
  `fitted_artifacts._ReconstructedModel.order_sha256`):
  **`b9ca520c640654ec8717bb39979bac912a926e003cc2d83a5abdfc7e9832657d`**.

### 4.2 `ridge_head`

```
impute(train-mean)                     # per-dim training mean of the 768 embedding columns
StandardScaler(with_mean=True, with_std=True)          # fit on 800 train rows
Ridge(alpha=3162.0, fit_intercept=True,
      solver="auto" -> resolved "cholesky", tol=1e-4)   # multi-output, 4,297 targets
```

- Feature space: **768** Geneformer CLS embedding dimensions, columns `0..767` of
  `geneformer_embeddings.csv`.
- Ordering fingerprint:
  **`c1dfea423cceb1cdc2b5d09a01bc6ab23858fc0a4599f2911e4a3376a944ccba`**.
- Geneformer provenance: checkpoint `Geneformer-V2-104M_CLcancer`, 4,096 context,
  `special_token=True`, CLS pooling from layer −1, `count_source =
  reconstructed_from_log_tpm`, Ensembl-mapping attrition zero (18,460/18,460). The
  embeddings are a committed data artifact and are **not** regenerable from a plain clone
  (produced in an interactive Kaggle GPU session) — this asymmetry is disclosed, not an
  oversight.

### 4.3 `mlp_head` (exploratory)

```
StandardScaler(with_mean=True, with_std=True)          # fit on 800 train rows
MLPRegressor(hidden_layer_sizes=(256,), alpha=0.1,     # alpha at grid ceiling (see §3)
             random_state=20260722)                    # multi-output, 4,297 targets
```

### 4.4 `lineage_mean` / `global_mean`

No fitted transform. `global_mean` predicts each target's training mean for every line.
`lineage_mean` predicts the training mean within each `OncotreeLineage`.

### 4.5 E1 random-projection control (validation-only)

```
impute(train-mean expression)
StandardScaler[fit on 800 train rows]
GaussianRandomProjection(n_components=768, random_state=20260722)[fit on standardised train]
train_head.run_ridge_head(...)         # its own StandardScaler + multi-output Ridge,
                                       #   alpha by patient-grouped inner 5-fold CV over
                                       #   HEAD_RIDGE_ALPHAS (13 pts, 1.0..1e6); selected 3162.0 (interior)
baseline.evaluate on the 170 val lines
```

Projection-component SHA-256 (`<f8`, C-order, no `.npy` header):
`d751f201d221c1b87048f9ef83fd93d91c810a98cbaabe2c9f14dd1c03828c38`. Byte-exact
regeneration needs the pinned stack (`requirements-e1.txt`, Python 3.14.6); portable
structural verification is `random_projection.py --validate-artifact` (24 checks). **No
`--split test` path — and none is to be added.**

---

## 5. Train / validation / test roles

| Split | n cell lines | Use | Status |
|---|---|---|---|
| train | **800** | fit every pipeline; inner CV for α selection | used |
| val | **170** | report Phase 1 headline + all exploratory controls; every model/hyperparameter decision | used |
| test | **170** | F1 only — run once, after this freeze, all F1-eligible models together | **untouched — never accessed for features, outcomes, predictions, metrics, or any decision** |

- Splits are **patient-grouped** (`PatientID`) and **lineage-stratified**
  (`OncotreeLineage`); `checks.py` fails the run if any patient group straddles a split.
- `RANDOM_SEED = 20260722`, `TEST_FRACTION = 0.15`, `VAL_FRACTION = 0.15`,
  `GROUP_SPLIT_BY_LINEAGE = False`.
- Train row set fingerprint (SHA-256 of the sorted train ModelID list, identical for both
  arms): **`8df915c92d9ab2c17a959555c2b313b4549b2ff840ec5f1d7f453f54f5580fad`**.
- Split assignment: `data/processed/splits.json`
  (`f1419abc7cbd31efc173a5857bab9eb318b53f8e535a17048bfcf0ea2f70aeef`). Only the
  assignment **labels** have ever been read (for split-integrity and role assertions); no
  test-split *data* has been loaded.

---

## 6. Inner cross-validation rule

- α is selected by **grouped 5-fold CV inside the 800 training rows** (`GroupKFold` on
  `PatientID`, no shuffle) — `baseline._select_alpha_inner_cv` /
  `train_head._select_mlp_alpha_inner_cv`.
- **The penalty is never tuned on the split being reported.** On `val`, α is selected by
  the inner train-set CV above. On `test` (F1), α uses the value selected against `val`
  — i.e. the frozen α in this document (`ridge_pca` 100000.0, `ridge_head` 3162.0,
  `mlp_head` 0.1). **No α is re-selected at F1.**
- Each feature space has its **own** α grid and the selected α must be interior
  (CLAUDE.md invariant 7): `ridge_pca` on `config.RIDGE_ALPHAS` (13 pts, 1.0–1e6, selected
  100000.0, interior); `ridge_head` and E1 on `train_head.HEAD_RIDGE_ALPHAS` (13 pts,
  1.0–1e6, selected 3162.0, interior); `mlp_head` on `MLP_ALPHAS` (3 pts, 1e-3–1e-1,
  selected 0.1, **at the ceiling** — see §3).

---

## 7. Feature and target counts, ordering, hashes

| Quantity | Value |
|---|---|
| Cell lines total | 1,140 (**800** train / **170** val / **170** test) |
| Expression features (`ridge_pca`) | **18,460** protein-coding genes, canonical `SYMBOL (ENTREZ)` order |
| Embedding dims (`ridge_head`, `mlp_head`, E1 output) | **768** |
| Selective CRISPR targets (all models) | **4,297**, frozen order = `selective_genes.json` `genes` filtered to `crispr_effect` columns, order preserved |
| PCA components (`ridge_pca`) | 200, 75.56 % cumulative explained variance |
| `ridge_pca` feature-order fingerprint | `b9ca520c640654ec8717bb39979bac912a926e003cc2d83a5abdfc7e9832657d` |
| `ridge_head` feature-order fingerprint | `c1dfea423cceb1cdc2b5d09a01bc6ab23858fc0a4599f2911e4a3376a944ccba` |
| Target-order fingerprint (all models) | `cc49219ee2f7eab7023fa5478f1e7dc7cce3641eaf96167d4f8369407fb26474` |

(Fingerprints = SHA-256 of the newline-joined label list, per
`fitted_artifacts._ReconstructedModel.order_sha256`; they are ordering checks, distinct
from the file SHA-256s in §9.)

---

## 8. Evaluation metric and aggregation

- **Metric:** per-target Spearman correlation between predicted and observed CRISPR
  GeneEffect, computed **across the held-out cell lines**, via
  `baseline.per_target_spearman`. **Imported from `baseline.py`, never reimplemented**
  (CLAUDE.md invariant 3). Any faster path must ship with an assertion that it matches
  `per_target_spearman` exactly on the unresampled data.
- **Aggregation:** mean over the 4,297 targets (also median, q25, q75, frac-positive).
  **R² mean is reported alongside** every model (required because `global_mean`'s Spearman
  is undefined).
- **Not reported:** pooled correlation over all (line, target) pairs — it is dominated by
  between-target mean-effect differences and flatters models that learned nothing about
  individual lines.
- **Uncertainty (Phase 1, val):** 1,000-resample cell-line bootstrap; paired Δ over the
  4,297 targets; Wilcoxon (reported with the "gene dependencies are correlated →
  optimistic" caveat); effective df = trace of the ridge hat matrix using raw singular
  values (`ridge_pca` 49.71, `ridge_head` 51.78). Bootstrap-with-replacement introduces
  ties into the rank correlation (~0.001 attenuation per model, cancels in the paired
  difference) — disclosed.

---

## 9. Relevant result-artifact hashes (SHA-256, computed 2026-09-03)

Frozen result artifacts:

```
b49169bd363a596f400b4faff8c21d354275b70404efe08b9109d38f1bdc0ffd          3,415  data/processed/baseline_results.json
1962206fa17646cbd1fec4b642a577cc2586c09c4cabd980541a7e11a8b6f894          3,312  data/processed/head_results.json
12431dad60d07f0bd2bea9a680367007c9e030e9f17c5c20ef0b0694dcb548f9          5,428  data/processed/analysis_results.json
4adfb78b24f613adf826e7202272bbe5d95fcb9001d2f46d80567f1af319d186          9,993  data/processed/random_projection_results.json
```

Frozen inputs the F1 run will consume:

```
3d5bfa0c3430584f8943fd2365be0eecf8b994b38bfc7d491d59d7b9ff251a2d     65,301,390  data/processed/expression.npz
d18005cc0aec3e4d5f0fd06c748ef66672256cd8c7a6f24ea4c441b0ca785983        334,679  data/processed/expression.labels.json
9214efa3ce172079e6ce4ca78853d8bf92fb8f6d4a55d0c6c71e4653b59e8826     75,224,750  data/processed/crispr_effect.npz
165906f07e61819c8fadb2bf3c95a73817e538a22a34f63c415a30222ac49b9f        334,743  data/processed/crispr_effect.labels.json
af8ee6d734bea11101d07884f1c72d2b4efaff9875506738a037102a712f1e46     10,387,948  data/processed/geneformer_embeddings.csv
f1419abc7cbd31efc173a5857bab9eb318b53f8e535a17048bfcf0ea2f70aeef         33,095  data/processed/splits.json
68c8fe39ae8965ce20b04f50870609cc21734386ceeff859f4d0bddd2e5bab35         95,170  data/processed/selective_genes.json
a4b8069cc93af48f01e745bb1a15f4eaf4a7b67c9f92ca44bef3bb9e44c6d0a1        932,973  data/processed/gene_columns.json
1c314197b57c1f8363eb44f8902b3733777e7c304d7f677c76d401e3cabe5180        166,502  data/processed/model_metadata.csv
```

(Consistent with `capstone/data-integrity-hashes.md`, which pins the original set to
commit `9381c5e…`; unchanged at `5fbf342`.)

---

## 10. Locked commands intended for F1 (do NOT run in this session or before the freeze is committed)

F1 is **the single irreversible step** in the project. Run it once, all F1-eligible
models together, only after this freeze is committed. Confirm each flag against `--help`
at F1 time.

```powershell
cd C:\Dev\Capstone
$env:Path = "C:\Program Files\Git\cmd;$env:Path"

py baseline.py   --split test          # global_mean, lineage_mean, ridge_pca on test
py train_head.py --split test          # ridge_head AND mlp_head on test
                                       #   -- do NOT pass --no-mlp (it rewrites head_results.json
                                       #      with "mlp": null; recover with git restore)
py analysis.py   --task crispr --split test   # bootstrap CI + paired comparison on test
py checks.py                           # integrity gate after
```

- **E1 (`random_projection.py`) is NOT part of F1** and has no `--split test` path. Do not
  add one.
- After F1, compare the test ranking of the models to the val ranking; the headline
  `ridge_pca` vs `ridge_head` direction and magnitude are the primary read.

---

## 11. Attestations

- **The test split has not been accessed.** No test-split expression feature has been
  loaded for inference; no test-split CRISPR outcome has been loaded for evaluation; no
  test-split prediction, ranking, metric, or performance number has been computed or
  reported; no model or hyperparameter decision has used test data. Only `splits.json`
  assignment **labels** have been read, for split-integrity and role assertions.
- **No further model selection may occur after this freeze.** The model set, every
  preprocessing pipeline, every feature and target set and ordering, and every selected α
  in this document are final. F1 runs this exact set on `test` once; it selects nothing.
- **The pre-specified headline comparison (`ridge_pca` vs `ridge_head`) was fixed before
  E1 and E3 existed** and stands regardless of the F1 outcome. E1 and E3 (`mlp_head`) are
  exploratory controls, reported regardless of outcome, and do not alter the headline.
- **Val-set optimism is acknowledged:** six models scored on 170 val lines; the mitigation
  is the pre-specified headline plus explicit exploratory labelling of E1 and E3 (§2.3).

---

## 12. Sign-off

Frozen 2026-09-03 at source commit `5fbf342f1536a125d80f5a3b0ed3c8f95dc58ba7`.
Next modelling action: **F1** (§10), in a separate session, after this document is
committed.
