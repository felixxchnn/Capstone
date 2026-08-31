# CLAUDE.md

Operating instructions for Claude Code in this repository. Read this before touching
anything. Last updated 2026-08-30. Phase 1 is validated, and the approved Phase 2
application layer is implemented and integrated into the repository-wide integrity checks:
`sample_profile.py`, `evidence.py`, `reconstruct_fitted.py`, `fitted_artifacts.py`,
`case_study.py`, `report.py`, and `checks.py` sections 9–12 are all committed. Nothing in
§7 or §9.7 is still pending as an uncommitted file.

---

## 1. What this project is

A high-school capstone testing one question: **do frozen Geneformer embeddings beat a
ridge-on-PCA-expression baseline at predicting CRISPR gene dependency across held-out
cancer cell lines, measured by per-target Spearman correlation?**

The answer is no, and the result is now quantified. This is a defensible negative finding
with identified mechanisms, not a failure. This comparison is **Phase 1** and is the
project's validated scientific core — none of it is discarded, rewritten, or reopened by
what follows.

**Phase 2** (approved 2026-08-25, §9.7): the October deliverable also includes an end-to-end
proof-of-concept — using Phase 1's frozen models to rank predicted dependencies for two real
samples, then connecting the top predictions to cited drug–gene evidence, shown in an offline
report. It demonstrates a possible future workflow; it does not predict patient treatment
response, and says so explicitly everywhere the distinction matters (§13). The presentation
is early October 2026.

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

Twenty-one modules today; Phase 2 (§9.7) adds more. Seventeen have been read in full.

| Module | Role | Read? |
|---|---|---|
| `config.py` | All constants, paths, file alias resolution, thresholds | yes |
| `io_utils.py` | `save_matrix` / `load_matrix` / `save_table` / `save_json`, parquet-or-npz fallback | yes |
| `baseline.py` | `prepare_task`, `impute_with_train_mean`, `per_target_spearman`, `evaluate`, `run_global_mean`, `run_lineage_mean`, `run_ridge_pca`, `_select_alpha_inner_cv`, `save_prediction_bundle`, `verify_prediction_bundle` | yes |
| `train_head.py` | `load_embeddings`, `prepare_task`, `run_ridge_head`, `run_mlp_head`, `HEAD_RIDGE_ALPHAS` | yes |
| `analysis.py` | A1–A4: bootstrap CI, Wilcoxon, per-target correlation, effective df; `--fast` vectorised path | yes (authored) |
| `random_projection.py` | E1 control (§9.3), **validation-only, exploratory**. `_fit_project_predict`: `baseline.impute_with_train_mean` → `StandardScaler`[train] → `GaussianRandomProjection(n_components=768, random_state=config.RANDOM_SEED)`[fit on standardized train] → `train_head.run_ridge_head` (its own `StandardScaler` + Ridge + patient-grouped inner-CV alpha on `HEAD_RIDGE_ALPHAS`) → `baseline.evaluate` on the 170 val lines. Writes `data/processed/random_projection_results.json` (schema `random-projection-control/1`, CRLF, deterministic within the recorded env). Gates every input against `EXPECTED_INPUT_SHA256`; hard-asserts 800/170/4297/18460/768, disjoint train/val, no patient group across the boundary, interior alpha, finite projection, projected-variance ≈ 18460/768. **No `--split test` path exists.** `--run` / `--run --save-predictions` (gitignored bundle, round-trip re-scored) / `--validate` (18 fail-closed checks incl. byte-identical recompute + protected-artifact hashes) / `--self-test` / `--check-determinism`. Carries its own validation — `checks.py` stays 55. | yes (authored) |
| `gene_ids.py` | `parse_gene_label`, `intersect_gene_spaces`, `canonical_labels`, `map_external_matrix` — Entrez is the join key throughout; no Ensembl handling here (see `prepare_geneformer_input.py`). `map_external_matrix`'s `fill_value` defaults to `0.0`, not NaN — callers reindexing a truly external sample must override it or "missing" and "measured zero" silently collide. Its symbol-matching pass uses `lookup.setdefault`, so duplicate external symbols silently keep whichever is seen first. For those two reasons `sample_profile.py` does its own canonical reindex rather than calling `map_external_matrix`. | yes |
| `build_dataset.py` | Joins DepMap CSVs into the processed matrices; `expression.npz` is `log2(TPM+1)`, confirmed both from this module's docstring and from the committed values (max 15.37) | yes (expression/CRISPR/metadata/PRISM loaders and the osteosarcoma section) |
| `prepare_geneformer_input.py` | Reconstructs pseudo-counts from log-TPM (`TPM = 2**log_tpm - 1`) since Geneformer wants raw counts and DepMap ships log-TPM; maps canonical `SYMBOL (ENTREZ)` to Ensembl via `ensembl_map.csv` (columns `entrez,ensembl_id`). The frozen Kaggle embeddings used a `mygene`-built map with 18,460/18,460 zero attrition, and that Kaggle cache was never carried back here. `data/processed/ensembl_map.csv` **now exists** — committed in `d78fbf8`, rebuilt from the static NCBI `gene2ensembl` reference (not `mygene`), coverage **18,459/18,460**; the one gap is Entrez `79400` (NOX5), which NCBI carries no Ensembl xref for. A local `prepare_geneformer_input` re-run would therefore drop NOX5 (`genes_dropped_no_ensembl == 1`); the frozen `geneformer_embeddings.csv` is unaffected. Full provenance: `capstone/data-integrity-hashes.md`. | yes |
| `sample_profile.py` | Phase 2 external-sample loader. `parse_gct` (GCT v1.2 schema validation) + `load_external_sample`: reads the committed `BG003082.gene_tpm.gct.gz` (linear TPM, versioned Ensembl IDs), strips Ensembl version suffixes, joins to canonical Entrez via `ensembl_map.csv` **only** (no symbol fallback — see the module docstring for why), sums linear TPM within a canonical gene, reindexes to `gene_columns.json` order, `log2(TPM+1)`, leaves unresolved canonical genes as explicit `NaN` (never 0), returns `(Series, provenance dict)`. Does not impute, call a model, or write to `data/processed/`. `--self-test` covers the schema/edge cases offline. | yes (authored) |
| `geneformer_sample_input.py` | Phase 2 repo-local builder (commit `0756a90`): turns the committed GCT + `ensembl_map.csv` into validated Geneformer-input frames (`X`, `var`, `obs`) — re-keys canonical Entrez → Ensembl, asserts no NaN/inf/negative, unique well-formed ENSG, no symbol fallback, extra-drop guard for a measured gene with no Ensembl id. No tokenisation/GPU/network. The GPU half (`capstone/kaggle_bg003082_embedding.py`) ran on Kaggle 2026-08-29 and produced the committed `geneformer_bg003082_embedding.csv` sidecar. `--self-test` green. | yes (authored) |
| `evidence.py` | Phase 2 DGIdb drug–gene interaction **evidence retrieval** (not treatment prediction; no model). `load_snapshot` / `get_evidence_for_gene(entrez_id, symbol=None, top_k=config.TOP_K_EVIDENCE_PER_GENE)` read the committed offline snapshot `data/external/dgidb/dgidb_2026-06b.*` — **no network**, missing snapshot raises rather than querying live. `build_snapshot` (`--refresh` / `--from-staging`) takes **5 pinned inputs** (3 DGIdb `2026-06b` TSVs + the HGNC `2026-06-02` monthly complete set + the DGIdb `2026-06b` SQL dump), gates each on exact size+SHA-256, filters to the 6 redistribution-verified interaction sources (`config.DGIDB_INCLUDED_SOURCES`; CC0/CC-BY/CC-BY-SA/US-Gov-PD only). **Gene identity is an identifier join**: `hgnc:<n>` → HGNC ID → Entrez via the HGNC file (verified 1:1; ambiguous = hard fail), kept iff canonical; symbols are a 3-way consistency check only, never a key. **PMIDs** are recovered by identifier join from the SQL dump (`interaction_claims`→`…_publications`→`publications`), `;`-joined, numerically sorted, deduped — never parsed from free text; empty where DGIdb has none (all ChEMBL/GtoP/FDA). Direction from DGIdb's `interactionClaimTypes.directionality` vocab → inhibitory/activating/unknown. `top_k` enforced per tier; every returned record carries `config.DGIDB_EVIDENCE_DISCLAIMER`. **Both committed files are byte-identical across rebuilds and `--refresh`≡`--from-staging`** — no wall-clock in any tracked output; `config.DGIDB_RETRIEVED_UTC` is the fixed provenance input, per-run times go to the git-ignored `build_runlog.jsonl`. Never commits the unfiltered TSVs, the HGNC file, or the SQL dump. `--self-test` (synthetic, offline; builds a tiny pg_dump `.sql.gz`) and `--validate` (34 checks on the committed snapshot). SQL linkage is reported both per `(gene,drug,source)` **group** (`structural_linkage`, 36,950/36,950) and per **snapshot row** (`row_level_linkage`, 37,343/37,343, 0 unlinked); the two counts differ because 373 groups hold >1 row. | yes (authored) |
| `reconstruct_fitted.py` | Phase 2 (2026-08-29). Builds **reconstructed** fitted state for the two frozen Phase 1 linear models: re-fits `impute(train-mean) → StandardScaler → PCA(200) → Ridge` (baseline) and `impute → StandardScaler → Ridge` (head) on **exactly the committed `train` split**, at the alpha **read from** `baseline_results.json` / `head_results.json` (100000.0 / 3162.0 — no selection re-run). Serialises plain `.npy` + `manifest.json` under `data/processed/reconstructed_fitted/{baseline_ridge_pca,head_ridge_head}/`. **Not the original Phase 1 objects** (never serialised) — a reproducibility convenience for `case_study.py`. `--build` / `--validate` (reproduces every committed val statistic exactly at 4 dp) / `--check-determinism` (byte-identical rebuild) / `--verify` (13/13 checklist) / `--self-test`. Ten committed inputs gated against their SHA-256 at base commit `12fab80`; refuses to build if any moved. | yes (authored) |
| `fitted_artifacts.py` | Phase 2 loader for the above. **numpy + stdlib only — no sklearn import, no `fit()` / `fit_transform()`.** `load_baseline_ridge_pca()` / `load_head_ridge()` → objects with `.predict(X)` (closed-form scale → PCA → ridge / scale → ridge), `.assert_feature_order()` / `.assert_target_order()`, `.alpha`. Every `.npy` and label file is SHA-256-verified against the artifact's own `manifest.json` on load; any mismatch / malformed manifest / shape-dtype disagreement hard-fails. | yes (authored) |
| `case_study.py` | Phase 2 (2026-08-29). Orchestrates the two reconstructed models over **ACH-000364** (U-2 OS, `val`, verification anchor — `held_out_prediction`/`measured_crispr`; observed CRISPR attached *after* ranking, verification only) and **BG003082** (osteosarcoma tumour, absent from every split — `exploratory_external_prediction`/`unavailable`, no observed data). Writes one deterministic `data/processed/case_study.json` (schema `case-study/1`): top-25 ranked predicted dependencies per model per sample (ascending predicted GeneEffect = stronger dependency, numeric-Entrez tie-break), drug-gene interaction evidence for the displayed genes (`evidence.get_evidence_for_gene`, retrieved **after** rankings freeze), and the locked five-line osteosarcoma descriptive aggregate. **Imports no sklearn, calls no `fit()`/`fit_transform()`, never imports `reconstruct_fitted`** (AST-checked); inference is `fitted_artifacts` only. `baseline.per_target_spearman` (pure numpy/scipy) is the sole `baseline` use, for the mandated osteosarcoma metric. `--build` / `--validate` (43/43, byte-identical regen ×2, protected artifacts unchanged) / `--self-test`. **Cross-platform reproducibility repair (2026-08-30, generator code only — `case_study.json` bytes and SHA-256 `a962c01a…` unchanged):** `_write_json_deterministic` now writes explicit CRLF bytes via `write_bytes` (was `Path.write_text`, which emitted LF on POSIX); the artifact's `environment` block is read from the reconstructed-fitted baseline manifest (Python 3.14.6 / numpy 2.5.0 / pandas 3.0.3 / scipy 1.18.0), not `platform.python_version()` / `np.__version__` at runtime; `_preflight` hard-fails unless the baseline and head reconstruction manifests carry identical `environment` data; `SCHEMA_VERSION` = `config.CASE_STUDY_SCHEMA_VERSION`. | yes (authored) |
| `report.py` | Phase 2 (2026-08-29). Renders the committed `case_study.json` into one self-contained offline **`phase2_report.html`** (schema `phase2-report/1`, repo root, `.gitattributes` `-text`). **No inference, no evidence lookup, no recomputation** — reads only hash-pinned committed artifacts (`case_study.json`, `baseline_results.json`/`head_results.json`/`analysis_results.json` for the section-B Phase 1 headline, DGIdb `manifest.json` for release/vintage/licence/~19%-coverage facts). CSS + JS embedded locally; no CDN / remote font / analytics; opens from `file://`; byte-identical rebuild; no wall-clock, no absolute path. Sections A–I: header + "Research demonstration — not clinical guidance", frozen Phase 1 result, sample/model selectors (ridge_pca and ridge_head never merged), 25-row ranked dependency tables (observed cols only for ACH-000364), expandable "Drug–gene interaction evidence" grouped under each gene, sample interpretation, osteosarcoma descriptive aggregate, collapsible methods/provenance/limitations. `--build` / `--validate` (25 structural checks + headless-Chrome interaction smoke: selectors/search/filter/expand) / `--self-test`. Claim-language gate hard-fails on prohibited framing in rendered prose. `REPORT_SCHEMA` / `CASE_STUDY_SCHEMA_EXPECTED` / `CASE_STUDY_SHA256` now reference `config.REPORT_SCHEMA_VERSION` / `config.CASE_STUDY_SCHEMA_VERSION` / `config.CASE_STUDY_JSON_SHA256` (2026-08-30; HTML bytes and SHA-256 `f4a093b0…` unchanged). | yes (authored) |
| `splits.py` | Patient-grouped, lineage-stratified split generation | no |
| `checks.py` | Integrity assertions; fails the run on group straddling. **12 sections, 55 fail-closed checks.** Sections 1–8 are the original dataset-integrity portion (**32/32**, unchanged — PRISM and osteosarcoma coverage included). Sections 9–12 (added 2026-08-30, commit after `6837121`) integrate the Phase 2 application layer: §9 committed-artifact identity + SHA-256 (`case_study.json` vs `config.CASE_STUDY_JSON_SHA256`, `phase2_report.html` vs `config.REPORT_HTML_SHA256`, strict JSON parse, schema version); §10 sample reconciliation + leakage prevention (`ACH-000364`/`BG003082` roles and split membership, a fresh `sample_profile.load_external_sample()` compared byte-for-byte against the committed reconciliation, the 18,427 + 33 = 18,460 / 1,407 + 17,020 = 18,427 counts, no id collisions, no symbol fallback); §11 rankings + evidence (`case_study.validate` 43/43, `evidence.validate_snapshot` 34/34, 2×25 ranked rows, derived 56-gene evidence union, coverage reconciliation, tier/source/disclaimer vocabulary, cross-checked DGIdb hashes); §12 offline-report integrity (`report.validate` 25/25, embedded JSON == committed `case_study.json`, no `<script src>` / external stylesheet / `@import url()` / `fetch(` / `XMLHttpRequest`, fixed non-efficacy disclaimer present). stdlib helpers: `_sha256_file`, `_load_strict_json` (rejects `NaN`/`Infinity`/`-Infinity`). A missing file or hash mismatch is a hard `[FAIL]`, never a skip. | yes (authored) |
| `make_fixture.py` | Synthetic fixture generation mirroring real DepMap file structure, including a synthetic Bone/Osteosarcoma subset and a synthetic PRISM matrix | yes |
| `run_geneformer_embeddings.py` | Embedding extraction; ran on Kaggle | no |
| `inspect_data.py` | Ad-hoc inspection, three hard-coded raw filenames read via bare `pd.read_csv`, not `config.resolve_file`/`RAW_DIR`. Unrunnable from a clone (raw CSVs gitignored) even before that. Ship-or-cut undecided. | yes |

**`data/processed`** holds 19 tracked files directly: the 15 Phase 1 artifacts, plus
`ensembl_map.csv` (added `d78fbf8`), plus the Phase 2 BG003082 Geneformer sidecar
`geneformer_bg003082_embedding.csv` + `geneformer_bg003082_embedding.provenance.json` (added
2026-08-29, Kaggle run), plus `case_study.json` (added 2026-08-29). The sidecar is a
`1 × 768` external-sample embedding, **not** part of the frozen `geneformer_embeddings.csv`
(1,140 × 768), which is untouched.
Plus the `reconstructed_fitted/` subtree (added 2026-08-29): `baseline_ridge_pca/` (14 files)
and `head_ridge_head/` (9 files), ~61 MiB of plain `.npy` + `manifest.json` — reconstructed
fitted state for `case_study.py`, **not** the original Phase 1 fitted objects. Hashes and the
reproduction check are in `capstone/data-integrity-hashes.md`.
`data/processed/predictions/` stays git-ignored (regenerable prediction matrices).
**`data/external/`** holds three tracked files:
`sid_osteosarc/BG003082.gene_tpm.gct.gz` (`d78fbf8`) and, added 2026-08-29,
`dgidb/dgidb_2026-06b.interactions.filtered.tsv` (**37,343** licence-filtered DGIdb
interaction records, ~14 MB, with `pmids` populated from the SQL dump) +
`dgidb/dgidb_2026-06b.manifest.json`. `.gitattributes` marks `data/external/** -text` so
those bytes survive a clone on any platform. The evidence build's 5 pinned inputs — the 3
unfiltered DGIdb release TSVs (`config.DGIDB_ASSETS`), the HGNC monthly complete set
(`config.DGIDB_HGNC_ASSET`) and the DGIdb SQL dump (`config.DGIDB_SQL_ASSET`) — are **never
committed**. `data/external/dgidb/build_runlog.jsonl` is git-ignored (the only place a
wall-clock build time is written). `data/processed/predictions/` is gitignored — 12 files,
~12 MB, regenerable in ~3 minutes from `--save-predictions` on both modules.
**PRISM (`prism_response`) is fully wired end to end in code — `config.py`'s three
`prism_*` aliases, `build_dataset.load_prism`, both `run_task` functions' `task="prism"`
branch, `checks.py` §7, `analysis.py --task prism` — but no raw PRISM file has ever been
downloaded; `join_report.json` records `"prism": {"file": None, "status": "absent"}`.** This
is a different thing from the Phase 2 drug-evidence layer in §9.7: PRISM is drug-response
*prediction* (a second ML target), the Phase 2 layer is gene–drug evidence *retrieval*, no
model attached.

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
`capstone/RESUME-HERE.md`, and `capstone/geneformer-provenance-findings.md` were referenced
here as if they existed and needed correcting. **They do not exist and never have** —
confirmed by a full `git log --all --diff-filter=A` search of this repository's entire
history. Only `capstone/CONTEXT-2026-08-24.md`, `capstone/RESUME-HERE-2026-08-23.md`, and the
Kaggle notebook snapshot have ever been committed. Treat any future reference to a capstone
doc not in that list as suspect until confirmed with `git log --all -- <path>`.
`capstone/data-integrity-hashes.md` is in the same position — referenced below as if it
already held a hash table; it does not exist yet and needs to be created, not corrected.
Corrections needed in the two docs that do exist:

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

Also: `data-integrity-hashes.md` needs to be created (it does not exist — see above), with
SHA-256 for the original 14 files plus `analysis_results.json` as the 15th, determinism
pinned to a specific commit.

The README rewrite is larger: it still describes the project as "Week 1–2," lists 8 of 13
modules, and predates `train_head.py`, `analysis.py` and the entire Geneformer arm.

### 9.2 Vectorised bootstrap — already complete

`analysis.py --fast` (commit `47fd91b`) does this: fully-observed target columns are ranked
with `rankdata(axis=0)` and correlated as one vectorised block; partially-missing columns
fall back to the original per-column loop, which is unchanged and still the default.
Verified against `baseline.per_target_spearman` (invariant 3) both on the unresampled data
and on one seeded resample with ties before being trusted; `--fast --bootstrap 1000`
reproduces the committed `analysis_results.json` byte-for-byte. Measured ~10.6× faster
(259.8s → 24.5s at 50 resamples). **Do not re-run this as new work.**

### 9.3 E1 — random projection control — DONE (committed; commit after `5db198b`)

Implemented as `random_projection.py` (§7). `GaussianRandomProjection(n_components=768,
random_state=config.RANDOM_SEED)` on train-standardised expression, then
`train_head.run_ridge_head` (identical ridge head, its own 13-point `HEAD_RIDGE_ALPHAS` grid,
patient-grouped inner-CV alpha selection on the 800 training lines), `baseline.evaluate` on
the 170 val lines. **Validation-only, exploratory** — no `--split test` path; a sixth model
on the same 170 val lines, so the val-set-optimism caveat (§11) applies.

**Result** (`data/processed/random_projection_results.json`, schema
`random-projection-control/1`, sha256 `915c4234ee6e54783d89a149fb8420d0d7fed3e00cf707855d24563cfe5ea6f7`,
9,436 B; run on Python 3.14.6 / numpy 2.5.0 / scipy 1.18.0 / scikit-learn 1.9.0 / pandas
3.0.3):

| | val Spearman mean |
|---|---|
| random-projection ridge (E1) | **0.2104** |
| Δ vs `ridge_pca` (0.2356) | **−0.0252** |
| Δ vs `ridge_head` (0.2047) | **+0.0057** |
| Δ vs `lineage_mean` (0.1500) | **+0.0604** |

Selected alpha **3162** — interior, identical to the Geneformer head; the inner-CV sweep is
unimodal (grid-min 1.0 → 0.0638, peak 3162 → 0.2066, grid-max 1e6 → 0.1767). Projected-column
variance mean **24.21**, consistent with the 18,460/768 ≈ 24.04 scale (ratio 1.007). Byte-identical
on rebuild (`--check-determinism`); `--validate` = 18/18.

**Reading (the outcome older documents treated as unlikely, and it is accepted as reported).**
A single fixed random 768-dim linear sketch of expression **matches — very slightly exceeds
(+0.0057, one seed, not bootstrapped) — the frozen Geneformer embedding** under a matched
linear head, and **both still lose to PCA-200** (−0.0252). So the embedding's deficit versus
PCA is *not* explained by the 768-dim bottleneck itself: a random bottleneck of the same
width does not pay that price. This does not beat PCA and is not framed as beating the head;
it is "no worse than the pretrained embedding at matched width." An unexpected result here
is a real finding, not a broken one, and is not to be re-run under another seed or altered
because of how it came out.

Model-set freeze (§9.6) and the one-time final test-split evaluation (F1) remain outstanding
and are unaffected by E1.

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

E1 (§9.3) is now done and committed; E3 (§9.5) is done. Freeze the model set in writing,
dated — `ridge_pca`, `ridge_head`, `lineage_mean`, `global_mean`, `mlp_head`, and the E1
random-projection control, with E1/E3 labelled exploratory. Then F1 (test split, once, all models together),
F2 (osteosarcoma — `analysis.py`'s saved-prediction machinery gives this nearly free by
resampling osteosarcoma rows instead of all 170; **5**, not ~4, val-split osteosarcoma lines,
state n explicitly), F3 (narrative and slides — the actual deliverable, draft complete two
weeks before).

### 9.7 Phase 2 — precision-oncology proof-of-concept demo

**Approved scope change, 2026-08-25** (recorded in `capstone/scope-decisions.md`): the
October deliverable is not just the Phase 1 comparison above. It adds an end-to-end demo —
dependency ranking plus a drug–gene evidence layer plus an offline UI — built on two samples,
each playing a different, explicitly labeled role:

- **ACH-000364 (U-2 OS, val split)** — internal verification anchor. Selected for
  osteosarcoma identity and complete artifacts (expression, CRISPR ground truth, embedding,
  all four saved val-split predictions), **not** for a favorable score. Shown alongside the
  aggregate result across all 5 val-split osteosarcoma lines so it doesn't read as
  cherry-picked. `prediction_status = held_out_prediction`, `outcome_status =
  measured_crispr`.
- **BG003082** — Sid Sijbrandij's real primary-tumor RNA-seq (`osteosarc.com`, CC0 1.0,
  resected 2022-12-16). Primary demo sample, exploratory only:
  `prediction_status = exploratory_external_prediction`, `outcome_status = unavailable` (no
  CRISPR screen exists for this tissue). Bulk tumor tissue is a real domain shift from the
  cultured cell lines the model was trained/validated on — every prediction on it carries
  that caveat, never phrased as measured performance.

This does not touch any Phase 1 invariant, split, or committed result. `E1`/`E3`/model-set
freeze/F1 (§9.3–9.6) proceed independently; `D1`-numbered items below are additive.
`E2` (§9.4) is deprioritized past October given the added workload — see the full plan for
the complete design (gene-ID reconciliation order, log-scale handling, the DGIdb evidence
schema, the interactive static-HTML report):
`C:\Users\Leo He\.claude\plans\moonlit-dazzling-dream.md`.

**Progress.** Committed inputs: `data/external/sid_osteosarc/BG003082.gene_tpm.gct.gz` and
`data/processed/ensembl_map.csv` (both `d78fbf8`). **`sample_profile.py` — the external-sample
loader — is implemented and validated** (`--self-test` green; real BG003082 load: 18,427 / 18,460
canonical genes resolved via Ensembl-ID join, 33 left as explicit `NaN`, 0 identifier collisions,
deterministic). It deliberately does **not** do the plan's step-4 symbol fallback (rescues only 4
of the 33; would reintroduce cross-provenance mixing the `ensembl_map.csv` single-provenance
decision rejected — see the module docstring).

**The BG003082 Geneformer embedding now exists.** `capstone/kaggle_bg003082_embedding.py` ran
on Kaggle 2026-08-29 (Tesla T4; `geneformer` pinned to `04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5`)
and produced the committed Phase 2 sidecar `data/processed/geneformer_bg003082_embedding.csv`
(`1 × 768`, SHA-256 `06a4ab9f85e5ac908975268ed502912317503ed277d28eeab1663d8305835080`) plus
`…embedding.provenance.json`. Tokenisation hard checks all passed (token length 4,096, in-vocab,
top-50 norm/rank replication). Frozen `geneformer_embeddings.csv` byte-unchanged (`af8ee6d7…`
before == after). Hashes in `capstone/data-integrity-hashes.md`. This is a **separate sidecar**,
not a row in the frozen matrix, and it changes no scientific claim: BG003082 stays
`exploratory_external_prediction` / `outcome_status = unavailable`, non-commensurable with the
Phase 1 embeddings (bulk-tumour input, NCBI-not-`mygene` map, fresh revision pin) — see
`capstone/geneformer-bg003082-feasibility.md` §3–§5.

**`evidence.py` and the committed DGIdb snapshot now exist** (2026-08-29; integrity repair
same day). Offline, licence-filtered DGIdb drug–gene interaction **evidence retrieval** —
retrieval only, no model, no efficacy/approval/indication/osteosarcoma inference. Committed
pair under `data/external/dgidb/`:
`dgidb_2026-06b.interactions.filtered.tsv` (**37,343** records, SHA-256
`f7d2089facc17ddac01e422cab8dc89d48aae463573094490f04bc42ef0a0bee`) +
`dgidb_2026-06b.manifest.json` (SHA-256
`9fb585c723cb2102a7cd335dbfac478b206d91cad04951f8ca7f70f495f6f912`; was
`d6b6f171…` before the 2026-08-29 linkage-metric repair — TSV byte-unchanged). **Both
regenerate byte-identically** across rebuilds and regardless of `--refresh` vs
`--from-staging`; no tracked output carries a wall-clock value. Records only from the 6 redistribution-verified
interaction sources (CIViC / ChEMBL / GuideToPharmacology / DoCM / NCI / FDA); the other 15
are excluded with the per-source decision recorded, and each retained record still names its
own source licence (compilation licence does not override it). **Gene identity = HGNC-ID →
Entrez join** (pinned HGNC `2026-06-02`; 1:1 verified, ambiguous = hard fail); symbols are a
consistency check only (5 disagreements recorded, e.g. `TMEM30A`/`CDC50A`). **`pmids`** are
recovered by identifier join from the pinned (temporary, uncommitted) DGIdb SQL dump —
`;`-joined/numeric-sorted/deduped, never from free text. SQL linkage is reported at two
granularities, kept distinct in the manifest: **per linkage group** (36,950 distinct
`(gene_concept,drug_concept,source)` triples, 36,950/36,950 = 100% linked;
`publications.structural_linkage`) and **per snapshot row** (37,343/37,343 rows routed
through the linkage path, 0 unlinked, categorised reasons; `publications.row_level_linkage`,
hard-asserted in `build_snapshot` + `--validate`). 36,950 ≠ the 37,343 row count because
373 groups hold >1 row. Source-skewed PMID coverage disclosed (CIViC/DoCM/NCI ~full,
ChEMBL/GtoP/FDA zero; 7,078 records overall). Three versions kept distinct: release tag
`2026-06b`, interaction data `Dec-2023`, app `v.5.0.11`/`v.5.0.12` — the records are **not**
June-2026 data. `--self-test`, `--validate` (34 checks),
`sample_profile`/`geneformer_sample_input` self-tests all green. Details: `capstone/data-integrity-hashes.md`,
`LICENSES.md`, and the manifest. (`py checks.py` is now **55/55** — §9 committed-artifact
identity + §11 `evidence.validate_snapshot` cover this layer end to end; the dataset-integrity
portion is still 32/32.)

**Reconstructed fitted artifacts now exist** (2026-08-29). The original Phase 1
`StandardScaler`/`PCA`/`Ridge` objects were never serialised, so `case_study.py` could not
run the frozen models without refitting. `reconstruct_fitted.py` + `fitted_artifacts.py`
(module map §7) close that gap: **"reconstructed fitted state at the frozen Phase 1 alpha
from the unchanged frozen training data"** — re-fit on exactly the committed `train` split
at the alpha read from the results JSONs, serialised as plain `.npy` under
`data/processed/reconstructed_fitted/`, loaded with **no `fit()`**. Reproduces every
committed Phase 1 val statistic exactly at 4 dp (`--validate`); byte-identical rebuild
(`--check-determinism`); `--verify` 13/13; `--self-test` green. They are **not** the
historical fitted objects and change no Phase 1 result.

**`case_study.py` and `data/processed/case_study.json` now exist** (2026-08-29, schema
`case-study/1`, sha256 `a962c01a…`, 169,235 B; ranking now sorts raw `float64`
predictions, rounding only for display after top-N freezes — audit-repaired 2026-08-29,
every ranked gene row byte-unchanged). It orchestrates the two reconstructed
models over ACH-000364 (`val` anchor) and BG003082 (external, absent from all splits),
producing per-model top-25 ranked predicted dependencies, drug-gene interaction evidence
for the 56 displayed genes (8 cited / 8 source-only / 40 none), and the locked five-line
osteosarcoma aggregate (mean per-target Spearman `ridge_pca` 0.119436 vs `ridge_head`
0.082773 over 4,255 common targets; Δ −0.036663; descriptive, unstable at n=5, **not** a
replacement for the frozen Phase 1 result). `--validate` 43/43, byte-identical regen.

**`report.py` and `phase2_report.html` now exist** (2026-08-29, sha256 `f4a093b0…`,
339,626 B). Self-contained offline HTML rendered from `case_study.json` — no inference /
evidence lookup / recomputation, no CDN or network, byte-identical rebuild, headless-Chrome
interaction smoke test green. **To view: open `phase2_report.html` in any modern browser
(double-click, or a `file://` URL) — no server, no build step, no network.**

**Phase 2 checks integration is complete** (2026-08-30, commit after `6837121`).
`checks.py` now runs **55 fail-closed checks in 12 sections**: the dataset-integrity
portion is unchanged at **32/32** (sections 1–8); sections 9–12 add the Phase 2
application-layer loop (committed-artifact identity + SHA-256, sample reconciliation and
leakage prevention, ranking/evidence consistency, offline-report integrity — see the
`checks.py` row in §7). `py checks.py` → `55/55 checks passed, 0 warnings, 0 failures`.

**Cross-platform deterministic `case_study.json` — repair (2026-08-30, generator code only).**
`case_study.json` was written with `Path.write_text`, which emits LF on POSIX and CRLF on
Windows, and stamped `environment` from the running interpreter — so a non-Windows rebuild
produced different bytes. Fix, in `case_study.py` / `config.py` / `report.py` only:
`_write_json_deterministic` writes explicit CRLF bytes; `environment` is read from the
reconstructed-fitted baseline manifest (the frozen build environment), not runtime;
`_preflight` hard-fails unless the two reconstruction manifests agree on `environment`; the
four schema/hash constants moved to `config.py` (`CASE_STUDY_SCHEMA_VERSION`,
`REPORT_SCHEMA_VERSION`, `CASE_STUDY_JSON_SHA256`, `REPORT_HTML_SHA256`). **The committed
`case_study.json` is byte-identical — SHA-256 still `a962c01a5b65a6ef579ea57dced67048bf9016ba0f66aab2355cf1f054796e8c`;
`phase2_report.html` still `f4a093b0…`.** No Phase 1 artifact, metric, or result changed.

**Phase 2 scope-drift check (2026-08-30) — no material change.** This task added integrity
checks and a cross-platform-reproducibility repair only. Unchanged: the dataset (DepMap
Public 26Q1 + the two named Phase 2 samples), the prediction target (per-target CRISPR gene
effect), the frozen models (`ridge_pca` α 100000, `ridge_head` α 3162 — reconstructed fitted
state, read not re-selected), the metric (per-target Spearman), the sample roles
(`ACH-000364` held_out_prediction / measured_crispr; `BG003082`
exploratory_external_prediction / unavailable), the offline-report deliverable, and every
Phase 1 artifact, result, and committed hash. `capstone/scope-decisions.md` is unchanged —
this task approves no scope change.

**Test-split discipline — precise statement.** `checks.py` reads `splits.json`
split-assignment labels for split-integrity assertions and to confirm the two Phase 2
sample roles. Beyond those label reads: no test-split expression features were loaded for
model inference; no test CRISPR outcomes were loaded for evaluation; no test predictions,
rankings, metrics, or performance results were computed or reported; and no model or
hyperparameter decision was made using test data.

---

## 10. Do not

- Do not run anything on the test split.
- Do not install `pyarrow`.
- Do not add scope beyond what's recorded in `capstone/scope-decisions.md`. Fine-tuning,
  alternative pooling, and a CPIC/PharmGKB pharmacogenomic *safety* layer (germline adverse-
  reaction prediction) are all still cut — the Phase 2 drug–gene evidence-retrieval layer in
  §9.7 is a different, narrower thing (cited interaction evidence, no model, no safety
  claims) and its approval does not reopen either of those.
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

## 13. Scope governance — mandatory

The approved research question is authoritative. Do not silently redefine it to match an available dataset, model, or implementation.

A proposed change is material if it alters any of the following:

- population or unit of analysis;
- input data;
- prediction target or label;
- evaluation metric;
- intended user or use case;
- scientific claim;
- distinction between experimental, patient-level, and clinical evidence.

Before making a material change:

1. State the current approved scope.
2. Describe the proposed change precisely.
3. Explain why it is necessary.
4. Identify which original claims would no longer be supported.
5. Present alternatives, including retaining the original scope.
6. Obtain Felix’s explicit approval.
7. Record the decision and date before implementation.

Never treat a shared motivation as proof that two research questions have the same scope. In particular, cancer-cell-line gene-dependency prediction is not equivalent to patient treatment-response prediction or personalized treatment recommendation.

At every major milestone, perform and record a scope-drift check:

> Does the current dataset, prediction target, evaluation, and deliverable still answer the explicitly approved research question?

If the answer is no or uncertain, stop implementation and ask Felix before proceeding.