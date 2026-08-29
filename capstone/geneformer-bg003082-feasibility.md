# BG003082 Geneformer feasibility — repo-local determination

Created 2026-08-29. Parent commit `38177d4` (the validated BG003082 external-sample
loader). This document is written from executed code, not from memory: every count below
was produced by running `py sample_profile.py --json` and
`py geneformer_sample_input.py --self-test` on the committed bytes on this machine.

Companion code added in the same change:

- `geneformer_sample_input.py` — deterministic repo-local builder that turns the committed
  GCT + `ensembl_map.csv` into validated Geneformer-input frames (`X`, `var`, `obs`). No
  tokenisation, no GPU, no network. `--self-test` is green.
- `capstone/kaggle_bg003082_embedding.py` — the GPU half. Executable, **not yet run**.
- One-line docstring fix in `run_geneformer_embeddings.py` (512 → 768 for Geneformer V2 104M).

---

## 1. Verdict

| Question | Answer |
|---|---|
| Can the committed BG003082 GCT + `ensembl_map.csv` produce a **schema-valid** Geneformer input? | **Yes.** Built and validated locally: `X` is `1 × 18,427`, no NaN/inf/negative, unique well-formed Ensembl IDs, `obs.n_counts` = exact retained-row sum (972,338.58), `var` carries an `ensembl_id` column. |
| Is that input **scientifically in-distribution**? | **No.** It is bulk primary-tumour TPM. Every training example was a cultured DepMap cell line, tokenised from reconstructed pseudo-counts. BG003082 mixes tumour, stromal, and immune transcriptional signal in proportions no cell line has. The input is *format-compatible*, not *distribution-compatible*. |
| Does a BG003082 **embedding** exist? | **No.** None exists unless `capstone/kaggle_bg003082_embedding.py` is actually run on Kaggle/Colab and produces the sidecar files. It has not been run. |
| Can the public repo **reproduce** the embedding? | **No.** It needs the `geneformer` package, the model's LFS `.pkl` dictionaries, and a GPU — none present here. Same reproducibility asymmetry the Phase 1 Geneformer arm already carries (CLAUDE.md §11), plus the map-provenance gap in §3. |

Nothing here touches a Phase 1 invariant, split, or committed artifact. `git status` after all
of the above shows only `run_geneformer_embeddings.py` (M) and two new files; the three
results JSONs and `geneformer_embeddings.csv` are byte-unchanged, asserted inside the
builder's self-test (`_frozen_hashes()` before/after).

---

## 2. Identifier mapping — verified

`sample_profile.load_external_sample()` on the committed bytes:

```
external GCT rows                : 74,628   (all versioned ENSG, 0 unexpected, 0 duplicate)
resolved via ensembl_map.csv     : 18,427
canonical genes mapped           : 18,427
canonical genes unresolved (NaN) : 33
canonical-Entrez collisions      : 0
retained linear-TPM mass         : 972,338.58 of 1,000,000  (97.23%; the 2.77% balance is
                                   on non-coding GCT rows outside the canonical space)
```

`geneformer_sample_input.build_geneformer_frames()` then re-keys canonical Entrez → Ensembl
through `ensembl_map.csv` and asserts:

- **18,460** canonical inputs (length of `gene_columns.json`);
- **18,427** mapped values → `X` has exactly 18,427 columns;
- **33** unresolved values, left out (never imputed here — imputation is a downstream
  model-fit step);
- `X` contains **no NaN, no inf, no negative** value;
- Ensembl IDs are **unique** and match `^ENSG\d+$`;
- **no symbol fallback** is used;
- source row order — of the GCT *and* of `ensembl_map.csv` — cannot change the output
  (semantic determinism; HDF5 bytes are deliberately not compared).

### NOX5 is not a 34th missing gene

Exactly one canonical Entrez has no Ensembl ID in the committed `ensembl_map.csv`: **79400
(NOX5)** — NCBI `gene2ensembl` carries no xref for it, and its true ID
(`ENSG00000255346`) was deliberately not patched in, to keep the file single-provenance
(`capstone/data-integrity-hashes.md`). NOX5 is **already one of the 33** genes `sample_profile`
leaves unresolved for BG003082 (it is in `symbol_fallback_candidates`). The builder asserts
the coincidence explicitly: any canonical gene with a *finite* sample value but *no* Ensembl
ID would be an additional drop beyond the 33, and that raises `GeneformerInputError`. On the
committed data that count is **0** — re-keying drops nothing new. Unresolved genes are
dropped exactly once.

---

## 3. Provenance gap: the training map is gone

The frozen 1,140 training embeddings (`geneformer_embeddings.csv`) were tokenised through an
Entrez→Ensembl map built by **`mygene`** on Kaggle (18,460/18,460, zero attrition — Kaggle
notebook cell 4). **That map was never carried back into the repo** and cannot be
regenerated offline.

The committed `data/processed/ensembl_map.csv` is a **different artifact**: built from the
static **NCBI `gene2ensembl`** reference, 18,459/18,460 (NOX5 the one gap). Coverage counts
are close; the *identity* of the chosen ENSG per Entrez is not guaranteed to agree
gene-for-gene. `data-integrity-hashes.md` records **18** Entrez IDs that NCBI `gene2ensembl`
maps to more than one Ensembl gene, which the NCBI build resolved by transcript count. That
**18 is the count of NCBI-side multi-mappings only — it is not a proven upper bound on how
many genes the two maps disagree on.** `mygene` and NCBI can also pick different primary
ENSGs for Entrez IDs that are one-to-one in both maps, and nothing in the repo lets us
enumerate those, because the `mygene` map was never carried back (first paragraph above).
Any gene where the two maps chose a different ENSG tokenises under a different Geneformer
token for BG003082 than the same gene did in training.

Magnitude: **the complete cross-source identifier-disagreement count is unknown, and is
unquantifiable without the original `mygene` map.** The 18 NCBI multi-mappings are the only
disagreements the repo can even point at; the true number of differing ENSG choices could be
higher or lower. This must be disclosed wherever a BG003082 Geneformer prediction is shown
next to the training-derived ones. It is not a reason to abandon the attempt; it is a reason
not to present BG003082's Geneformer output as commensurable with the Phase 1 numbers at
face value.

### The training-run Geneformer code revision is also unrecorded

The 2026-08-06 Kaggle run cloned `ctheodoris/Geneformer` with no revision argument and
called `snapshot_download` with no `revision=` (notebook cells 1 and 21). The commit that
produced the frozen 1,140 embeddings was never captured and **cannot be recovered** from
this repo or from Kaggle. `capstone/kaggle_bg003082_embedding.py` therefore does **not**
claim to reproduce it — it pins a **new** verified revision,
`04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5`, for the BG003082 run only. That commit is the
current HEAD of `main` (an "Update README.md" commit dated 2026-05-26; repo not gated) and
predates the 2026-08-06 training run, so it is *probably* the same code Phase 1 used — but
that is an inference from commit dates, not a verified match, and it does not lift the
commensurability caveat above.

---

## 4. Pseudo-count basis (disclosed)

`geneformer_sample_input` feeds **linear TPM directly** as the pseudo-count row
(`sample_profile.load_external_sample(log2_transform=False)`). Geneformer expects raw
counts; this is an approximation. It is the *same* approximation the training embeddings
used: `logtpm_to_pseudocounts` computes `2**(log2(TPM+1)) − 1`, which round-trips to linear
TPM (float-exact to ~1e-9). Using linear TPM here is therefore the consistency-preserving
choice. Real RSEM `expected_count` is downloadable for BG003082 and is **deliberately not
used** — doing so would stack a second confound (real vs reconstructed counts) on top of the
domain shift. This matches the plan's provenance-consistency decision
(`moonlit-dazzling-dream.md`).

`n_counts` for the retained 18,427 genes is **972,338.58**, comfortably inside the training
distribution (Kaggle-reported median 782,190). Geneformer's rank encoding is
median-normalised and scale-invariant to a global factor, so the absolute value is not
load-bearing; it is recorded for provenance and because the tokeniser reads it.

---

## 5. What the Kaggle script does (and has not yet done)

`capstone/kaggle_bg003082_embedding.py` is executable, not prose. When run on Kaggle/Colab
it:

1. records `sys.version`, platform, and versions of `numpy/pandas/torch/transformers/
   anndata/scanpy/geneformer/huggingface_hub/datasets/scipy`, plus the GPU name;
2. pins the Geneformer/HF revision via `GENEFORMER_REVISION`, set in the script to the
   verified immutable commit `04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5` (current HEAD of
   `ctheodoris/Geneformer` on the Hub, "Update README.md", 2026-05-26, repo not gated) —
   the *same* revision for the `geneformer` package source, the LFS `.pkl` dictionaries, and
   the `Geneformer-V2-104M_CLcancer` weight download — and records both the pinned commit
   and the commit resolved at run time, warning if a cloned checkout's HEAD differs. This is
   a fresh pin for BG003082, **not** the Phase 1 training revision, which was never recorded
   (§3);
3. repairs the LFS `.pkl` stubs and applies the pandas-2 `.iloc[coding_miRNA_loc]` tokenizer
   patch (notebook cells 7–8), recording what was changed;
4. builds the AnnData through `geneformer_sample_input.build_bg003082_input()` — the same
   validated frames as the local run — with `model_input_size=4096`, `special_token=True`,
   `model_version="V2"`;
5. states the exact normalisation — `norm[g] = raw[g] / n_counts * 10_000 / gene_median[g]`,
   then rank-descending, keep non-zero, truncate to `MODEL_INPUT_SIZE − 2`, prepend `<cls>`,
   append `<eos>` — and **re-implements it independently** to confirm the top-50 token
   ordering matches what the tokeniser emitted;
6. verifies, before extracting anything: `ModelID == "BG003082"` survived tokenisation; the
   sequence starts with `<cls>` and ends with `<eos>`; length ≤ 4096; every token id is in
   the model vocabulary; and records a deterministic SHA-256 of the token-id list;
7. extracts one **finite `1 × 768`** CLS embedding (layer −1);
8. writes it to a **separate Phase 2 sidecar** —
   `data/processed/geneformer_bg003082_embedding.csv` and
   `…embedding.provenance.json` — with an explicit guard that the output path is not
   `geneformer_embeddings.csv`, and a before/after SHA-256 check that the frozen matrix is
   untouched;
9. prints the sidecar SHA-256 and provenance block, ready to append to
   `capstone/data-integrity-hashes.md` pinned to the commit that adds the sidecar.

**It has not been run.** No claim of a working BG003082 tokenisation or embedding may be made
in this repository until it has been, and the sidecar files exist.

---

## 6. Baseline-only fallback (the current repository state)

If no embedding is generated, the repository stays exactly as it is now:
`geneformer_sample_input` produces a validated input, and the BG003082 demo arm runs
**baseline-only**.

Concrete path: `sample_profile.load_external_sample(log2_transform=True)` → canonical
`log2(TPM+1)` row → `baseline.impute_with_train_mean` for the 33 unresolved genes → the
persisted `ridge_pca` alpha (100,000, read from `baseline_results.json`, **not** reselected)
refit on Phase 1's unchanged `X_train`/`Y_train` → ranked predicted dependencies (ascending,
most-negative-first) → DGIdb evidence retrieval.

### Claims that survive the fallback

- The end-to-end **baseline** workflow on a real external primary-tumour sample: expression
  → canonical reindex → held-alpha ridge → ranked predicted dependencies → cited drug–gene
  evidence.
- The entire Phase 1 comparison — untouched.
- The **Geneformer-vs-baseline** contrast on **ACH-000364**, whose embedding already exists
  in the frozen 1,140-row matrix.
- The 5-line osteosarcoma-aggregate statistic (computed from saved predictions, no refit).

### Claims that do NOT survive the fallback

- Any **Geneformer-vs-baseline** contrast **on BG003082**.
- Any statement about the pretrained transformer's behaviour on bulk-tumour input.

The fallback costs nothing in scientific integrity: Geneformer was already Phase 1's negative
result, and BG003082 is scoped as `exploratory_external_prediction` with
`outcome_status = unavailable`. Shipping baseline-only for BG003082, with "Geneformer
embedding unavailable for this sample" stated as a limitation, is a legitimate final state.

---

## 7. Reproduce the local determination

```powershell
cd C:\Dev\Capstone
py sample_profile.py --json                    # 18,427 mapped / 33 missing
py geneformer_sample_input.py --self-test      # synthetic + real-data gates, all green
py geneformer_sample_input.py                  # build + validate, print summary (exit 0)
py geneformer_sample_input.py --write-h5ad     # exit 1 here: clear "anndata absent" message,
                                               # no placeholder written
git status --short                             # results JSONs / embeddings unchanged
```
