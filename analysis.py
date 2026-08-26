"""
analysis.py
===========
Uncertainty and capacity analysis for the baseline-vs-Geneformer comparison.

    py analysis.py                      # val split, 1000 bootstrap resamples
    py analysis.py --bootstrap 200      # quicker pass while iterating
    py analysis.py --split test         # final; only after the model set is frozen

Output
------
    analysis_results.json   every number below, with its inputs recorded

What this closes
----------------
The headline of this project is a difference of two Spearman means:

    ridge on Geneformer embeddings   0.2047
    ridge on PCA of expression       0.2356
                                   ---------
                                    -0.0309

and until now that difference carried no standard error, no confidence
interval, and no significance test. "Geneformer loses" was an unquantified
claim. This module quantifies it, and reports two further things the record
was missing.

A1 -- Bootstrap confidence interval over held-out cell lines
    Resample the held-out lines with replacement, recompute per-target Spearman
    for every model on the resampled lines, take the difference of means, and
    repeat. The 2.5th and 97.5th percentiles of that distribution are the
    interval.

    The resampling unit is the cell line, not the target, because the claim
    being tested is a generalisation claim: would this difference survive a
    different draw of held-out cell lines? Resampling targets would answer a
    weaker question. With ~170 lines the interval will be wide, and that width
    is the finding, not a defect in it.

    Resampling with replacement duplicates cell lines, which introduces ties
    into the rank correlation. scipy's spearmanr handles ties by averaging
    ranks, so the computation is correct -- but the tie behaviour is a
    consequence of the method and must be stated in the write-up, not left for
    a reader to discover.

A2 -- Wilcoxon signed-rank test over targets
    Both models are scored on the same held-out lines and the same targets, so
    the per-target correlations are paired. This tests whether the median
    paired difference is zero.

    Caveat, to be stated wherever the p-value is: gene dependencies are
    correlated -- co-essential genes move together -- so the effective number
    of independent targets is well below the nominal count and the p-value is
    optimistic. It is reported because a reader will ask for it, not because it
    is the strongest evidence here. The bootstrap in A1 is.

A3 -- Correlation between the two models' per-target performance
    If the two per-target correlation vectors are themselves strongly
    correlated, the embeddings are a lossier compression of the same signal:
    both models do well and badly on the same targets. If they are weakly
    correlated, the embeddings carry different information rather than less of
    it. This settles a reading the project has so far marked "suggestive, not
    proven".

A4 -- Effective degrees of freedom, both models, one convention
    Effective df for ridge is the trace of the hat matrix,

        df(alpha) = sum_i  s_i^2 / (s_i^2 + alpha)

    where s_i are the singular values of the *centred design matrix the Ridge
    estimator actually sees*, on the training rows it was actually fitted on.

    The convention matters more than it looks. sklearn's PCA exposes
    `explained_variance_`, which is s_i^2 / (n - 1). Substituting that into the
    formula above rescales alpha by a factor of (n - 1) -- roughly 800 here --
    and produces a plausible-looking number that is wrong by two orders of
    magnitude. Both models are therefore measured by one function, on the same
    convention, from the same kind of matrix.

    Without both numbers, "the two models were given comparable capacity" --
    the premise of the whole comparison -- is unsupported.

Also reported: the share of the expression-over-lineage gain that the
embeddings recover. The baseline beats the tissue-identity control by some
margin; the head beats it by a smaller one; the ratio of those two margins is
a far more interpretable statement of the result than a raw difference of
correlations, and it gets its own bootstrap interval.

How predictions are obtained
----------------------------
Preferred: read them from <processed-dir>/predictions/, written by
`baseline.py --save-predictions` and `train_head.py --save-predictions`.

Fallback: refit both models here, using the same prepare/impute/fit functions
the two pipelines use, and assert that the refit reproduces the published
Spearman means before any analysis runs. If it does not, the run stops -- a
confidence interval around a number that is not the reported number would be
worse than no interval at all.

This module never modifies the pipeline's own outputs. It reads
baseline_results.json and head_results.json, and writes only
analysis_results.json.

Memory note: the refit path loads the full expression matrix and converts it to
float64, which peaks in the region of half a gigabyte. The effective-df step
needs the training design matrices regardless of where the predictions came
from, so the expression matrix is loaded either way and released as soon as it
is no longer needed.
"""

from __future__ import annotations

import sys
import time
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import config
import io_utils
import baseline
import train_head


# --------------------------------------------------------------------------
# Effective degrees of freedom
# --------------------------------------------------------------------------

def squared_singular_values(design: np.ndarray) -> np.ndarray:
    """
    Squared singular values of a centred design matrix.

    Ridge with fit_intercept=True centres its input, so the hat matrix is built
    from the centred design. Centring here rather than assuming it lets the same
    function take PCA scores (already centred, so this is a no-op) and
    standardised embeddings (also centred) without either caller needing to know.

    Zero and non-finite values are dropped: they contribute nothing to the trace
    and would only introduce division noise.
    """
    design = np.asarray(design, dtype=float)
    centred = design - design.mean(axis=0, keepdims=True)
    singular = np.linalg.svd(centred, compute_uv=False)
    squared = np.asarray(singular, dtype=float) ** 2
    return squared[np.isfinite(squared) & (squared > 0.0)]


def effective_df(squared: np.ndarray, alpha: float) -> float:
    """
    Trace of the ridge hat matrix: sum s_i^2 / (s_i^2 + alpha).

    Ranges from min(n_samples, n_features) at alpha -> 0 down towards zero as
    alpha grows. At large alpha it converges to trace(Z'Z) / alpha, which is a
    useful sanity check on which matrix a reported figure was measured from.
    """
    alpha = float(alpha)
    return float(np.sum(squared / (squared + alpha)))


def effective_df_sweep(squared: np.ndarray, alphas: list[float]) -> dict:
    """Effective df at every alpha in a grid, from one SVD."""
    return {str(alpha): round(effective_df(squared, alpha), 2) for alpha in alphas}


# --------------------------------------------------------------------------
# Task assembly
# --------------------------------------------------------------------------

def assemble(
    prepared: tuple,
    metadata: pd.DataFrame,
    assignment: pd.Series,
    eval_split: str,
) -> dict:
    """
    Turn a prepare_task result into the arrays the fitting functions expect.

    This reproduces, line for line, what baseline.run_task and
    train_head.run_task do between prepare_task and the model call: the
    train/eval slicing, the train-mean imputation, the val arrays used for
    alpha selection on the test split, and the patient grouping vector.

    It is duplicated rather than imported because both run_task functions score
    their predictions and then discard them, so there is no way to reach this
    intermediate state through the public path. The duplication is guarded: the
    refit is checked against the published Spearman means before it is used,
    and a mis-copy here would fail that check immediately.

    Note on the grouping vector: `.astype(str)` converts a missing PatientID to
    the string "nan" before `.fillna` is reached, so the fillna never fires and
    lines with no patient ID are pooled into a single group. That is what both
    pipelines do. It is conservative rather than leaky -- pooling forces those
    lines into the same CV fold -- and it is copied verbatim, because "fixing"
    it here would make the refit diverge from the published run.
    """
    X, Y, train_index, eval_index = prepared

    X_train_raw = X.loc[train_index].to_numpy(dtype=float)
    X_eval_raw = X.loc[eval_index].to_numpy(dtype=float)
    Y_train_raw = Y.loc[train_index].to_numpy(dtype=float)
    Y_eval_raw = Y.loc[eval_index].to_numpy(dtype=float)

    X_train, X_eval = baseline.impute_with_train_mean(X_train_raw, X_eval_raw)
    Y_train_filled, = baseline.impute_with_train_mean(Y_train_raw)

    X_val = Y_val = None
    if eval_split == "test":
        val_index = X.index[(assignment.reindex(X.index) == "val").to_numpy()]
        if len(val_index) >= 5:
            X_val_raw = X.loc[val_index].to_numpy(dtype=float)
            X_val, = baseline.impute_with_train_mean(X_train_raw, X_val_raw)[1:]
            Y_val = Y.loc[val_index].to_numpy(dtype=float)

    train_groups = None
    if config.GROUP_COL in metadata.columns:
        train_groups = (
            metadata.loc[train_index, config.GROUP_COL]
            .astype(str)
            .fillna(pd.Series(train_index, index=train_index).astype(str))
            .to_numpy()
        )

    return {
        "X_train": X_train,
        "X_eval": X_eval,
        "Y_train_filled": Y_train_filled,
        "Y_train_raw": Y_train_raw,
        "Y_eval_raw": Y_eval_raw,
        "X_val": X_val,
        "Y_val": Y_val,
        "train_groups": train_groups,
        "train_index": train_index,
        "eval_index": eval_index,
        "target_columns": [str(c) for c in Y.columns],
    }


def baseline_design_matrix(X_train: np.ndarray, seed: int) -> np.ndarray:
    """
    Rebuild the matrix baseline.run_ridge_pca hands to Ridge.

    Standardise, then PCA with the same component clamp and the same seed.
    PCA's randomised solver is seeded, so this is deterministic and reproduces
    the fitted model's design exactly rather than approximately.
    """
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)

    n_components = int(min(
        config.PCA_COMPONENTS, X_train_s.shape[0] - 1, X_train_s.shape[1]
    ))
    n_components = max(n_components, 1)

    pca = PCA(n_components=n_components, random_state=seed)
    return pca.fit_transform(X_train_s)


def head_design_matrix(X_train: np.ndarray) -> np.ndarray:
    """Rebuild the matrix train_head.run_ridge_head hands to Ridge."""
    scaler = StandardScaler()
    return scaler.fit_transform(X_train)


# --------------------------------------------------------------------------
# Prediction acquisition
# --------------------------------------------------------------------------

def discover_saved_predictions(
    processed_dir: Path,
    source: str,
    task: str,
    eval_split: str,
) -> dict:
    """
    Find prediction matrices written by --save-predictions.

    Returns a mapping of model name -> stem path, including "y_true". Empty if
    the directory or the files are absent, which is the signal to refit.
    """
    pred_dir = Path(processed_dir) / "predictions"
    if not pred_dir.is_dir():
        return {}

    prefix = f"{source}_{task}_{eval_split}_"
    found: dict[str, Path] = {}
    for path in sorted(pred_dir.iterdir()):
        if not path.is_file() or path.suffix not in (".npz", ".parquet"):
            continue
        stem = path.stem
        if not stem.startswith(prefix):
            continue
        found[stem[len(prefix):]] = pred_dir / stem
    return found


def load_saved_predictions(stems: dict) -> tuple[np.ndarray, dict, pd.Index, list]:
    """
    Load a saved bundle.

    Returns (y_true, {model: array}, eval_index, target_columns). The row index
    and column labels come from the truth matrix and every model matrix is
    checked against them, so a bundle assembled from mismatched runs cannot pass
    silently.
    """
    if "y_true" not in stems:
        raise FileNotFoundError(
            "Found saved prediction matrices but no y_true matrix alongside "
            "them. The bundle is incomplete; re-run with --save-predictions "
            "or pass --force-refit."
        )

    truth_frame = io_utils.load_matrix(stems["y_true"])
    eval_index = truth_frame.index
    target_columns = [str(c) for c in truth_frame.columns]

    models: dict[str, np.ndarray] = {}
    for name, stem in stems.items():
        if name == "y_true":
            continue
        frame = io_utils.load_matrix(stem)
        if not frame.index.equals(eval_index):
            raise ValueError(
                f"Saved predictions for {name!r} have a different cell-line "
                f"index from the truth matrix. The bundle is inconsistent."
            )
        if [str(c) for c in frame.columns] != target_columns:
            raise ValueError(
                f"Saved predictions for {name!r} have different target columns "
                f"from the truth matrix. The bundle is inconsistent."
            )
        models[name] = frame.to_numpy(dtype=float)

    return truth_frame.to_numpy(dtype=float), models, eval_index, target_columns


# --------------------------------------------------------------------------
# Bootstrap
# --------------------------------------------------------------------------

def bootstrap_over_cell_lines(
    y_true: np.ndarray,
    predictions: dict,
    n_boot: int,
    seed: int,
    progress_every: int = 50,
) -> dict:
    """
    Resample held-out cell lines with replacement and rescore every model.

    Each resample draws n row indices with replacement from the n held-out
    lines, slices the truth matrix and every prediction matrix by the same
    indices -- so the models stay paired within a resample, which is what makes
    the difference of means meaningful -- and recomputes the mean per-target
    Spearman.

    The metric is baseline.per_target_spearman, not a faster reimplementation
    of it, so the bootstrap distribution is centred on the same measurement the
    project reports everywhere else. That costs runtime: the metric loops over
    targets in Python, so expect a few minutes per thousand resamples per model.

    Returns {model_name: array of length n_boot}.
    """
    rng = np.random.default_rng(seed)
    n_lines = y_true.shape[0]
    names = list(predictions)
    draws = {name: np.full(int(n_boot), np.nan, dtype=float) for name in names}

    for b in range(int(n_boot)):
        idx = rng.integers(0, n_lines, size=n_lines)
        y_resampled = y_true[idx]
        for name in names:
            rho = baseline.per_target_spearman(y_resampled, predictions[name][idx])
            if np.isfinite(rho).any():
                draws[name][b] = float(np.nanmean(rho))
        if progress_every and (b + 1) % progress_every == 0:
            print(f"        {b + 1}/{int(n_boot)} resamples", flush=True)

    return draws


def _colwise_pearson(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Pearson correlation between corresponding columns of two (n, k) arrays.

    Applied to ranks, this is exactly what scipy.stats.spearmanr computes for
    an untied-enough column: the rank correlation. A constant column (all
    values tied, hence all ranks tied) gives a zero-variance column and a 0/0
    division, which numpy turns into NaN -- the same "undefined, not zero"
    convention baseline.per_target_spearman uses for constant input.
    """
    a = a - a.mean(axis=0, keepdims=True)
    b = b - b.mean(axis=0, keepdims=True)
    num = np.sum(a * b, axis=0)
    denom = np.sqrt(np.sum(a * a, axis=0) * np.sum(b * b, axis=0))
    with np.errstate(divide="ignore", invalid="ignore"):
        return num / denom


def fast_per_target_spearman(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    min_samples: int = 5,
    fully_observed_mask: np.ndarray | None = None,
) -> np.ndarray:
    """
    Vectorised Spearman per target column, section 9.2 of CLAUDE.md.

    Columns with no missing values in either y_true or y_pred are ranked and
    correlated as one block with scipy.stats.rankdata(axis=0) and a vectorised
    Pearson correlation of the ranks -- no per-column Python loop. Columns with
    any NaN fall back to baseline.per_target_spearman, unchanged, so this
    function's output is defined to match it column for column.

    fully_observed_mask lets a caller who resamples rows (the bootstrap) supply
    a mask computed once from the unresampled data, since resampling with
    replacement never introduces a NaN into a column that had none, or removes
    one from a column that did -- missingness is a property of which rows
    exist, not of which are drawn.

    This function is never used for anything until
    verify_fast_bootstrap_matches has confirmed it reproduces
    baseline.per_target_spearman exactly, on the unresampled data, for every
    model it will be run on. That check, not this docstring, is the guarantee.
    """
    n_targets = y_true.shape[1]
    n_rows = y_true.shape[0]

    if fully_observed_mask is None:
        fully_observed_mask = (
            np.all(np.isfinite(y_true), axis=0) & np.all(np.isfinite(y_pred), axis=0)
        )
    if n_rows < min_samples:
        fully_observed_mask = np.zeros(n_targets, dtype=bool)

    rho = np.full(n_targets, np.nan, dtype=float)

    if fully_observed_mask.any():
        t_ranks = stats.rankdata(y_true[:, fully_observed_mask], axis=0)
        p_ranks = stats.rankdata(y_pred[:, fully_observed_mask], axis=0)
        rho[fully_observed_mask] = _colwise_pearson(t_ranks, p_ranks)

    remainder = ~fully_observed_mask
    if remainder.any():
        rho[remainder] = baseline.per_target_spearman(
            y_true[:, remainder], y_pred[:, remainder], min_samples=min_samples
        )

    return rho


def self_test_fast_fallback_branch(seed: int = 0, tol: float = 1e-9) -> None:
    """
    Exercise fast_per_target_spearman's fallback branch on synthetic data.

    Whether the fallback (the per-column baseline.per_target_spearman path for
    columns with a NaN in y_true or y_pred) ever runs on real data depends on
    whether the truth matrix happens to have a missing value in some target
    column -- not guaranteed, and not something this module controls. This
    builds a small matrix with NaNs planted in some columns and none in
    others, so both the vectorised block and the fallback loop run in the same
    call, and asserts fast_per_target_spearman reproduces
    baseline.per_target_spearman exactly on it. Independent of whatever main()
    finds when it checks the real val truth matrix.
    """
    rng = np.random.default_rng(seed)
    n_rows, n_targets = 40, 6
    y_true = rng.normal(size=(n_rows, n_targets))
    y_pred = rng.normal(size=(n_rows, n_targets))

    # Columns 1 and 3 get a NaN in y_true; column 4 gets one in y_pred.
    # Columns 0, 2, 5 stay fully observed, so the vectorised block also runs.
    y_true[0:3, 1] = np.nan
    y_true[5, 3] = np.nan
    y_pred[2, 4] = np.nan

    fully_observed = (
        np.all(np.isfinite(y_true), axis=0) & np.all(np.isfinite(y_pred), axis=0)
    )
    if fully_observed.all() or not fully_observed.any():
        raise AssertionError(
            "self_test_fast_fallback_branch's synthetic data does not mix "
            "fully-observed and partially-missing columns; it would not "
            "exercise both branches. Fix the test."
        )

    rho_reference = baseline.per_target_spearman(y_true, y_pred)
    rho_fast = fast_per_target_spearman(y_true, y_pred)

    nan_reference = np.isnan(rho_reference)
    nan_fast = np.isnan(rho_fast)
    if not np.array_equal(nan_reference, nan_fast):
        raise AssertionError(
            "self_test_fast_fallback_branch: fast_per_target_spearman's "
            "fallback branch disagrees with baseline.per_target_spearman on "
            "which targets are defined."
        )
    finite = ~nan_reference
    if finite.any() and not np.allclose(
        rho_reference[finite], rho_fast[finite], atol=tol, rtol=0.0
    ):
        max_diff = float(np.max(np.abs(rho_reference[finite] - rho_fast[finite])))
        raise AssertionError(
            f"self_test_fast_fallback_branch: fast_per_target_spearman's "
            f"fallback branch diverges from baseline.per_target_spearman by "
            f"up to {max_diff:.3e}."
        )


def _assert_fast_matches_reference(
    y_true: np.ndarray,
    predictions: dict,
    case_label: str,
    tol: float,
) -> None:
    """One case (unresampled, or one resample) of the --fast gate."""
    for name, pred in predictions.items():
        rho_reference = baseline.per_target_spearman(y_true, pred)
        rho_fast = fast_per_target_spearman(y_true, pred)

        nan_reference = np.isnan(rho_reference)
        nan_fast = np.isnan(rho_fast)
        if not np.array_equal(nan_reference, nan_fast):
            mismatched = int(np.sum(nan_reference != nan_fast))
            raise AssertionError(
                f"--fast bootstrap verification failed on {case_label} data "
                f"for model {name!r}: the vectorised path and "
                f"baseline.per_target_spearman disagree on which of "
                f"{rho_reference.size} targets are defined ({mismatched} "
                f"targets differ). Refusing to run the fast path."
            )

        finite = ~nan_reference
        if finite.any() and not np.allclose(
            rho_reference[finite], rho_fast[finite], atol=tol, rtol=0.0
        ):
            max_diff = float(np.max(np.abs(rho_reference[finite] - rho_fast[finite])))
            raise AssertionError(
                f"--fast bootstrap verification failed on {case_label} data "
                f"for model {name!r}: vectorised Spearman diverges from "
                f"baseline.per_target_spearman by up to {max_diff:.3e} "
                f"(tolerance {tol:.1e}). Refusing to run the fast path."
            )


def verify_fast_bootstrap_matches(
    y_true: np.ndarray,
    predictions: dict,
    seed: int,
    tol: float = 1e-9,
) -> None:
    """
    Hard gate on the --fast path (CLAUDE.md invariant 3).

    Checks two cases, for every model about to be bootstrapped:

    1. The unresampled data -- no duplicate rows, so no ties.
    2. One resample drawn from np.random.default_rng(seed) then
       rng.integers(0, n_lines, size=n_lines) -- exactly the first draw
       bootstrap_over_cell_lines_fast will make with this seed. This is the
       case most likely to catch a real divergence: resampling with
       replacement duplicates rows and therefore introduces ties, and
       Spearman-with-ties equals Pearson-of-ranks only when tie-averaged ranks
       are used throughout. The unresampled data alone cannot exercise that.

    Raises AssertionError on any mismatch in either case rather than falling
    back silently: a bootstrap distribution built on a metric that quietly
    diverged from the one every other number in this project uses would be
    worse than not having one.
    """
    _assert_fast_matches_reference(y_true, predictions, "unresampled", tol)

    n_lines = y_true.shape[0]
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n_lines, size=n_lines)
    resampled_predictions = {name: pred[idx] for name, pred in predictions.items()}
    _assert_fast_matches_reference(
        y_true[idx], resampled_predictions, "one seeded resample (with ties)", tol
    )


def bootstrap_over_cell_lines_fast(
    y_true: np.ndarray,
    predictions: dict,
    n_boot: int,
    seed: int,
    progress_every: int = 50,
) -> dict:
    """
    Vectorised counterpart to bootstrap_over_cell_lines.

    Draws the identical resample sequence -- same seed, same rng calls in the
    same order, one draw of n_lines indices per resample, models scored in the
    same order within it -- so its output is bit-for-bit interchangeable with
    the loop version's. The only difference is how each resample is scored:
    fully-observed columns (precomputed once, before the loop, since
    missingness does not change under resampling) go through
    fast_per_target_spearman's vectorised block instead of a per-column loop.

    Caller must have already run verify_fast_bootstrap_matches on these same
    (y_true, predictions) before calling this -- this function does not check.
    """
    rng = np.random.default_rng(seed)
    n_lines = y_true.shape[0]
    names = list(predictions)
    draws = {name: np.full(int(n_boot), np.nan, dtype=float) for name in names}

    fully_observed = {
        name: (
            np.all(np.isfinite(y_true), axis=0)
            & np.all(np.isfinite(predictions[name]), axis=0)
        )
        for name in names
    }

    for b in range(int(n_boot)):
        idx = rng.integers(0, n_lines, size=n_lines)
        y_resampled = y_true[idx]
        for name in names:
            rho = fast_per_target_spearman(
                y_resampled,
                predictions[name][idx],
                fully_observed_mask=fully_observed[name],
            )
            if np.isfinite(rho).any():
                draws[name][b] = float(np.nanmean(rho))
        if progress_every and (b + 1) % progress_every == 0:
            print(f"        {b + 1}/{int(n_boot)} resamples [fast]", flush=True)

    return draws


def percentile_ci(draws: np.ndarray, level: float = 95.0) -> tuple:
    """Percentile interval of a bootstrap distribution, ignoring NaN draws."""
    finite = np.asarray(draws, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return None, None
    tail = (100.0 - float(level)) / 2.0
    return (
        float(np.percentile(finite, tail)),
        float(np.percentile(finite, 100.0 - tail)),
    )


def summarise_draws(draws: np.ndarray, level: float = 95.0) -> dict:
    """Point summary of a bootstrap distribution."""
    finite = np.asarray(draws, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return {"n_draws": 0}
    low, high = percentile_ci(finite, level)
    return {
        "n_draws": int(finite.size),
        "mean": round(float(np.mean(finite)), 4),
        "std_error": round(float(np.std(finite, ddof=1)), 4) if finite.size > 1 else None,
        "ci_level": float(level),
        "ci_low": round(low, 4),
        "ci_high": round(high, 4),
    }


# --------------------------------------------------------------------------
# Paired analysis over targets
# --------------------------------------------------------------------------

def paired_target_analysis(
    y_true: np.ndarray,
    pred_reference: np.ndarray,
    pred_candidate: np.ndarray,
    name_reference: str,
    name_candidate: str,
) -> dict:
    """
    A2 and A3 in one pass: both use the same two per-target correlation vectors.

    Targets are kept only where both models produced a defined correlation, so
    the pairing is exact. A2 is the Wilcoxon signed-rank test on the paired
    differences; A3 is the correlation between the two vectors themselves.
    """
    rho_reference = baseline.per_target_spearman(y_true, pred_reference)
    rho_candidate = baseline.per_target_spearman(y_true, pred_candidate)

    mask = np.isfinite(rho_reference) & np.isfinite(rho_candidate)
    ref = rho_reference[mask]
    cand = rho_candidate[mask]
    diff = cand - ref

    out: dict = {
        "reference_model": name_reference,
        "candidate_model": name_candidate,
        "n_targets_total": int(rho_reference.size),
        "n_targets_paired": int(mask.sum()),
        "n_targets_dropped": int(rho_reference.size - mask.sum()),
    }

    if mask.sum() == 0:
        out["status"] = "no targets scored by both models"
        return out

    out.update({
        "mean_reference_rho": round(float(np.mean(ref)), 4),
        "mean_candidate_rho": round(float(np.mean(cand)), 4),
        "mean_difference": round(float(np.mean(diff)), 4),
        "median_difference": round(float(np.median(diff)), 4),
        "frac_targets_candidate_wins": round(float((diff > 0).mean()), 4),
        "frac_targets_tied": round(float((diff == 0).mean()), 4),
    })

    # ---- A2: Wilcoxon signed-rank on the paired differences ----------------
    try:
        statistic, p_value = stats.wilcoxon(cand, ref)
        out["wilcoxon"] = {
            "statistic": float(statistic),
            "p_value": float(p_value),
            "caveat": (
                "Gene dependencies are correlated -- co-essential genes move "
                "together -- so the effective number of independent targets is "
                "well below n_targets_paired and this p-value is optimistic. "
                "Report it with this sentence attached. The bootstrap over "
                "held-out cell lines is the stronger evidence."
            ),
        }
    except ValueError as exc:
        out["wilcoxon"] = {"status": f"not computed: {exc}"}

    # ---- A3: correlation between the two per-target rho vectors ------------
    pearson_r, pearson_p = stats.pearsonr(ref, cand)
    spearman_r, spearman_p = stats.spearmanr(ref, cand)
    out["per_target_rho_correlation"] = {
        "pearson_r": round(float(pearson_r), 4),
        "pearson_p": float(pearson_p),
        "spearman_r": round(float(spearman_r), 4),
        "spearman_p": float(spearman_p),
        "reading": (
            "High correlation means the two models succeed and fail on the "
            "same targets, which supports reading the embeddings as a lossier "
            "compression of the same signal. Low correlation means they carry "
            "different information rather than less of it."
        ),
    }

    return out


# --------------------------------------------------------------------------
# Verification
# --------------------------------------------------------------------------

def published_spearman(results: dict | None, task: str, model: str):
    """Pull one model's reported spearman_mean out of a results JSON."""
    if not results:
        return None
    return (
        results.get("tasks", {})
        .get(task, {})
        .get("models", {})
        .get(model, {})
        .get("spearman_mean")
    )


def check_against_published(
    y_true: np.ndarray,
    predictions: dict,
    published: dict,
    label: str,
) -> bool:
    """
    Score every prediction matrix and compare to the published number.

    This is the guard on everything downstream. Whether the predictions were
    loaded from disk or refitted here, they must reproduce the values the
    project reports, or the analysis is describing a different model than the
    one on the slides.
    """
    all_ok = True
    for model, preds in predictions.items():
        expected = published.get(model)
        actual = baseline.evaluate(y_true, preds)["spearman_mean"]

        if expected is None:
            print(f"        {label}/{model:<14}: {actual}   "
                  f"[no published value to compare]")
            continue

        ok = actual is not None and abs(float(expected) - float(actual)) < 1e-9
        all_ok = all_ok and ok
        status = "ok" if ok else "MISMATCH"
        print(f"        {label}/{model:<14}: published {expected} "
              f"-> computed {actual}   [{status}]")

    return all_ok


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Bootstrap confidence interval, paired significance test, "
            "per-target correlation, and matched effective degrees of freedom "
            "for the baseline-vs-embeddings comparison."
        )
    )
    parser.add_argument(
        "--split", choices=["val", "test"], default="val",
        help="Which split to analyse. Default 'val'.",
    )
    parser.add_argument(
        "--task", choices=["crispr", "prism"], default="crispr",
        help="Which prediction task to analyse. Default 'crispr'.",
    )
    parser.add_argument(
        "--bootstrap", type=int, default=1000,
        help="Number of bootstrap resamples over held-out cell lines. "
             "Default 1000. Use a smaller number while iterating.",
    )
    parser.add_argument(
        "--seed", type=int, default=config.RANDOM_SEED,
        help=f"Seed for the bootstrap resampler. Default {config.RANDOM_SEED}.",
    )
    parser.add_argument(
        "--processed-dir", default=str(config.PROCESSED_DIR),
        help="Directory holding the processed artifacts and embeddings.",
    )
    parser.add_argument(
        "--force-refit", action="store_true",
        help="Ignore any saved prediction matrices and refit both models here.",
    )
    parser.add_argument(
        "--fast", action="store_true",
        help=(
            "Use the vectorised bootstrap (section 9.2 of CLAUDE.md) instead of "
            "the per-column loop. Verified against baseline.per_target_spearman "
            "on the unresampled data before it is used; aborts rather than "
            "falling back silently if that check fails. Default off -- the "
            "loop is what produced every committed analysis_results.json."
        ),
    )
    args = parser.parse_args()

    out = Path(args.processed_dir)
    task = args.task
    eval_split = args.split

    print("=" * 74)
    print("UNCERTAINTY AND CAPACITY ANALYSIS")
    print("=" * 74)

    if eval_split == "test":
        print("\n  *** Analysing the TEST split. ***")
        print("  Only do this once the model set is frozen and the test run has")
        print("  already happened. If you are still choosing between models,")
        print("  stop and use --split val.\n")

    # ---- load ----------------------------------------------------------
    try:
        expression = io_utils.load_matrix(out / "expression")
        crispr = io_utils.load_matrix(out / "crispr_effect")
        metadata = io_utils.load_table(out / "model_metadata")
        selective = io_utils.load_json(out / "selective_genes.json")
        splits_payload = io_utils.load_json(out / "splits.json")
        embeddings = train_head.load_embeddings(out)
    except FileNotFoundError as exc:
        print(f"\n{exc}")
        print("\nRun build_dataset.py, splits.py and run_geneformer_embeddings.py first.")
        return 1

    metadata.index.name = config.MODEL_ID
    assignment = pd.Series(splits_payload["assignment"], name="split")
    selective_genes = [g for g in selective["genes"] if g in crispr.columns]

    try:
        prism = io_utils.load_matrix(out / "prism_response")
    except FileNotFoundError:
        prism = None

    baseline_results = None
    head_results = None
    try:
        baseline_results = io_utils.load_json(out / "baseline_results.json")
    except FileNotFoundError:
        print("\n  note: baseline_results.json not found -- the refit cannot be "
              "checked against the published baseline number.")
    try:
        head_results = io_utils.load_json(out / "head_results.json")
    except FileNotFoundError:
        print("\n  note: head_results.json not found -- the refit cannot be "
              "checked against the published head number.")

    # ---- assemble both tasks -------------------------------------------
    baseline_prepared = baseline.prepare_task(
        task, expression, crispr, prism, selective_genes, assignment, eval_split
    )
    head_prepared = train_head.prepare_task(
        task, embeddings, crispr, prism, selective_genes, assignment, eval_split
    )

    if baseline_prepared is None or head_prepared is None:
        print(f"\n  Task {task!r} is unavailable for one or both models. "
              f"Nothing to analyse.")
        return 1

    baseline_arrays = assemble(baseline_prepared, metadata, assignment, eval_split)
    head_arrays = assemble(head_prepared, metadata, assignment, eval_split)

    # ---- the alignment guard -------------------------------------------
    # baseline.prepare_task intersects the CRISPR targets with the *expression*
    # index; train_head.prepare_task intersects them with the *embedding* index.
    # Nothing in either pipeline asserts those two intersections are the same
    # set of cell lines. A paired bootstrap and a paired signed-rank test are
    # both meaningless if they are not, and the failure is silent -- the numbers
    # would still compute. So it is checked here, before anything else runs.
    print(f"\n{'=' * 74}")
    print("ALIGNMENT")
    print("=" * 74)

    base_eval = pd.Index([str(i) for i in baseline_arrays["eval_index"]])
    head_eval = pd.Index([str(i) for i in head_arrays["eval_index"]])

    print(f"\n  baseline held-out lines : {len(base_eval)}")
    print(f"  head held-out lines     : {len(head_eval)}")

    if not base_eval.equals(head_eval):
        only_baseline = sorted(set(base_eval) - set(head_eval))
        only_head = sorted(set(head_eval) - set(base_eval))
        print("\n  *** The two models are not evaluated on the same cell lines.")
        print(f"      In the baseline only : {len(only_baseline)} "
              f"{only_baseline[:5]}")
        print(f"      In the head only     : {len(only_head)} {only_head[:5]}")
        print("\n      A paired test across models requires identical held-out")
        print("      rows in identical order. Reconcile the expression and")
        print("      embedding indices before running this analysis.")
        return 1

    if baseline_arrays["target_columns"] != head_arrays["target_columns"]:
        print("\n  *** The two models do not share target columns. "
              "Reconcile before pairing.")
        return 1

    print("  Held-out cell lines and target columns match exactly.")

    n_targets = len(baseline_arrays["target_columns"])
    n_lines = len(base_eval)
    n_train = len(baseline_arrays["train_index"])
    print(f"  targets                 : {n_targets}")
    print(f"  training lines          : {n_train}")

    # ---- A4: effective degrees of freedom, one convention, both models ---
    print(f"\n{'=' * 74}")
    print("A4  EFFECTIVE DEGREES OF FREEDOM")
    print("=" * 74)
    print("\n  Building the design matrices Ridge actually sees...")

    baseline_design = baseline_design_matrix(baseline_arrays["X_train"], config.RANDOM_SEED)
    head_design = head_design_matrix(head_arrays["X_train"])

    baseline_s2 = squared_singular_values(baseline_design)
    head_s2 = squared_singular_values(head_design)

    baseline_alpha = published_alpha(baseline_results, task, "ridge_pca")
    head_alpha = published_alpha(head_results, task, "ridge_head")

    df_block: dict = {
        "convention": (
            "trace of the ridge hat matrix, sum s_i^2 / (s_i^2 + alpha), where "
            "s_i are the singular values of the centred design matrix the Ridge "
            "estimator receives, computed on the training rows only. Not "
            "sklearn's PCA.explained_variance_, which is s_i^2 / (n - 1) and "
            "would rescale alpha by a factor of n - 1."
        ),
        "n_train_rows": int(n_train),
        "models": {},
    }

    for name, squared, alpha, design, alphas in (
        ("ridge_pca", baseline_s2, baseline_alpha, baseline_design, config.RIDGE_ALPHAS),
        ("ridge_head", head_s2, head_alpha, head_design, train_head.HEAD_RIDGE_ALPHAS),
    ):
        entry = {
            "design_shape": [int(design.shape[0]), int(design.shape[1])],
            "n_nonzero_singular_values": int(squared.size),
            "trace_ZtZ": round(float(np.sum(squared)), 2),
            "effective_df_sweep": effective_df_sweep(squared, alphas),
        }
        if alpha is not None:
            entry["selected_alpha"] = float(alpha)
            entry["effective_df_at_selected_alpha"] = round(
                effective_df(squared, alpha), 2
            )
        df_block["models"][name] = entry

        print(f"\n  {name}")
        print(f"    design            : {design.shape[0]} x {design.shape[1]}")
        print(f"    trace(Z'Z)        : {entry['trace_ZtZ']}")
        if alpha is not None:
            print(f"    selected alpha    : {alpha:g}")
            print(f"    effective df      : "
                  f"{entry['effective_df_at_selected_alpha']} "
                  f"of {design.shape[1]}")
        else:
            print("    selected alpha    : unknown (results JSON not found)")

    both = df_block["models"]
    if ("effective_df_at_selected_alpha" in both["ridge_pca"]
            and "effective_df_at_selected_alpha" in both["ridge_head"]):
        print(f"\n  Capacity on one axis: "
              f"{both['ridge_pca']['effective_df_at_selected_alpha']} "
              f"(baseline, {both['ridge_pca']['design_shape'][1]} dims) vs "
              f"{both['ridge_head']['effective_df_at_selected_alpha']} "
              f"(head, {both['ridge_head']['design_shape'][1]} dims)")

    del baseline_design, head_design

    # ---- predictions: load if saved, otherwise refit --------------------
    print(f"\n{'=' * 74}")
    print("PREDICTIONS")
    print("=" * 74)

    baseline_published = {
        "ridge_pca": published_spearman(baseline_results, task, "ridge_pca"),
        "lineage_mean": published_spearman(baseline_results, task, "lineage_mean"),
        "global_mean": published_spearman(baseline_results, task, "global_mean"),
    }
    head_published = {
        "ridge_head": published_spearman(head_results, task, "ridge_head"),
        "mlp_head": published_spearman(head_results, task, "mlp_head"),
    }

    baseline_stems = {} if args.force_refit else discover_saved_predictions(
        out, "baseline", task, eval_split
    )
    head_stems = {} if args.force_refit else discover_saved_predictions(
        out, "head", task, eval_split
    )

    prediction_source = "saved matrices"

    if baseline_stems and head_stems:
        print(f"\n  Loading saved matrices from {out / 'predictions'}")
        y_true, baseline_preds, saved_index, saved_columns = load_saved_predictions(
            baseline_stems
        )
        y_true_head, head_preds, head_index, head_columns = load_saved_predictions(
            head_stems
        )

        if not saved_index.astype(str).equals(head_index.astype(str)):
            print("\n  *** The saved baseline and head bundles cover different "
                  "cell lines. Re-save both from the same split.")
            return 1
        if not np.allclose(y_true, y_true_head, equal_nan=True):
            print("\n  *** The saved baseline and head truth matrices differ. "
                  "Re-save both from the same run.")
            return 1
        if not pd.Index([str(i) for i in saved_index]).equals(base_eval):
            print("\n  *** The saved matrices do not match the cell lines this "
                  "run assembled. Re-save, or pass --force-refit.")
            return 1
        if saved_columns != baseline_arrays["target_columns"]:
            print("\n  *** The saved matrices do not match the target columns "
                  "this run assembled. Re-save, or pass --force-refit.")
            return 1
    else:
        prediction_source = "refitted in analysis.py"
        if args.force_refit:
            print("\n  --force-refit given; refitting both models here.")
        else:
            print("\n  No saved prediction matrices found; refitting both models here.")
            print("  Run baseline.py and train_head.py with --save-predictions to")
            print("  skip this step in future.")

        y_true = baseline_arrays["Y_eval_raw"]
        baseline_preds = {}
        head_preds = {}

        print("\n  Fitting ridge on PCA of expression...")
        preds, _info = baseline.run_ridge_pca(
            baseline_arrays["X_train"],
            baseline_arrays["Y_train_filled"],
            baseline_arrays["X_eval"],
            baseline_arrays["X_val"],
            baseline_arrays["Y_val"],
            n_components=config.PCA_COMPONENTS,
            alphas=config.RIDGE_ALPHAS,
            seed=config.RANDOM_SEED,
            train_groups=baseline_arrays["train_groups"],
        )
        baseline_preds["ridge_pca"] = preds

        lineage_col = config.STRATIFY_COL
        if lineage_col in metadata.columns:
            print("  Fitting the lineage-mean control...")
            baseline_preds["lineage_mean"] = baseline.run_lineage_mean(
                baseline_arrays["Y_train_raw"],
                metadata.loc[baseline_arrays["train_index"], lineage_col]
                .astype(str).to_numpy(),
                metadata.loc[baseline_arrays["eval_index"], lineage_col]
                .astype(str).to_numpy(),
            )

        print("  Fitting ridge on Geneformer embeddings...")
        preds, _info = train_head.run_ridge_head(
            head_arrays["X_train"],
            head_arrays["Y_train_filled"],
            head_arrays["X_eval"],
            head_arrays["X_val"],
            head_arrays["Y_val"],
            alphas=train_head.HEAD_RIDGE_ALPHAS,
            seed=config.RANDOM_SEED,
            train_groups=head_arrays["train_groups"],
        )
        head_preds["ridge_head"] = preds

    # Expression and embeddings are no longer needed; the analysis runs on the
    # held-out matrices alone.
    del expression, embeddings, baseline_arrays, head_arrays

    # ---- verify against the published numbers ---------------------------
    print("\n  Checking every prediction matrix against the published values:")
    ok_baseline = check_against_published(
        y_true, baseline_preds, baseline_published, "baseline"
    )
    ok_head = check_against_published(
        y_true, head_preds, head_published, "head"
    )

    if not (ok_baseline and ok_head):
        print("\n  *** At least one model does not reproduce its published "
              "Spearman mean.")
        print("      Stopping. A confidence interval around a number that is "
              "not the reported")
        print("      number would be worse than no interval at all.")
        return 1
    print("\n  Every model reproduces its published value.")

    # ---- A2 and A3 ------------------------------------------------------
    print(f"\n{'=' * 74}")
    print("A2  PAIRED TEST OVER TARGETS   /   A3  PER-TARGET CORRELATION")
    print("=" * 74)

    paired = paired_target_analysis(
        y_true,
        baseline_preds["ridge_pca"],
        head_preds["ridge_head"],
        "ridge_pca",
        "ridge_head",
    )

    print(f"\n  targets paired            : {paired['n_targets_paired']} "
          f"of {paired['n_targets_total']}")
    print(f"  mean rho, baseline        : {paired.get('mean_reference_rho')}")
    print(f"  mean rho, head            : {paired.get('mean_candidate_rho')}")
    print(f"  mean difference           : {paired.get('mean_difference')}")
    print(f"  median difference         : {paired.get('median_difference')}")
    print(f"  targets where head wins   : "
          f"{paired.get('frac_targets_candidate_wins')}")
    wilcoxon_block = paired.get("wilcoxon", {})
    if "p_value" in wilcoxon_block:
        print(f"  Wilcoxon signed-rank p    : {wilcoxon_block['p_value']:.3e}")
        print("    (optimistic -- gene dependencies are correlated, so the")
        print("     effective number of independent targets is far below "
              f"{paired['n_targets_paired']})")
    correlation_block = paired.get("per_target_rho_correlation", {})
    if correlation_block:
        print(f"  per-target rho, Pearson r : {correlation_block['pearson_r']}")
        print(f"  per-target rho, Spearman  : {correlation_block['spearman_r']}")

    # ---- A1 -------------------------------------------------------------
    print(f"\n{'=' * 74}")
    print("A1  BOOTSTRAP OVER HELD-OUT CELL LINES")
    print("=" * 74)

    bootstrap_inputs = {"ridge_pca": baseline_preds["ridge_pca"],
                        "ridge_head": head_preds["ridge_head"]}
    if "lineage_mean" in baseline_preds:
        bootstrap_inputs["lineage_mean"] = baseline_preds["lineage_mean"]

    print(f"\n  {args.bootstrap} resamples of {n_lines} cell lines, "
          f"{len(bootstrap_inputs)} models, seed {args.seed}")
    print("  Resampling with replacement duplicates cell lines, which puts ties")
    print("  into the rank correlation. spearmanr averages tied ranks, so this")
    print("  is correct -- but state it in the methods section.\n")

    cols_with_nan = np.any(~np.isfinite(y_true), axis=0)
    print(f"  val truth matrix    : {int(cols_with_nan.sum())} of {y_true.shape[1]} "
          f"target columns contain a NaN ({int(np.sum(~np.isfinite(y_true)))} "
          f"NaN cells total)")
    if cols_with_nan.any():
        print("    fast_per_target_spearman's per-column fallback branch runs "
              "on real data for these columns.")
    else:
        print("    every target column is fully observed here, so "
              "fast_per_target_spearman's fallback branch never runs on real "
              "data -- see self_test_fast_fallback_branch for synthetic "
              "coverage of it.")

    if args.fast:
        print("\n  --fast given. Running the synthetic fallback-branch self-test...")
        try:
            self_test_fast_fallback_branch()
        except AssertionError as exc:
            print(f"\n  *** {exc}")
            print("\n      Stopping.")
            return 1
        print("  Self-test passed.")

        print("\n  Verifying the vectorised path against "
              "baseline.per_target_spearman,")
        print("  on the unresampled data and on one seeded resample (with ties)...")
        try:
            verify_fast_bootstrap_matches(y_true, bootstrap_inputs, seed=args.seed)
        except AssertionError as exc:
            print(f"\n  *** {exc}")
            print("\n      Stopping. Not falling back to the loop silently -- ")
            print("      re-run without --fast, or fix the mismatch.")
            return 1
        print("  Verified: exact match in both cases. Running the vectorised "
              "bootstrap.\n")
        bootstrap_start = time.perf_counter()
        draws = bootstrap_over_cell_lines_fast(
            y_true, bootstrap_inputs, args.bootstrap, args.seed
        )
        bootstrap_elapsed = time.perf_counter() - bootstrap_start
    else:
        bootstrap_start = time.perf_counter()
        draws = bootstrap_over_cell_lines(
            y_true, bootstrap_inputs, args.bootstrap, args.seed
        )
        bootstrap_elapsed = time.perf_counter() - bootstrap_start

    delta_draws = draws["ridge_head"] - draws["ridge_pca"]
    delta_summary = summarise_draws(delta_draws)

    bootstrap_block: dict = {
        "n_resamples": int(args.bootstrap),
        "seed": int(args.seed),
        "resampling_unit": "held-out cell line",
        "n_cell_lines": int(n_lines),
        "tie_note": (
            "Resampling with replacement duplicates cell lines and therefore "
            "introduces ties into the per-target rank correlation. scipy's "
            "spearmanr handles ties by averaging ranks, so the computation is "
            "correct, but the tie behaviour is a property of the method and "
            "must be disclosed in the methods section."
        ),
        "per_model": {
            name: summarise_draws(values) for name, values in draws.items()
        },
        "delta_head_minus_baseline": delta_summary,
    }

    print(f"\n  ridge_pca   mean rho : {bootstrap_block['per_model']['ridge_pca']['mean']}"
          f"  95% CI [{bootstrap_block['per_model']['ridge_pca']['ci_low']}, "
          f"{bootstrap_block['per_model']['ridge_pca']['ci_high']}]")
    print(f"  ridge_head  mean rho : {bootstrap_block['per_model']['ridge_head']['mean']}"
          f"  95% CI [{bootstrap_block['per_model']['ridge_head']['ci_low']}, "
          f"{bootstrap_block['per_model']['ridge_head']['ci_high']}]")

    # ---- share of the expression-over-lineage gain ----------------------
    if "lineage_mean" in draws:
        baseline_gain = draws["ridge_pca"] - draws["lineage_mean"]
        head_gain = draws["ridge_head"] - draws["lineage_mean"]
        with np.errstate(divide="ignore", invalid="ignore"):
            share = np.where(
                np.abs(baseline_gain) > 1e-12, head_gain / baseline_gain, np.nan
            )
        share_summary = summarise_draws(share)
        bootstrap_block["share_of_expression_gain"] = {
            **share_summary,
            "definition": (
                "(head rho - lineage rho) / (baseline rho - lineage rho), "
                "recomputed inside every bootstrap resample. The fraction of "
                "the gain that raw expression buys over tissue identity alone "
                "which the embeddings recover. More interpretable than a raw "
                "difference of correlations, and unlike that difference it is "
                "expressed on a scale a non-specialist audience can read."
            ),
            "caveat": (
                "lineage_mean is not a zero-information floor, so this is a "
                "share of a gain over tissue identity, not a share of all "
                "available signal. It also inherits the noise of three "
                "estimates rather than two, which is why it is reported with "
                "its own interval."
            ),
        }
        print(f"  lineage_mean mean rho: "
              f"{bootstrap_block['per_model']['lineage_mean']['mean']}")
        print(f"\n  Share of the expression-over-lineage gain recovered by the")
        print(f"  embeddings           : {share_summary['mean']}  "
              f"95% CI [{share_summary['ci_low']}, {share_summary['ci_high']}]")

    print(f"\n  {'-' * 70}")
    print(f"  DELTA (head - baseline) : {delta_summary['mean']}  "
          f"95% CI [{delta_summary['ci_low']}, {delta_summary['ci_high']}]")
    print(f"  bootstrap standard error: {delta_summary['std_error']}")
    crosses_zero = (
        delta_summary["ci_low"] is not None
        and delta_summary["ci_high"] is not None
        and delta_summary["ci_low"] <= 0.0 <= delta_summary["ci_high"]
    )
    if crosses_zero:
        print("\n  The interval contains zero: on this evidence the difference")
        print("  is not distinguishable from no difference at the 95% level.")
        print("  That is a legitimate result and a more honest headline than an")
        print("  unqualified point estimate. State the interval, not just -0.03.")
    else:
        print("\n  The interval excludes zero: the deficit survives resampling")
        print("  of the held-out cell lines at the 95% level.")
    print(f"  {'-' * 70}")
    print(f"  bootstrap elapsed       : {bootstrap_elapsed:.2f}s "
          f"({'vectorised --fast' if args.fast else 'per-column loop'} path, "
          f"{args.bootstrap} resamples)")
    print("  Not written to analysis_results.json -- wall-clock time is not a")
    print("  reproducible artifact and would break the byte-identical guarantee.")

    # ---- write ----------------------------------------------------------
    payload = {
        "eval_split": eval_split,
        "task": task,
        "seed": int(config.RANDOM_SEED),
        "bootstrap_seed": int(args.seed),
        "prediction_source": prediction_source,
        "n_cell_lines": int(n_lines),
        "n_targets": int(n_targets),
        "n_train_lines": int(n_train),
        "models_analysed": sorted(set(baseline_preds) | set(head_preds)),
        "A1_bootstrap": bootstrap_block,
        "A2_A3_paired_over_targets": paired,
        "A4_effective_degrees_of_freedom": df_block,
    }
    path = io_utils.save_json(payload, out / "analysis_results.json")

    print(f"\n  Written to: {path}")
    print("\n  Report the delta as an interval, not a point estimate, and carry")
    print("  the two caveats with their numbers: the Wilcoxon p-value is")
    print("  optimistic because targets are correlated, and the bootstrap")
    print("  introduces ties by duplicating cell lines.")
    return 0


def published_alpha(results: dict | None, task: str, model: str):
    """Pull one model's selected alpha out of a results JSON."""
    if not results:
        return None
    return (
        results.get("tasks", {})
        .get(task, {})
        .get("models", {})
        .get(model, {})
        .get("alpha")
    )


if __name__ == "__main__":
    sys.exit(main())
