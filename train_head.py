"""
train_head.py
=============
Trains prediction heads on the Geneformer cell embeddings and reports the
result as a delta against the ridge-on-expression baseline, on the *same*
train/val/test split, with the *same* per-target Spearman metric.

    python train_head.py                 # evaluate on the validation split
    python train_head.py --split test    # final evaluation; use ONCE, at the end
    python train_head.py --no-mlp        # ridge head only (fast)

Output
------
    head_results.json   metrics per head, per task, plus explicit deltas
                        against baseline_results.json

Where this sits in the pipeline
-------------------------------
    build_dataset.py            -> expression, crispr_effect, metadata, selective_genes
    splits.py                   -> splits.json         (frozen, leakage-checked)
    baseline.py                 -> baseline_results.json  ("the number to beat")
    prepare_geneformer_input.py -> geneformer_input.h5ad
    run_geneformer_embeddings.py-> geneformer_embeddings.csv   (ModelID x dims)
    train_head.py  (this file)  -> head_results.json

Why two heads, not one
----------------------
This module runs a ladder, exactly as baseline.py does (global mean ->
lineage mean -> ridge). Here the rungs are:

1. **Ridge head.** A linear model on the embeddings. This is the control that
   matters. baseline.py already established what a linear model achieves on raw
   *expression* (ridge_pca). Running a linear model on the *embeddings* isolates
   one thing and one thing only: does the Geneformer representation carry more
   linearly-accessible signal than raw expression? Same estimator family, same
   metric, same split -- only the features change.

2. **MLP head.** A small learned neural head on the embeddings. Reported as a
   further delta over the ridge head. If the MLP beats ridge-on-expression but
   not ridge-on-embeddings, the gain is the embedding, not the nonlinearity --
   and saying so is the honest finding. Presenting the MLP without the linear
   control on the same features would make that attribution impossible.

The discipline this module inherits from baseline.py
----------------------------------------------------
* The metric is imported from baseline.py, not re-implemented, so "0.24 here"
  and "0.2313 there" are the same measurement.
* The regularisation strength is never chosen on the split being reported.
  On --split val it is chosen by patient-grouped cross-validation *inside the
  training set*. On --split test it is chosen on the held-out val split. This
  is the same rule run_ridge_pca follows.
* The test split is read only when you ask for it, and reading it is announced.

HONESTY NOTE
------------
The heads and the whole evaluation path here were written and tested against a
correctly-formatted *synthetic* embeddings file, because real Geneformer
embeddings require the GPU step in run_geneformer_embeddings.py. The code paths
(loading, alignment, split handling, metric, hyperparameter selection, delta
reporting) are exercised. The *numbers* only mean something once
geneformer_embeddings.csv is the real thing.
"""

from __future__ import annotations

import sys
import argparse
import warnings

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GroupKFold

import config
import io_utils

# The metric and the imputation discipline come straight from baseline.py so
# the comparison is literally the same measurement, not a re-implementation of
# it. _select_alpha_inner_cv is reused for the ridge head for the same reason.
from baseline import (
    per_target_spearman,
    evaluate,
    impute_with_train_mean,
    _select_alpha_inner_cv,
    save_prediction_bundle,
    verify_prediction_bundle,
)


# --------------------------------------------------------------------------
# MLP configuration
# --------------------------------------------------------------------------
# A modest, well-regularised head. The point is not to win an architecture
# search; it is to give the embeddings a fair non-linear model while keeping
# the overfitting surface small on ~1,000 cell lines. The L2 penalty (alpha) is
# selected with the same split discipline as the ridge penalty.

MLP_HIDDEN_DEFAULT = 256
MLP_ALPHAS = [1e-3, 1e-2, 1e-1]
MLP_MAX_ITER = 300
MLP_N_ITER_NO_CHANGE = 10


# --------------------------------------------------------------------------
# Ridge head alpha grid
# --------------------------------------------------------------------------
# The baseline regularises 200 PCA components of expression; the head
# regularises 768 raw embedding dimensions with no PCA. The same alpha does
# not mean the same shrinkage on two different spectra, so the head gets its
# own grid.
#
# Effective degrees of freedom below are the trace of the ridge hat matrix,
#
#     df(alpha) = sum_i  s_i^2 / (s_i^2 + alpha)
#
# where s_i are the singular values of the centred design matrix Ridge
# actually receives: the standardised embeddings, on the 800 training rows
# only. These are computed by analysis.py, not by hand. Regenerate with
# `py analysis.py --split val` and read A4_effective_degrees_of_freedom in
# analysis_results.json.
#
#     alpha=1     -> 611.44     alpha=3.16e3 ->  51.78
#     alpha=10    -> 450.46     alpha=1e4    ->  25.89
#     alpha=1e2   -> 251.77     alpha=1e5    ->   4.87
#     alpha=1e3   ->  95.87     alpha=1e6    ->   0.60
#
# trace(Z'Z) = 614400.0, which is exactly 800 rows x 768 unit-variance
# columns. That identity is the cheapest available check that a df figure was
# measured on the right matrix: at large alpha the sum converges to
# trace/alpha, so alpha=1e6 cannot exceed 0.614 on this design.
#
# A hand-computed table previously in this file reported 676.2 / 513.0 /
# 297.0 / 117.1 / 64.3 / 32.7 / 6.5 / 0.8 across the same grid. Every entry
# was high: by 11% at alpha=1, rising to 33% at alpha=1e6. Its large-alpha
# limit implies trace(Z'Z) of roughly 819,000 rather than 614,400, so it was
# measured on some other matrix -- more rows than the 800 the model was
# fitted on, or an unstandardised design, or both. Which of those is not
# recoverable. That it was not the fitted design is certain. The 64.2 / 64.3
# figure quoted throughout the project documents comes from that table and is
# superseded by 51.78.
#
# The grid does not extend below alpha=1 because 611 of 768 dimensions are
# already live there. The ceiling is 1e6 because the model is shrunk to 0.60
# effective degrees of freedom, below a single parameter, so nothing above it
# can be optimal.
#
# Where the two models land is the reason for measuring this at all. Inner CV
# chose alpha=1e5 for the baseline and alpha=3162 for the head -- a factor of
# 32 apart, on feature spaces of 200 and 768 dimensions -- and the two arrive
# at 49.71 and 51.78 effective parameters, within 4% of each other. Nothing
# was tuned to produce that. It means the capacity available to the two
# models was matched by measurement rather than assumed, and it suggests that
# roughly 50 effective parameters is what 800 training lines support on this
# task, largely independent of the representation.

HEAD_RIDGE_ALPHAS = [
    1.0, 3.16, 10.0, 31.6, 100.0, 316.0,
    1000.0, 3162.0, 10000.0, 31623.0,
    100000.0, 316228.0, 1000000.0,
]
MLP_VALIDATION_FRACTION = 0.1   # internal early-stopping split, taken from train only


def make_mlp(alpha: float, hidden: int, seed: int) -> MLPRegressor:
    """Construct the MLP head. early_stopping carves its validation set out of
    the training data only, so it never sees val or test."""
    return MLPRegressor(
        hidden_layer_sizes=(hidden,),
        activation="relu",
        solver="adam",
        alpha=alpha,
        max_iter=MLP_MAX_ITER,
        early_stopping=True,
        n_iter_no_change=MLP_N_ITER_NO_CHANGE,
        validation_fraction=MLP_VALIDATION_FRACTION,
        random_state=seed,
    )


# --------------------------------------------------------------------------
# Embedding loading
# --------------------------------------------------------------------------

def load_embeddings(processed_dir, filename_stem: str = "geneformer_embeddings") -> pd.DataFrame:
    """
    Load the Geneformer embeddings written by run_geneformer_embeddings.py.

    That script writes a CSV (and optionally a parquet) indexed by ModelID with
    one numeric column per embedding dimension. It does *not* go through
    io_utils.save_matrix, so this loader reads the file directly rather than
    assuming io_utils' matrix layout.
    """
    parquet = processed_dir / f"{filename_stem}.parquet"
    csv = processed_dir / f"{filename_stem}.csv"

    if parquet.is_file():
        frame = pd.read_parquet(parquet)
    elif csv.is_file():
        frame = pd.read_csv(csv, index_col=0)
    else:
        raise FileNotFoundError(
            f"\nNo embeddings found in {processed_dir}.\n"
            f"  Looked for : {parquet.name}, {csv.name}\n"
            f"  Fix        : run run_geneformer_embeddings.py (on a GPU box) "
            f"first; it writes geneformer_embeddings.csv here.\n"
        )

    frame.index.name = config.MODEL_ID

    # Keep only numeric embedding dimensions; guard against a stray label column.
    numeric = frame.select_dtypes(include=[np.number])
    if numeric.shape[1] == 0:
        raise ValueError(
            "Embeddings file has no numeric columns. Expected one column per "
            "embedding dimension."
        )
    if numeric.shape[1] != frame.shape[1]:
        dropped = [c for c in frame.columns if c not in numeric.columns]
        print(f"  note: dropped {len(dropped)} non-numeric column(s) from "
              f"embeddings: {dropped[:5]}")
    return numeric


# --------------------------------------------------------------------------
# Task assembly (parallel to baseline.prepare_task, but embeddings -> targets)
# --------------------------------------------------------------------------

def prepare_task(
    task: str,
    embeddings: pd.DataFrame,
    crispr: pd.DataFrame,
    prism: pd.DataFrame | None,
    selective_genes: list[str],
    assignment: pd.Series,
    eval_split: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Index, pd.Index] | None:
    """
    Assemble (X=embeddings, Y=targets) for one task on the cell lines that have
    an embedding, a target row, and a split assignment.

    Returns (X, Y, train_index, eval_index) or None if the task is unavailable.
    """
    if task == "crispr":
        targets = crispr[selective_genes] if selective_genes else None
        if targets is None or targets.shape[1] == 0:
            return None
    elif task == "prism":
        if prism is None or prism.shape[0] == 0:
            return None
        targets = prism
    else:
        raise ValueError(f"Unknown task {task!r}")

    lines = embeddings.index.intersection(targets.index)
    lines = pd.Index(sorted(lines), name=config.MODEL_ID)
    if len(lines) == 0:
        return None

    X = embeddings.loc[lines]
    Y = targets.loc[lines]

    split_here = assignment.reindex(lines)
    train_index = lines[(split_here == "train").to_numpy()]
    eval_index = lines[(split_here == eval_split).to_numpy()]

    if len(train_index) < 10 or len(eval_index) < 5:
        return None

    return X, Y, train_index, eval_index


# --------------------------------------------------------------------------
# Heads
# --------------------------------------------------------------------------

def run_ridge_head(
    X_train: np.ndarray,
    Y_train: np.ndarray,
    X_eval: np.ndarray,
    X_val: np.ndarray | None,
    Y_val: np.ndarray | None,
    alphas: list[float],
    seed: int,
    train_groups: np.ndarray | None,
) -> tuple[np.ndarray, dict]:
    """
    StandardScaler -> Ridge on the embeddings.

    No PCA: the embedding is already a dense, low-dimensional learned
    representation, so projecting it further would discard signal rather than
    denoise features (which is what PCA does for the 18k-wide raw expression
    matrix in baseline.py). Everything else -- the alpha grid, the selection
    rule, the metric -- matches the baseline ridge.
    """
    scaler = StandardScaler()
    Z_train = scaler.fit_transform(X_train)
    Z_eval = scaler.transform(X_eval)

    info: dict = {"n_features": int(Z_train.shape[1])}

    if X_val is not None and Y_val is not None and len(X_val) >= 5:
        Z_val = scaler.transform(X_val)
        best_score = -np.inf
        chosen = alphas[len(alphas) // 2]
        sweep = {}
        for alpha in alphas:
            model = Ridge(alpha=alpha)
            model.fit(Z_train, Y_train)
            rho = per_target_spearman(Y_val, model.predict(Z_val))
            score = float(np.nanmean(rho)) if np.isfinite(rho).any() else -np.inf
            sweep[str(alpha)] = round(score, 4) if np.isfinite(score) else None
            if score > best_score:
                best_score, chosen = score, alpha
        info["alpha_sweep"] = sweep
        info["alpha_selected_on"] = "held-out validation split"
    else:
        chosen, cv_info = _select_alpha_inner_cv(
            Z_train, Y_train, train_groups, alphas, seed
        )
        info["alpha_sweep"] = cv_info.get("sweep")
        info["alpha_selected_on"] = (
            f"grouped {cv_info.get('n_folds', '?')}-fold CV inside the training set"
            if cv_info.get("status") == "ok"
            else f"default ({cv_info.get('status')})"
        )

    info["alpha"] = chosen
    if len(alphas) > 1 and chosen in (min(alphas), max(alphas)):
        info["alpha_at_grid_boundary"] = True
        info["alpha_warning"] = (
            f"alpha={chosen:g} is at the edge of the grid "
            f"[{min(alphas):g}, {max(alphas):g}]; widen MLP/ridge alphas and re-run."
        )
    else:
        info["alpha_at_grid_boundary"] = False

    model = Ridge(alpha=chosen)
    model.fit(Z_train, Y_train)
    return model.predict(Z_eval), info


def _select_mlp_alpha_inner_cv(
    X_train: np.ndarray,
    Y_train: np.ndarray,
    groups: np.ndarray | None,
    alphas: list[float],
    hidden: int,
    seed: int,
    n_folds: int = 5,
) -> tuple[float, dict]:
    """
    Choose the MLP L2 penalty by patient-grouped CV inside the training set.

    Directly parallel to baseline._select_alpha_inner_cv, but fits the MLP head
    instead of Ridge and standardises within each fold so the scaler never sees
    the inner-test fold.
    """
    n_samples = X_train.shape[0]
    if groups is None:
        groups = np.arange(n_samples)

    unique_groups = np.unique(groups)
    n_folds = int(min(n_folds, len(unique_groups)))
    if n_folds < 2:
        return alphas[len(alphas) // 2], {"status": "too few groups for inner CV"}

    splitter = GroupKFold(n_splits=n_folds)
    scores: dict[float, list[float]] = {alpha: [] for alpha in alphas}

    for inner_train, inner_test in splitter.split(X_train, groups=groups):
        if len(inner_test) < 3:
            continue
        scaler = StandardScaler()
        Z_tr = scaler.fit_transform(X_train[inner_train])
        Z_te = scaler.transform(X_train[inner_test])
        for alpha in alphas:
            model = make_mlp(alpha, hidden, seed)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")  # convergence chatter in inner CV
                model.fit(Z_tr, Y_train[inner_train])
            rho = per_target_spearman(Y_train[inner_test], model.predict(Z_te))
            if np.isfinite(rho).any():
                scores[alpha].append(float(np.nanmean(rho)))

    mean_scores = {
        alpha: (float(np.mean(vals)) if vals else -np.inf)
        for alpha, vals in scores.items()
    }
    best = max(mean_scores, key=lambda a: mean_scores[a])
    sweep = {
        str(alpha): (round(s, 4) if np.isfinite(s) else None)
        for alpha, s in mean_scores.items()
    }
    return best, {"status": "ok", "n_folds": n_folds, "sweep": sweep}


def run_mlp_head(
    X_train: np.ndarray,
    Y_train: np.ndarray,
    X_eval: np.ndarray,
    X_val: np.ndarray | None,
    Y_val: np.ndarray | None,
    alphas: list[float],
    hidden: int,
    seed: int,
    train_groups: np.ndarray | None,
) -> tuple[np.ndarray, dict]:
    """StandardScaler -> MLPRegressor on the embeddings, same selection rule."""
    scaler = StandardScaler()
    Z_train = scaler.fit_transform(X_train)
    Z_eval = scaler.transform(X_eval)

    info: dict = {"n_features": int(Z_train.shape[1]), "hidden_layer_sizes": [hidden]}

    if X_val is not None and Y_val is not None and len(X_val) >= 5:
        Z_val = scaler.transform(X_val)
        best_score = -np.inf
        chosen = alphas[len(alphas) // 2]
        sweep = {}
        for alpha in alphas:
            model = make_mlp(alpha, hidden, seed)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model.fit(Z_train, Y_train)
            rho = per_target_spearman(Y_val, model.predict(Z_val))
            score = float(np.nanmean(rho)) if np.isfinite(rho).any() else -np.inf
            sweep[str(alpha)] = round(score, 4) if np.isfinite(score) else None
            if score > best_score:
                best_score, chosen = score, alpha
        info["alpha_sweep"] = sweep
        info["alpha_selected_on"] = "held-out validation split"
    else:
        chosen, cv_info = _select_mlp_alpha_inner_cv(
            X_train, Y_train, train_groups, alphas, hidden, seed
        )
        info["alpha_sweep"] = cv_info.get("sweep")
        info["alpha_selected_on"] = (
            f"grouped {cv_info.get('n_folds', '?')}-fold CV inside the training set"
            if cv_info.get("status") == "ok"
            else f"default ({cv_info.get('status')})"
        )

    info["alpha"] = chosen
    info["alpha_at_grid_boundary"] = bool(
        len(alphas) > 1 and chosen in (min(alphas), max(alphas))
    )

    model = make_mlp(chosen, hidden, seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model.fit(Z_train, Y_train)
    return model.predict(Z_eval), info


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------

def run_task(
    task: str,
    embeddings: pd.DataFrame,
    crispr: pd.DataFrame,
    prism: pd.DataFrame | None,
    selective_genes: list[str],
    assignment: pd.Series,
    metadata: pd.DataFrame,
    eval_split: str,
    hidden: int,
    run_mlp: bool,
    baseline_ref: dict | None,
    predictions_out: dict | None = None,
) -> dict | None:
    """
    Run the head ladder for one prediction task and compare to baseline.

    If `predictions_out` is a dict it is filled in place with each head's
    held-out predictions and the matching truth matrix, in the same shape
    baseline.run_task produces, so baseline.save_prediction_bundle can write
    both. Passing None -- the default -- leaves behaviour unchanged.
    """
    prepared = prepare_task(
        task, embeddings, crispr, prism, selective_genes, assignment, eval_split
    )
    if prepared is None:
        print(f"\n  Task {task!r}: insufficient overlap between embeddings, "
              f"targets and splits -- skipping.")
        return None

    X, Y, train_index, eval_index = prepared

    # ---- alignment / integrity guards (checks.py is not importable here, so
    # ---- the guards train_head must satisfy are asserted inline) -----------
    assert X.index.is_unique, "Embedding index has duplicate ModelIDs"
    assert Y.index.equals(X.index), "Embeddings and targets are misaligned"
    assert set(train_index).isdisjoint(set(eval_index)), (
        "train and eval indices overlap -- leakage"
    )

    print(f"\n{'=' * 74}")
    print(f"TASK: {task}")
    print("=" * 74)
    print(f"  embedding dims     : {X.shape[1]}")
    print(f"  targets            : {Y.shape[1]}")
    print(f"  training lines     : {len(train_index)}")
    print(f"  {eval_split + ' lines':<19}: {len(eval_index)}")

    X_train_raw = X.loc[train_index].to_numpy(dtype=float)
    X_eval_raw = X.loc[eval_index].to_numpy(dtype=float)
    Y_train_raw = Y.loc[train_index].to_numpy(dtype=float)
    Y_eval_raw = Y.loc[eval_index].to_numpy(dtype=float)

    # Impute exactly as baseline does: fill values from the training set only.
    # Embeddings should have no NaNs, but the same defensive fill is applied so
    # a stray missing value cannot crash the head or leak eval statistics.
    X_train, X_eval = impute_with_train_mean(X_train_raw, X_eval_raw)
    Y_train_filled, = impute_with_train_mean(Y_train_raw)

    # Validation data for hyperparameter selection when reporting on test.
    X_val = Y_val = None
    if eval_split == "test":
        val_index = X.index[(assignment.reindex(X.index) == "val").to_numpy()]
        if len(val_index) >= 5:
            X_val_raw = X.loc[val_index].to_numpy(dtype=float)
            X_val, = impute_with_train_mean(X_train_raw, X_val_raw)[1:]
            Y_val = Y.loc[val_index].to_numpy(dtype=float)

    train_groups = None
    if config.GROUP_COL in metadata.columns:
        train_groups = (
            metadata.loc[train_index, config.GROUP_COL]
            .astype(str)
            .fillna(pd.Series(train_index, index=train_index).astype(str))
            .to_numpy()
        )

    results: dict = {
        "task": task,
        "eval_split": eval_split,
        "embedding_dims": int(X.shape[1]),
        "n_targets": int(Y.shape[1]),
        "n_train_lines": int(len(train_index)),
        "n_eval_lines": int(len(eval_index)),
        "models": {},
    }

    # Truth matrix and labels recorded once; each head appends its own array.
    # Y_eval_raw is stored unfilled, NaNs intact, because per_target_spearman
    # masks non-finite entries per column and imputing them here would silently
    # change the metric for anything computed downstream.
    if predictions_out is not None:
        predictions_out[task] = {
            "eval_index": [str(i) for i in eval_index],
            "target_columns": [str(c) for c in Y.columns],
            "y_true": Y_eval_raw,
            "models": {},
        }

    # ---- ridge head (linear control on embeddings) ---------------------
    print("\n  [1/2] ridge head (linear control on embeddings)")
    preds, info = run_ridge_head(
        X_train, Y_train_filled, X_eval, X_val, Y_val,
        alphas=HEAD_RIDGE_ALPHAS, seed=config.RANDOM_SEED,
        train_groups=train_groups,
    )
    metrics = evaluate(Y_eval_raw, preds)
    if predictions_out is not None:
        predictions_out[task]["models"]["ridge_head"] = preds
    results["models"]["ridge_head"] = {
        "description": "StandardScaler -> Ridge on Geneformer embeddings.",
        **info, **metrics,
    }
    print(f"        alpha            : {info['alpha']} ({info['alpha_selected_on']})")
    if info.get("alpha_at_grid_boundary"):
        print(f"        WARNING          : {info.get('alpha_warning', 'alpha at grid edge')}")
    print(f"        spearman mean    : {metrics['spearman_mean']}")
    print(f"        spearman median  : {metrics['spearman_median']}")
    print(f"        targets positive : {metrics['spearman_frac_positive']}")
    print(f"        r2 mean          : {metrics['r2_mean']}")

    # ---- MLP head (learned nonlinear head) -----------------------------
    if run_mlp:
        print("\n  [2/2] MLP head (learned nonlinear head on embeddings)")
        preds, info = run_mlp_head(
            X_train, Y_train_filled, X_eval, X_val, Y_val,
            alphas=MLP_ALPHAS, hidden=hidden, seed=config.RANDOM_SEED,
            train_groups=train_groups,
        )
        metrics = evaluate(Y_eval_raw, preds)
        if predictions_out is not None:
            predictions_out[task]["models"]["mlp_head"] = preds
        results["models"]["mlp_head"] = {
            "description": "StandardScaler -> MLPRegressor on Geneformer embeddings.",
            **info, **metrics,
        }
        print(f"        alpha            : {info['alpha']} ({info['alpha_selected_on']})")
        print(f"        spearman mean    : {metrics['spearman_mean']}")
        print(f"        spearman median  : {metrics['spearman_median']}")
        print(f"        targets positive : {metrics['spearman_frac_positive']}")
        print(f"        r2 mean          : {metrics['r2_mean']}")
    else:
        print("\n  [2/2] MLP head -- skipped (--no-mlp)")

    # ---- deltas vs the baseline number to beat -------------------------
    ref = _baseline_reference_for(baseline_ref, task)
    results["baseline_reference"] = ref
    if ref:
        results["deltas"] = _compute_deltas(results["models"], ref)
    return results


def _baseline_reference_for(baseline_ref: dict | None, task: str) -> dict:
    """Pull the ridge_pca and lineage_mean spearman means for this task from
    baseline_results.json, if present."""
    if not baseline_ref:
        return {}
    task_block = baseline_ref.get("tasks", {}).get(task, {})
    models = task_block.get("models", {})
    ref = {}
    if "ridge_pca" in models:
        ref["ridge_pca_spearman_mean"] = models["ridge_pca"].get("spearman_mean")
    if "lineage_mean" in models:
        ref["lineage_mean_spearman_mean"] = models["lineage_mean"].get("spearman_mean")
    return ref


def _compute_deltas(head_models: dict, ref: dict) -> dict:
    """Delta of each head's spearman mean against the baseline references."""
    deltas: dict = {}
    ridge_pca = ref.get("ridge_pca_spearman_mean")
    lineage = ref.get("lineage_mean_spearman_mean")
    for head_name, block in head_models.items():
        val = block.get("spearman_mean")
        if val is None:
            continue
        entry = {}
        if ridge_pca is not None:
            entry["vs_ridge_pca"] = round(val - ridge_pca, 4)
        if lineage is not None:
            entry["vs_lineage_mean"] = round(val - lineage, 4)
        deltas[head_name] = entry
    return deltas


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Train prediction heads on Geneformer embeddings and "
                    "compare to the ridge baseline on the same split."
    )
    parser.add_argument(
        "--split", choices=["val", "test"], default="val",
        help="Which split to evaluate on. Default 'val'. Use 'test' only for "
             "the final reported number.",
    )
    parser.add_argument(
        "--no-mlp", action="store_true",
        help="Run only the linear ridge head (fast). The MLP head is on by default.",
    )
    parser.add_argument(
        "--mlp-hidden", type=int, default=MLP_HIDDEN_DEFAULT,
        help=f"Hidden layer width for the MLP head (default {MLP_HIDDEN_DEFAULT}).",
    )
    parser.add_argument(
        "--processed-dir", default=str(config.PROCESSED_DIR),
        help="Directory holding the processed artifacts and embeddings.",
    )
    parser.add_argument(
        "--save-predictions", action="store_true",
        help=(
            "Also write each head's held-out prediction matrix and the "
            "matching truth matrix to <processed-dir>/predictions/, then "
            "re-score them from disk to prove the files match what was "
            "reported. Off by default, so the plain command remains the one "
            "that produced the published head_results.json."
        ),
    )
    args = parser.parse_args()

    from pathlib import Path
    out = Path(args.processed_dir)

    print("=" * 74)
    print("PREDICTION HEADS ON GENEFORMER EMBEDDINGS")
    print("=" * 74)

    if args.split == "test":
        print("\n  *** Evaluating on the TEST split. ***")
        print("  This should happen once, at the end of the project.")
        print("  If you are still tuning anything, stop and use --split val.\n")

    # ---- load everything, with the same failure messages as baseline ----
    try:
        embeddings = load_embeddings(out)
        crispr = io_utils.load_matrix(out / "crispr_effect")
        metadata = io_utils.load_table(out / "model_metadata")
        selective = io_utils.load_json(out / "selective_genes.json")
        splits_payload = io_utils.load_json(out / "splits.json")
    except FileNotFoundError as exc:
        print(f"\n{exc}")
        print("\nRun build_dataset.py, splits.py, and run_geneformer_embeddings.py first.")
        return 1

    metadata.index.name = config.MODEL_ID
    assignment = pd.Series(splits_payload["assignment"], name="split")
    selective_genes = [g for g in selective["genes"] if g in crispr.columns]

    try:
        prism = io_utils.load_matrix(out / "prism_response")
    except FileNotFoundError:
        prism = None

    # baseline_results.json is optional but expected: it is the number to beat.
    baseline_ref = None
    try:
        baseline_ref = io_utils.load_json(out / "baseline_results.json")
    except FileNotFoundError:
        print("\n  note: baseline_results.json not found -- deltas will be "
              "omitted. Run baseline.py to produce the number to beat.")

    # ---- coverage report: how many cell lines actually have an embedding ---
    covered = embeddings.index.intersection(crispr.index)
    print(f"\n  embeddings : {embeddings.shape[0]} cell lines x "
          f"{embeddings.shape[1]} dims")
    print(f"  crispr     : {len(selective_genes)} selective targets")
    print(f"  overlap    : {len(covered)} cell lines have both an embedding "
          f"and CRISPR data")
    if len(covered) < 0.5 * len(crispr.index):
        print(f"  WARNING    : embeddings cover under half of the CRISPR cell "
              f"lines ({len(covered)}/{len(crispr.index)}). Check that "
              f"run_geneformer_embeddings.py ran on the full panel.")
    print(f"  prism      : {'absent' if prism is None else f'{prism.shape[1]} compounds'}")
    print(f"  eval split : {args.split}")

    all_results: dict = {
        "eval_split": args.split,
        "seed": config.RANDOM_SEED,
        "embedding_source": "geneformer_embeddings",
        "embedding_dims": int(embeddings.shape[1]),
        "mlp": (None if args.no_mlp else {"hidden": args.mlp_hidden, "alphas": MLP_ALPHAS}),
        "metric_note": (
            "Per-target Spearman across held-out cell lines, averaged. Computed "
            "by baseline.evaluate, the same function baseline.py uses, so head "
            "numbers and baseline numbers are the same measurement."
        ),
        "tasks": {},
    }

    predictions_out: dict | None = {} if args.save_predictions else None

    for task in ("crispr", "prism"):
        result = run_task(
            task, embeddings, crispr, prism, selective_genes,
            assignment, metadata, args.split, args.mlp_hidden,
            run_mlp=not args.no_mlp, baseline_ref=baseline_ref,
            predictions_out=predictions_out,
        )
        if result is not None:
            all_results["tasks"][task] = result

    if not all_results["tasks"]:
        print("\nNo task could be run. Check that embeddings overlap the "
              "CRISPR/PRISM cell lines and the splits.")
        return 1

    path = io_utils.save_json(all_results, out / "head_results.json")

    # ---- optional: persist held-out predictions, then prove they round-trip
    if predictions_out:
        print(f"\n{'=' * 74}")
        print("SAVING HELD-OUT PREDICTIONS")
        print("=" * 74)
        written = save_prediction_bundle(predictions_out, out, "head", args.split)
        print(f"\n  {len(written)} file(s) written to {out / 'predictions'}")
        print("\n  Round-trip check -- re-scoring each matrix from disk:")
        ok = verify_prediction_bundle(
            predictions_out, all_results, out, "head", args.split
        )
        if ok:
            print("\n  All saved matrices re-score to the reported values.")
        else:
            print("\n  *** MISMATCH: at least one saved matrix does not re-score")
            print("      to the value reported above. Do not use these files.")
            return 1

    # ------------------------------------------------------------ summary
    print(f"\n{'=' * 74}")
    print("HEAD vs BASELINE")
    print("=" * 74)
    for task, result in all_results["tasks"].items():
        ref = result.get("baseline_reference", {})
        deltas = result.get("deltas", {})
        print(f"\n  {task}  (evaluated on {result['eval_split']}, "
              f"{result['n_eval_lines']} held-out cell lines)")
        if "ridge_pca_spearman_mean" in ref:
            print(f"    baseline ridge on expression : "
                  f"spearman {ref['ridge_pca_spearman_mean']}")
        if "lineage_mean_spearman_mean" in ref:
            print(f"    lineage-mean control         : "
                  f"spearman {ref['lineage_mean_spearman_mean']}")
        for head_name in ("ridge_head", "mlp_head"):
            block = result["models"].get(head_name)
            if not block:
                continue
            line = f"    {head_name:<28}: spearman {block.get('spearman_mean')}"
            d = deltas.get(head_name, {})
            if "vs_ridge_pca" in d:
                line += f"   (vs ridge on expression: {d['vs_ridge_pca']:+})"
            print(line)
        # the honest verdict, stated plainly
        rh = result["models"].get("ridge_head", {}).get("spearman_mean")
        rp = ref.get("ridge_pca_spearman_mean")
        if rh is not None and rp is not None:
            if rh - rp <= 0.01:
                print("      -> embeddings add little beyond raw expression under a")
                print("         linear model. A real, reportable result.")
            else:
                print("      -> embeddings carry linearly-accessible signal beyond")
                print("         raw expression. The gain is the representation.")
        mh = result["models"].get("mlp_head", {}).get("spearman_mean")
        if mh is not None and rh is not None and mh - rh <= 0.01:
            print("      -> the MLP does not beat the linear head on the same")
            print("         features: no evidence the nonlinearity helps here.")

    print(f"\n  Written to: {path}")
    print("\n  Report every head number as a delta against the baseline ridge.")
    print("  A head that does not beat ridge-on-expression is a finding, not a")
    print("  failure -- state it plainly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
