"""
baseline.py
===========
Establishes the baseline that any later model must beat.

    python baseline.py                # evaluate on the validation split
    python baseline.py --split test   # final evaluation; use ONCE, at the end

Output
------
    baseline_results.json   all metrics, per model, per task

Why this module matters more than the model that follows it
-----------------------------------------------------------
A fine-tuned transformer that scores 0.34 means nothing on its own. It means
something only next to the number a linear model on raw expression already
achieved. If ridge regression scores 0.33, the transformer has learned almost
nothing and the honest conclusion is that a single-cell foundation model does
not transfer to this task -- which is a real finding, publishable in a capstone,
and considerably more interesting than a vague claim of success.

So this runs first, the number goes on the whiteboard, and every later result
is reported as a delta against it.

The baseline ladder
-------------------
Three models of increasing information, so a gain can be attributed:

1. **Global mean.** Predicts each gene's training-set mean for every cell line.
   The null model. Its per-target rank correlation is undefined (a constant
   prediction has no variance to correlate), which is precisely the point:
   it demonstrates why per-target correlation is the honest metric and pooled
   correlation is not.

2. **Lineage mean.** Predicts the training-set mean *within each tissue
   lineage*. This is the control that matters. Much of what looks like learned
   biology is really "this is a bone tumour, and bone tumours behave like
   this." Any model must beat this to have learned anything context-specific
   beyond tissue identity.

3. **Ridge on PCA of expression.** A real, if simple, learned model.

Metric
------
Per-target Spearman correlation, computed across held-out cell lines for each
target separately, then summarised. Pooled correlation over all
(cell line, gene) pairs is not used: it is dominated by between-gene
differences in mean effect and looks impressive even for a model that has
learned nothing about individual cell lines.
"""

from __future__ import annotations

import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

import config
import io_utils


# --------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------

def per_target_spearman(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    min_samples: int = 5,
) -> np.ndarray:
    """
    Spearman correlation for each target column, across samples.

    Returns an array of length n_targets. Entries are NaN where the target has
    too few finite observations, or where either side is constant (correlation
    is undefined, not zero -- conflating the two would flatter a constant
    predictor).
    """
    n_targets = y_true.shape[1]
    out = np.full(n_targets, np.nan, dtype=float)

    for j in range(n_targets):
        true_col = y_true[:, j]
        pred_col = y_pred[:, j]
        mask = np.isfinite(true_col) & np.isfinite(pred_col)
        if mask.sum() < min_samples:
            continue
        t = true_col[mask]
        p = pred_col[mask]
        if np.all(t == t[0]) or np.all(p == p[0]):
            continue
        rho, _ = stats.spearmanr(t, p)
        out[j] = rho

    return out


def per_target_r2(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
    """
    Coefficient of determination per target, against the held-out target mean.

    Unlike Spearman this is defined for a constant predictor, so it is the
    metric that lets the global-mean baseline be scored at all. Values below
    zero mean the model is worse than predicting the held-out mean.
    """
    n_targets = y_true.shape[1]
    out = np.full(n_targets, np.nan, dtype=float)

    for j in range(n_targets):
        mask = np.isfinite(y_true[:, j]) & np.isfinite(y_pred[:, j])
        if mask.sum() < 2:
            continue
        t = y_true[mask, j]
        p = y_pred[mask, j]
        ss_res = float(np.sum((t - p) ** 2))
        ss_tot = float(np.sum((t - t.mean()) ** 2))
        if ss_tot <= 0:
            continue
        out[j] = 1.0 - ss_res / ss_tot

    return out


def summarise_metric(values: np.ndarray, name: str) -> dict:
    """Summary statistics for a per-target metric array."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {
            f"{name}_mean": None,
            f"{name}_median": None,
            f"{name}_q25": None,
            f"{name}_q75": None,
            f"{name}_frac_positive": None,
            f"{name}_n_targets_scored": 0,
            f"{name}_n_targets_undefined": int(values.size),
        }
    return {
        f"{name}_mean": round(float(np.mean(finite)), 4),
        f"{name}_median": round(float(np.median(finite)), 4),
        f"{name}_q25": round(float(np.percentile(finite, 25)), 4),
        f"{name}_q75": round(float(np.percentile(finite, 75)), 4),
        f"{name}_frac_positive": round(float((finite > 0).mean()), 4),
        f"{name}_n_targets_scored": int(finite.size),
        f"{name}_n_targets_undefined": int(values.size - finite.size),
    }


def evaluate(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    """Full metric bundle for one model on one task."""
    rho = per_target_spearman(y_true, y_pred)
    r2 = per_target_r2(y_true, y_pred)
    metrics = {}
    metrics.update(summarise_metric(rho, "spearman"))
    metrics.update(summarise_metric(r2, "r2"))
    return metrics


# --------------------------------------------------------------------------
# Data preparation
# --------------------------------------------------------------------------

def impute_with_train_mean(
    train: np.ndarray,
    *others: np.ndarray,
) -> tuple[np.ndarray, ...]:
    """
    Replace NaNs using column means computed on the training set only.

    Computing the fill values on the full dataset would leak held-out
    information into training. It is a small leak, and it is still a leak.
    """
    with np.errstate(invalid="ignore"):
        means = np.nanmean(train, axis=0)
    means = np.where(np.isfinite(means), means, 0.0)

    def _fill(arr: np.ndarray) -> np.ndarray:
        filled = arr.copy()
        idx = np.where(~np.isfinite(filled))
        if idx[0].size:
            filled[idx] = np.take(means, idx[1])
        return filled

    return (_fill(train), *(_fill(o) for o in others))


def prepare_task(
    task: str,
    expression: pd.DataFrame,
    crispr: pd.DataFrame,
    prism: pd.DataFrame | None,
    selective_genes: list[str],
    assignment: pd.Series,
    eval_split: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Index, pd.Index] | None:
    """
    Assemble features and targets for one task.

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

    lines = expression.index.intersection(targets.index)
    lines = pd.Index(sorted(lines), name=config.MODEL_ID)
    if len(lines) == 0:
        return None

    X = expression.loc[lines]
    Y = targets.loc[lines]

    split_here = assignment.reindex(lines)
    train_index = lines[(split_here == "train").to_numpy()]
    eval_index = lines[(split_here == eval_split).to_numpy()]

    if len(train_index) < 10 or len(eval_index) < 5:
        return None

    return X, Y, train_index, eval_index


# --------------------------------------------------------------------------
# Baseline models
# --------------------------------------------------------------------------

def run_global_mean(Y_train: np.ndarray, n_eval: int) -> np.ndarray:
    """Predict each target's training mean for every held-out cell line."""
    with np.errstate(invalid="ignore"):
        means = np.nanmean(Y_train, axis=0)
    means = np.where(np.isfinite(means), means, 0.0)
    return np.tile(means, (n_eval, 1))


def run_lineage_mean(
    Y_train: np.ndarray,
    lineage_train: np.ndarray,
    lineage_eval: np.ndarray,
) -> np.ndarray:
    """
    Predict the training mean within each tissue lineage.

    Cell lines whose lineage is absent from training fall back to the global
    training mean.
    """
    with np.errstate(invalid="ignore"):
        global_means = np.nanmean(Y_train, axis=0)
    global_means = np.where(np.isfinite(global_means), global_means, 0.0)

    per_lineage: dict[str, np.ndarray] = {}
    for lineage in np.unique(lineage_train):
        mask = lineage_train == lineage
        if mask.sum() == 0:
            continue
        with np.errstate(invalid="ignore"):
            means = np.nanmean(Y_train[mask], axis=0)
        per_lineage[lineage] = np.where(np.isfinite(means), means, global_means)

    predictions = np.vstack([
        per_lineage.get(lineage, global_means) for lineage in lineage_eval
    ])
    return predictions


def _select_alpha_inner_cv(
    Z_train: np.ndarray,
    Y_train: np.ndarray,
    groups: np.ndarray | None,
    alphas: list[float],
    seed: int,
    n_folds: int = 5,
) -> tuple[float, dict]:
    """
    Choose the ridge penalty by cross-validation *within the training set*.

    Used when the evaluation split is `val`, because tuning on the split you
    are about to report would make that number meaningless. Folds are grouped
    by patient so the inner evaluation carries the same leakage guarantee as
    the outer split.
    """
    n_samples = Z_train.shape[0]
    if groups is None:
        groups = np.arange(n_samples)

    unique_groups = np.unique(groups)
    n_folds = int(min(n_folds, len(unique_groups)))
    if n_folds < 2:
        return alphas[len(alphas) // 2], {"status": "too few groups for inner CV"}

    splitter = GroupKFold(n_splits=n_folds)
    scores: dict[float, list[float]] = {alpha: [] for alpha in alphas}

    for inner_train, inner_test in splitter.split(Z_train, groups=groups):
        if len(inner_test) < 3:
            continue
        for alpha in alphas:
            model = Ridge(alpha=alpha)
            model.fit(Z_train[inner_train], Y_train[inner_train])
            preds = model.predict(Z_train[inner_test])
            rho = per_target_spearman(Y_train[inner_test], preds)
            if np.isfinite(rho).any():
                scores[alpha].append(float(np.nanmean(rho)))

    mean_scores = {
        alpha: (float(np.mean(vals)) if vals else -np.inf)
        for alpha, vals in scores.items()
    }
    best = max(mean_scores, key=lambda a: mean_scores[a])

    sweep = {
        str(alpha): (round(score, 4) if np.isfinite(score) else None)
        for alpha, score in mean_scores.items()
    }
    return best, {"status": "ok", "n_folds": n_folds, "sweep": sweep}


def run_ridge_pca(
    X_train: np.ndarray,
    Y_train: np.ndarray,
    X_eval: np.ndarray,
    X_val: np.ndarray | None,
    Y_val: np.ndarray | None,
    n_components: int,
    alphas: list[float],
    seed: int,
    train_groups: np.ndarray | None = None,
) -> tuple[np.ndarray, dict]:
    """
    Standardise -> PCA -> multi-output ridge.

    The regularisation strength is never chosen on the split being reported.
    When a separate validation split is available (i.e. the final test-set
    evaluation) it is used. Otherwise the penalty is selected by grouped
    cross-validation inside the training set.
    """
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_eval_s = scaler.transform(X_eval)

    n_components = int(min(n_components, X_train_s.shape[0] - 1, X_train_s.shape[1]))
    n_components = max(n_components, 1)

    pca = PCA(n_components=n_components, random_state=seed)
    Z_train = pca.fit_transform(X_train_s)
    Z_eval = pca.transform(X_eval_s)

    info = {
        "n_components": n_components,
        "explained_variance_ratio": round(
            float(np.sum(pca.explained_variance_ratio_)), 4
        ),
    }

    if X_val is not None and Y_val is not None and len(X_val) >= 5:
        Z_val = pca.transform(scaler.transform(X_val))
        best_score = -np.inf
        chosen_alpha = alphas[len(alphas) // 2]
        sweep = {}
        for alpha in alphas:
            model = Ridge(alpha=alpha)
            model.fit(Z_train, Y_train)
            preds = model.predict(Z_val)
            rho = per_target_spearman(Y_val, preds)
            score = float(np.nanmean(rho)) if np.isfinite(rho).any() else -np.inf
            sweep[str(alpha)] = round(score, 4) if np.isfinite(score) else None
            if score > best_score:
                best_score = score
                chosen_alpha = alpha
        info["alpha_sweep"] = sweep
        info["alpha_selected_on"] = "held-out validation split"
    else:
        chosen_alpha, cv_info = _select_alpha_inner_cv(
            Z_train, Y_train, train_groups, alphas, seed
        )
        info["alpha_sweep"] = cv_info.get("sweep")
        info["alpha_selected_on"] = (
            f"grouped {cv_info.get('n_folds', '?')}-fold CV inside the training set"
            if cv_info.get("status") == "ok"
            else f"default ({cv_info.get('status')})"
        )

    info["alpha"] = chosen_alpha

    model = Ridge(alpha=chosen_alpha)
    model.fit(Z_train, Y_train)
    return model.predict(Z_eval), info


# --------------------------------------------------------------------------
# Held-out prediction persistence
# --------------------------------------------------------------------------
# `evaluate` returns summary statistics and discards the per-target correlation
# vectors that produced them. That is the right shape for a results file, but it
# means every downstream question about uncertainty -- a bootstrap over held-out
# cell lines, a paired test across targets, the correlation between two models'
# per-target performance -- requires refitting both models just to recover
# arrays that already existed in memory.
#
# These helpers write the held-out predictions and the matching truth matrix to
# disk so that work becomes free. They are opt-in (`--save-predictions`): the
# default run is unchanged, so the command that produced the published
# baseline_results.json still produces exactly that file.

def save_prediction_bundle(
    predictions: dict,
    processed_dir,
    source: str,
    eval_split: str,
) -> list:
    """
    Write held-out predictions and the matching truth matrix to disk.

    One matrix per (task, model) plus one truth matrix per task, each indexed by
    ModelID with one column per target, written through io_utils.save_matrix so
    they carry the same format and labelling discipline as every other matrix in
    the project.

    Files land in <processed_dir>/predictions/ rather than <processed_dir>
    itself, so the fourteen-file integrity record for data/processed is
    unaffected by turning this on.

    Parameters
    ----------
    predictions
        Bundle built by run_task when it is handed a `predictions_out` dict.
    source
        "baseline" or "head" -- keeps the two modules' outputs distinguishable.
    eval_split
        The split the predictions were made on. Part of the filename so a val
        run and a test run cannot overwrite each other.

    Returns
    -------
    List of paths actually written.
    """
    out_dir = Path(processed_dir) / "predictions"
    written = []

    for task, bundle in predictions.items():
        index = pd.Index(bundle["eval_index"], name=config.MODEL_ID)
        columns = bundle["target_columns"]

        truth = pd.DataFrame(bundle["y_true"], index=index, columns=columns)
        written.append(io_utils.save_matrix(
            truth, out_dir / f"{source}_{task}_{eval_split}_y_true"
        ))

        for model_name, preds in bundle["models"].items():
            frame = pd.DataFrame(preds, index=index, columns=columns)
            written.append(io_utils.save_matrix(
                frame, out_dir / f"{source}_{task}_{eval_split}_{model_name}"
            ))

    return written


def verify_prediction_bundle(
    predictions: dict,
    results: dict,
    processed_dir,
    source: str,
    eval_split: str,
) -> bool:
    """
    Re-score the saved matrices from disk and compare to the in-memory metrics.

    Without this the failure mode is silent. Predictions could be written with
    the wrong row order, or a model's array captured from the wrong branch, and
    the printed metrics would still be correct -- because they were computed
    from the in-memory array, not from the file. Round-tripping is the only
    check that the thing on disk is the thing that was scored.

    Returns True if every saved model matrix re-scores to the recorded
    spearman_mean, False otherwise.
    """
    out_dir = Path(processed_dir) / "predictions"
    all_ok = True

    for task, bundle in predictions.items():
        task_block = results.get("tasks", {}).get(task)
        if task_block is None:
            continue

        truth = io_utils.load_matrix(
            out_dir / f"{source}_{task}_{eval_split}_y_true"
        ).to_numpy(dtype=float)

        for model_name in bundle["models"]:
            recorded = task_block["models"].get(model_name, {}).get("spearman_mean")
            reloaded = io_utils.load_matrix(
                out_dir / f"{source}_{task}_{eval_split}_{model_name}"
            ).to_numpy(dtype=float)
            round_tripped = evaluate(truth, reloaded)["spearman_mean"]

            if recorded is None and round_tripped is None:
                ok = True
            elif recorded is None or round_tripped is None:
                ok = False
            else:
                ok = abs(float(recorded) - float(round_tripped)) < 1e-9

            all_ok = all_ok and ok
            status = "ok" if ok else "MISMATCH"
            print(f"        {task}/{model_name:<14}: recorded {recorded} "
                  f"-> reloaded {round_tripped}   [{status}]")

    return all_ok


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------

def run_task(
    task: str,
    expression: pd.DataFrame,
    crispr: pd.DataFrame,
    prism: pd.DataFrame | None,
    selective_genes: list[str],
    assignment: pd.Series,
    metadata: pd.DataFrame,
    eval_split: str,
    predictions_out: dict | None = None,
) -> dict | None:
    """
    Run the full baseline ladder for one prediction task.

    If `predictions_out` is a dict it is filled in place with the held-out
    predictions of every model, the matching truth matrix, and the row and
    column labels needed to write them out. Passing None -- the default --
    leaves behaviour byte-for-byte as it was before this parameter existed.
    """
    prepared = prepare_task(
        task, expression, crispr, prism, selective_genes, assignment, eval_split
    )
    if prepared is None:
        print(f"\n  Task {task!r}: insufficient data -- skipping.")
        return None

    X, Y, train_index, eval_index = prepared

    print(f"\n{'=' * 74}")
    print(f"TASK: {task}")
    print("=" * 74)
    print(f"  targets            : {Y.shape[1]}")
    print(f"  training lines     : {len(train_index)}")
    print(f"  {eval_split + ' lines':<19}: {len(eval_index)}")

    X_train_raw = X.loc[train_index].to_numpy(dtype=float)
    X_eval_raw = X.loc[eval_index].to_numpy(dtype=float)
    Y_train_raw = Y.loc[train_index].to_numpy(dtype=float)
    Y_eval_raw = Y.loc[eval_index].to_numpy(dtype=float)

    X_train, X_eval = impute_with_train_mean(X_train_raw, X_eval_raw)
    Y_train_filled, = impute_with_train_mean(Y_train_raw)

    # Validation data for hyperparameter selection (only when evaluating on test).
    X_val = Y_val = None
    if eval_split == "test":
        val_index = X.index[(assignment.reindex(X.index) == "val").to_numpy()]
        if len(val_index) >= 5:
            X_val_raw = X.loc[val_index].to_numpy(dtype=float)
            X_val, = impute_with_train_mean(X_train_raw, X_val_raw)[1:]
            Y_val = Y.loc[val_index].to_numpy(dtype=float)

    results: dict = {
        "task": task,
        "eval_split": eval_split,
        "n_targets": int(Y.shape[1]),
        "n_train_lines": int(len(train_index)),
        "n_eval_lines": int(len(eval_index)),
        "models": {},
    }

    # Truth matrix and labels are recorded once; each model appends its own
    # prediction array below. Y_eval_raw is stored *unfilled*, with its NaNs
    # intact, because per_target_spearman masks non-finite entries per column
    # and imputing them here would silently change the metric downstream.
    if predictions_out is not None:
        predictions_out[task] = {
            "eval_index": [str(i) for i in eval_index],
            "target_columns": [str(c) for c in Y.columns],
            "y_true": Y_eval_raw,
            "models": {},
        }

    # ---- 1. global mean ------------------------------------------------
    print("\n  [1/3] global mean (null model)")
    preds = run_global_mean(Y_train_raw, len(eval_index))
    metrics = evaluate(Y_eval_raw, preds)
    if predictions_out is not None:
        predictions_out[task]["models"]["global_mean"] = preds
    results["models"]["global_mean"] = {
        "description": "Predicts each target's training mean for every line.",
        **metrics,
    }
    print(f"        spearman : {metrics['spearman_mean']} "
          f"(undefined for {metrics['spearman_n_targets_undefined']} targets "
          f"-- a constant prediction has no variance to rank)")
    print(f"        r2       : {metrics['r2_mean']}")

    # ---- 2. lineage mean -----------------------------------------------
    lineage_col = config.STRATIFY_COL
    if lineage_col in metadata.columns:
        print("\n  [2/3] lineage mean (tissue-identity control)")
        lineage_train = (
            metadata.loc[train_index, lineage_col].astype(str).to_numpy()
        )
        lineage_eval = (
            metadata.loc[eval_index, lineage_col].astype(str).to_numpy()
        )
        preds = run_lineage_mean(Y_train_raw, lineage_train, lineage_eval)
        metrics = evaluate(Y_eval_raw, preds)
        if predictions_out is not None:
            predictions_out[task]["models"]["lineage_mean"] = preds
        results["models"]["lineage_mean"] = {
            "description": (
                "Predicts the training mean within each tissue lineage. Any "
                "model must beat this to have learned something beyond "
                "tissue identity."
            ),
            **metrics,
        }
        print(f"        spearman : {metrics['spearman_mean']}")
        print(f"        r2       : {metrics['r2_mean']}")
    else:
        print(f"\n  [2/3] lineage mean -- skipped ({lineage_col} not in metadata)")

    # ---- 3. ridge on PCA ------------------------------------------------
    print("\n  [3/3] ridge regression on PCA of expression")
    if config.GROUP_COL in metadata.columns:
        train_groups = (
            metadata.loc[train_index, config.GROUP_COL]
            .astype(str)
            .fillna(pd.Series(train_index, index=train_index).astype(str))
            .to_numpy()
        )
    else:
        train_groups = None

    preds, info = run_ridge_pca(
        X_train,
        Y_train_filled,
        X_eval,
        X_val,
        Y_val,
        n_components=config.PCA_COMPONENTS,
        alphas=config.RIDGE_ALPHAS,
        seed=config.RANDOM_SEED,
        train_groups=train_groups,
    )
    metrics = evaluate(Y_eval_raw, preds)
    if predictions_out is not None:
        predictions_out[task]["models"]["ridge_pca"] = preds
    results["models"]["ridge_pca"] = {
        "description": "StandardScaler -> PCA -> multi-output Ridge.",
        **info,
        **metrics,
    }
    print(f"        components        : {info['n_components']} "
          f"({info['explained_variance_ratio']:.1%} of variance)")
    print(f"        alpha             : {info['alpha']} "
          f"({info['alpha_selected_on']})")
    print(f"        spearman mean     : {metrics['spearman_mean']}")
    print(f"        spearman median   : {metrics['spearman_median']}")
    print(f"        targets positive  : {metrics['spearman_frac_positive']}")
    print(f"        r2 mean           : {metrics['r2_mean']}")

    return results


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Baseline models for DepMap perturbation response."
    )
    parser.add_argument(
        "--split",
        choices=["val", "test"],
        default="val",
        help=(
            "Which split to evaluate on. Default 'val'. Use 'test' only for "
            "the final reported number -- every extra look at the test set "
            "quietly turns it into a validation set."
        ),
    )
    parser.add_argument(
        "--save-predictions",
        action="store_true",
        help=(
            "Also write the held-out prediction matrices and the matching "
            "truth matrix to data/processed/predictions/, then re-score them "
            "from disk to prove the files match what was reported. Off by "
            "default, so the plain command remains the one that produced the "
            "published baseline_results.json."
        ),
    )
    args = parser.parse_args()

    out = config.PROCESSED_DIR

    print("=" * 74)
    print("BASELINE MODELS")
    print("=" * 74)

    if args.split == "test":
        print("\n  *** Evaluating on the TEST split. ***")
        print("  This should happen once, at the end of the project.")
        print("  If you are still tuning anything, stop and use --split val.\n")

    try:
        expression = io_utils.load_matrix(out / "expression")
        crispr = io_utils.load_matrix(out / "crispr_effect")
        metadata = io_utils.load_table(out / "model_metadata")
        selective = io_utils.load_json(out / "selective_genes.json")
        splits_payload = io_utils.load_json(out / "splits.json")
    except FileNotFoundError as exc:
        print(f"\n{exc}")
        print("\nRun build_dataset.py and splits.py first.")
        return 1

    metadata.index.name = config.MODEL_ID
    assignment = pd.Series(splits_payload["assignment"], name="split")
    selective_genes = [g for g in selective["genes"] if g in crispr.columns]

    try:
        prism = io_utils.load_matrix(out / "prism_response")
    except FileNotFoundError:
        prism = None

    print(f"\n  expression : {expression.shape[0]} lines x {expression.shape[1]} genes")
    print(f"  crispr     : {len(selective_genes)} selective targets")
    print(f"  prism      : {'absent' if prism is None else f'{prism.shape[1]} compounds'}")
    print(f"  eval split : {args.split}")

    all_results: dict = {
        "eval_split": args.split,
        "seed": config.RANDOM_SEED,
        "pca_components": config.PCA_COMPONENTS,
        "ridge_alphas": config.RIDGE_ALPHAS,
        "metric_note": (
            "Spearman is computed per target across held-out cell lines, then "
            "averaged. Pooled correlation over all (line, target) pairs is not "
            "reported: it is dominated by between-target differences in mean "
            "effect and flatters models that have learned nothing about "
            "individual cell lines."
        ),
        "tasks": {},
    }

    predictions_out: dict | None = {} if args.save_predictions else None

    for task in ("crispr", "prism"):
        result = run_task(
            task, expression, crispr, prism, selective_genes,
            assignment, metadata, args.split,
            predictions_out=predictions_out,
        )
        if result is not None:
            all_results["tasks"][task] = result

    if not all_results["tasks"]:
        print("\nNo task could be run. Check the dataset.")
        return 1

    path = io_utils.save_json(all_results, out / "baseline_results.json")

    # ---- optional: persist held-out predictions, then prove they round-trip
    if predictions_out:
        print(f"\n{'=' * 74}")
        print("SAVING HELD-OUT PREDICTIONS")
        print("=" * 74)
        written = save_prediction_bundle(
            predictions_out, out, "baseline", args.split
        )
        print(f"\n  {len(written)} file(s) written to {out / 'predictions'}")
        print("\n  Round-trip check -- re-scoring each matrix from disk:")
        ok = verify_prediction_bundle(
            predictions_out, all_results, out, "baseline", args.split
        )
        if ok:
            print("\n  All saved matrices re-score to the reported values.")
        else:
            print("\n  *** MISMATCH: at least one saved matrix does not re-score")
            print("      to the value reported above. Do not use these files.")
            return 1

    # ---------------------------------------------------------- summary
    print(f"\n{'=' * 74}")
    print("THE NUMBER TO BEAT")
    print("=" * 74)
    for task, result in all_results["tasks"].items():
        ridge = result["models"].get("ridge_pca", {})
        lineage = result["models"].get("lineage_mean", {})
        print(f"\n  {task}  (evaluated on {result['eval_split']}, "
              f"{result['n_eval_lines']} held-out cell lines)")
        print(f"    ridge on expression : spearman {ridge.get('spearman_mean')}")
        if lineage:
            print(f"    lineage mean control: spearman {lineage.get('spearman_mean')}")
            r_val = ridge.get("spearman_mean")
            l_val = lineage.get("spearman_mean")
            if r_val is not None and l_val is not None:
                delta = round(r_val - l_val, 4)
                print(f"    gain from expression: {delta:+}")
                if delta <= 0.01:
                    print("      -> expression adds little beyond tissue identity.")
                    print("         Worth reporting plainly; it is a real result.")

    print(f"\n  Written to: {path}")
    print("\n  Put the ridge number on a whiteboard. Every later model gets")
    print("  reported as a delta against it.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
