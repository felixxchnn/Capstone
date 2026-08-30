"""
reconstruct_fitted.py
=====================
Generate (and validate) **reconstructed** fitted state for the two frozen
Phase 1 linear models:

* ``baseline_ridge_pca`` -- impute(train-mean) -> StandardScaler -> PCA(200)
  -> multi-output Ridge(alpha = 100000.0)
* ``head_ridge_head``   -- impute(train-mean) -> StandardScaler
  -> multi-output Ridge(alpha = 3162.0)

    py reconstruct_fitted.py --build              # write the artifacts
    py reconstruct_fitted.py --validate           # reload + reproduce Phase 1 val metrics
    py reconstruct_fitted.py --check-determinism  # build twice, require byte-identical
    py reconstruct_fitted.py --self-test          # synthetic end-to-end, offline

Honesty
-------
The original Phase 1 ``StandardScaler`` / ``PCA`` / ``Ridge`` objects were never
serialised. What this module writes is

    "reconstructed fitted state at the frozen Phase 1 alpha from the unchanged
     frozen training data"

The pipeline, estimator options, random seed, feature order and target order are
identical to ``baseline.run_ridge_pca`` / ``train_head.run_ridge_head``; the
alpha is **read from** ``baseline_results.json`` / ``head_results.json`` (never
re-selected). Fitting touches **only** the committed ``train`` split. The result
reproduces the committed Phase 1 validation metrics exactly at the precision
those metrics were recorded (``--validate`` enforces this and hard-stops on any
disagreement).

Serialisation
-------------
Plain ``.npy`` arrays (no pickle / joblib), original ``float64`` dtype preserved,
plus a ``manifest.json``. ``.npy`` headers carry no timestamp, so a rebuild is
byte-identical (``--check-determinism``). Loaded, with no ``fit`` and no training
data, by ``fitted_artifacts.py``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

import config
import io_utils
import fitted_artifacts
from fitted_artifacts import (
    SCHEMA_VERSION,
    RECONSTRUCTION_STATUS,
    ReconstructedArtifactError,
)

import baseline
import train_head

# The repository commit whose committed training data + frozen results these
# artifacts reconstruct. Not `git rev-parse HEAD` at build time: the build must
# stay byte-reproducible after the artifacts themselves are committed. Integrity
# is enforced instead by EXPECTED_INPUT_SHA256 below.
BASE_COMMIT = "12fab80a705d0adf473ca07dd9b455f1b807fc35"

# Exact SHA-256 of every committed input this reconstruction consumes, as of
# BASE_COMMIT. If any of these no longer matches, the frozen training data or
# frozen results have moved and a *new* reconstruction base must be established
# deliberately -- the build hard-fails rather than silently mislabel its output.
EXPECTED_INPUT_SHA256 = {
    "data/processed/expression.npz":
        "3d5bfa0c3430584f8943fd2365be0eecf8b994b38bfc7d491d59d7b9ff251a2d",
    "data/processed/expression.labels.json":
        "d18005cc0aec3e4d5f0fd06c748ef66672256cd8c7a6f24ea4c441b0ca785983",
    "data/processed/crispr_effect.npz":
        "9214efa3ce172079e6ce4ca78853d8bf92fb8f6d4a55d0c6c71e4653b59e8826",
    "data/processed/crispr_effect.labels.json":
        "165906f07e61819c8fadb2bf3c95a73817e538a22a34f63c415a30222ac49b9f",
    "data/processed/selective_genes.json":
        "68c8fe39ae8965ce20b04f50870609cc21734386ceeff859f4d0bddd2e5bab35",
    "data/processed/splits.json":
        "f1419abc7cbd31efc173a5857bab9eb318b53f8e535a17048bfcf0ea2f70aeef",
    "data/processed/model_metadata.csv":
        "1c314197b57c1f8363eb44f8902b3733777e7c304d7f677c76d401e3cabe5180",
    "data/processed/baseline_results.json":
        "b49169bd363a596f400b4faff8c21d354275b70404efe08b9109d38f1bdc0ffd",
    "data/processed/head_results.json":
        "1962206fa17646cbd1fec4b642a577cc2586c09c4cabd980541a7e11a8b6f894",
    "data/processed/geneformer_embeddings.csv":
        "af8ee6d734bea11101d07884f1c72d2b4efaff9875506738a037102a712f1e46",
}

# Metric-reproduction tolerance. The reconstruction re-executes the identical
# fit, so agreement is expected to full float64 precision; Phase 1 recorded its
# summary statistics rounded to 4 decimals (baseline.summarise_metric), so the
# committed check is "equal after rounding to 4 dp".
METRIC_ROUND_DP = 4
_METRIC_KEYS = (
    "spearman_mean", "spearman_median", "spearman_q25", "spearman_q75",
    "spearman_frac_positive", "spearman_n_targets_scored",
    "spearman_n_targets_undefined",
    "r2_mean", "r2_median", "r2_q25", "r2_q75", "r2_frac_positive",
    "r2_n_targets_scored", "r2_n_targets_undefined",
)

_HASH_CHUNK = 1 << 20


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(_HASH_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _order_sha256(names) -> str:
    return _sha256_bytes("\n".join(str(n) for n in names).encode("utf-8"))


def _train_column_means(raw: np.ndarray) -> np.ndarray:
    """
    The per-column fill vector used by ``baseline.impute_with_train_mean``:
    ``np.nanmean`` over the training rows, with all-NaN columns guarded to 0.0.
    Replicated here (the function inlines it and returns filled arrays, not the
    vector) so the vector can be serialised for artifact-only imputation.
    """
    with np.errstate(invalid="ignore"):
        means = np.nanmean(raw, axis=0)
    return np.where(np.isfinite(means), means, 0.0)


def _save_npy(arr: np.ndarray, path: Path) -> dict:
    """Write one C-contiguous float64 .npy (deterministic; no timestamp)."""
    arr = np.ascontiguousarray(arr)
    if arr.dtype != np.float64:
        raise ValueError(f"{path.name}: dtype {arr.dtype} != float64 "
                         f"(refusing to quantize)")
    with open(path, "wb") as handle:
        np.save(handle, arr, allow_pickle=False)
    return {
        "file": path.name,
        "shape": list(arr.shape),
        "dtype": str(arr.dtype),
        "sha256": _sha256_file(path),
    }


def _write_json_deterministic(obj, path: Path) -> None:
    path.write_text(
        json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _verify_inputs(processed_dir: Path) -> dict:
    """Gate every committed input against EXPECTED_INPUT_SHA256. Returns the map."""
    root = processed_dir.parent.parent  # .../data/processed -> repo root
    got: dict[str, str] = {}
    for rel, expected in EXPECTED_INPUT_SHA256.items():
        p = root / rel
        if not p.is_file():
            raise ReconstructedArtifactError(f"missing committed input: {rel}")
        actual = _sha256_file(p)
        got[rel] = actual
        if actual != expected:
            raise ReconstructedArtifactError(
                f"{rel}: sha256 {actual} != expected {expected} for base commit "
                f"{BASE_COMMIT}. The frozen training data or results have "
                f"changed; establish a new reconstruction base deliberately."
            )
    return got


# --------------------------------------------------------------------------
# data assembly (mirrors baseline.py / train_head.py main() exactly)
# --------------------------------------------------------------------------

def _load_common(processed_dir: Path):
    expression = io_utils.load_matrix(processed_dir / "expression")
    crispr = io_utils.load_matrix(processed_dir / "crispr_effect")
    metadata = io_utils.load_table(processed_dir / "model_metadata")
    metadata.index.name = config.MODEL_ID
    selective = io_utils.load_json(processed_dir / "selective_genes.json")
    splits_payload = io_utils.load_json(processed_dir / "splits.json")
    assignment = pd.Series(splits_payload["assignment"], name="split")
    selective_genes = [g for g in selective["genes"] if g in crispr.columns]
    return expression, crispr, metadata, assignment, selective_genes


def _split_indices(feature_index: pd.Index, target_index: pd.Index,
                   assignment: pd.Series):
    """Reproduce prepare_task's row selection: sorted feature∩target, then split."""
    lines = pd.Index(sorted(feature_index.intersection(target_index)),
                     name=config.MODEL_ID)
    s = assignment.reindex(lines)
    train = lines[(s == "train").to_numpy()]
    val = lines[(s == "val").to_numpy()]
    test = lines[(s == "test").to_numpy()]
    return lines, train, val, test


# --------------------------------------------------------------------------
# build one model
# --------------------------------------------------------------------------

def _estimator_params(est) -> dict:
    return {k: (v if isinstance(v, (int, float, str, bool, type(None))) else str(v))
            for k, v in sorted(est.get_params().items())}


def _build_baseline(processed_dir: Path, input_hashes: dict) -> tuple[dict, dict]:
    """Returns (arrays_by_name, manifest_core) for baseline_ridge_pca."""
    expression, crispr, metadata, assignment, selective_genes = _load_common(processed_dir)
    targets = crispr[selective_genes]
    lines, train_idx, val_idx, test_idx = _split_indices(
        expression.index, targets.index, assignment)

    # cross-check against the real prepare_task
    X, Y, pt_train, pt_eval = baseline.prepare_task(
        "crispr", expression, crispr, None, selective_genes, assignment, "val")
    assert list(pt_train) == list(train_idx), "train index disagrees with prepare_task"
    assert list(pt_eval) == list(val_idx), "val index disagrees with prepare_task"
    assert set(train_idx).isdisjoint(val_idx) and set(train_idx).isdisjoint(test_idx), \
        "train overlaps val/test"

    feature_names = [str(c) for c in X.columns]
    target_names = [str(c) for c in Y.columns]
    assert len(feature_names) == 18460, len(feature_names)
    assert len(target_names) == 4297, len(target_names)

    X_train_raw = X.loc[train_idx].to_numpy(dtype=float)
    Y_train_raw = Y.loc[train_idx].to_numpy(dtype=float)
    X_train, = baseline.impute_with_train_mean(X_train_raw)
    Y_train_filled, = baseline.impute_with_train_mean(Y_train_raw)
    impute_mean = _train_column_means(X_train_raw)

    alpha, alpha_block = _read_frozen_alpha(
        processed_dir / "baseline_results.json",
        ("tasks", "crispr", "models", "ridge_pca"))

    # ---- fit, exactly as baseline.run_ridge_pca (X_val=None path) ----
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    n_components = int(min(config.PCA_COMPONENTS, X_train_s.shape[0] - 1,
                           X_train_s.shape[1]))
    n_components = max(n_components, 1)
    pca = PCA(n_components=n_components, random_state=config.RANDOM_SEED)
    Z_train = pca.fit_transform(X_train_s)
    ridge = Ridge(alpha=alpha)
    ridge.fit(Z_train, Y_train_filled)

    evr_sum = round(float(np.sum(pca.explained_variance_ratio_)), 4)
    committed_evr = alpha_block.get("explained_variance_ratio")
    assert committed_evr is None or abs(evr_sum - committed_evr) < 1e-9, \
        f"explained_variance_ratio {evr_sum} != committed {committed_evr}"

    arrays = {
        "impute_mean": impute_mean,
        "scaler_mean": scaler.mean_,
        "scaler_scale": scaler.scale_,
        "scaler_var": scaler.var_,
        "pca_mean": pca.mean_,
        "pca_components": pca.components_,
        "pca_explained_variance": pca.explained_variance_,
        "pca_explained_variance_ratio": pca.explained_variance_ratio_,
        "pca_singular_values": pca.singular_values_,
        "ridge_coef": ridge.coef_,
        "ridge_intercept": np.atleast_1d(ridge.intercept_).astype(float),
    }
    core = {
        "model": "baseline_ridge_pca",
        "pipeline": ("impute(train-mean) -> StandardScaler -> PCA -> "
                     "multi-output Ridge"),
        "estimator_params": {
            "StandardScaler": _estimator_params(scaler),
            "PCA": _estimator_params(pca),
            "Ridge": _estimator_params(ridge),
        },
        "resolved_solvers": {
            "PCA_fit_svd_solver": str(getattr(pca, "_fit_svd_solver", "unknown")),
            "Ridge_solver": "auto -> cholesky (dense multi-output, no sample_weight)",
        },
        "pca_n_components": int(n_components),
        "pca_noise_variance": float(pca.noise_variance_),
        "scaler_n_samples_seen": int(np.atleast_1d(scaler.n_samples_seen_)[0]),
        "explained_variance_ratio_sum": evr_sum,
        "frozen_alpha": {
            "value": float(alpha),
            "source_file": "data/processed/baseline_results.json",
            "source_json_path": "tasks.crispr.models.ridge_pca.alpha",
            "selection": "read verbatim from the committed results; NOT re-selected",
            "recorded_alpha_selected_on": alpha_block.get("alpha_selected_on"),
            "recorded_alpha_at_grid_boundary": alpha_block.get("alpha_at_grid_boundary"),
        },
        "train_split": {
            "n_train_lines": int(len(train_idx)),
            "n_val_lines": int(len(val_idx)),
            "n_test_lines": int(len(test_idx)),
            "row_source": ("sorted(expression.index intersect "
                           "crispr[selective].index), then splits.json "
                           "assignment == 'train'"),
            "train_model_ids_sha256": _order_sha256(list(train_idx)),
            "disjoint_from_val_test": True,
            "no_val_or_test_row_used_for_fitting": True,
        },
        "feature_order": {
            "n_features": len(feature_names),
            "source": "expression matrix columns (expression.labels.json)",
            "sha256": _order_sha256(feature_names),
            "file": "feature_names.json",
        },
        "target_order": {
            "n_targets": len(target_names),
            "source": ("selective_genes.json 'genes', filtered to crispr_effect "
                       "columns, order preserved"),
            "sha256": _order_sha256(target_names),
            "file": "target_names.json",
        },
        "_feature_names": feature_names,
        "_target_names": target_names,
        "_input_hashes": input_hashes,
    }
    return arrays, core


def _build_head(processed_dir: Path, input_hashes: dict) -> tuple[dict, dict]:
    """Returns (arrays_by_name, manifest_core) for head_ridge_head."""
    _expr, crispr, metadata, assignment, selective_genes = _load_common(processed_dir)
    embeddings = train_head.load_embeddings(processed_dir)
    targets = crispr[selective_genes]
    lines, train_idx, val_idx, test_idx = _split_indices(
        embeddings.index, targets.index, assignment)

    X, Y, pt_train, pt_eval = train_head.prepare_task(
        "crispr", embeddings, crispr, None, selective_genes, assignment, "val")
    assert list(pt_train) == list(train_idx), "train index disagrees with prepare_task"
    assert list(pt_eval) == list(val_idx), "val index disagrees with prepare_task"
    assert set(train_idx).isdisjoint(val_idx) and set(train_idx).isdisjoint(test_idx), \
        "train overlaps val/test"

    feature_names = [str(c) for c in X.columns]
    target_names = [str(c) for c in Y.columns]
    assert len(feature_names) == 768, len(feature_names)
    assert len(target_names) == 4297, len(target_names)

    X_train_raw = X.loc[train_idx].to_numpy(dtype=float)
    Y_train_raw = Y.loc[train_idx].to_numpy(dtype=float)
    X_train, = baseline.impute_with_train_mean(X_train_raw)
    Y_train_filled, = baseline.impute_with_train_mean(Y_train_raw)
    impute_mean = _train_column_means(X_train_raw)

    alpha, alpha_block = _read_frozen_alpha(
        processed_dir / "head_results.json",
        ("tasks", "crispr", "models", "ridge_head"))

    # ---- fit, exactly as train_head.run_ridge_head (X_val=None path) ----
    scaler = StandardScaler()
    Z_train = scaler.fit_transform(X_train)
    ridge = Ridge(alpha=alpha)
    ridge.fit(Z_train, Y_train_filled)

    arrays = {
        "impute_mean": impute_mean,
        "scaler_mean": scaler.mean_,
        "scaler_scale": scaler.scale_,
        "scaler_var": scaler.var_,
        "ridge_coef": ridge.coef_,
        "ridge_intercept": np.atleast_1d(ridge.intercept_).astype(float),
    }
    core = {
        "model": "head_ridge_head",
        "pipeline": "impute(train-mean) -> StandardScaler -> multi-output Ridge",
        "estimator_params": {
            "StandardScaler": _estimator_params(scaler),
            "Ridge": _estimator_params(ridge),
        },
        "resolved_solvers": {
            "Ridge_solver": "auto -> cholesky (dense multi-output, no sample_weight)",
        },
        "scaler_n_samples_seen": int(np.atleast_1d(scaler.n_samples_seen_)[0]),
        "frozen_alpha": {
            "value": float(alpha),
            "source_file": "data/processed/head_results.json",
            "source_json_path": "tasks.crispr.models.ridge_head.alpha",
            "selection": "read verbatim from the committed results; NOT re-selected",
            "recorded_alpha_selected_on": alpha_block.get("alpha_selected_on"),
            "recorded_alpha_at_grid_boundary": alpha_block.get("alpha_at_grid_boundary"),
        },
        "train_split": {
            "n_train_lines": int(len(train_idx)),
            "n_val_lines": int(len(val_idx)),
            "n_test_lines": int(len(test_idx)),
            "row_source": ("sorted(geneformer_embeddings.index intersect "
                           "crispr[selective].index), then splits.json "
                           "assignment == 'train'"),
            "train_model_ids_sha256": _order_sha256(list(train_idx)),
            "disjoint_from_val_test": True,
            "no_val_or_test_row_used_for_fitting": True,
        },
        "feature_order": {
            "n_features": len(feature_names),
            "source": "geneformer_embeddings.csv columns (0..767)",
            "sha256": _order_sha256(feature_names),
            "file": "feature_names.json",
        },
        "target_order": {
            "n_targets": len(target_names),
            "source": ("selective_genes.json 'genes', filtered to crispr_effect "
                       "columns, order preserved"),
            "sha256": _order_sha256(target_names),
            "file": "target_names.json",
        },
        "_feature_names": feature_names,
        "_target_names": target_names,
        "_input_hashes": input_hashes,
    }
    return arrays, core


def _read_frozen_alpha(results_path: Path, json_path: tuple) -> tuple[float, dict]:
    obj = io_utils.load_json(results_path)
    block = obj
    for key in json_path:
        block = block[key]
    alpha = block["alpha"]
    grid = (config.RIDGE_ALPHAS if "baseline" in results_path.name
            else train_head.HEAD_RIDGE_ALPHAS)
    if float(alpha) not in {float(a) for a in grid}:
        raise ReconstructedArtifactError(
            f"{results_path.name}: frozen alpha {alpha} is not on the recorded "
            f"grid {grid}")
    return float(alpha), block


# --------------------------------------------------------------------------
# serialise a built model to a directory
# --------------------------------------------------------------------------

def _write_model_dir(model_dir: Path, arrays: dict, core: dict) -> dict:
    if model_dir.exists():
        shutil.rmtree(model_dir)
    model_dir.mkdir(parents=True)

    feature_names = core.pop("_feature_names")
    target_names = core.pop("_target_names")
    input_hashes = core.pop("_input_hashes")

    _write_json_deterministic(feature_names, model_dir / "feature_names.json")
    _write_json_deterministic(target_names, model_dir / "target_names.json")

    array_specs = {}
    for name, arr in arrays.items():
        array_specs[name] = _save_npy(np.asarray(arr, dtype=float),
                                      model_dir / f"{name}.npy")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "reconstruction_status": RECONSTRUCTION_STATUS,
        "not_original_note": (
            "These are NOT the original Phase 1 fitted objects -- those were "
            "never serialised. Scaler/PCA/Ridge were re-fit on exactly the "
            "committed train split, at the alpha read from the committed "
            "results file (no hyper-parameter selection). Verified to reproduce "
            "the committed Phase 1 validation metrics at the recorded precision."
        ),
        "base_commit": BASE_COMMIT,
        "generated_by": "reconstruct_fitted.py",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "scikit_learn": __import__("sklearn").__version__,
            "pandas": pd.__version__,
        },
        "input_artifact_sha256": input_hashes,
        "arrays": array_specs,
        **core,
    }

    _write_json_deterministic(manifest, model_dir / "manifest.json")

    # output_sha256 covers every file except manifest.json (which cannot hash
    # itself); recorded into the manifest with a second, final write.
    outputs = {}
    for p in sorted(model_dir.iterdir()):
        if p.name == "manifest.json":
            continue
        outputs[p.name] = _sha256_file(p)
    manifest["output_sha256"] = outputs
    _write_json_deterministic(manifest, model_dir / "manifest.json")
    return manifest


# --------------------------------------------------------------------------
# public build
# --------------------------------------------------------------------------

def build(out_root: Path | str = config.RECONSTRUCTED_FITTED_DIR,
          processed_dir: Path | str = config.PROCESSED_DIR,
          *, verbose: bool = True) -> dict:
    out_root = Path(out_root)
    processed_dir = Path(processed_dir)
    input_hashes = _verify_inputs(processed_dir)

    out_root.mkdir(parents=True, exist_ok=True)
    manifests = {}
    for name, builder in (("baseline_ridge_pca", _build_baseline),
                          ("head_ridge_head", _build_head)):
        arrays, core = builder(processed_dir, input_hashes)
        manifests[name] = _write_model_dir(out_root / name, arrays, core)
        if verbose:
            total = sum(s["shape"] and int(np.prod(s["shape"])) or 1
                        for s in manifests[name]["arrays"].values())
            print(f"  [built] {name}: {len(manifests[name]['arrays'])} arrays, "
                  f"{total:,} float64 values")
    return manifests


# --------------------------------------------------------------------------
# validation: reload (no fit) and reproduce the committed Phase 1 val metrics
# --------------------------------------------------------------------------

def _rebuild_val_arrays(processed_dir: Path):
    """Raw (unimputed) val X and unfilled val Y for both arms, as Phase 1 scores them."""
    expression, crispr, metadata, assignment, selective_genes = _load_common(processed_dir)

    Xb, Yb, trb, evb = baseline.prepare_task(
        "crispr", expression, crispr, None, selective_genes, assignment, "val")
    Xb_eval_raw = Xb.loc[evb].to_numpy(dtype=float)
    Yb_eval_raw = Yb.loc[evb].to_numpy(dtype=float)

    embeddings = train_head.load_embeddings(processed_dir)
    Xh, Yh, trh, evh = train_head.prepare_task(
        "crispr", embeddings, crispr, None, selective_genes, assignment, "val")
    Xh_eval_raw = Xh.loc[evh].to_numpy(dtype=float)
    Yh_eval_raw = Yh.loc[evh].to_numpy(dtype=float)

    return {
        "baseline_ridge_pca": {
            "X_eval_raw": Xb_eval_raw, "Y_eval_raw": Yb_eval_raw,
            "feature_names": [str(c) for c in Xb.columns],
            "target_names": [str(c) for c in Yb.columns],
            "eval_ids": [str(i) for i in evb],
        },
        "head_ridge_head": {
            "X_eval_raw": Xh_eval_raw, "Y_eval_raw": Yh_eval_raw,
            "feature_names": [str(c) for c in Xh.columns],
            "target_names": [str(c) for c in Yh.columns],
            "eval_ids": [str(i) for i in evh],
        },
    }


def _committed_metric_block(processed_dir: Path, model: str) -> dict:
    if model == "baseline_ridge_pca":
        obj = io_utils.load_json(processed_dir / "baseline_results.json")
        return obj["tasks"]["crispr"]["models"]["ridge_pca"]
    obj = io_utils.load_json(processed_dir / "head_results.json")
    return obj["tasks"]["crispr"]["models"]["ridge_head"]


def validate(out_root: Path | str = config.RECONSTRUCTED_FITTED_DIR,
             processed_dir: Path | str = config.PROCESSED_DIR,
             *, verbose: bool = True) -> dict:
    out_root = Path(out_root)
    processed_dir = Path(processed_dir)

    loaders = {
        "baseline_ridge_pca": fitted_artifacts.load_baseline_ridge_pca(out_root),
        "head_ridge_head": fitted_artifacts.load_head_ridge(out_root),
    }
    val = _rebuild_val_arrays(processed_dir)
    gi_pred_dir = processed_dir / "predictions"
    report: dict = {"models": {}, "all_ok": True}

    for model, loader in loaders.items():
        v = val[model]
        loader.assert_feature_order(v["feature_names"])
        loader.assert_target_order(v["target_names"])

        pred = loader.predict(v["X_eval_raw"])
        metrics = baseline.evaluate(v["Y_eval_raw"], pred)
        committed = _committed_metric_block(processed_dir, model)

        diffs = {}
        ok = True
        for k in _METRIC_KEYS:
            got = metrics.get(k)
            exp = committed.get(k)
            if isinstance(got, float):
                got_r = round(got, METRIC_ROUND_DP)
            else:
                got_r = got
            match = (got_r == exp)
            if not match:
                ok = False
            diffs[k] = {"reconstructed": got_r, "committed": exp, "match": match}

        # optional, non-authoritative: compare to the gitignored prediction matrix
        gi = {"status": "absent"}
        stem = ("baseline_crispr_val_ridge_pca" if model == "baseline_ridge_pca"
                else "head_crispr_val_ridge_head")
        if (gi_pred_dir / f"{stem}.npz").is_file():
            saved = io_utils.load_matrix(gi_pred_dir / stem).to_numpy(dtype=float)
            gi = {
                "status": "present",
                "authoritative": False,
                "shape_match": list(saved.shape) == list(pred.shape),
                "max_abs_diff": float(np.max(np.abs(saved - pred)))
                if saved.shape == pred.shape else None,
            }

        report["models"][model] = {
            "eval_split": "val",
            "n_eval_lines": len(v["eval_ids"]),
            "tolerance": (f"equal after rounding to {METRIC_ROUND_DP} dp (the "
                          f"precision Phase 1 recorded; underlying agreement is "
                          f"full float64 -- the identical fit is re-executed)"),
            "committed_source": (
                "baseline_results.json::tasks.crispr.models.ridge_pca"
                if model == "baseline_ridge_pca"
                else "head_results.json::tasks.crispr.models.ridge_head"),
            "metrics": diffs,
            "all_within_tolerance": ok,
            "gitignored_prediction_matrix_compare": gi,
        }
        report["all_ok"] = report["all_ok"] and ok
        if verbose:
            flag = "ok" if ok else "MISMATCH"
            print(f"  [{flag}] {model}: "
                  f"spearman_mean recon={diffs['spearman_mean']['reconstructed']} "
                  f"committed={diffs['spearman_mean']['committed']}"
                  + ("" if gi["status"] == "absent"
                     else f"  | vs gitignored matrix max-abs-diff={gi.get('max_abs_diff')}"))
    return report


# --------------------------------------------------------------------------
# determinism
# --------------------------------------------------------------------------

def check_determinism(processed_dir: Path | str = config.PROCESSED_DIR,
                      *, verbose: bool = True) -> dict:
    processed_dir = Path(processed_dir)
    tmp = Path(tempfile.mkdtemp(prefix="recon_det_"))
    try:
        a, b = tmp / "a", tmp / "b"
        build(a, processed_dir, verbose=False)
        build(b, processed_dir, verbose=False)
        mismatches = []
        files_a = sorted(p.relative_to(a) for p in a.rglob("*") if p.is_file())
        files_b = sorted(p.relative_to(b) for p in b.rglob("*") if p.is_file())
        if files_a != files_b:
            mismatches.append(f"file set differs: {files_a} vs {files_b}")
        for rel in files_a:
            if (a / rel).read_bytes() != (b / rel).read_bytes():
                mismatches.append(str(rel))
        ok = not mismatches
        if verbose:
            print(f"  determinism: {'byte-identical' if ok else 'MISMATCH'} "
                  f"across {len(files_a)} files"
                  + ("" if ok else f"  -> {mismatches}"))
        return {"ok": ok, "n_files": len(files_a), "mismatches": mismatches}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------
# verify: run the full step-7 checklist against the committed artifacts
# --------------------------------------------------------------------------

def verify(out_root: Path | str = config.RECONSTRUCTED_FITTED_DIR,
           processed_dir: Path | str = config.PROCESSED_DIR,
           *, verbose: bool = True) -> dict:
    out_root = Path(out_root)
    processed_dir = Path(processed_dir)
    checks: list[tuple[str, bool, str]] = []

    def _c(name: str, ok: bool, detail: str = "") -> None:
        checks.append((name, bool(ok), detail))
        if verbose:
            print(f"  [{'ok' if ok else 'FAIL'}] {name}"
                  + (f"  -- {detail}" if detail else ""))

    # 1. committed inputs unchanged
    try:
        _verify_inputs(processed_dir)
        _c("protected Phase 1 inputs unchanged (EXPECTED_INPUT_SHA256)", True)
    except ReconstructedArtifactError as exc:
        _c("protected Phase 1 inputs unchanged (EXPECTED_INPUT_SHA256)", False, str(exc))

    # 2. on-disk artifacts == a fresh deterministic rebuild
    tmp = Path(tempfile.mkdtemp(prefix="recon_verify_"))
    try:
        build(tmp, processed_dir, verbose=False)
        rebuilt = sorted(p.relative_to(tmp) for p in tmp.rglob("*") if p.is_file())
        on_disk = sorted(p.relative_to(out_root) for p in out_root.rglob("*") if p.is_file())
        same_set = rebuilt == on_disk
        byte_same = same_set and all(
            (tmp / rel).read_bytes() == (out_root / rel).read_bytes() for rel in rebuilt)
        _c("committed artifacts are byte-identical to a fresh rebuild",
           same_set and byte_same,
           "" if byte_same else f"set_match={same_set}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # 3. loaders (hash-verify every file on load) + no-fit source guarantee
    try:
        bl = fitted_artifacts.load_baseline_ridge_pca(out_root)
        hd = fitted_artifacts.load_head_ridge(out_root)
        _c("both artifacts load with every file hash-verified", True)
    except ReconstructedArtifactError as exc:
        _c("both artifacts load with every file hash-verified", False, str(exc))
        return {"checks": checks, "n_fail": sum(1 for _, ok, _ in checks if not ok)}

    src = Path(fitted_artifacts.__file__).read_text(encoding="utf-8")
    no_fit = not any(b in src for b in (".fit(", ".fit_transform(",
                                        "import sklearn", "from sklearn"))
    _c("fitted_artifacts.py: no sklearn import, no fit()/fit_transform()", no_fit)

    # 4. alpha in each manifest == alpha in the frozen results file (re-read)
    a_bl, _ = _read_frozen_alpha(processed_dir / "baseline_results.json",
                                 ("tasks", "crispr", "models", "ridge_pca"))
    a_hd, _ = _read_frozen_alpha(processed_dir / "head_results.json",
                                 ("tasks", "crispr", "models", "ridge_head"))
    _c("baseline manifest alpha == frozen baseline_results.json alpha",
       bl.alpha == a_bl == 100000.0, f"{bl.alpha} / {a_bl}")
    _c("head manifest alpha == frozen head_results.json alpha",
       hd.alpha == a_hd == 3162.0, f"{hd.alpha} / {a_hd}")
    _c("manifests record alpha as read-not-selected",
       bl.manifest["frozen_alpha"]["selection"].startswith("read verbatim")
       and hd.manifest["frozen_alpha"]["selection"].startswith("read verbatim"))

    # 5. training rows exactly equal the committed train split; disjoint from val/test
    expression, crispr, _md, assignment, selective_genes = _load_common(processed_dir)
    targets = crispr[selective_genes]
    _, tr_b, val_b, test_b = _split_indices(expression.index, targets.index, assignment)
    emb = train_head.load_embeddings(processed_dir)
    _, tr_h, val_h, test_h = _split_indices(emb.index, targets.index, assignment)
    _c("baseline train rows == splits.json 'train' (independent recompute)",
       _order_sha256(list(tr_b)) == bl.manifest["train_split"]["train_model_ids_sha256"]
       and len(tr_b) == 800)
    _c("head train rows == splits.json 'train' (independent recompute)",
       _order_sha256(list(tr_h)) == hd.manifest["train_split"]["train_model_ids_sha256"]
       and len(tr_h) == 800)
    _c("no val/test row is in either training set",
       set(tr_b).isdisjoint(val_b) and set(tr_b).isdisjoint(test_b)
       and set(tr_h).isdisjoint(val_h) and set(tr_h).isdisjoint(test_h)
       and len(val_b) == len(test_b) == 170)

    # 6. exact feature + target ordering, checked against the committed sources
    try:
        bl.assert_feature_order([str(c) for c in expression.columns])
        bl.assert_target_order([g for g in
                                io_utils.load_json(processed_dir / "selective_genes.json")["genes"]
                                if g in crispr.columns])
        hd.assert_feature_order([str(c) for c in emb.columns])
        hd.assert_target_order([g for g in
                                io_utils.load_json(processed_dir / "selective_genes.json")["genes"]
                                if g in crispr.columns])
        _c("feature + target ordering match the committed source files", True)
    except ReconstructedArtifactError as exc:
        _c("feature + target ordering match the committed source files", False, str(exc))

    # 7. a mismatched artifact hash hard-fails the loader (on a throwaway copy)
    tmp2 = Path(tempfile.mkdtemp(prefix="recon_tamper_"))
    try:
        dst = tmp2 / "baseline_ridge_pca"
        shutil.copytree(out_root / "baseline_ridge_pca", dst)
        blob = bytearray((dst / "ridge_coef.npy").read_bytes())
        blob[-1] ^= 0x01
        (dst / "ridge_coef.npy").write_bytes(bytes(blob))
        try:
            fitted_artifacts.ReconstructedBaselineRidgePCA(dst)
            _c("tampered artifact hard-fails the loader", False, "no error raised")
        except ReconstructedArtifactError:
            _c("tampered artifact hard-fails the loader", True)
    finally:
        shutil.rmtree(tmp2, ignore_errors=True)

    # 8. reconstructed metrics reproduce the frozen Phase 1 results
    rep = validate(out_root, processed_dir, verbose=False)
    _c("reconstructed val metrics reproduce baseline_results.json / head_results.json",
       rep["all_ok"])

    n_fail = sum(1 for _, ok, _ in checks if not ok)
    return {"checks": checks, "n_fail": n_fail, "n_pass": len(checks) - n_fail,
            "metric_report": rep}


# --------------------------------------------------------------------------
# self-test (synthetic, offline)
# --------------------------------------------------------------------------

def _self_test() -> int:  # noqa: C901 -- one linear scenario
    print("Running reconstruct_fitted.py self-test...")
    tmp = Path(tempfile.mkdtemp(prefix="recon_selftest_"))
    try:
        rng = np.random.default_rng(0)
        n, p, d, t = 60, 40, 12, 7
        # --- synthetic loader round-trip: fit tiny sklearn, serialise, reload ---
        Xtr = rng.normal(size=(n, p))
        Ytr = rng.normal(size=(n, t))
        sc = StandardScaler().fit(Xtr)
        pca = PCA(n_components=d, random_state=1).fit(sc.transform(Xtr))
        rg = Ridge(alpha=123.0).fit(pca.transform(sc.transform(Xtr)), Ytr)

        mdir = tmp / "baseline_ridge_pca"
        arrays = {
            "impute_mean": _train_column_means(Xtr),
            "scaler_mean": sc.mean_, "scaler_scale": sc.scale_, "scaler_var": sc.var_,
            "pca_mean": pca.mean_, "pca_components": pca.components_,
            "pca_explained_variance": pca.explained_variance_,
            "pca_explained_variance_ratio": pca.explained_variance_ratio_,
            "pca_singular_values": pca.singular_values_,
            "ridge_coef": rg.coef_, "ridge_intercept": np.atleast_1d(rg.intercept_).astype(float),
        }
        core = {
            "model": "baseline_ridge_pca", "pipeline": "test",
            "estimator_params": {"Ridge": _estimator_params(rg)},
            "frozen_alpha": {"value": 123.0, "source_file": "x", "source_json_path": "x",
                             "selection": "test", "recorded_alpha_selected_on": "x",
                             "recorded_alpha_at_grid_boundary": False},
            "train_split": {"n_train_lines": n, "n_val_lines": 0, "n_test_lines": 0,
                            "row_source": "test", "train_model_ids_sha256": "x",
                            "disjoint_from_val_test": True,
                            "no_val_or_test_row_used_for_fitting": True},
            "feature_order": {"n_features": p, "source": "test",
                              "sha256": _order_sha256([f"f{i}" for i in range(p)]),
                              "file": "feature_names.json"},
            "target_order": {"n_targets": t, "source": "test",
                             "sha256": _order_sha256([f"t{i}" for i in range(t)]),
                             "file": "target_names.json"},
            "_feature_names": [f"f{i}" for i in range(p)],
            "_target_names": [f"t{i}" for i in range(t)],
            "_input_hashes": {"synthetic": "0" * 64},
        }
        _write_model_dir(mdir, dict(arrays), dict(core))

        m1 = _write_model_dir(tmp / "a", dict(arrays), dict(core))
        m2 = _write_model_dir(tmp / "b", dict(arrays), dict(core))
        for p_ in (tmp / "a").iterdir():
            assert p_.read_bytes() == (tmp / "b" / p_.name).read_bytes(), \
                f"non-deterministic file: {p_.name}"
        print("  [ok] repeated serialisation is byte-identical")

        loader = fitted_artifacts.ReconstructedBaselineRidgePCA(mdir)
        Xnew = rng.normal(size=(5, p))
        got = loader.predict(Xnew)
        want = rg.predict(pca.transform(sc.transform(Xnew)))
        assert np.allclose(got, want, rtol=0, atol=1e-9), np.max(np.abs(got - want))
        print("  [ok] artifact-only inference matches the fitted sklearn pipeline "
              "(no fit/fit_transform in the loader)")

        # source-level guarantee: the loader module never fits
        src = Path(fitted_artifacts.__file__).read_text(encoding="utf-8")
        for banned in (".fit(", ".fit_transform(", "import sklearn", "from sklearn"):
            assert banned not in src, f"fitted_artifacts.py contains {banned!r}"
        print("  [ok] fitted_artifacts.py imports no sklearn and calls no fit()")

        # hash tamper -> hard fail
        comp = mdir / "pca_components.npy"
        raw = bytearray(comp.read_bytes())
        raw[-1] ^= 0x01
        comp.write_bytes(bytes(raw))
        try:
            fitted_artifacts.ReconstructedBaselineRidgePCA(mdir)
            raise AssertionError("tampered array should have hard-failed")
        except ReconstructedArtifactError as exc:
            assert "sha256" in str(exc)
        print("  [ok] a mismatched artifact hash hard-fails the loader")

        # feature-order mismatch -> hard fail
        _write_model_dir(mdir, dict(arrays), dict(core))
        good = fitted_artifacts.ReconstructedBaselineRidgePCA(mdir)
        try:
            good.assert_feature_order([f"f{i}" for i in range(p - 1)] + ["WRONG"])
            raise AssertionError("feature-order mismatch should raise")
        except ReconstructedArtifactError:
            pass
        print("  [ok] feature/target order assertions catch a mismatch")

    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("\nSelf-test passed.")
    return 0


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def _print_sizes(out_root: Path) -> None:
    import gzip
    print(f"\n  artifact sizes under {out_root}:")
    grand_raw = grand_gz = 0
    for model_dir in sorted(p for p in out_root.iterdir() if p.is_dir()):
        raw = gz = 0
        for f in sorted(model_dir.iterdir()):
            b = f.read_bytes()
            raw += len(b)
            gz += len(gzip.compress(b, 9, mtime=0))
        grand_raw += raw
        grand_gz += gz
        print(f"    {model_dir.name:22s} {raw:>12,} B  ({gz:>12,} B gzip-9)")
    print(f"    {'TOTAL':22s} {grand_raw:>12,} B  ({grand_gz:>12,} B gzip-9)")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Generate/validate reconstructed fitted state for the two "
                    "frozen Phase 1 linear models.")
    ap.add_argument("--build", action="store_true",
                    help="Fit and serialise the artifacts.")
    ap.add_argument("--validate", action="store_true",
                    help="Reload (no fit) and reproduce the committed Phase 1 "
                         "val metrics; non-zero exit on any disagreement.")
    ap.add_argument("--check-determinism", action="store_true",
                    help="Build twice into temp dirs; require byte-identical output.")
    ap.add_argument("--verify", action="store_true",
                    help="Run the full checklist against the committed artifacts "
                         "(inputs unchanged, byte-identical rebuild, no-fit loader, "
                         "alpha read-not-selected, train split, ordering, tamper "
                         "hard-fail, metric reproduction).")
    ap.add_argument("--self-test", action="store_true",
                    help="Synthetic end-to-end check, offline.")
    ap.add_argument("--out-dir", default=str(config.RECONSTRUCTED_FITTED_DIR))
    args = ap.parse_args()

    if args.self_test:
        return _self_test()

    out_root = Path(args.out_dir)
    did_something = False
    rc = 0

    if args.build:
        did_something = True
        print("=" * 74)
        print("BUILD reconstructed fitted state")
        print("=" * 74)
        build(out_root)
        _print_sizes(out_root)

    if args.check_determinism:
        did_something = True
        print("\n" + "=" * 74)
        print("DETERMINISM CHECK")
        print("=" * 74)
        det = check_determinism()
        rc = rc or (0 if det["ok"] else 1)

    if args.validate:
        did_something = True
        print("\n" + "=" * 74)
        print("VALIDATE against committed Phase 1 metrics")
        print("=" * 74)
        rep = validate(out_root)
        print()
        print("  ALL METRICS REPRODUCE" if rep["all_ok"]
              else "  *** RECONSTRUCTION DOES NOT REPRODUCE THE FROZEN RESULTS ***")
        rc = rc or (0 if rep["all_ok"] else 1)

    if args.verify:
        did_something = True
        print("\n" + "=" * 74)
        print("VERIFY checklist (committed artifacts)")
        print("=" * 74)
        vr = verify(out_root)
        print()
        print(f"  {'ALL CHECKS PASSED' if vr['n_fail'] == 0 else str(vr['n_fail']) + ' FAILED'}"
              f"  ({vr.get('n_pass', 0)}/{vr.get('n_pass', 0) + vr['n_fail']})")
        rc = rc or (0 if vr["n_fail"] == 0 else 1)

    if not did_something:
        ap.print_help()
        return 2
    return rc


if __name__ == "__main__":
    sys.exit(main())
