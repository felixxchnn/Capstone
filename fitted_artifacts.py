"""
fitted_artifacts.py
===================
Loader for the **reconstructed** fitted state of the two frozen Phase 1 linear
models (baseline ``ridge_pca``; Geneformer ``ridge_head``).

What these artifacts are -- and are not
--------------------------------------
The original Phase 1 fitted objects (the ``StandardScaler`` / ``PCA`` / ``Ridge``
instances that produced ``baseline_results.json`` and ``head_results.json``) were
**never serialised** and cannot be recovered. The artifacts this module loads are

    "reconstructed fitted state at the frozen Phase 1 alpha from the unchanged
     frozen training data"

i.e. the identical scikit-learn pipeline re-fit on exactly the committed training
split, at the alpha **read from** the committed results files (no hyper-parameter
selection was re-run). They are byte-for-byte reproducible from the committed
inputs, and they reproduce the committed Phase 1 validation metrics exactly at
the precision those metrics were recorded (see ``reconstruct_fitted.py`` and
``capstone/data-integrity-hashes.md``). They are still a *reconstruction*, not
the historical objects.

This module
----------
* imports **numpy + stdlib only** -- no scikit-learn, no ``baseline`` /
  ``train_head`` import, no access to any training array;
* performs inference with plain array arithmetic -- **no** ``fit`` /
  ``fit_transform`` anywhere;
* verifies every ``.npy`` and label file against the SHA-256 recorded in the
  artifact's own ``manifest.json`` on load, and hard-fails on any mismatch,
  malformed manifest, or shape/dtype disagreement.

Transform math (exact equivalents of the fitted sklearn objects, ``whiten=False``,
``with_mean=with_std=True``, multi-output ``Ridge``):

    impute : x[~finite] <- impute_mean[col]           (mirrors baseline.impute_with_train_mean)
    scale  : Xs = (X - scaler_mean) / scaler_scale
    pca    : Z  = (Xs - pca_mean) @ pca_components.T   (baseline only)
    ridge  : Y_hat = Z @ ridge_coef.T + ridge_intercept
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

import config

SCHEMA_VERSION = "reconstructed-fitted/1"

RECONSTRUCTION_STATUS = (
    "reconstructed fitted state at the frozen Phase 1 alpha from the unchanged "
    "frozen training data"
)

_HASH_CHUNK = 1 << 20


class ReconstructedArtifactError(RuntimeError):
    """Raised on a missing / malformed / hash-mismatched reconstructed artifact."""


# --------------------------------------------------------------------------
# low-level helpers
# --------------------------------------------------------------------------

def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(_HASH_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_npy_verified(path: Path, spec: dict) -> np.ndarray:
    if not path.is_file():
        raise ReconstructedArtifactError(f"missing array file: {path}")
    actual_sha = _sha256_file(path)
    if actual_sha != spec["sha256"]:
        raise ReconstructedArtifactError(
            f"{path.name}: sha256 {actual_sha} != manifest {spec['sha256']}"
        )
    arr = np.load(path, allow_pickle=False)
    if list(arr.shape) != list(spec["shape"]):
        raise ReconstructedArtifactError(
            f"{path.name}: shape {list(arr.shape)} != manifest {spec['shape']}"
        )
    if str(arr.dtype) != spec["dtype"]:
        raise ReconstructedArtifactError(
            f"{path.name}: dtype {arr.dtype} != manifest {spec['dtype']}"
        )
    return arr


def _impute(X: np.ndarray, impute_mean: np.ndarray) -> np.ndarray:
    """
    Replace non-finite entries with the per-column training mean.

    Byte-for-byte the same operation as ``baseline.impute_with_train_mean``'s
    inner ``_fill`` (which itself guards all-NaN training columns to 0.0 before
    this point -- that guard is baked into ``impute_mean`` at build time).
    """
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X.reshape(1, -1)
    out = X.copy()
    bad = ~np.isfinite(out)
    if bad.any():
        cols = np.where(bad)[1]
        out[bad] = np.take(impute_mean, cols)
    return out


# --------------------------------------------------------------------------
# model wrappers
# --------------------------------------------------------------------------

class _ReconstructedModel:
    """Common manifest handling for the two reconstructed models."""

    model_name: str = ""

    def __init__(self, artifact_dir: str | Path):
        self.artifact_dir = Path(artifact_dir)
        manifest_path = self.artifact_dir / "manifest.json"
        if not manifest_path.is_file():
            raise ReconstructedArtifactError(f"no manifest.json in {self.artifact_dir}")
        self.manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

        if self.manifest.get("schema_version") != SCHEMA_VERSION:
            raise ReconstructedArtifactError(
                f"{manifest_path}: schema_version "
                f"{self.manifest.get('schema_version')!r} != {SCHEMA_VERSION!r}"
            )
        if self.manifest.get("reconstruction_status") != RECONSTRUCTION_STATUS:
            raise ReconstructedArtifactError(
                f"{manifest_path}: reconstruction_status wording does not match "
                f"the mandated phrasing"
            )
        if self.manifest.get("model") != self.model_name:
            raise ReconstructedArtifactError(
                f"{manifest_path}: model {self.manifest.get('model')!r} != "
                f"{self.model_name!r}"
            )

        # every declared output file (except the manifest) must be present and
        # match its recorded hash -- covers the .npy arrays *and* the label JSON.
        for rel, sha in sorted(self.manifest["output_sha256"].items()):
            fpath = self.artifact_dir / rel
            if not fpath.is_file():
                raise ReconstructedArtifactError(f"missing output file: {fpath}")
            actual = _sha256_file(fpath)
            if actual != sha:
                raise ReconstructedArtifactError(
                    f"{rel}: sha256 {actual} != manifest {sha}"
                )

        self.feature_names: list[str] = json.loads(
            (self.artifact_dir / "feature_names.json").read_text(encoding="utf-8")
        )
        self.target_names: list[str] = json.loads(
            (self.artifact_dir / "target_names.json").read_text(encoding="utf-8")
        )
        arrays = self.manifest["arrays"]
        self._arr = {
            name: _load_npy_verified(self.artifact_dir / spec["file"], spec)
            for name, spec in arrays.items()
        }

        self.alpha: float = float(self.manifest["frozen_alpha"]["value"])
        self.n_features: int = int(self.manifest["feature_order"]["n_features"])
        self.n_targets: int = int(self.manifest["target_order"]["n_targets"])

    # -- ordering hashes, for callers that must assert alignment ----------
    @property
    def feature_order_sha256(self) -> str:
        return self.manifest["feature_order"]["sha256"]

    @property
    def target_order_sha256(self) -> str:
        return self.manifest["target_order"]["sha256"]

    @staticmethod
    def order_sha256(names) -> str:
        """SHA-256 of the newline-joined name list -- the ordering fingerprint."""
        return _sha256_bytes("\n".join(str(n) for n in names).encode("utf-8"))

    def assert_feature_order(self, names) -> None:
        got = self.order_sha256(names)
        if got != self.feature_order_sha256:
            raise ReconstructedArtifactError(
                f"feature order mismatch: caller {got} != artifact "
                f"{self.feature_order_sha256}"
            )

    def assert_target_order(self, names) -> None:
        got = self.order_sha256(names)
        if got != self.target_order_sha256:
            raise ReconstructedArtifactError(
                f"target order mismatch: caller {got} != artifact "
                f"{self.target_order_sha256}"
            )


class ReconstructedBaselineRidgePCA(_ReconstructedModel):
    """impute(train-mean) -> StandardScaler -> PCA(200) -> multi-output Ridge."""

    model_name = "baseline_ridge_pca"

    def transform(self, X: np.ndarray) -> np.ndarray:
        """Raw expression rows (n, n_features) -> PCA scores (n, n_components)."""
        X = np.asarray(X, dtype=float)
        if X.ndim == 1:
            X = X.reshape(1, -1)
        if X.shape[1] != self.n_features:
            raise ReconstructedArtifactError(
                f"expected {self.n_features} expression features, got {X.shape[1]}"
            )
        Xi = _impute(X, self._arr["impute_mean"])
        Xs = (Xi - self._arr["scaler_mean"]) / self._arr["scaler_scale"]
        return (Xs - self._arr["pca_mean"]) @ self._arr["pca_components"].T

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Raw expression rows (n, n_features) -> predicted GeneEffect (n, n_targets)."""
        Z = self.transform(X)
        return Z @ self._arr["ridge_coef"].T + self._arr["ridge_intercept"]


class ReconstructedHeadRidge(_ReconstructedModel):
    """impute(train-mean) -> StandardScaler -> multi-output Ridge on embeddings."""

    model_name = "head_ridge_head"

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Geneformer embedding rows (n, 768) -> predicted GeneEffect (n, n_targets)."""
        X = np.asarray(X, dtype=float)
        if X.ndim == 1:
            X = X.reshape(1, -1)
        if X.shape[1] != self.n_features:
            raise ReconstructedArtifactError(
                f"expected {self.n_features} embedding dims, got {X.shape[1]}"
            )
        Xi = _impute(X, self._arr["impute_mean"])
        Xs = (Xi - self._arr["scaler_mean"]) / self._arr["scaler_scale"]
        return Xs @ self._arr["ridge_coef"].T + self._arr["ridge_intercept"]


# --------------------------------------------------------------------------
# public entry points
# --------------------------------------------------------------------------

def load_baseline_ridge_pca(
    root: str | Path = config.RECONSTRUCTED_FITTED_DIR,
) -> ReconstructedBaselineRidgePCA:
    return ReconstructedBaselineRidgePCA(Path(root) / "baseline_ridge_pca")


def load_head_ridge(
    root: str | Path = config.RECONSTRUCTED_FITTED_DIR,
) -> ReconstructedHeadRidge:
    return ReconstructedHeadRidge(Path(root) / "head_ridge_head")
