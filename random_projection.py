"""
random_projection.py
====================
E1 -- the pre-planned random-projection control (CLAUDE.md section 9.3).

    py random_projection.py --run                 # compute + write the val result
    py random_projection.py --run --save-predictions
    py random_projection.py --validate-artifact    # portable: structure + provenance, no recompute
    py random_projection.py --validate            # full numerical regen; needs the pinned env
    py random_projection.py --self-test           # synthetic, offline, portable
    py random_projection.py --check-determinism    # build twice, require identical JSON

Reproducibility environment
---------------------------
A fixed integer seed does **not** make sklearn's ``GaussianRandomProjection``
components byte-identical across NumPy / scikit-learn versions.  The committed
projection-component hash reproduces only under the exact stack recorded in
``EXPECTED_E1_ENVIRONMENT`` (Python 3.14.6 + the pins in ``requirements-e1.txt``:
numpy 2.5.0, pandas 3.0.3, scipy 1.18.0, scikit-learn 1.9.0).  ``--run``, full
``--validate`` and ``--check-determinism`` refuse to recompute outside that
stack and point the user at ``requirements-e1.txt``.  ``--validate-artifact``
and ``--self-test`` are portable and never recompute the real projection.

Scientific question
-------------------
Does Geneformer pretraining provide more useful predictive information than a
*fixed random linear compression of expression at the same 768-dimensional
width*?  If a seeded Gaussian random projection of expression, fed through the
identical ridge head, matches or beats the frozen Geneformer ridge head, then
"a random 768-dim linear sketch of expression" is the more honest baseline for
the embeddings -- a stronger result, not a broken one.  It is reported here
regardless of which way it comes out.

This is an **exploratory validation control**.  It is a sixth model scored on
the same 170 val cell lines; the pre-specified headline comparison
(``ridge_pca`` vs ``ridge_head``) was fixed before E1 existed and is not
reopened.  No test-split feature, outcome, prediction, metric, or performance
number is produced anywhere in this module -- there is deliberately no
``--split test`` path.

Pipeline (exactly, and every fit touches the 800 training rows only)
-------------------------------------------------------------------
    impute(train-mean expression)
      -> StandardScaler                     [fit on train]
      -> GaussianRandomProjection(n_components=768,
                                  random_state=config.RANDOM_SEED)
                                            [fit on standardized train]
      -> StandardScaler                     [fit on train, inside
                                             train_head.run_ridge_head]
      -> multi-output Ridge, alpha by patient-grouped inner 5-fold CV on the
         800 training lines (train_head.HEAD_RIDGE_ALPHAS, 13 points)
      -> per-target Spearman across the 170 val cell lines
         (baseline.evaluate / baseline.per_target_spearman)

Alpha-grid discipline (CLAUDE.md invariant 7)
--------------------------------------------
The projected columns carry variance ~= 18460/768 ~= 24 each -- a third
spectrum, between the head's unit-variance embeddings and the baseline's PCA
eigenvalues.  Alpha is selected on its own copy of the 13-point grid and **must
be interior**.  If it lands on the grid minimum or maximum this module writes
the result (so the sweep is on record) and exits non-zero: the grid decision is
then a separate, approved step -- the grid is never widened automatically and
the seed is never changed to chase an interior alpha.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.random_projection import GaussianRandomProjection

import config
import io_utils
import baseline
import train_head


# --------------------------------------------------------------------------
# constants
# --------------------------------------------------------------------------

SCHEMA_VERSION = "random-projection-control/1"
RESULT_FILE = config.PROCESSED_DIR / "random_projection_results.json"

N_COMPONENTS = 768                       # matched to the Geneformer embedding width
N_EXPRESSION_FEATURES = 18460            # canonical expression columns, pre-projection
N_TARGETS = 4297                         # selective CRISPR targets
N_TRAIN_LINES = 800
N_VAL_LINES = 170

# Theoretical projected-column variance for a unit-variance input under a
# Gaussian random projection with entries ~ N(0, 1/n_components):
#   var(proj_j) ~= n_features * (1/n_components) * var(standardized column)
#              ~= 18460 / 768
THEORETICAL_PROJECTED_VARIANCE = N_EXPRESSION_FEATURES / N_COMPONENTS   # ~= 24.036458

# The observed mean projected variance must sit within this multiplicative band
# of the theoretical value or the run aborts for investigation rather than
# reporting a number produced by a mis-wired pipeline.
_VARIANCE_CONSISTENCY_BAND = (0.5, 2.0)

# Reference values E1 is reported as deltas against. Read from the committed
# results files at run time AND asserted equal to these literals (fail closed).
REFERENCE_EXPECTED = {
    "ridge_pca_spearman_mean": 0.2356,
    "ridge_head_spearman_mean": 0.2047,
    "lineage_mean_spearman_mean": 0.15,
}

# The one authoritative environment in which the committed E1 numbers regenerate
# byte-for-byte. A fixed seed is NOT enough: sklearn's GaussianRandomProjection
# draws its matrix through NumPy's Generator/BitGenerator machinery, whose exact
# byte stream has changed across NumPy releases, so the seeded component matrix
# -- and therefore its SHA-256 -- differs between versions. An independent
# checkout on Python 3.12.13 / numpy 2.3.5 / pandas 2.2.3 / scipy 1.17.0 /
# scikit-learn 1.8.0 produced a different component hash and 16/18 on --validate.
# requirements-e1.txt pins the library half; Python 3.14.6 is the interpreter.
EXPECTED_E1_ENVIRONMENT = {
    "python": "3.14.6",
    "numpy": "2.5.0",
    "pandas": "3.0.3",
    "scipy": "1.18.0",
    "scikit_learn": "1.9.0",
}

# Approved, committed E1 artifact fingerprints. --validate-artifact checks the
# committed random_projection_results.json against these WITHOUT recomputing the
# projection, so it runs anywhere. These describe the LF-terminated artifact
# regenerated under EXPECTED_E1_ENVIRONMENT; the numbers below never change
# without a deliberate, reviewed re-approval.
RESULT_SHA256 = (
    "4adfb78b24f613adf826e7202272bbe5d95fcb9001d2f46d80567f1af319d186"
)
RESULT_SIZE_BYTES = 9993
COMPONENT_SHA256 = (
    "d751f201d221c1b87048f9ef83fd93d91c810a98cbaabe2c9f14dd1c03828c38"
)
RESULT_SPEARMAN_MEAN = 0.2104
RESULT_SELECTED_ALPHA = 3162.0
RESULT_DELTAS = {
    "random_projection_minus_ridge_pca": -0.0252,
    "random_projection_minus_ridge_head": 0.0057,
    "random_projection_minus_lineage_mean": 0.0604,
}
# Softened interpretation (one seed, no paired uncertainty analysis on the
# +0.0057). Stored verbatim in the artifact and checked by --validate-artifact.
RESULT_INTERPRETATION = (
    "This result makes bottleneck width alone an unlikely explanation for "
    "Geneformer's deficit and is consistent with the learned representation "
    "discarding or distorting useful expression signal. Because E1 uses one "
    "projection seed and the +0.0057 difference was not given its own paired "
    "uncertainty analysis, it does not establish a causal mechanism or prove "
    "that random projection is superior."
)

# Exact SHA-256 of every committed input this control consumes. If any no longer
# matches, a frozen Phase 1 artifact has moved and the build hard-fails rather
# than quietly producing a control against a changed baseline.
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
    "data/processed/gene_columns.json":
        "a4b8069cc93af48f01e745bb1a15f4eaf4a7b67c9f92ca44bef3bb9e44c6d0a1",
    "data/processed/baseline_results.json":
        "b49169bd363a596f400b4faff8c21d354275b70404efe08b9109d38f1bdc0ffd",
    "data/processed/head_results.json":
        "1962206fa17646cbd1fec4b642a577cc2586c09c4cabd980541a7e11a8b6f894",
}

# Protected artifacts this experiment must never touch -- verified unchanged by
# --validate.  (baseline_results.json / head_results.json are already gated by
# EXPECTED_INPUT_SHA256 above; these are the ones this module never reads.)
PROTECTED_ARTIFACT_SHA256 = {
    "data/processed/analysis_results.json":
        "12431dad60d07f0bd2bea9a680367007c9e030e9f17c5c20ef0b0694dcb548f9",
    "data/processed/case_study.json": config.CASE_STUDY_JSON_SHA256,
    "phase2_report.html": config.REPORT_HTML_SHA256,
}

_HASH_CHUNK = 1 << 20
_ROUND_DP = 4          # metrics / deltas: match baseline.summarise_metric
_VAR_ROUND_DP = 6      # projected-variance diagnostics


class RandomProjectionControlError(RuntimeError):
    """Raised on any integrity failure in this module (fail closed, never skip)."""


class E1EnvironmentError(RandomProjectionControlError):
    """The runtime library stack differs from EXPECTED_E1_ENVIRONMENT, so a
    real-data recompute cannot reproduce the committed component hash."""


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------

def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(_HASH_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


COMPONENT_SHA256_BYTE_CONVENTION = (
    "hashlib.sha256 over "
    "numpy.ascontiguousarray(components_, dtype='<f8').tobytes(order='C') -- "
    "little-endian IEEE-754 float64, C-order, no .npy header; "
    "768x18460 row-major matrix"
)


def _component_sha256(components: np.ndarray) -> str:
    """
    Deterministic SHA-256 of the projection component matrix.

    Byte convention (documented, stable): little-endian IEEE-754 float64,
    C-order, no ``.npy`` header. ``dtype="<f8"`` forces little-endian
    explicitly so the digest is identical on a big-endian host; the array is
    made C-contiguous and hashed over ``ndarray.tobytes(order="C")``.
    """
    arr = np.ascontiguousarray(np.asarray(components, dtype="<f8"))
    return hashlib.sha256(arr.tobytes(order="C")).hexdigest()


def _reject_json_constant(token: str):
    raise RandomProjectionControlError(
        f"non-standard JSON constant {token!r} is not allowed in "
        f"random_projection_results.json"
    )


def _load_strict_json(path: Path):
    text = Path(path).read_text(encoding="utf-8")
    return json.loads(text, parse_constant=_reject_json_constant)


def _result_text(payload: dict) -> str:
    # allow_nan=False: a NaN / Infinity anywhere in the payload is a hard error
    # at serialisation time, not a silently written non-standard literal.
    return json.dumps(
        payload, indent=2, sort_keys=False, ensure_ascii=False, allow_nan=False
    ) + "\n"


def _result_bytes(payload: dict) -> bytes:
    # Explicit LF bytes. data/processed/** is marked -text in .gitattributes, so
    # exactly these bytes are committed on every platform; LF (not CRLF) keeps a
    # POSIX `git diff --check` from flagging every added line as a trailing CR.
    # (case_study.json keeps its own deliberate CRLF convention -- not this file.)
    return _result_text(payload).encode("utf-8")


def _write_result(payload: dict, path: Path) -> None:
    Path(path).write_bytes(_result_bytes(payload))


def _environment() -> dict:
    import scipy
    import sklearn
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "pandas": pd.__version__,
    }


def _assert_e1_environment() -> None:
    """
    Gate every real-data numerical recompute (--run, full --validate,
    --check-determinism). A fixed seed does NOT make GaussianRandomProjection
    byte-stable across NumPy / scikit-learn versions, so the committed component
    hash only reproduces under EXPECTED_E1_ENVIRONMENT. Fails BEFORE any
    projection or ridge fit. Portable checks (--validate-artifact, --self-test)
    never call this.
    """
    observed = _environment()
    mismatch = [
        k for k in EXPECTED_E1_ENVIRONMENT
        if observed.get(k) != EXPECTED_E1_ENVIRONMENT[k]
    ]
    if not mismatch:
        return
    lines = [
        "E1 numerical regeneration requires the exact pinned environment.",
        "Seeded sklearn GaussianRandomProjection components are NOT byte-stable",
        "across NumPy / scikit-learn versions, so the committed projection hash",
        "(and full --validate) only reproduce on the stack below.",
        "",
        f"  {'component':<14} {'expected':<12} observed",
        f"  {'-' * 14} {'-' * 12} {'-' * 20}",
    ]
    for k in EXPECTED_E1_ENVIRONMENT:
        exp = EXPECTED_E1_ENVIRONMENT[k]
        got = observed.get(k)
        flag = "" if exp == got else "   <- mismatch"
        lines.append(f"  {k:<14} {exp:<12} {got}{flag}")
    lines += [
        "",
        "Fix: install the pinned library stack from requirements-e1.txt on",
        "Python 3.14.6, then retry -- or run `--validate-artifact` for a portable",
        "structural / provenance check that does not recompute the projection.",
    ]
    raise E1EnvironmentError("\n".join(lines))


def _estimator_params(est) -> dict:
    return {
        k: (v if isinstance(v, (int, float, str, bool, type(None))) else str(v))
        for k, v in sorted(est.get_params().items())
    }


def _verify_inputs(processed_dir: Path) -> dict:
    """Gate every committed input against EXPECTED_INPUT_SHA256. Returns rel->hash."""
    root = processed_dir.parent.parent           # data/processed -> repo root
    got: dict[str, str] = {}
    for rel, expected in EXPECTED_INPUT_SHA256.items():
        p = root / rel
        if not p.is_file():
            raise RandomProjectionControlError(f"missing committed input: {rel}")
        actual = _sha256_file(p)
        got[rel] = actual
        if actual != expected:
            raise RandomProjectionControlError(
                f"{rel}: sha256 {actual} != expected {expected}. A frozen "
                f"Phase 1 artifact has changed; do not run E1 against it."
            )
    return got


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


def _patient_group_vector(metadata: pd.DataFrame, index: pd.Index) -> pd.Series:
    """PatientID per line; a missing / 'nan' id falls back to the line's own id
    (same convention as checks.py section 6)."""
    g = metadata[config.GROUP_COL].astype(str).reindex(index)
    g = g.where(g.notna() & (g != "nan"),
                pd.Series(index, index=index).astype(str))
    return g


def _assert_train_eval_disjoint(train_index, eval_index) -> None:
    """Hard-fail if any line is in both the training and the evaluation set."""
    if not set(train_index).isdisjoint(set(eval_index)):
        overlap = sorted(set(train_index) & set(eval_index))
        raise RandomProjectionControlError(
            f"train and val indices overlap -- leakage ({len(overlap)} shared: {overlap[:5]})"
        )


def _groups_crossing_train_val(metadata, train_index, eval_index) -> list[str]:
    used = pd.Index(list(train_index) + list(eval_index))
    split_of = pd.Series(
        ["train"] * len(train_index) + ["val"] * len(eval_index),
        index=used, name="split",
    )
    groups = _patient_group_vector(metadata, used)
    return [str(gp) for gp, sub in split_of.groupby(groups) if sub.nunique() > 1]


# --------------------------------------------------------------------------
# the numeric core (shared by run() and the self-test)
# --------------------------------------------------------------------------

def _fit_project_predict(
    X_train_raw: np.ndarray,
    X_eval_raw: np.ndarray,
    Y_train_raw: np.ndarray,
    Y_eval_raw: np.ndarray,
    train_groups: np.ndarray | None,
    *,
    seed: int,
    alphas: list[float],
    n_components: int,
) -> dict:
    """
    impute(train-mean) -> StandardScaler[train] -> GaussianRandomProjection[train]
    -> train_head.run_ridge_head (its own StandardScaler + Ridge + inner-CV alpha).

    Every fit here sees training rows only; the eval arrays are transform-only.
    Returns the fitted objects, the projected matrices, the predictions and the
    train_head info dict.
    """
    X_train, X_eval = baseline.impute_with_train_mean(X_train_raw, X_eval_raw)
    Y_train_filled, = baseline.impute_with_train_mean(Y_train_raw)

    scaler_expr = StandardScaler()
    Xs_train = scaler_expr.fit_transform(X_train)
    Xs_eval = scaler_expr.transform(X_eval)

    projector = GaussianRandomProjection(
        n_components=n_components, random_state=seed
    )
    proj_train = projector.fit_transform(Xs_train)
    proj_eval = projector.transform(Xs_eval)

    # train_head.run_ridge_head: StandardScaler(fit on proj_train) -> Ridge,
    # alpha chosen by baseline._select_alpha_inner_cv (patient-grouped 5-fold
    # CV inside the training rows) because X_val is None.
    preds, info = train_head.run_ridge_head(
        proj_train, Y_train_filled, proj_eval,
        None, None,
        alphas=list(alphas), seed=seed, train_groups=train_groups,
    )
    return {
        "scaler_expr": scaler_expr,
        "projector": projector,
        "proj_train": proj_train,
        "proj_eval": proj_eval,
        "predictions": preds,
        "info": info,
        "Y_eval_raw": Y_eval_raw,
    }


# --------------------------------------------------------------------------
# run: compute the real val result (does NOT write)
# --------------------------------------------------------------------------

def run(processed_dir: Path | str = config.PROCESSED_DIR,
        *, verbose: bool = True) -> tuple[dict, dict]:
    """
    Compute the E1 val-split result.

    Returns ``(payload, artifacts)`` where ``payload`` is the JSON-ready result
    dict and ``artifacts`` carries the arrays needed for --save-predictions.
    Writing to disk is the caller's job.
    """
    _assert_e1_environment()          # fail before any expensive projection / fit
    processed_dir = Path(processed_dir)
    input_hashes = _verify_inputs(processed_dir)

    expression, crispr, metadata, assignment, selective_genes = _load_common(processed_dir)
    gene_columns = io_utils.load_json(processed_dir / "gene_columns.json")

    prepared = baseline.prepare_task(
        "crispr", expression, crispr, None, selective_genes, assignment, "val"
    )
    if prepared is None:
        raise RandomProjectionControlError("baseline.prepare_task returned None for crispr/val")
    X, Y, train_index, eval_index = prepared

    # ---- hard asserts: shapes, ordering, split discipline --------------------
    canonical = list(gene_columns["canonical_columns"])
    if list(X.columns) != canonical:
        raise RandomProjectionControlError("expression feature order != gene_columns.json canonical order")
    target_names = [str(c) for c in Y.columns]
    if target_names != [g for g in selective_genes if g in crispr.columns]:
        raise RandomProjectionControlError("target order != selective_genes.json (filtered) order")
    if len(train_index) != N_TRAIN_LINES:
        raise RandomProjectionControlError(f"train rows = {len(train_index)} (expected {N_TRAIN_LINES})")
    if len(eval_index) != N_VAL_LINES:
        raise RandomProjectionControlError(f"val rows = {len(eval_index)} (expected {N_VAL_LINES})")
    if Y.shape[1] != N_TARGETS:
        raise RandomProjectionControlError(f"targets = {Y.shape[1]} (expected {N_TARGETS})")
    if X.shape[1] != N_EXPRESSION_FEATURES:
        raise RandomProjectionControlError(
            f"expression features = {X.shape[1]} (expected {N_EXPRESSION_FEATURES})"
        )
    _assert_train_eval_disjoint(train_index, eval_index)

    # every fitted / evaluated row is train or val; never test
    used_split = {i: assignment.get(i) for i in list(train_index) + list(eval_index)}
    if any(used_split[i] != "train" for i in train_index):
        raise RandomProjectionControlError("a non-train row entered the training set")
    if any(used_split[i] != "val" for i in eval_index):
        raise RandomProjectionControlError("a non-val row entered the evaluation set")
    if any(v not in ("train", "val") for v in used_split.values()):
        raise RandomProjectionControlError("a test-split row was selected for E1")

    crossing = _groups_crossing_train_val(metadata, train_index, eval_index)
    if crossing:
        raise RandomProjectionControlError(
            f"{len(crossing)} patient group(s) straddle the train/val boundary: {crossing[:5]}"
        )

    # ---- arrays ------------------------------------------------------------
    X_train_raw = X.loc[train_index].to_numpy(dtype=float)
    X_eval_raw = X.loc[eval_index].to_numpy(dtype=float)
    Y_train_raw = Y.loc[train_index].to_numpy(dtype=float)
    Y_eval_raw = Y.loc[eval_index].to_numpy(dtype=float)

    # patient grouping vector, built exactly as train_head.run_task builds it
    # (the .astype(str) -> "nan" fillna no-op is preserved deliberately; see
    # CLAUDE.md "Known footguns").
    if config.GROUP_COL in metadata.columns:
        train_groups = (
            metadata.loc[train_index, config.GROUP_COL]
            .astype(str)
            .fillna(pd.Series(train_index, index=train_index).astype(str))
            .to_numpy()
        )
    else:
        train_groups = None

    fitted = _fit_project_predict(
        X_train_raw, X_eval_raw, Y_train_raw, Y_eval_raw, train_groups,
        seed=config.RANDOM_SEED, alphas=train_head.HEAD_RIDGE_ALPHAS,
        n_components=N_COMPONENTS,
    )
    projector = fitted["projector"]
    proj_train = fitted["proj_train"]
    proj_eval = fitted["proj_eval"]
    preds = fitted["predictions"]
    info = fitted["info"]

    components = np.ascontiguousarray(np.asarray(projector.components_, dtype=np.float64))
    if components.shape != (N_COMPONENTS, N_EXPRESSION_FEATURES):
        raise RandomProjectionControlError(f"components_ shape {components.shape}")
    if proj_train.shape != (N_TRAIN_LINES, N_COMPONENTS):
        raise RandomProjectionControlError(f"proj_train shape {proj_train.shape}")
    if proj_eval.shape != (N_VAL_LINES, N_COMPONENTS):
        raise RandomProjectionControlError(f"proj_eval shape {proj_eval.shape}")
    if not (np.isfinite(proj_train).all() and np.isfinite(proj_eval).all()):
        raise RandomProjectionControlError("projection produced non-finite values")

    # ---- projected-feature variance diagnostic (pre ridge-head StandardScaler)
    col_var = proj_train.var(axis=0, ddof=1)
    var_mean = float(np.mean(col_var))
    ratio = var_mean / THEORETICAL_PROJECTED_VARIANCE
    consistent = _VARIANCE_CONSISTENCY_BAND[0] <= ratio <= _VARIANCE_CONSISTENCY_BAND[1]
    if not consistent:
        raise RandomProjectionControlError(
            f"projected-feature variance mean {var_mean:.4f} is grossly "
            f"inconsistent with the theoretical scale "
            f"{THEORETICAL_PROJECTED_VARIANCE:.4f} (ratio {ratio:.3f}); "
            f"investigate the pipeline before accepting this result."
        )

    metrics = baseline.evaluate(Y_eval_raw, preds)
    if not isinstance(metrics.get("spearman_mean"), float) or \
            not isinstance(metrics.get("r2_mean"), float):
        raise RandomProjectionControlError("baseline.evaluate produced a non-finite summary metric")

    rp_spearman = metrics["spearman_mean"]

    # ---- reference values: read from disk AND assert equal to the literals ---
    b_obj = io_utils.load_json(processed_dir / "baseline_results.json")
    h_obj = io_utils.load_json(processed_dir / "head_results.json")
    ref_read = {
        "ridge_pca_spearman_mean":
            b_obj["tasks"]["crispr"]["models"]["ridge_pca"]["spearman_mean"],
        "ridge_head_spearman_mean":
            h_obj["tasks"]["crispr"]["models"]["ridge_head"]["spearman_mean"],
        "lineage_mean_spearman_mean":
            b_obj["tasks"]["crispr"]["models"]["lineage_mean"]["spearman_mean"],
    }
    for k, expect in REFERENCE_EXPECTED.items():
        if ref_read[k] != expect:
            raise RandomProjectionControlError(
                f"reference {k} on disk = {ref_read[k]} != expected {expect}"
            )

    selected_alpha = float(info["alpha"])
    grid = [float(a) for a in train_head.HEAD_RIDGE_ALPHAS]
    at_min = selected_alpha == min(grid)
    at_max = selected_alpha == max(grid)
    at_boundary = bool(info.get("alpha_at_grid_boundary", at_min or at_max))

    payload = {
        "schema": SCHEMA_VERSION,
        "experiment": "E1",
        "title": (
            "Random-projection control: a fixed 768-dim Gaussian random "
            "compression of expression through the Geneformer ridge head"
        ),
        "scientific_question": (
            "Does Geneformer pretraining provide more useful predictive "
            "information than a fixed random linear compression of expression "
            "at the same 768-dimensional width?"
        ),
        "status": {
            "kind": "exploratory validation control",
            "split_evaluated": "val",
            "test_split_touched": False,
            "test_evaluation_performed": False,
            "preplanned": (
                "CLAUDE.md section 9.3; covered by the 2026-08-25 Phase 2 scope "
                "decision as already-open work -- not a new material scope change"
            ),
            "val_set_optimism": (
                "This is a sixth model scored on the same 170 val cell lines. "
                "Selecting among several models on 170 lines inflates reported "
                "numbers (CLAUDE.md section 11). The pre-specified headline "
                "comparison (ridge_pca vs ridge_head) was fixed before E1 "
                "existed; E1 is exploratory and is reported regardless of "
                "outcome."
            ),
        },
        "seed": int(config.RANDOM_SEED),
        "pipeline": (
            "impute(train-mean expression) -> StandardScaler[fit on train] -> "
            "GaussianRandomProjection(n_components=768, "
            "random_state=config.RANDOM_SEED)[fit on standardized train] -> "
            "StandardScaler[fit on train, inside train_head.run_ridge_head] -> "
            "multi-output Ridge (alpha by patient-grouped inner 5-fold CV on "
            "the 800 training lines) -> per-target Spearman across the 170 val "
            "cell lines"
        ),
        "reused_components": [
            "baseline.impute_with_train_mean",
            "baseline.evaluate / baseline.per_target_spearman",
            "baseline._select_alpha_inner_cv (via train_head.run_ridge_head)",
            "train_head.run_ridge_head",
            "train_head.HEAD_RIDGE_ALPHAS",
            "baseline.save_prediction_bundle / baseline.verify_prediction_bundle "
            "(only with --save-predictions)",
        ],
        "input_artifact_sha256": input_hashes,
        "counts": {
            "n_train_lines": int(len(train_index)),
            "n_val_lines": int(len(eval_index)),
            "n_test_lines_used": 0,
            "n_targets": int(Y.shape[1]),
            "n_expression_features_before_projection": int(X.shape[1]),
            "n_projected_features": int(components.shape[0]),
            "train_val_indices_disjoint": True,
            "patient_groups_crossing_train_val_boundary": 0,
        },
        "leakage_guards": {
            "imputation_fit_on": "train rows only (baseline.impute_with_train_mean)",
            "standardscaler_expression_fit_on": "train rows only",
            "random_projection_fit_on": (
                "standardized train rows only; GaussianRandomProjection draws "
                "its matrix from the seed and n_features -- no data values "
                "enter the component matrix"
            ),
            "alpha_selection": (
                "patient-grouped 5-fold CV inside the 800 training lines "
                "(X_val=None); no val row is scored during selection"
            ),
            "ridge_fit_on": "train rows only",
            "feature_order_matches": "gene_columns.json canonical_columns (18460, in order)",
            "target_order_matches": (
                "selective_genes.json 'genes' filtered to crispr_effect columns "
                "(4297, in order)"
            ),
            "no_validation_statistic_enters": (
                "preprocessing, projection construction, imputation, alpha "
                "selection, or model fitting"
            ),
        },
        "projection": {
            "class": "sklearn.random_projection.GaussianRandomProjection",
            "params": _estimator_params(projector),
            "random_state": int(config.RANDOM_SEED),
            "n_features_in": int(components.shape[1]),
            "n_components_out": int(components.shape[0]),
            "components_shape": [int(components.shape[0]), int(components.shape[1])],
            "components_dtype": str(components.dtype),
            "components_sha256": _component_sha256(components),
            "components_sha256_byte_convention": COMPONENT_SHA256_BYTE_CONVENTION,
            "distribution": (
                "sklearn draws entries i.i.d. from N(0, 1/n_components) = N(0, 1/768)"
            ),
        },
        "projected_feature_variance_train": {
            "note": (
                "column variance of the 768 projected TRAINING features, before "
                "the ridge-head StandardScaler"
            ),
            "ddof": 1,
            "mean": round(var_mean, _VAR_ROUND_DP),
            "median": round(float(np.median(col_var)), _VAR_ROUND_DP),
            "min": round(float(np.min(col_var)), _VAR_ROUND_DP),
            "max": round(float(np.max(col_var)), _VAR_ROUND_DP),
            "q25": round(float(np.percentile(col_var, 25)), _VAR_ROUND_DP),
            "q75": round(float(np.percentile(col_var, 75)), _VAR_ROUND_DP),
            "theoretical_approx": round(THEORETICAL_PROJECTED_VARIANCE, _VAR_ROUND_DP),
            "theoretical_formula": "n_expression_features / n_components = 18460 / 768",
            "observed_mean_over_theoretical": round(ratio, _VAR_ROUND_DP),
            "consistency_band": list(_VARIANCE_CONSISTENCY_BAND),
            "consistent_with_theoretical_scale": bool(consistent),
        },
        "alpha_grid": grid,
        "alpha_grid_source": "train_head.HEAD_RIDGE_ALPHAS (13 points, 1.0 .. 1e6)",
        "alpha_grid_n_points": len(grid),
        "alpha_selection_method": info.get("alpha_selected_on"),
        "alpha_sweep": info.get("alpha_sweep"),
        "selected_alpha": selected_alpha,
        "selected_alpha_is_grid_min": bool(at_min),
        "selected_alpha_is_grid_max": bool(at_max),
        "selected_alpha_at_grid_boundary": at_boundary,
        "selected_alpha_interior": (not at_boundary),
        "metrics": {"eval_split": "val", **metrics},
        "reference_values": {
            "ridge_pca_spearman_mean": {
                "value": ref_read["ridge_pca_spearman_mean"],
                "source": "data/processed/baseline_results.json::tasks.crispr.models.ridge_pca.spearman_mean",
                "matches_expected_literal": True,
            },
            "ridge_head_spearman_mean": {
                "value": ref_read["ridge_head_spearman_mean"],
                "source": "data/processed/head_results.json::tasks.crispr.models.ridge_head.spearman_mean",
                "matches_expected_literal": True,
            },
            "lineage_mean_spearman_mean": {
                "value": ref_read["lineage_mean_spearman_mean"],
                "source": "data/processed/baseline_results.json::tasks.crispr.models.lineage_mean.spearman_mean",
                "matches_expected_literal": True,
            },
        },
        "deltas": {
            "random_projection_minus_ridge_pca":
                round(rp_spearman - ref_read["ridge_pca_spearman_mean"], _ROUND_DP),
            "random_projection_minus_ridge_head":
                round(rp_spearman - ref_read["ridge_head_spearman_mean"], _ROUND_DP),
            "random_projection_minus_lineage_mean":
                round(rp_spearman - ref_read["lineage_mean_spearman_mean"], _ROUND_DP),
        },
        "interpretation": RESULT_INTERPRETATION,
        "environment": _environment(),
        "determinism": (
            "Byte-identical on rebuild within the recorded environment: "
            "StandardScaler, GaussianRandomProjection (seeded), GroupKFold "
            "(no shuffle) and Ridge (cholesky) are all deterministic. Verified "
            "by --check-determinism; --validate re-runs and byte-compares to "
            "the committed file."
        ),
        "no_test_evaluation": (
            "No test-split expression feature was loaded for inference, no "
            "test-split CRISPR outcome was loaded for evaluation, and no test "
            "prediction, ranking, metric, or performance number was computed or "
            "reported by this module. There is no --split test code path."
        ),
        "limitations": [
            "This is ONE fixed random projection (seed 20260722). A different "
            "seed gives a different projection and could give a different "
            "score; the result is not averaged over projections.",
            "Exploratory validation control, not a final result. It adds a "
            "sixth model scored on the same 170 val lines (val-set optimism, "
            "CLAUDE.md section 11).",
            "No test-split evaluation was performed. The one-time test-split "
            "run happens only after the model set is frozen in writing.",
            "Reproducibility is asymmetric: the baseline arm and this control "
            "reproduce end to end from the public repo; the Geneformer arm's "
            "embeddings do not (Kaggle GPU artifact).",
            "Byte-exact regeneration of the projection needs the pinned stack in "
            "requirements-e1.txt on Python 3.14.6: a fixed seed alone does not "
            "make sklearn's GaussianRandomProjection components identical across "
            "NumPy / scikit-learn versions. Portable structural / provenance "
            "verification is `random_projection.py --validate-artifact`.",
            "The projected-feature spectrum (variance ~24 per column) is a "
            "third spectrum, between the head's unit-variance embeddings and "
            "the baseline's PCA eigenvalues; alpha is selected on its own "
            "13-point grid and is required to be interior (CLAUDE.md invariant 7).",
        ],
    }

    artifacts = {
        "predictions": preds,
        "y_true": Y_eval_raw,
        "eval_index": [str(i) for i in eval_index],
        "target_names": target_names,
        "spearman_mean": rp_spearman,
    }

    if verbose:
        print(f"  random projection ridge  spearman mean : {rp_spearman}")
        print(f"  selected alpha                         : {selected_alpha:g} "
              f"({'INTERIOR' if not at_boundary else 'GRID BOUNDARY'})")
        print(f"  projected-feature variance (train) mean: {var_mean:.4f} "
              f"(theory ~= {THEORETICAL_PROJECTED_VARIANCE:.4f})")
        print(f"  vs ridge_pca  : {payload['deltas']['random_projection_minus_ridge_pca']:+}")
        print(f"  vs ridge_head : {payload['deltas']['random_projection_minus_ridge_head']:+}")
        print(f"  vs lineage    : {payload['deltas']['random_projection_minus_lineage_mean']:+}")

    return payload, artifacts


# --------------------------------------------------------------------------
# optional: persist the held-out prediction bundle (gitignored)
# --------------------------------------------------------------------------

def save_predictions(artifacts: dict,
                     processed_dir: Path | str = config.PROCESSED_DIR) -> bool:
    """
    Write val truth + random_projection_ridge predictions under
    data/processed/predictions/ (gitignored) and re-score them from disk.

    Returns True iff the reloaded matrix re-scores to the recorded Spearman mean.
    """
    processed_dir = Path(processed_dir)
    predictions_out = {
        "crispr": {
            "eval_index": list(artifacts["eval_index"]),
            "target_columns": list(artifacts["target_names"]),
            "y_true": np.asarray(artifacts["y_true"], dtype=float),
            "models": {"random_projection_ridge": np.asarray(artifacts["predictions"], dtype=float)},
        }
    }
    written = baseline.save_prediction_bundle(
        predictions_out, processed_dir, "random_projection", "val"
    )
    print(f"\n  {len(written)} file(s) written to {processed_dir / 'predictions'}")
    for p in written:
        print(f"    {p}")

    results_for_verify = {
        "tasks": {"crispr": {"models": {
            "random_projection_ridge": {"spearman_mean": artifacts["spearman_mean"]}
        }}}
    }
    print("\n  Round-trip check -- re-scoring the saved matrix from disk:")
    ok = baseline.verify_prediction_bundle(
        predictions_out, results_for_verify, processed_dir,
        "random_projection", "val",
    )
    print("\n  " + ("Saved matrix re-scores to the reported value."
                    if ok else "*** MISMATCH: saved matrix does not re-score."))
    return ok


# --------------------------------------------------------------------------
# validate-artifact: portable structural / provenance checks, NO recompute
# --------------------------------------------------------------------------

def validate_artifact(processed_dir: Path | str = config.PROCESSED_DIR,
                      result_file: Path | str = RESULT_FILE,
                      *, verbose: bool = True) -> dict:
    """
    Validate the committed random_projection_results.json without recomputing
    the projection or the ridge fit -- so it runs on ANY environment, not just
    EXPECTED_E1_ENVIRONMENT. Structural, provenance and pinned-fingerprint
    checks only. The full numerical regeneration + byte-compare is ``validate``.
    """
    processed_dir = Path(processed_dir)
    result_file = Path(result_file)
    root = processed_dir.parent.parent
    checks: list[tuple[str, bool, str]] = []

    def _c(name: str, ok: bool, detail: str = "") -> None:
        checks.append((name, bool(ok), detail))
        if verbose:
            print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))

    def _done() -> dict:
        n_fail = sum(1 for _, ok, _ in checks if not ok)
        return {"checks": checks, "n_fail": n_fail, "n_pass": len(checks) - n_fail}

    # 1. present
    if not result_file.is_file():
        _c("random_projection_results.json present", False, f"missing: {result_file}")
        return _done()
    _c("random_projection_results.json present", True)

    # 2. strict JSON
    try:
        committed = _load_strict_json(result_file)
        _c("parses under strict JSON (no NaN / Infinity)", True)
    except RandomProjectionControlError as exc:
        _c("parses under strict JSON (no NaN / Infinity)", False, str(exc))
        return _done()

    # 3. schema
    _c("schema == random-projection-control/1",
       committed.get("schema") == SCHEMA_VERSION, f"schema={committed.get('schema')!r}")

    # 4. committed result SHA-256 + size (pinned; LF artifact)
    raw = result_file.read_bytes()
    got_sha = hashlib.sha256(raw).hexdigest()
    _c("committed result SHA-256 == RESULT_SHA256",
       got_sha == RESULT_SHA256, f"{got_sha} vs pinned {RESULT_SHA256}")
    _c("committed result size == RESULT_SIZE_BYTES",
       len(raw) == RESULT_SIZE_BYTES, f"{len(raw)} vs pinned {RESULT_SIZE_BYTES}")
    _c("committed result is LF-terminated (no CRLF)",
       b"\r\n" not in raw, "" if b"\r\n" not in raw else "file contains CRLF bytes")

    # 5. all committed input hashes still match
    try:
        _verify_inputs(processed_dir)
        _c("committed Phase 1 inputs unchanged (EXPECTED_INPUT_SHA256)", True)
    except RandomProjectionControlError as exc:
        _c("committed Phase 1 inputs unchanged (EXPECTED_INPUT_SHA256)", False, str(exc))
    ia = committed.get("input_artifact_sha256", {})
    _c("recorded input_artifact_sha256 == EXPECTED_INPUT_SHA256",
       ia == EXPECTED_INPUT_SHA256, f"{sorted(set(ia) ^ set(EXPECTED_INPUT_SHA256))}")

    # 6. protected Phase 1 + Phase 2 artifacts unchanged
    bad = []
    for rel, want in PROTECTED_ARTIFACT_SHA256.items():
        p = root / rel
        if not p.is_file():
            bad.append(f"{rel}: missing")
        elif _sha256_file(p) != want:
            bad.append(f"{rel}: changed")
    _c("protected artifacts unchanged (analysis_results.json, case_study.json, phase2_report.html)",
       not bad, f"{bad}")

    # 7. recorded environment == the required E1 environment
    _c("recorded environment == EXPECTED_E1_ENVIRONMENT",
       committed.get("environment") == EXPECTED_E1_ENVIRONMENT,
       f"{committed.get('environment')}")

    # 8. counts + split labels (split labels re-derived from splits.json, cheap)
    c = committed.get("counts", {})
    counts_ok = (c.get("n_train_lines") == 800 and c.get("n_val_lines") == 170
                 and c.get("n_test_lines_used") == 0 and c.get("n_targets") == 4297
                 and c.get("n_expression_features_before_projection") == 18460
                 and c.get("n_projected_features") == 768
                 and c.get("train_val_indices_disjoint") is True
                 and c.get("patient_groups_crossing_train_val_boundary") == 0)
    _c("counts: 800/170/0 lines, 4297 targets, 18460->768 features, disjoint, 0 crossing",
       counts_ok, f"{c}")
    try:
        assign = _load_strict_json(processed_dir / "splits.json")["assignment"]
        sizes = {s: sum(1 for v in assign.values() if v == s) for s in ("train", "val", "test")}
        _c("splits.json labels: 800 train / 170 val / 170 test",
           sizes == {"train": 800, "val": 170, "test": 170}, f"{sizes}")
    except Exception as exc:  # noqa: BLE001
        _c("splits.json labels: 800 train / 170 val / 170 test", False, f"{type(exc).__name__}: {exc}")

    # 9. seed + projection dimensions
    pj = committed.get("projection", {})
    pj_dims = {k: pj.get(k) for k in
               ("random_state", "n_components_out", "n_features_in", "components_shape")}
    _c("seed == config.RANDOM_SEED (20260722); projection dims 768 x 18460",
       (committed.get("seed") == config.RANDOM_SEED == 20260722
        and pj.get("random_state") == 20260722
        and pj.get("n_components_out") == 768
        and pj.get("n_features_in") == 18460
        and pj.get("components_shape") == [768, 18460]),
       f"seed={committed.get('seed')} {pj_dims}")

    # 10. committed projection-component hash == pinned COMPONENT_SHA256
    _c("projection components_sha256 == COMPONENT_SHA256 (pinned)",
       pj.get("components_sha256") == COMPONENT_SHA256,
       f"{pj.get('components_sha256')} vs pinned {COMPONENT_SHA256}")
    _c("component-hash byte convention recorded",
       pj.get("components_sha256_byte_convention") == COMPONENT_SHA256_BYTE_CONVENTION)

    # 11. full 13-point alpha grid + sweep
    _c("alpha_grid == train_head.HEAD_RIDGE_ALPHAS (13 points)",
       committed.get("alpha_grid") == [float(a) for a in train_head.HEAD_RIDGE_ALPHAS]
       and committed.get("alpha_grid_n_points") == 13)
    sweep = committed.get("alpha_sweep")
    _c("alpha_sweep is a dict of all 13 grid points",
       isinstance(sweep, dict) and len(sweep) == 13
       and sorted(float(k) for k in sweep) == sorted(float(a) for a in train_head.HEAD_RIDGE_ALPHAS),
       f"{len(sweep) if isinstance(sweep, dict) else sweep}")

    # 12. selected alpha == 3162 and interior
    _c("selected_alpha == 3162.0 and INTERIOR",
       (committed.get("selected_alpha") == RESULT_SELECTED_ALPHA == 3162.0
        and committed.get("selected_alpha_interior") is True
        and committed.get("selected_alpha_at_grid_boundary") is False
        and committed.get("selected_alpha_is_grid_min") is False
        and committed.get("selected_alpha_is_grid_max") is False),
       f"selected_alpha={committed.get('selected_alpha')}")

    # 13. metrics + reference values
    m = committed.get("metrics", {})
    _c("metrics.spearman_mean == 0.2104 (pinned); headline metrics finite",
       (m.get("spearman_mean") == RESULT_SPEARMAN_MEAN == 0.2104
        and all(isinstance(m.get(k), (int, float)) and np.isfinite(m.get(k))
                for k in ("spearman_median", "spearman_q25", "spearman_q75", "r2_mean"))),
       f"spearman_mean={m.get('spearman_mean')}")
    rv = committed.get("reference_values", {})
    b_obj = io_utils.load_json(processed_dir / "baseline_results.json")
    h_obj = io_utils.load_json(processed_dir / "head_results.json")
    live = {
        "ridge_pca_spearman_mean": b_obj["tasks"]["crispr"]["models"]["ridge_pca"]["spearman_mean"],
        "ridge_head_spearman_mean": h_obj["tasks"]["crispr"]["models"]["ridge_head"]["spearman_mean"],
        "lineage_mean_spearman_mean": b_obj["tasks"]["crispr"]["models"]["lineage_mean"]["spearman_mean"],
    }
    rv_values = {k: rv.get(k, {}).get("value") for k in REFERENCE_EXPECTED}
    _c("reference values == expected literals == current baseline/head results",
       all(rv_values[k] == REFERENCE_EXPECTED[k] == live[k] for k in REFERENCE_EXPECTED),
       f"committed={rv_values} live={live}")

    # 14. deltas == pinned, and == spearman_mean minus each reference
    recomputed = {
        "random_projection_minus_ridge_pca": round(m.get("spearman_mean", float("nan")) - REFERENCE_EXPECTED["ridge_pca_spearman_mean"], _ROUND_DP),
        "random_projection_minus_ridge_head": round(m.get("spearman_mean", float("nan")) - REFERENCE_EXPECTED["ridge_head_spearman_mean"], _ROUND_DP),
        "random_projection_minus_lineage_mean": round(m.get("spearman_mean", float("nan")) - REFERENCE_EXPECTED["lineage_mean_spearman_mean"], _ROUND_DP),
    }
    _c("deltas == pinned RESULT_DELTAS == spearman_mean minus each reference",
       committed.get("deltas") == RESULT_DELTAS == recomputed,
       f"committed={committed.get('deltas')} pinned={RESULT_DELTAS} recomputed={recomputed}")

    # 15. no-test-evaluation declarations
    st = committed.get("status", {})
    nte = committed.get("no_test_evaluation", "")
    _c("no-test-evaluation declarations present and val-only",
       (st.get("test_split_touched") is False
        and st.get("test_evaluation_performed") is False
        and st.get("split_evaluated") == "val"
        and isinstance(nte, str)
        and "no --split test" in nte.lower()),
       f"status={st} no_test_evaluation={nte[:60]!r}")

    # 16. limitations + exploratory status + softened interpretation
    _c("status.kind == 'exploratory validation control'; limitations non-empty list",
       (st.get("kind") == "exploratory validation control"
        and isinstance(committed.get("limitations"), list)
        and len(committed["limitations"]) >= 3))
    _c("interpretation == softened RESULT_INTERPRETATION (no 'rules out' claim)",
       committed.get("interpretation") == RESULT_INTERPRETATION
       and "rules out" not in committed.get("interpretation", "").lower())

    return _done()


# --------------------------------------------------------------------------
# validate: re-run, byte-compare, fail closed on the full list
#           (requires EXPECTED_E1_ENVIRONMENT -- calls run())
# --------------------------------------------------------------------------

def validate(processed_dir: Path | str = config.PROCESSED_DIR,
             result_file: Path | str = RESULT_FILE,
             *, verbose: bool = True) -> dict:
    processed_dir = Path(processed_dir)
    result_file = Path(result_file)
    checks: list[tuple[str, bool, str]] = []

    def _c(name: str, ok: bool, detail: str = "") -> None:
        checks.append((name, bool(ok), detail))
        if verbose:
            print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))

    # 1. committed inputs unchanged
    try:
        _verify_inputs(processed_dir)
        _c("committed Phase 1 inputs unchanged (EXPECTED_INPUT_SHA256)", True)
    except RandomProjectionControlError as exc:
        _c("committed Phase 1 inputs unchanged (EXPECTED_INPUT_SHA256)", False, str(exc))
        return {"checks": checks, "n_fail": sum(1 for _, ok, _ in checks if not ok),
                "n_pass": sum(1 for _, ok, _ in checks if ok)}

    # 2. result file present + strict JSON (no NaN / Infinity)
    if not result_file.is_file():
        _c("random_projection_results.json present", False, f"missing: {result_file}")
        return {"checks": checks, "n_fail": sum(1 for _, ok, _ in checks if not ok),
                "n_pass": sum(1 for _, ok, _ in checks if ok)}
    _c("random_projection_results.json present", True)
    try:
        committed = _load_strict_json(result_file)
        _c("result file parses under strict JSON (no NaN / Infinity)", True)
    except RandomProjectionControlError as exc:
        _c("result file parses under strict JSON (no NaN / Infinity)", False, str(exc))
        return {"checks": checks, "n_fail": sum(1 for _, ok, _ in checks if not ok),
                "n_pass": sum(1 for _, ok, _ in checks if ok)}

    # 3. schema
    _c("schema == random-projection-control/1",
       committed.get("schema") == SCHEMA_VERSION, f"schema={committed.get('schema')!r}")

    # 4. recompute the whole result and byte-compare to the committed file
    fresh_payload, _artifacts = run(processed_dir, verbose=False)
    same_bytes = _result_bytes(fresh_payload) == result_file.read_bytes()
    _c("committed result is byte-identical to a fresh recompute", same_bytes,
       "" if same_bytes else "re-run produced different bytes")

    # 5. split / sample / feature / target counts
    c = committed.get("counts", {})
    _c("counts: 800 train, 170 val, 0 test, 4297 targets, 18460 features, 768 projected",
       (c.get("n_train_lines") == 800 and c.get("n_val_lines") == 170
        and c.get("n_test_lines_used") == 0 and c.get("n_targets") == 4297
        and c.get("n_expression_features_before_projection") == 18460
        and c.get("n_projected_features") == 768),
       f"{c}")
    _c("train/val indices disjoint; no patient group crosses the boundary",
       (c.get("train_val_indices_disjoint") is True
        and c.get("patient_groups_crossing_train_val_boundary") == 0))

    # 6. projection: seed / dimensions / component hash
    pj = committed.get("projection", {})
    _c("projection seed == config.RANDOM_SEED and dims 768 x 18460",
       (pj.get("random_state") == config.RANDOM_SEED
        and pj.get("n_components_out") == 768
        and pj.get("n_features_in") == 18460
        and pj.get("components_shape") == [768, 18460]),
       f"{ {k: pj.get(k) for k in ('random_state','n_components_out','n_features_in','components_shape')} }")
    fresh_hash = fresh_payload["projection"]["components_sha256"]
    _c("projection component SHA-256 matches a fresh projection",
       pj.get("components_sha256") == fresh_hash,
       "" if pj.get("components_sha256") == fresh_hash
       else f"committed {pj.get('components_sha256')} != fresh {fresh_hash}")

    # 7. alpha grid + interior alpha
    _c("alpha grid == train_head.HEAD_RIDGE_ALPHAS (13 points)",
       committed.get("alpha_grid") == [float(a) for a in train_head.HEAD_RIDGE_ALPHAS]
       and committed.get("alpha_grid_n_points") == 13)
    _c("alpha sweep has all 13 grid points",
       isinstance(committed.get("alpha_sweep"), dict)
       and len(committed["alpha_sweep"]) == 13,
       f"{len(committed.get('alpha_sweep') or {})} points")
    interior = (committed.get("selected_alpha_interior") is True
                and committed.get("selected_alpha_at_grid_boundary") is False
                and committed.get("selected_alpha_is_grid_min") is False
                and committed.get("selected_alpha_is_grid_max") is False)
    _c("selected alpha is INTERIOR (not grid min or max)", interior,
       f"selected_alpha={committed.get('selected_alpha')}")

    # 8. metrics finite
    m = committed.get("metrics", {})
    finite_ok = all(
        isinstance(m.get(k), (int, float)) and np.isfinite(m.get(k))
        for k in ("spearman_mean", "spearman_median", "r2_mean")
    )
    _c("headline metrics are finite numbers", finite_ok, f"spearman_mean={m.get('spearman_mean')}")

    # 9. projected-variance diagnostic consistent with 18460/768
    v = committed.get("projected_feature_variance_train", {})
    _c("projected-feature variance mean is consistent with the theoretical scale",
       v.get("consistent_with_theoretical_scale") is True
       and abs(v.get("theoretical_approx", 0) - N_EXPRESSION_FEATURES / N_COMPONENTS) < 1e-6,
       f"mean={v.get('mean')} theory={v.get('theoretical_approx')} ratio={v.get('observed_mean_over_theoretical')}")

    # 10. reference values: match the literals AND the current results files
    b_obj = io_utils.load_json(processed_dir / "baseline_results.json")
    h_obj = io_utils.load_json(processed_dir / "head_results.json")
    live = {
        "ridge_pca_spearman_mean": b_obj["tasks"]["crispr"]["models"]["ridge_pca"]["spearman_mean"],
        "ridge_head_spearman_mean": h_obj["tasks"]["crispr"]["models"]["ridge_head"]["spearman_mean"],
        "lineage_mean_spearman_mean": b_obj["tasks"]["crispr"]["models"]["lineage_mean"]["spearman_mean"],
    }
    rv = committed.get("reference_values", {})
    ref_ok = all(
        rv.get(k, {}).get("value") == REFERENCE_EXPECTED[k] == live[k]
        for k in REFERENCE_EXPECTED
    )
    rv_values = {k: rv.get(k, {}).get("value") for k in REFERENCE_EXPECTED}
    _c("reference values == expected literals == current baseline/head results",
       ref_ok, f"committed={rv_values} live={live}")

    # 11. deltas recompute correctly from the recorded metric + references
    rp = m.get("spearman_mean")
    exp_deltas = {
        "random_projection_minus_ridge_pca": round(rp - REFERENCE_EXPECTED["ridge_pca_spearman_mean"], _ROUND_DP),
        "random_projection_minus_ridge_head": round(rp - REFERENCE_EXPECTED["ridge_head_spearman_mean"], _ROUND_DP),
        "random_projection_minus_lineage_mean": round(rp - REFERENCE_EXPECTED["lineage_mean_spearman_mean"], _ROUND_DP),
    }
    _c("deltas == random_projection spearman_mean minus each reference",
       committed.get("deltas") == exp_deltas, f"committed={committed.get('deltas')} expected={exp_deltas}")

    # 12. no evidence of test evaluation
    st = committed.get("status", {})
    _c("result asserts no test-split evaluation",
       st.get("test_split_touched") is False
       and st.get("test_evaluation_performed") is False
       and isinstance(committed.get("no_test_evaluation"), str)
       and st.get("split_evaluated") == "val")

    # 13. protected artifacts unchanged
    root = processed_dir.parent.parent
    bad = []
    for rel, want in PROTECTED_ARTIFACT_SHA256.items():
        p = root / rel
        if not p.is_file():
            bad.append(f"{rel}: missing")
        elif _sha256_file(p) != want:
            bad.append(f"{rel}: changed")
    _c("protected artifacts unchanged (analysis_results.json, case_study.json, phase2_report.html)",
       not bad, f"{bad}")

    n_fail = sum(1 for _, ok, _ in checks if not ok)
    return {"checks": checks, "n_fail": n_fail, "n_pass": len(checks) - n_fail}


# --------------------------------------------------------------------------
# determinism
# --------------------------------------------------------------------------

def check_determinism(processed_dir: Path | str = config.PROCESSED_DIR,
                      *, verbose: bool = True) -> dict:
    processed_dir = Path(processed_dir)
    tmp = Path(tempfile.mkdtemp(prefix="rp_det_"))
    try:
        p1, _ = run(processed_dir, verbose=False)
        p2, _ = run(processed_dir, verbose=False)
        b1, b2 = _result_bytes(p1), _result_bytes(p2)
        (tmp / "a.json").write_bytes(b1)
        (tmp / "b.json").write_bytes(b2)
        ok = (tmp / "a.json").read_bytes() == (tmp / "b.json").read_bytes()
        if verbose:
            print(f"  determinism: {'byte-identical' if ok else 'MISMATCH'} "
                  f"({len(b1)} vs {len(b2)} bytes)")
        return {"ok": ok, "n_bytes": len(b1)}
    finally:
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------
# self-test (synthetic, offline)
# --------------------------------------------------------------------------

def _self_test() -> int:
    print("Running random_projection.py self-test...")
    import shutil

    rng = np.random.default_rng(0)
    n_tr, n_va, p, t, k = 60, 20, 200, 8, 16
    groups = np.repeat(np.arange(12), 5)[:n_tr]     # 12 patient groups over 60 rows

    X_tr = rng.normal(size=(n_tr, p))
    X_va = rng.normal(size=(n_va, p))
    Y_tr = rng.normal(size=(n_tr, t))
    Y_va = rng.normal(size=(n_va, t))
    alphas = list(train_head.HEAD_RIDGE_ALPHAS)

    a = _fit_project_predict(X_tr, X_va, Y_tr, Y_va, groups,
                             seed=123, alphas=alphas, n_components=k)
    b = _fit_project_predict(X_tr, X_va, Y_tr, Y_va, groups,
                             seed=123, alphas=alphas, n_components=k)
    assert np.array_equal(a["projector"].components_, b["projector"].components_), \
        "same seed produced different projection components"
    assert np.array_equal(a["predictions"], b["predictions"]), \
        "same seed produced different predictions"
    print("  [ok] same seed -> identical projection components and predictions")

    c = _fit_project_predict(X_tr, X_va, Y_tr, Y_va, groups,
                             seed=999, alphas=alphas, n_components=k)
    assert not np.array_equal(a["projector"].components_, c["projector"].components_), \
        "different seed produced the same projection"
    print("  [ok] different seed -> different projection")

    assert a["proj_train"].shape == (n_tr, k) and a["proj_eval"].shape == (n_va, k), \
        f"projection shape wrong: {a['proj_train'].shape} / {a['proj_eval'].shape}"
    print("  [ok] projection shape is (n_rows, n_components)")

    # preprocessing fit on training rows only: perturb the eval arrays wildly and
    # confirm the fitted scaler, the projection, and the selected alpha are all
    # unchanged.
    X_va_bad = X_va + 1e6
    Y_va_bad = Y_va + 1e6
    d = _fit_project_predict(X_tr, X_va_bad, Y_tr, Y_va_bad, groups,
                             seed=123, alphas=alphas, n_components=k)
    assert np.array_equal(a["scaler_expr"].mean_, d["scaler_expr"].mean_), "scaler mean moved with eval data"
    assert np.array_equal(a["scaler_expr"].scale_, d["scaler_expr"].scale_), "scaler scale moved with eval data"
    assert np.array_equal(a["projector"].components_, d["projector"].components_), "projection moved with eval data"
    assert np.array_equal(a["proj_train"], d["proj_train"]), "projected TRAIN features moved with eval data"
    assert a["info"]["alpha"] == d["info"]["alpha"], "selected alpha moved with eval data"
    assert a["info"]["alpha_sweep"] == d["info"]["alpha_sweep"], "alpha sweep moved with eval data"
    print("  [ok] preprocessing + alpha selection depend on training rows only "
          "(validation perturbation changes nothing fitted)")

    # train/val overlap must hard-fail; disjoint indices must pass
    idx = [f"L{i}" for i in range(10)]
    _assert_train_eval_disjoint(idx[:6], idx[6:])          # disjoint -> no raise
    raised = False
    try:
        _assert_train_eval_disjoint(idx[:7], idx[5:])      # 2 shared -> raise
    except RandomProjectionControlError:
        raised = True
    assert raised, "_assert_train_eval_disjoint did not flag an overlap"
    print("  [ok] train/val overlap is a hard failure")

    # no --split test path exists
    parser = _build_parser()
    assert not any(getattr(act, "dest", None) == "split" for act in parser._actions), \
        "a --split option exists; E1 must not offer a test path"
    import contextlib
    import io as _io
    try:
        with contextlib.redirect_stderr(_io.StringIO()):
            parser.parse_args(["--split", "test"])
        raise AssertionError("--split test was accepted")
    except SystemExit:
        pass
    print("  [ok] there is no --split test code path (argparse rejects it)")

    # environment gate: silent in the pinned env, hard-fails on a version mismatch
    _assert_e1_environment()                         # we ARE in EXPECTED_E1_ENVIRONMENT
    _mod = sys.modules[__name__]
    _saved_env = _mod.EXPECTED_E1_ENVIRONMENT
    try:
        _mod.EXPECTED_E1_ENVIRONMENT = {**_saved_env, "numpy": "0.0.0-not-real"}
        raised = False
        try:
            _assert_e1_environment()
        except E1EnvironmentError:
            raised = True
        assert raised, "_assert_e1_environment did not flag a version mismatch"
    finally:
        _mod.EXPECTED_E1_ENVIRONMENT = _saved_env
    print("  [ok] environment gate: silent in the pinned env, raises on a mismatch")

    # LF (not CRLF) result bytes, so a POSIX `git diff --check` stays quiet
    assert b"\r\n" not in _result_bytes({"k": 1}), "_result_bytes still emits CRLF"
    print("  [ok] result writer emits LF line endings")

    # strict JSON rejects NaN / Infinity, both on write and on read
    for bad_val in (float("nan"), float("inf"), float("-inf")):
        try:
            _result_bytes({"x": bad_val})
            raise AssertionError(f"_result_bytes accepted {bad_val}")
        except ValueError:
            pass
    tmp = Path(tempfile.mkdtemp(prefix="rp_selftest_"))
    try:
        bad_path = tmp / "bad.json"
        bad_path.write_text('{"x": NaN}', encoding="utf-8")
        raised = False
        try:
            _load_strict_json(bad_path)
        except RandomProjectionControlError:
            raised = True
        assert raised, "_load_strict_json accepted NaN"
        print("  [ok] strict JSON rejects NaN / Infinity on write and on read")

        # prediction round-trip preserves row + target order
        eval_index = [f"ACH-{i:06d}" for i in range(n_va)]
        target_names = [f"G{i} ({i})" for i in range(t)]
        artifacts = {
            "predictions": a["predictions"], "y_true": Y_va,
            "eval_index": eval_index, "target_names": target_names,
            "spearman_mean": baseline.evaluate(Y_va, a["predictions"])["spearman_mean"],
        }
        predictions_out = {"crispr": {
            "eval_index": eval_index, "target_columns": target_names,
            "y_true": Y_va, "models": {"random_projection_ridge": a["predictions"]},
        }}
        baseline.save_prediction_bundle(predictions_out, tmp, "random_projection", "val")
        truth = io_utils.load_matrix(tmp / "predictions" / "random_projection_crispr_val_y_true")
        reloaded = io_utils.load_matrix(
            tmp / "predictions" / "random_projection_crispr_val_random_projection_ridge")
        assert [str(i) for i in truth.index] == eval_index, "row order not preserved"
        assert [str(c) for c in reloaded.columns] == target_names, "target order not preserved"
        assert np.allclose(reloaded.to_numpy(dtype=float), a["predictions"], rtol=0, atol=0), \
            "prediction values not preserved on round-trip"
        rescored = baseline.evaluate(truth.to_numpy(dtype=float),
                                     reloaded.to_numpy(dtype=float))["spearman_mean"]
        assert abs(rescored - artifacts["spearman_mean"]) < 1e-9, "round-trip Spearman drifted"
        print("  [ok] prediction bundle round-trips with row + target order intact")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print("\nSelf-test passed.")
    return 0


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="E1 random-projection control (validation-only, exploratory)."
    )
    ap.add_argument("--run", action="store_true",
                    help="Compute the val-split result and write "
                         "data/processed/random_projection_results.json.")
    ap.add_argument("--save-predictions", action="store_true",
                    help="With --run: also write the val truth + "
                         "random_projection_ridge prediction matrix under "
                         "data/processed/predictions/ (gitignored) and re-score "
                         "them from disk.")
    ap.add_argument("--validate-artifact", action="store_true",
                    help="Portable structural / provenance validation of the "
                         "committed result: strict JSON, schema, pinned result + "
                         "component hashes, input + protected-artifact hashes, "
                         "recorded environment, counts, split labels, alpha grid "
                         "/ sweep, interior alpha, metrics, references, deltas, "
                         "no-test-evaluation, softened interpretation. Does NOT "
                         "recompute the projection -- runs on any environment.")
    ap.add_argument("--validate", action="store_true",
                    help="Full numerical regeneration: re-run the projection + "
                         "ridge, byte-compare to the committed result, fail "
                         "closed on any integrity problem. Requires the pinned "
                         "EXPECTED_E1_ENVIRONMENT (see requirements-e1.txt).")
    ap.add_argument("--check-determinism", action="store_true",
                    help="Compute the real result twice into temp files; "
                         "require byte-identical JSON. Does not touch the "
                         "committed result.")
    ap.add_argument("--self-test", action="store_true",
                    help="Synthetic, offline end-to-end check.")
    ap.add_argument("--out-file", default=str(RESULT_FILE),
                    help="Override the result path (testing only).")
    return ap


def main(argv: list[str] | None = None) -> int:
    ap = _build_parser()
    args = ap.parse_args(argv)

    if args.self_test:
        return _self_test()

    result_file = Path(args.out_file)
    did_something = False
    rc = 0

    # ---- portable: no recompute, no environment gate ----------------------
    if args.validate_artifact:
        did_something = True
        print("=" * 74)
        print("VALIDATE-ARTIFACT  (portable: structure + provenance, no recompute)")
        print("=" * 74)
        rep = validate_artifact(config.PROCESSED_DIR, result_file)
        print()
        print(f"  {rep['n_pass']}/{rep['n_pass'] + rep['n_fail']} checks passed"
              + ("" if rep["n_fail"] == 0 else f"  ({rep['n_fail']} FAILED)"))
        rc = rc or (0 if rep["n_fail"] == 0 else 1)

    # ---- everything below recomputes real data: gated on the pinned env ---
    try:
        if args.run:
            did_something = True
            print("=" * 74)
            print("E1  RANDOM-PROJECTION CONTROL  (val split, exploratory)")
            print("=" * 74)
            payload, artifacts = run(config.PROCESSED_DIR)
            _write_result(payload, result_file)
            print(f"\n  Written to: {result_file}")

            if args.save_predictions:
                print(f"\n{'=' * 74}")
                print("SAVING HELD-OUT PREDICTIONS (gitignored)")
                print("=" * 74)
                if not save_predictions(artifacts):
                    rc = rc or 1

            if payload["selected_alpha_at_grid_boundary"]:
                print("\n" + "!" * 74)
                print("  BOUNDARY ALPHA: the inner CV selected the grid "
                      f"{'minimum' if payload['selected_alpha_is_grid_min'] else 'maximum'} "
                      f"({payload['selected_alpha']:g}).")
                print("  Full alpha sweep:")
                for k_alpha, v_alpha in payload["alpha_sweep"].items():
                    print(f"    alpha={k_alpha:>12}  mean inner-CV spearman = {v_alpha}")
                print("  Per CLAUDE.md section 9.3 / invariant 7: NOT widening the "
                      "grid, NOT changing the seed, NOT committing this as a final")
                print("  E1 result. The grid decision is a separate, approved step.")
                print("!" * 74)
                rc = rc or 1
            else:
                print(f"\n  Selected alpha {payload['selected_alpha']:g} is interior. "
                      f"E1 result is on record.")

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
            print("VALIDATE random_projection_results.json  (full numerical regen)")
            print("=" * 74)
            rep = validate(config.PROCESSED_DIR, result_file)
            print()
            print(f"  {rep['n_pass']}/{rep['n_pass'] + rep['n_fail']} checks passed"
                  + ("" if rep["n_fail"] == 0 else f"  ({rep['n_fail']} FAILED)"))
            rc = rc or (0 if rep["n_fail"] == 0 else 1)

    except E1EnvironmentError as exc:
        print("\n" + "!" * 74)
        print(str(exc))
        print("!" * 74)
        return 1

    if not did_something:
        ap.print_help()
        return 2
    return rc


if __name__ == "__main__":
    sys.exit(main())
