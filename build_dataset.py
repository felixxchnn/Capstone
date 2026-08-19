"""
build_dataset.py
================
Builds the aligned modelling dataset from raw DepMap files.

Run this once, before anything else:

    python build_dataset.py

Outputs (written to config.PROCESSED_DIR)
-----------------------------------------
    expression.parquet        ModelID x canonical genes (log2 TPM + 1)
    crispr_effect.parquet     ModelID x canonical genes (Chronos gene effect)
    prism_response.parquet    ModelID x compounds (optional; only if PRISM present)
    model_metadata.parquet    ModelID x curated metadata columns
    gene_columns.json         canonical gene order + Entrez IDs  <-- the socket
    gene_id_map.parquet       entrez / symbol / per-source labels
    selective_genes.json      CRISPR targets worth predicting
    join_report.txt           human-readable audit of every filtering step
    join_report.json          the same numbers, machine-readable

The two structural traps this module exists to defuse
-----------------------------------------------------
1. The expression matrix is *not* indexed by ModelID. Its first columns are
   sequencing metadata. Slicing it naively drags string columns into the
   feature matrix.

2. The expression matrix contains multiple sequencing profiles per cell line.
   Without filtering on `IsDefaultEntryForModel`, ModelID is not unique, and
   merging on a duplicated key multiplies rows. The dataset then looks larger
   than it is, near-identical rows straddle the train/test boundary, and the
   reported score comes out inflated. No error is raised at any point. This is
   the single most likely way for this project to produce a confident wrong
   answer, so the filter is applied here and uniqueness is asserted after.
"""

from __future__ import annotations

import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd

import config
import gene_ids
import io_utils


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------

def _log(message: str = "") -> None:
    print(message, flush=True)


def _rule(title: str = "") -> None:
    if title:
        _log(f"\n{'=' * 74}\n{title}\n{'=' * 74}")
    else:
        _log("-" * 74)


def _drop_unnamed(df: pd.DataFrame) -> pd.DataFrame:
    """Remove pandas' auto-generated 'Unnamed: N' columns."""
    unnamed = [c for c in df.columns if str(c).startswith("Unnamed:")]
    return df.drop(columns=unnamed) if unnamed else df


def _coerce_bool(series: pd.Series) -> pd.Series:
    """
    Coerce a column to boolean regardless of whether pandas read it as bool,
    as the strings 'True'/'False', or as 1/0.
    """
    if series.dtype == bool:
        return series
    as_str = series.astype(str).str.strip().str.lower()
    return as_str.isin({"true", "1", "1.0", "yes", "t"})


# --------------------------------------------------------------------------
# Loaders
# --------------------------------------------------------------------------

def load_expression() -> tuple[pd.DataFrame, dict]:
    """
    Load the expression matrix, filtered to one canonical profile per cell line.

    Returns
    -------
    (frame, report)
        `frame` is indexed by ModelID, columns are gene labels only.
    """
    path = config.resolve_file("expression")
    _log(f"Reading expression matrix: {path.name}")
    _log("  (this is a large file; expect this step to take a minute or two)")

    frame = pd.read_csv(path, nrows=config.DEBUG_NROWS, low_memory=False)
    frame = _drop_unnamed(frame)

    report: dict = {
        "file": path.name,
        "rows_raw": int(frame.shape[0]),
        "columns_raw": int(frame.shape[1]),
    }

    if config.MODEL_ID not in frame.columns:
        raise KeyError(
            f"Expression matrix has no {config.MODEL_ID!r} column. "
            f"First 10 columns are: {list(frame.columns[:10])}. "
            f"Update config.EXPRESSION_META_COLS if this release differs."
        )

    # --- Trap 2: collapse to one profile per cell line -------------------
    flag = config.EXPRESSION_DEFAULT_FLAG
    if flag in frame.columns:
        is_default = _coerce_bool(frame[flag])
        report["rows_before_default_filter"] = int(frame.shape[0])
        frame = frame.loc[is_default].copy()
        report["rows_after_default_filter"] = int(frame.shape[0])
        report["rows_dropped_non_default"] = (
            report["rows_before_default_filter"] - report["rows_after_default_filter"]
        )
        _log(f"  {flag} filter: kept {frame.shape[0]} of "
             f"{report['rows_before_default_filter']} rows")
    else:
        report["rows_before_default_filter"] = int(frame.shape[0])
        report["rows_after_default_filter"] = int(frame.shape[0])
        report["rows_dropped_non_default"] = 0
        _log(f"  WARNING: no {flag!r} column found; skipping the default-entry "
             f"filter. Verify ModelID uniqueness below.")

    # --- Deduplicate defensively -----------------------------------------
    duplicated = frame[config.MODEL_ID].duplicated(keep="first")
    report["duplicate_model_ids_dropped"] = int(duplicated.sum())
    if duplicated.any():
        _log(f"  WARNING: {int(duplicated.sum())} duplicate ModelIDs remained "
             f"after filtering; keeping the first occurrence of each.")
        frame = frame.loc[~duplicated].copy()

    frame = frame.set_index(config.MODEL_ID)
    frame.index.name = config.MODEL_ID

    # --- Trap 1: strip metadata columns ----------------------------------
    meta_present = [
        c for c in config.EXPRESSION_META_COLS
        if c in frame.columns and c != config.MODEL_ID
    ]
    frame = frame.drop(columns=meta_present)
    report["metadata_columns_dropped"] = meta_present

    # Drop anything non-numeric that survived.
    non_numeric = [
        c for c in frame.columns
        if not pd.api.types.is_numeric_dtype(frame[c])
    ]
    if non_numeric:
        _log(f"  Dropping {len(non_numeric)} non-numeric columns "
             f"(e.g. {non_numeric[:3]})")
        frame = frame.drop(columns=non_numeric)
    report["non_numeric_columns_dropped"] = len(non_numeric)

    frame = frame.astype(config.FLOAT_DTYPE)

    report["cell_lines"] = int(frame.shape[0])
    report["gene_columns"] = int(frame.shape[1])
    _log(f"  -> {frame.shape[0]} cell lines x {frame.shape[1]} gene columns")
    return frame, report


def load_crispr() -> tuple[pd.DataFrame, dict]:
    """Load the CRISPR (Chronos) gene effect matrix, indexed by ModelID."""
    path = config.resolve_file("crispr")
    _log(f"Reading CRISPR gene effect matrix: {path.name}")
    _log("  (also large; please be patient)")

    frame = pd.read_csv(path, index_col=0, nrows=config.DEBUG_NROWS,
                        low_memory=False)
    frame.index.name = config.MODEL_ID

    report: dict = {
        "file": path.name,
        "rows_raw": int(frame.shape[0]),
        "columns_raw": int(frame.shape[1]),
    }

    duplicated = frame.index.duplicated(keep="first")
    report["duplicate_model_ids_dropped"] = int(duplicated.sum())
    if duplicated.any():
        _log(f"  WARNING: dropping {int(duplicated.sum())} duplicate ModelIDs")
        frame = frame.loc[~duplicated]

    non_numeric = [
        c for c in frame.columns
        if not pd.api.types.is_numeric_dtype(frame[c])
    ]
    if non_numeric:
        frame = frame.drop(columns=non_numeric)
    report["non_numeric_columns_dropped"] = len(non_numeric)

    frame = frame.astype(config.FLOAT_DTYPE)

    report["cell_lines"] = int(frame.shape[0])
    report["gene_columns"] = int(frame.shape[1])
    _log(f"  -> {frame.shape[0]} cell lines x {frame.shape[1]} gene columns")
    return frame, report


def load_model_metadata() -> tuple[pd.DataFrame, dict]:
    """Load Model.csv, keeping the curated subset of columns."""
    path = config.resolve_file("model")
    _log(f"Reading cell line metadata: {path.name}")

    frame = pd.read_csv(path, index_col=0, low_memory=False)
    frame.index.name = config.MODEL_ID

    report: dict = {
        "file": path.name,
        "rows_raw": int(frame.shape[0]),
        "columns_raw": int(frame.shape[1]),
    }

    duplicated = frame.index.duplicated(keep="first")
    if duplicated.any():
        frame = frame.loc[~duplicated]
    report["duplicate_model_ids_dropped"] = int(duplicated.sum())

    keep = [c for c in config.MODEL_KEEP_COLS if c in frame.columns]
    missing = [c for c in config.MODEL_KEEP_COLS if c not in frame.columns]
    if missing:
        _log(f"  Note: these expected metadata columns are absent: {missing}")
    report["metadata_columns_missing"] = missing

    frame = frame[keep].copy()
    report["cell_lines"] = int(frame.shape[0])
    _log(f"  -> {frame.shape[0]} cell lines x {frame.shape[1]} metadata columns")
    return frame, report


def load_prism() -> tuple[pd.DataFrame | None, dict]:
    """
    Load the PRISM drug-response matrix if present.

    Optional by design: the pipeline is fully functional without it, and the
    drug-response head is a later phase. Handles both the older
    `primary-screen-replicate-collapsed-logfold-change.csv` layout and the
    newer `Repurposing_Public_*_Data_Matrix.csv` layout, and transposes
    automatically if compounds are on the rows.
    """
    path = config.resolve_file("prism_matrix", required=False)
    if path is None:
        _log("PRISM drug-response matrix not found -- skipping.")
        _log("  (This is expected for now. The drug head is a later phase;")
        _log("   drop the file into the data directory and re-run to add it.)")
        return None, {"file": None, "status": "absent"}

    _log(f"Reading PRISM drug response matrix: {path.name}")
    frame = pd.read_csv(path, index_col=0, low_memory=False)

    report: dict = {
        "file": path.name,
        "rows_raw": int(frame.shape[0]),
        "columns_raw": int(frame.shape[1]),
        "status": "present",
    }

    def _looks_like_model_ids(values) -> bool:
        sample = [str(v) for v in list(values)[:50]]
        if not sample:
            return False
        hits = sum(1 for v in sample if v.startswith("ACH-"))
        return hits >= max(1, len(sample) // 2)

    if not _looks_like_model_ids(frame.index):
        if _looks_like_model_ids(frame.columns):
            _log("  Cell lines are on the columns; transposing.")
            frame = frame.T
            report["transposed"] = True
        else:
            _log("  WARNING: could not find ACH- style cell line IDs on either "
                 "axis. PRISM will be skipped.")
            report["status"] = "unrecognised_layout"
            return None, report
    else:
        report["transposed"] = False

    frame.index = frame.index.astype(str)
    frame.index.name = config.MODEL_ID

    duplicated = frame.index.duplicated(keep="first")
    report["duplicate_model_ids_dropped"] = int(duplicated.sum())
    if duplicated.any():
        frame = frame.loc[~duplicated]

    non_numeric = [
        c for c in frame.columns
        if not pd.api.types.is_numeric_dtype(frame[c])
    ]
    if non_numeric:
        frame = frame.drop(columns=non_numeric)
    report["non_numeric_columns_dropped"] = len(non_numeric)

    # Deduplicate compound columns (a known quirk of the 24Q2 compound list:
    # at least one BRD ID is reported twice, another is missing).
    dup_cols = frame.columns.duplicated(keep="first")
    report["duplicate_compound_columns_dropped"] = int(dup_cols.sum())
    if dup_cols.any():
        _log(f"  Dropping {int(dup_cols.sum())} duplicate compound columns")
        frame = frame.loc[:, ~dup_cols]

    frame = frame.astype(config.FLOAT_DTYPE)
    report["cell_lines"] = int(frame.shape[0])
    report["compounds"] = int(frame.shape[1])
    _log(f"  -> {frame.shape[0]} cell lines x {frame.shape[1]} compounds")
    return frame, report


# --------------------------------------------------------------------------
# Target selection
# --------------------------------------------------------------------------

def select_selective_genes(crispr: pd.DataFrame) -> tuple[list[str], dict]:
    """
    Choose which CRISPR genes are worth predicting.

    Chronos gene effect is scaled so that 0 is no effect and -1 is the median
    of pan-essential genes. Two classes of gene are excluded:

    * Pan-essential genes (nearly every line depends on them). A model scores
      well on these by learning a constant, which says nothing about biology
      but inflates aggregate metrics.
    * Never-essential genes (no line depends on them). These are noise.

    What remains are *selective* dependencies -- genes essential in some
    genetic contexts and not others. This is where context-specific signal
    lives, and it is the only honest place to measure whether an expression
    model has learned anything.
    """
    n_lines = crispr.shape[0]

    nan_fraction = crispr.isna().mean(axis=0)
    dependent_counts = (crispr < config.DEPENDENCY_THRESHOLD).sum(axis=0)
    dependent_fraction = dependent_counts / max(n_lines, 1)
    gene_std = crispr.std(axis=0, skipna=True)

    keep_mask = (
        (nan_fraction <= config.MAX_GENE_NAN_FRACTION)
        & (dependent_counts >= config.MIN_DEPENDENT_LINES)
        & (dependent_fraction <= config.MAX_DEPENDENT_FRACTION)
        & (gene_std >= config.MIN_GENE_EFFECT_STD)
    )

    selected = crispr.columns[keep_mask].tolist()

    report = {
        "cell_lines_considered": int(n_lines),
        "genes_considered": int(crispr.shape[1]),
        "dependency_threshold": config.DEPENDENCY_THRESHOLD,
        "min_dependent_lines": config.MIN_DEPENDENT_LINES,
        "max_dependent_fraction": config.MAX_DEPENDENT_FRACTION,
        "min_gene_effect_std": config.MIN_GENE_EFFECT_STD,
        "max_gene_nan_fraction": config.MAX_GENE_NAN_FRACTION,
        "excluded_too_many_nans": int((nan_fraction > config.MAX_GENE_NAN_FRACTION).sum()),
        "excluded_never_essential": int((dependent_counts < config.MIN_DEPENDENT_LINES).sum()),
        "excluded_pan_essential": int((dependent_fraction > config.MAX_DEPENDENT_FRACTION).sum()),
        "excluded_low_variance": int((gene_std < config.MIN_GENE_EFFECT_STD).sum()),
        "selective_genes_kept": len(selected),
    }
    return selected, report


def filter_expression_features(expression: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Drop uninformative expression features (all-zero or flat genes)."""
    report: dict = {"genes_before": int(expression.shape[1])}

    keep = pd.Series(True, index=expression.columns)

    if config.DROP_ALL_ZERO_EXPRESSION:
        all_zero = (expression.fillna(0.0) == 0.0).all(axis=0)
        report["dropped_all_zero"] = int(all_zero.sum())
        keep &= ~all_zero
    else:
        report["dropped_all_zero"] = 0

    if config.MIN_EXPRESSION_STD > 0:
        low_var = expression.std(axis=0, skipna=True) < config.MIN_EXPRESSION_STD
        report["dropped_low_variance"] = int((low_var & keep).sum())
        keep &= ~low_var
    else:
        report["dropped_low_variance"] = 0

    filtered = expression.loc[:, keep]
    report["genes_after"] = int(filtered.shape[1])
    return filtered, report


# --------------------------------------------------------------------------
# Report writing
# --------------------------------------------------------------------------

def write_report(report: dict, processed_dir: Path) -> tuple[Path, Path]:
    """Write the audit report in both human and machine readable form."""
    json_path = processed_dir / "join_report.json"
    io_utils.save_json(report, json_path)

    lines: list[str] = []
    lines.append("=" * 74)
    lines.append("DEPMAP DATASET BUILD REPORT")
    lines.append("=" * 74)
    lines.append("")
    lines.append("This file records every filtering and intersection step that")
    lines.append("produced the modelling dataset. Cite these numbers in your")
    lines.append("methods section; they are the provenance of everything")
    lines.append("downstream.")
    lines.append("")

    def _section(title: str, payload: dict) -> None:
        lines.append("-" * 74)
        lines.append(title)
        lines.append("-" * 74)
        for key, value in payload.items():
            if isinstance(value, list) and len(value) > 6:
                shown = f"{value[:6]} ... ({len(value)} total)"
            else:
                shown = value
            lines.append(f"  {key:<38}: {shown}")
        lines.append("")

    for title, key in [
        ("1. EXPRESSION MATRIX", "expression"),
        ("2. CRISPR GENE EFFECT MATRIX", "crispr"),
        ("3. CELL LINE METADATA", "model_metadata"),
        ("4. PRISM DRUG RESPONSE (optional)", "prism"),
        ("5. GENE SPACE RECONCILIATION", "gene_space"),
        ("6. EXPRESSION FEATURE FILTER", "expression_features"),
        ("7. SELECTIVE CRISPR TARGET SELECTION", "selective_genes"),
        ("8. CELL LINE INTERSECTION", "intersection"),
        ("9. OSTEOSARCOMA COVERAGE", "osteosarcoma"),
        ("10. FINAL DATASET", "final"),
    ]:
        if key in report:
            _section(title, report[key])

    text = "\n".join(lines)
    txt_path = processed_dir / "join_report.txt"
    txt_path.write_text(text, encoding="utf-8")
    return txt_path, json_path


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main() -> int:
    config.ensure_dirs()

    _rule("CONFIGURATION")
    _log(config.describe())

    report: dict = {}

    # ---------------------------------------------------------------- load
    _rule("LOADING RAW FILES")
    try:
        expression, report["expression"] = load_expression()
        crispr, report["crispr"] = load_crispr()
        metadata, report["model_metadata"] = load_model_metadata()
        prism, report["prism"] = load_prism()
    except config.MissingDataFile as exc:
        _log(str(exc))
        return 1

    # ---------------------------------------------------------- gene space
    _rule("RECONCILING GENE SPACE")
    space = gene_ids.intersect_gene_spaces(
        expression.columns,
        crispr.columns,
        source_a="expression",
        source_b="crispr",
    )
    report["gene_space"] = space.report
    _log(f"  expression genes parsed : {space.report['expression_parsed']}")
    _log(f"  crispr genes parsed     : {space.report['crispr_parsed']}")
    _log(f"  shared (by Entrez ID)   : {space.report['shared_genes']}")
    if space.report["symbol_conflicts"]:
        _log(f"  symbol disagreements    : {space.report['symbol_conflicts']} "
             f"(expression symbol preferred)")

    if not space.entrez_ids:
        _log("\nFATAL: the two matrices share no genes. Check that both files "
             "use the 'SYMBOL (ENTREZ)' column format.")
        return 1

    canonical_cols = gene_ids.canonical_labels(space.entrez_ids, space.symbols)

    expression = expression[space.labels_a]
    expression.columns = canonical_cols
    crispr = crispr[space.labels_b]
    crispr.columns = canonical_cols

    # -------------------------------------------------- expression filter
    _rule("FILTERING EXPRESSION FEATURES")
    expression, report["expression_features"] = filter_expression_features(expression)
    _log(f"  {report['expression_features']['genes_before']} -> "
         f"{report['expression_features']['genes_after']} genes")

    # Keep the canonical gene record aligned to the surviving features.
    surviving = set(expression.columns)
    kept_pairs = [
        (ent, sym, lab)
        for ent, sym, lab in zip(space.entrez_ids, space.symbols, canonical_cols)
        if lab in surviving
    ]
    final_entrez = [p[0] for p in kept_pairs]
    final_symbols = [p[1] for p in kept_pairs]
    final_labels = [p[2] for p in kept_pairs]

    # ------------------------------------------------ cell line alignment
    _rule("INTERSECTING CELL LINES")
    shared_lines = (
        expression.index
        .intersection(crispr.index)
        .intersection(metadata.index)
    )
    shared_lines = pd.Index(sorted(shared_lines), name=config.MODEL_ID)

    intersection_report = {
        "expression_lines": int(expression.shape[0]),
        "crispr_lines": int(crispr.shape[0]),
        "metadata_lines": int(metadata.shape[0]),
        "shared_expression_crispr_metadata": int(len(shared_lines)),
        "lost_from_expression": int(expression.shape[0] - len(shared_lines)),
        "lost_from_crispr": int(crispr.shape[0] - len(shared_lines)),
    }
    _log(f"  expression : {intersection_report['expression_lines']}")
    _log(f"  crispr     : {intersection_report['crispr_lines']}")
    _log(f"  metadata   : {intersection_report['metadata_lines']}")
    _log(f"  -> shared  : {intersection_report['shared_expression_crispr_metadata']}")

    if len(shared_lines) == 0:
        _log("\nFATAL: no cell lines are shared across the three files. "
             "Check that all files come from the same DepMap release.")
        return 1

    expression = expression.loc[shared_lines]
    crispr = crispr.loc[shared_lines]
    metadata = metadata.loc[shared_lines]

    if prism is not None:
        prism_shared = prism.index.intersection(shared_lines)
        intersection_report["prism_lines"] = int(prism.shape[0])
        intersection_report["prism_overlap_with_core"] = int(len(prism_shared))
        _log(f"  prism      : {intersection_report['prism_lines']} "
             f"({len(prism_shared)} overlap with the core set)")
        prism = prism.loc[sorted(prism_shared)]

    report["intersection"] = intersection_report

    # ------------------------------------------------- selective targets
    _rule("SELECTING CRISPR TARGETS")
    selective, report["selective_genes"] = select_selective_genes(crispr)
    sel = report["selective_genes"]
    _log(f"  considered            : {sel['genes_considered']}")
    _log(f"  excluded, too many NA : {sel['excluded_too_many_nans']}")
    _log(f"  excluded, never ess.  : {sel['excluded_never_essential']}")
    _log(f"  excluded, pan-ess.    : {sel['excluded_pan_essential']}")
    _log(f"  excluded, low variance: {sel['excluded_low_variance']}")
    _log(f"  -> selective targets  : {sel['selective_genes_kept']}")

    if not selective:
        _log("\nWARNING: no selective genes survived the filter. The "
             "thresholds in config.py are probably too strict for this "
             "dataset (or DEBUG_NROWS is set very low).")

    # ------------------------------------------- osteosarcoma coverage
    _rule("OSTEOSARCOMA COVERAGE")
    os_report: dict = {}
    if "OncotreeLineage" in metadata.columns:
        bone_mask = metadata["OncotreeLineage"] == config.OSTEOSARCOMA_LINEAGE
        os_report["bone_lineage_lines"] = int(bone_mask.sum())
    os_mask = config.osteosarcoma_mask(metadata)
    os_report["osteosarcoma_lines"] = int(os_mask.sum())
    os_report["osteosarcoma_model_ids"] = metadata.index[os_mask].tolist()
    if prism is not None:
        os_report["osteosarcoma_lines_in_prism"] = int(
            len(prism.index.intersection(metadata.index[os_mask]))
        )
    report["osteosarcoma"] = os_report
    for key, value in os_report.items():
        if key != "osteosarcoma_model_ids":
            _log(f"  {key:<34}: {value}")

    # ------------------------------------------------------------- save
    _rule("WRITING OUTPUTS")
    out = config.PROCESSED_DIR

    paths = {
        "expression": io_utils.save_matrix(expression, out / "expression"),
        "crispr_effect": io_utils.save_matrix(crispr, out / "crispr_effect"),
    }
    if prism is not None and prism.shape[0] > 0:
        paths["prism_response"] = io_utils.save_matrix(prism, out / "prism_response")

    paths["model_metadata"] = io_utils.save_table(metadata, out / "model_metadata")

    gene_map = pd.DataFrame(
        {
            "entrez": final_entrez,
            "symbol": final_symbols,
            "canonical_label": final_labels,
        }
    )
    paths["gene_id_map"] = io_utils.save_table(
        gene_map.set_index("entrez"), out / "gene_id_map"
    )

    io_utils.save_json(
        {
            "description": (
                "Canonical gene feature space. Any external expression matrix "
                "scored by the trained model MUST be reindexed to exactly "
                "these columns, in exactly this order. Use "
                "gene_ids.map_external_matrix()."
            ),
            "n_genes": len(final_labels),
            "entrez_ids": final_entrez,
            "symbols": final_symbols,
            "canonical_columns": final_labels,
        },
        out / "gene_columns.json",
    )
    paths["gene_columns"] = out / "gene_columns.json"

    io_utils.save_json(
        {
            "description": (
                "CRISPR genes retained as prediction targets: selective "
                "dependencies only. Pan-essential and never-essential genes "
                "are excluded because predicting them is trivial and inflates "
                "aggregate metrics."
            ),
            "criteria": report["selective_genes"],
            "n_genes": len(selective),
            "genes": selective,
        },
        out / "selective_genes.json",
    )
    paths["selective_genes"] = out / "selective_genes.json"

    report["final"] = {
        "cell_lines": int(len(shared_lines)),
        "expression_features": int(expression.shape[1]),
        "crispr_targets_all": int(crispr.shape[1]),
        "crispr_targets_selective": len(selective),
        "prism_compounds": int(prism.shape[1]) if prism is not None else 0,
        "prism_cell_lines": int(prism.shape[0]) if prism is not None else 0,
    }

    txt_path, json_path = write_report(report, out)
    paths["join_report"] = txt_path

    for name, path in paths.items():
        _log(f"  {name:<22}-> {Path(path).name}")

    # ------------------------------------------------------------ summary
    _rule("SUMMARY")
    final = report["final"]
    _log(f"  cell lines               : {final['cell_lines']}")
    _log(f"  expression features      : {final['expression_features']}")
    _log(f"  CRISPR targets (all)     : {final['crispr_targets_all']}")
    _log(f"  CRISPR targets (selective): {final['crispr_targets_selective']}")
    _log(f"  PRISM compounds          : {final['prism_compounds']}")
    _log("")
    _log(f"  Full audit written to    : {txt_path}")
    _log("")
    _log("  Next: python splits.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
