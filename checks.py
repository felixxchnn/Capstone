"""
checks.py
=========
Standalone integrity checks for the processed dataset.

Run any time, and especially after changing anything upstream:

    python checks.py

Exits 0 if everything passes, 1 if any check fails.

Purpose
-------
Every check here corresponds to a specific way this dataset could be silently
wrong -- wrong in a way that produces no exception, no warning, and a plausible
number at the end. Run it after every change to the pipeline, and run it once
more before you report any result.

The checks are deliberately independent of the code that built the dataset. If
`build_dataset.py` has a bug, a check written from the same assumptions would
share it, so these are written against the *properties the data must have*
rather than against the steps that produced it.
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

import config
import io_utils


class CheckResult:
    """Accumulates pass/fail outcomes and prints a readable summary."""

    def __init__(self) -> None:
        self.passed: list[str] = []
        self.failed: list[tuple[str, str]] = []
        self.warnings: list[tuple[str, str]] = []

    def check(self, name: str, condition: bool, detail: str = "") -> bool:
        if condition:
            self.passed.append(name)
            print(f"  [PASS] {name}")
        else:
            self.failed.append((name, detail))
            print(f"  [FAIL] {name}")
            if detail:
                print(f"         {detail}")
        return condition

    def warn(self, name: str, condition: bool, detail: str = "") -> bool:
        """A soft check: worth knowing about, not fatal."""
        if condition:
            self.passed.append(name)
            print(f"  [PASS] {name}")
        else:
            self.warnings.append((name, detail))
            print(f"  [WARN] {name}")
            if detail:
                print(f"         {detail}")
        return condition

    def summary(self) -> int:
        print()
        print("=" * 74)
        print(f"  passed   : {len(self.passed)}")
        print(f"  warnings : {len(self.warnings)}")
        print(f"  failed   : {len(self.failed)}")
        print("=" * 74)

        if self.warnings:
            print("\nWarnings:")
            for name, detail in self.warnings:
                print(f"  - {name}: {detail}")

        if self.failed:
            print("\nFailures:")
            for name, detail in self.failed:
                print(f"  - {name}: {detail}")
            print("\nDo not train on this dataset until these are resolved.")
            return 1

        print("\nDataset integrity verified.")
        return 0


def _section(title: str) -> None:
    print(f"\n{title}")
    print("-" * 74)


def main() -> int:
    out = config.PROCESSED_DIR
    result = CheckResult()

    print("=" * 74)
    print("DATASET INTEGRITY CHECKS")
    print("=" * 74)
    print(f"Reading from: {out}")

    # ------------------------------------------------------------- load
    try:
        expression = io_utils.load_matrix(out / "expression")
        crispr = io_utils.load_matrix(out / "crispr_effect")
        metadata = io_utils.load_table(out / "model_metadata")
        gene_columns = io_utils.load_json(out / "gene_columns.json")
        selective = io_utils.load_json(out / "selective_genes.json")
    except FileNotFoundError as exc:
        print(f"\nCould not load the processed dataset:\n  {exc}")
        return 1

    metadata.index.name = config.MODEL_ID

    try:
        splits_payload = io_utils.load_json(out / "splits.json")
        assignment = pd.Series(splits_payload["assignment"], name="split")
    except FileNotFoundError:
        assignment = None
        print("\nNote: splits.json not found; split checks will be skipped.")

    try:
        prism = io_utils.load_matrix(out / "prism_response")
    except FileNotFoundError:
        prism = None

    # ------------------------------------------------------- shape sanity
    _section("1. Shapes and non-emptiness")
    result.check(
        "expression matrix is non-empty",
        expression.shape[0] > 0 and expression.shape[1] > 0,
        f"shape={expression.shape}",
    )
    result.check(
        "crispr matrix is non-empty",
        crispr.shape[0] > 0 and crispr.shape[1] > 0,
        f"shape={crispr.shape}",
    )
    result.check(
        "metadata is non-empty",
        metadata.shape[0] > 0,
        f"shape={metadata.shape}",
    )
    print(f"         expression: {expression.shape[0]} lines x "
          f"{expression.shape[1]} genes")
    print(f"         crispr    : {crispr.shape[0]} lines x "
          f"{crispr.shape[1]} genes")

    # ------------------------------------------------- index integrity
    _section("2. Cell line index integrity")
    result.check(
        "expression ModelIDs are unique",
        expression.index.is_unique,
        f"{int(expression.index.duplicated().sum())} duplicates. This is the "
        f"IsDefaultEntryForModel trap -- see build_dataset.py.",
    )
    result.check(
        "crispr ModelIDs are unique",
        crispr.index.is_unique,
        f"{int(crispr.index.duplicated().sum())} duplicates",
    )
    result.check(
        "metadata ModelIDs are unique",
        metadata.index.is_unique,
        f"{int(metadata.index.duplicated().sum())} duplicates",
    )
    result.check(
        "expression and crispr share an identical, aligned index",
        expression.index.equals(crispr.index),
        "Indexes differ in content or order. Row i of one matrix does not "
        "correspond to row i of the other.",
    )
    result.check(
        "metadata covers every cell line in the matrices",
        expression.index.isin(metadata.index).all(),
        f"{int((~expression.index.isin(metadata.index)).sum())} cell lines "
        f"have no metadata row",
    )
    result.check(
        "all ModelIDs look like DepMap identifiers",
        pd.Series(expression.index.astype(str)).str.startswith("ACH-").all(),
        "Some index values are not ACH- prefixed; the wrong column may have "
        "been used as the index.",
    )

    # ------------------------------------------------- column integrity
    _section("3. Gene column integrity")
    canonical = gene_columns["canonical_columns"]
    result.check(
        "expression columns are unique",
        expression.columns.is_unique,
        f"{int(pd.Index(expression.columns).duplicated().sum())} duplicates",
    )
    result.check(
        "crispr columns are unique",
        crispr.columns.is_unique,
        f"{int(pd.Index(crispr.columns).duplicated().sum())} duplicates",
    )
    result.check(
        "expression columns match gene_columns.json exactly, in order",
        list(expression.columns) == list(canonical),
        "The saved canonical order does not match the expression matrix. "
        "Any external data mapped through gene_columns.json would be "
        "misaligned against the trained model.",
    )
    result.check(
        "crispr columns are a subset of the canonical space",
        set(crispr.columns) <= set(canonical),
        f"{len(set(crispr.columns) - set(canonical))} crispr columns are "
        f"outside the canonical gene space",
    )
    result.check(
        "gene_columns.json is internally consistent",
        len(gene_columns["entrez_ids"]) == len(gene_columns["symbols"])
        == len(canonical) == gene_columns["n_genes"],
        "entrez_ids, symbols, canonical_columns and n_genes disagree",
    )
    result.check(
        "Entrez IDs are unique",
        len(set(gene_columns["entrez_ids"])) == len(gene_columns["entrez_ids"]),
        "Duplicate Entrez IDs in the canonical space",
    )

    # ----------------------------------------------------- value sanity
    _section("4. Value sanity")
    expr_values = expression.to_numpy()
    crispr_values = crispr.to_numpy()

    result.check(
        "expression contains no infinities",
        not np.isinf(expr_values).any(),
        f"{int(np.isinf(expr_values).sum())} infinite values",
    )
    result.check(
        "crispr contains no infinities",
        not np.isinf(crispr_values).any(),
        f"{int(np.isinf(crispr_values).sum())} infinite values",
    )
    result.check(
        "expression has no all-NaN rows",
        not np.isnan(expr_values).all(axis=1).any(),
        f"{int(np.isnan(expr_values).all(axis=1).sum())} cell lines have no "
        f"expression data at all",
    )

    expr_nan_fraction = float(np.isnan(expr_values).mean())
    crispr_nan_fraction = float(np.isnan(crispr_values).mean())
    result.warn(
        "expression missingness is low",
        expr_nan_fraction < 0.01,
        f"{expr_nan_fraction:.2%} of expression values are NaN",
    )
    result.warn(
        "crispr missingness is low",
        crispr_nan_fraction < 0.05,
        f"{crispr_nan_fraction:.2%} of gene effect values are NaN",
    )

    finite_expr = expr_values[np.isfinite(expr_values)]
    if finite_expr.size:
        result.warn(
            "expression values are on a log scale (non-negative, modest range)",
            bool(finite_expr.min() >= -0.001 and finite_expr.max() < 100),
            f"range [{finite_expr.min():.2f}, {finite_expr.max():.2f}] -- "
            f"log2(TPM+1) should be non-negative and typically under ~20. "
            f"Raw TPM may have been loaded by mistake.",
        )

    finite_crispr = crispr_values[np.isfinite(crispr_values)]
    if finite_crispr.size:
        result.warn(
            "crispr values are on the Chronos scale",
            bool(finite_crispr.min() > -10 and finite_crispr.max() < 10),
            f"range [{finite_crispr.min():.2f}, {finite_crispr.max():.2f}] -- "
            f"Chronos gene effect is centred near 0 with pan-essentials "
            f"near -1.",
        )
        median = float(np.median(finite_crispr))
        result.warn(
            "crispr median is near zero",
            abs(median) < 0.25,
            f"median gene effect is {median:.3f}; expected close to 0",
        )

    # ------------------------------------------------ selective targets
    _section("5. Selective CRISPR targets")
    selective_genes = selective["genes"]
    result.check(
        "at least some selective genes were retained",
        len(selective_genes) > 0,
        "No genes passed the selectivity filter. Thresholds in config.py are "
        "too strict for this dataset.",
    )
    result.check(
        "every selective gene exists in the crispr matrix",
        set(selective_genes) <= set(crispr.columns),
        f"{len(set(selective_genes) - set(crispr.columns))} selective genes "
        f"are missing from the matrix",
    )
    if selective_genes:
        subset = crispr[selective_genes].to_numpy()
        finite = subset[np.isfinite(subset)]
        if finite.size:
            result.warn(
                "selective genes actually vary across cell lines",
                float(np.nanstd(subset)) > 0.05,
                f"std={float(np.nanstd(subset)):.4f}; near-constant targets "
                f"cannot be predicted meaningfully",
            )
        print(f"         {len(selective_genes)} selective targets retained")

    # -------------------------------------------------------- splits
    if assignment is not None:
        _section("6. Split integrity")

        aligned = assignment.reindex(expression.index)
        result.check(
            "every cell line has a split assignment",
            aligned.notna().all(),
            f"{int(aligned.isna().sum())} cell lines are unassigned",
        )
        result.check(
            "split labels are valid",
            set(aligned.dropna().unique()) <= {"train", "val", "test"},
            f"unexpected labels: {sorted(set(aligned.dropna().unique()))}",
        )

        sizes = aligned.value_counts().to_dict()
        result.check(
            "no split is empty",
            all(sizes.get(name, 0) > 0 for name in ("train", "val", "test")),
            f"sizes: {sizes}",
        )

        # The decisive one: no patient group straddles a split.
        if config.GROUP_COL in metadata.columns:
            groups = metadata[config.GROUP_COL].astype(str).reindex(aligned.index)
            groups = groups.where(groups.notna() & (groups != "nan"),
                                  pd.Series(aligned.index, index=aligned.index))
            straddling = [
                str(g) for g, sub in aligned.groupby(groups) if sub.nunique() > 1
            ]
            result.check(
                f"no {config.GROUP_COL} straddles a split boundary",
                len(straddling) == 0,
                f"{len(straddling)} groups appear in more than one split, "
                f"e.g. {straddling[:5]}. Sibling cell lines from one donor on "
                f"both sides of the split is leakage.",
            )

        for name in ("train", "val", "test"):
            print(f"         {name:<6}: {sizes.get(name, 0)} cell lines")

    # --------------------------------------------------------- prism
    if prism is not None:
        _section("7. PRISM drug response")
        result.check(
            "prism ModelIDs are unique",
            prism.index.is_unique,
            f"{int(prism.index.duplicated().sum())} duplicates",
        )
        result.check(
            "prism cell lines are a subset of the core dataset",
            set(prism.index) <= set(expression.index),
            f"{len(set(prism.index) - set(expression.index))} PRISM cell "
            f"lines are absent from the expression matrix",
        )
        result.warn(
            "prism overlap is large enough to model",
            prism.shape[0] >= 50,
            f"only {prism.shape[0]} cell lines overlap; a drug-response head "
            f"trained on this will be underpowered",
        )
        print(f"         {prism.shape[0]} lines x {prism.shape[1]} compounds")

    # ---------------------------------------------- osteosarcoma coverage
    _section("8. Osteosarcoma coverage")
    os_mask = config.osteosarcoma_mask(metadata)
    os_lines = metadata.index[os_mask]
    in_dataset = expression.index.intersection(os_lines)
    result.warn(
        "osteosarcoma lines survive into the modelling dataset",
        len(in_dataset) >= 5,
        f"only {len(in_dataset)} osteosarcoma lines remain; the "
        f"tumour-type sanity check will be weak",
    )
    print(f"         {len(in_dataset)} osteosarcoma cell lines in the "
          f"final dataset")
    if assignment is not None and len(in_dataset):
        os_splits = assignment.reindex(in_dataset).value_counts().to_dict()
        print(f"         by split: {os_splits}")

    return result.summary()


if __name__ == "__main__":
    sys.exit(main())
