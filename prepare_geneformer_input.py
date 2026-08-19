"""
prepare_geneformer_input.py
===========================
Turns the processed DepMap expression matrix into a Geneformer-ready `.h5ad`.

    python prepare_geneformer_input.py

Output (to config.PROCESSED_DIR)
--------------------------------
    geneformer_input.h5ad        cells x genes, raw-count-style X, ready to tokenize
    ensembl_map.csv              entrez -> ensembl_id (cached; reused on re-run)
    geneformer_prep_report.json  coverage and provenance

Why this module is needed, and why it is fiddly
-----------------------------------------------
Geneformer's tokenizer does not accept what you have. Two hard requirements:

1. **Raw counts, not TPM.** Geneformer expects raw read counts in `X`; it does
   its own median normalisation internally. Your DepMap file is log2(TPM+1) --
   already normalised and logged. The definitive fix is to download DepMap's
   gene-level *expected counts* expression file and point this module at it
   (see `config.FILE_ALIASES["expression_counts"]`). If that file is absent,
   this module reconstructs an approximate count-like matrix from the log-TPM
   so the pipeline runs today -- but that is an approximation, flagged loudly,
   and the counts file is the correct input for any reported result.

2. **Ensembl IDs, not symbols.** Geneformer keys its vocabulary on Ensembl gene
   IDs. Your columns are `SYMBOL (ENTREZ)`. This module maps Entrez -> Ensembl
   (via a cached CSV, or `mygene` if the CSV is absent) and drops genes with no
   mapping. Geneformer then filters again to its own vocabulary, so genes it
   does not recognise fall away at tokenisation regardless.

Testability note
----------------
`anndata` is not installed in every environment, and the final `.h5ad` write
depends on it. Everything upstream of that write -- count handling, Ensembl
mapping, frame construction, coverage reporting -- is plain pandas/numpy and is
covered by `_self_test()` (run `python prepare_geneformer_input.py --self-test`).
The AnnData assembly itself is a thin wrapper. If `anndata` is missing, the
module still writes the three component frames to disk and tells you how to
finish the write where anndata is available.
"""

from __future__ import annotations

import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import config
import io_utils
import gene_ids


# --------------------------------------------------------------------------
# Count handling
# --------------------------------------------------------------------------

def logtpm_to_pseudocounts(log_tpm: pd.DataFrame) -> pd.DataFrame:
    """
    Reconstruct an approximate count-like matrix from log2(TPM+1).

    TPM = 2**log_tpm - 1. TPM is already depth-normalised (each sample sums to
    ~1e6), so the row sum used as `n_counts` is ~1e6 for every cell line. This
    preserves within-sample relative abundances -- which is what Geneformer's
    rank encoding uses -- but carries TPM's gene-length normalisation, which
    true counts do not. It is a stand-in, not the real thing.
    """
    tpm = np.power(2.0, log_tpm.to_numpy(dtype=np.float64)) - 1.0
    tpm = np.clip(tpm, 0.0, None)
    return pd.DataFrame(tpm, index=log_tpm.index, columns=log_tpm.columns)


def load_expression_counts() -> tuple[pd.DataFrame, dict]:
    """
    Load an expression matrix suitable for Geneformer, in the canonical gene
    space produced by build_dataset.py.

    Prefers a real expected-counts file if one is configured and present;
    otherwise reconstructs pseudo-counts from the processed log-TPM matrix.

    Returns
    -------
    (counts_df, report)
        counts_df: cell lines x canonical genes, non-negative, count-like.
    """
    processed_expr = io_utils.load_matrix(config.PROCESSED_DIR / "expression")
    canonical_cols = list(processed_expr.columns)

    counts_path = config.resolve_file("expression_counts", required=False)
    if counts_path is not None:
        print(f"Using real expected-counts file: {counts_path.name}")
        raw = pd.read_csv(counts_path, index_col=0, low_memory=False)

        # The counts file uses the same SYMBOL (ENTREZ) columns and ModelID
        # index as the TPM file. Reindex it into the canonical gene space and
        # to the exact cell lines already in the processed dataset.
        raw = raw.loc[raw.index.intersection(processed_expr.index)]
        raw = raw.reindex(index=processed_expr.index)

        available = [c for c in canonical_cols if c in raw.columns]
        missing = [c for c in canonical_cols if c not in raw.columns]
        counts = pd.DataFrame(
            0.0, index=processed_expr.index, columns=canonical_cols
        )
        counts[available] = raw[available].to_numpy(dtype=np.float64)
        counts = counts.clip(lower=0.0)

        report = {
            "count_source": "depmap_expected_counts",
            "count_file": counts_path.name,
            "genes_in_counts_file": int(len(available)),
            "genes_missing_from_counts_file": int(len(missing)),
        }
    else:
        print("No expected-counts file found -- reconstructing pseudo-counts "
              "from log-TPM.")
        print("  For a reported result, download DepMap's gene-level expected-"
              "counts expression file and re-run; see README.")
        counts = logtpm_to_pseudocounts(processed_expr)
        report = {
            "count_source": "reconstructed_from_log_tpm",
            "count_file": None,
            "note": (
                "Approximate. Geneformer expects raw counts; the definitive "
                "input is DepMap's expected-counts file."
            ),
        }

    report["cell_lines"] = int(counts.shape[0])
    report["canonical_genes"] = int(counts.shape[1])
    return counts, report


# --------------------------------------------------------------------------
# Ensembl mapping
# --------------------------------------------------------------------------

def build_ensembl_map_via_mygene(entrez_ids: list[str]) -> pd.DataFrame:
    """
    Query mygene.info for Entrez -> Ensembl gene ID.

    Requires the `mygene` package and network access (both available on Kaggle
    and Colab). Returns a DataFrame with columns [entrez, ensembl_id]; genes
    with no mapping are omitted. Where mygene returns several Ensembl IDs for
    one Entrez ID, the first is taken.
    """
    import mygene  # imported lazily so the module loads without it

    mg = mygene.MyGeneInfo()
    print(f"Querying mygene.info for {len(entrez_ids)} Entrez IDs...")
    results = mg.querymany(
        entrez_ids,
        scopes="entrezgene",
        fields="ensembl.gene",
        species="human",
        as_dataframe=False,
        returnall=False,
    )

    rows: list[dict] = []
    for hit in results:
        if hit.get("notfound"):
            continue
        query = str(hit.get("query"))
        ensembl = hit.get("ensembl")
        gene_id = None
        if isinstance(ensembl, dict):
            gene_id = ensembl.get("gene")
        elif isinstance(ensembl, list) and ensembl:
            first = ensembl[0]
            gene_id = first.get("gene") if isinstance(first, dict) else None
        if isinstance(gene_id, list):
            gene_id = gene_id[0] if gene_id else None
        if gene_id:
            rows.append({"entrez": query, "ensembl_id": gene_id})

    return pd.DataFrame(rows).drop_duplicates(subset="entrez")


def load_ensembl_map(entrez_ids: list[str]) -> tuple[pd.DataFrame, dict]:
    """
    Obtain an Entrez -> Ensembl mapping, cached to disk.

    Order of preference:
      1. ``ensembl_map.csv`` in the processed directory (from a previous run,
         or supplied manually) -- fully reproducible, no network needed.
      2. A live mygene.info query, whose result is then cached for next time.

    Returns
    -------
    (map_df, report)
        map_df: columns [entrez, ensembl_id], one row per mapped gene.
    """
    cache_path = config.PROCESSED_DIR / "ensembl_map.csv"

    if cache_path.is_file():
        print(f"Loading cached Ensembl map: {cache_path.name}")
        map_df = pd.read_csv(cache_path, dtype={"entrez": str, "ensembl_id": str})
        source = "cache"
    else:
        try:
            map_df = build_ensembl_map_via_mygene(entrez_ids)
            map_df.to_csv(cache_path, index=False)
            source = "mygene"
            print(f"Cached Ensembl map to {cache_path.name}")
        except ImportError:
            raise SystemExit(
                "\nNo ensembl_map.csv found and the `mygene` package is not "
                "installed.\nEither:\n"
                "  - pip install mygene   (needs network; easiest on Kaggle), or\n"
                f"  - place a CSV at {cache_path} with columns "
                "'entrez,ensembl_id'.\n"
            )

    map_df = map_df.dropna().drop_duplicates(subset="entrez")
    requested = set(entrez_ids)
    map_df = map_df[map_df["entrez"].astype(str).isin(requested)]

    report = {
        "ensembl_map_source": source,
        "entrez_requested": len(entrez_ids),
        "entrez_mapped": int(map_df.shape[0]),
        "entrez_unmapped": int(len(entrez_ids) - map_df.shape[0]),
    }
    return map_df, report


# --------------------------------------------------------------------------
# Frame construction (the testable core)
# --------------------------------------------------------------------------

def build_frames(
    counts: pd.DataFrame,
    gene_map: pd.DataFrame,
    ensembl_map: pd.DataFrame,
    metadata: pd.DataFrame,
    assignment: pd.Series | None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """
    Assemble the three AnnData components, keyed on Ensembl IDs.

    Parameters
    ----------
    counts
        cell lines x canonical `SYMBOL (ENTREZ)` columns.
    gene_map
        canonical gene record: columns include [entrez, canonical_label].
    ensembl_map
        columns [entrez, ensembl_id].
    metadata
        cell line metadata, indexed by ModelID.
    assignment
        split labels indexed by ModelID, or None.

    Returns
    -------
    (X_df, var_df, obs_df, report)
        X_df   : cells x genes, columns are Ensembl IDs (the count matrix).
        var_df : indexed by Ensembl ID, with 'ensembl_id' column.
        obs_df : indexed by ModelID, with n_counts, lineage, split, etc.
    """
    label_to_entrez = dict(zip(gene_map["canonical_label"], gene_map["entrez"].astype(str)))
    entrez_to_ensembl = dict(zip(ensembl_map["entrez"].astype(str), ensembl_map["ensembl_id"]))

    # Keep only genes that have an Ensembl mapping, preserving column order.
    kept_labels: list[str] = []
    kept_ensembl: list[str] = []
    seen_ensembl: set[str] = set()
    for label in counts.columns:
        entrez = label_to_entrez.get(label)
        if entrez is None:
            continue
        ensembl = entrez_to_ensembl.get(entrez)
        if ensembl is None or ensembl in seen_ensembl:
            continue  # unmapped, or a duplicate Ensembl collision
        kept_labels.append(label)
        kept_ensembl.append(ensembl)
        seen_ensembl.add(ensembl)

    X = counts[kept_labels].copy()
    X.columns = kept_ensembl

    var_df = pd.DataFrame({"ensembl_id": kept_ensembl}, index=kept_ensembl)
    var_df.index.name = "ensembl_id"

    obs_df = pd.DataFrame(index=counts.index)
    obs_df.index.name = config.MODEL_ID
    obs_df["ModelID"] = counts.index.astype(str)
    obs_df["n_counts"] = X.to_numpy().sum(axis=1)

    for col in ("OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype"):
        if col in metadata.columns:
            obs_df[col] = metadata.reindex(counts.index)[col].astype(str).values

    if assignment is not None:
        obs_df["split"] = assignment.reindex(counts.index).astype(str).values

    report = {
        "genes_in": int(counts.shape[1]),
        "genes_mapped_to_ensembl": int(len(kept_labels)),
        "genes_dropped_no_ensembl": int(counts.shape[1] - len(kept_labels)),
        "ensembl_collisions_dropped": int(
            len(kept_labels) - len(set(kept_ensembl))
        ),
        "cells": int(X.shape[0]),
        "median_n_counts": float(np.median(obs_df["n_counts"])),
        "cells_with_zero_counts": int((obs_df["n_counts"] <= 0).sum()),
    }
    return X, var_df, obs_df, report


def write_h5ad(
    X: pd.DataFrame,
    var_df: pd.DataFrame,
    obs_df: pd.DataFrame,
    path: Path,
) -> bool:
    """
    Assemble and write the AnnData. Thin wrapper over `anndata`.

    Returns True on success, False if anndata is unavailable (in which case the
    component frames are written separately so nothing is lost).
    """
    try:
        import anndata as ad
        from scipy import sparse
    except ImportError:
        print("\n  `anndata` is not installed, so the .h5ad was not written.")
        print("  The three component frames were saved instead:")
        io_utils.save_table(X, path.parent / "geneformer_X")
        io_utils.save_table(var_df, path.parent / "geneformer_var")
        io_utils.save_table(obs_df, path.parent / "geneformer_obs")
        print("    geneformer_X, geneformer_var, geneformer_obs")
        print("  Install anndata (pip install anndata) and re-run to produce "
              "the .h5ad, or assemble it from these frames on Kaggle.")
        return False

    adata = ad.AnnData(
        X=sparse.csr_matrix(X.to_numpy(dtype=np.float32)),
        obs=obs_df.copy(),
        var=var_df.copy(),
    )
    adata.write_h5ad(path)
    return True


# --------------------------------------------------------------------------
# Self test (runs without anndata / mygene / network)
# --------------------------------------------------------------------------

def _self_test() -> int:
    """Exercise the testable core on tiny synthetic inputs."""
    print("Running self-test of the testable core...")

    idx = pd.Index([f"ACH-{i:06d}" for i in range(6)], name=config.MODEL_ID)
    labels = ["AAA (1)", "BBB (2)", "CCC (3)", "DDD (4)"]
    log_tpm = pd.DataFrame(
        np.array([
            [0.0, 1.0, 2.0, 3.0],
            [3.0, 2.0, 1.0, 0.0],
            [1.0, 1.0, 1.0, 1.0],
            [0.0, 0.0, 5.0, 0.0],
            [2.0, 2.0, 2.0, 2.0],
            [4.0, 0.0, 0.0, 1.0],
        ]),
        index=idx, columns=labels,
    )

    # 1. count reconstruction is a monotone function of log-TPM
    counts = logtpm_to_pseudocounts(log_tpm)
    assert (counts.to_numpy() >= 0).all(), "counts must be non-negative"
    assert np.isclose(counts.iloc[0, 0], 0.0), "log2(0+... ) -> 2^0-1 = 0"
    assert np.isclose(counts.iloc[0, 3], 2**3 - 1), "2^3 - 1 = 7"
    print("  [ok] pseudo-count reconstruction")

    # 2. frame construction maps to Ensembl, drops unmapped, sets n_counts
    gene_map = pd.DataFrame({
        "entrez": ["1", "2", "3", "4"],
        "symbol": ["AAA", "BBB", "CCC", "DDD"],
        "canonical_label": labels,
    })
    ensembl_map = pd.DataFrame({
        "entrez": ["1", "2", "3"],          # gene 4 deliberately unmapped
        "ensembl_id": ["ENSG001", "ENSG002", "ENSG003"],
    })
    metadata = pd.DataFrame(
        {"OncotreeLineage": ["Bone"] * 6,
         "OncotreePrimaryDisease": ["Osteosarcoma"] * 6},
        index=idx,
    )
    assignment = pd.Series(
        ["train", "train", "val", "test", "train", "val"], index=idx
    )

    X, var_df, obs_df, rep = build_frames(
        counts, gene_map, ensembl_map, metadata, assignment
    )
    assert list(X.columns) == ["ENSG001", "ENSG002", "ENSG003"], \
        f"unexpected Ensembl columns: {list(X.columns)}"
    assert rep["genes_dropped_no_ensembl"] == 1, "gene 4 should have dropped"
    assert "n_counts" in obs_df.columns and "split" in obs_df.columns
    assert np.allclose(
        obs_df["n_counts"].to_numpy(),
        X.to_numpy().sum(axis=1),
    ), "n_counts must equal the row sum of the count matrix"
    assert var_df.index.name == "ensembl_id"
    print("  [ok] frame construction, Ensembl mapping, n_counts")

    # 3. duplicate-Ensembl collision is dropped, not duplicated
    ensembl_dup = pd.DataFrame({
        "entrez": ["1", "2", "3"],
        "ensembl_id": ["ENSG001", "ENSG001", "ENSG003"],  # 1 and 2 collide
    })
    X2, _, _, rep2 = build_frames(
        counts, gene_map, ensembl_dup, metadata, assignment
    )
    assert list(X2.columns) == ["ENSG001", "ENSG003"], \
        f"collision not handled: {list(X2.columns)}"
    print("  [ok] Ensembl collision handling")

    print("\nSelf-test passed.")
    return 0


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare Geneformer input from processed DepMap data."
    )
    parser.add_argument("--self-test", action="store_true",
                        help="Run the offline self-test and exit.")
    args = parser.parse_args()

    if args.self_test:
        return _self_test()

    out = config.PROCESSED_DIR
    print("=" * 74)
    print("PREPARING GENEFORMER INPUT")
    print("=" * 74)

    # ---- load the canonical gene record and split -----------------------
    try:
        gene_map = io_utils.load_table(out / "gene_id_map")
        metadata = io_utils.load_table(out / "model_metadata")
        gene_columns = io_utils.load_json(out / "gene_columns.json")
    except FileNotFoundError as exc:
        print(f"\n{exc}\nRun build_dataset.py first.")
        return 1

    metadata.index.name = config.MODEL_ID
    gene_map = gene_map.reset_index()
    if "entrez" not in gene_map.columns and "index" in gene_map.columns:
        gene_map = gene_map.rename(columns={"index": "entrez"})
    gene_map["entrez"] = gene_map["entrez"].astype(str)

    try:
        splits_payload = io_utils.load_json(out / "splits.json")
        assignment = pd.Series(splits_payload["assignment"], name="split")
    except FileNotFoundError:
        assignment = None
        print("Note: splits.json not found; the split column will be omitted.")

    report: dict = {}

    # ---- counts ---------------------------------------------------------
    print("\n[1/4] Expression counts")
    counts, report["counts"] = load_expression_counts()
    print(f"      source     : {report['counts']['count_source']}")
    print(f"      cell lines : {report['counts']['cell_lines']}")

    # ---- ensembl mapping ------------------------------------------------
    print("\n[2/4] Ensembl mapping")
    entrez_ids = [str(e) for e in gene_columns["entrez_ids"]]
    ensembl_map, report["ensembl"] = load_ensembl_map(entrez_ids)
    print(f"      mapped   : {report['ensembl']['entrez_mapped']} / "
          f"{report['ensembl']['entrez_requested']}")
    print(f"      unmapped : {report['ensembl']['entrez_unmapped']}")

    # ---- frames ---------------------------------------------------------
    print("\n[3/4] Building AnnData frames")
    X, var_df, obs_df, report["frames"] = build_frames(
        counts, gene_map, ensembl_map, metadata, assignment
    )
    fr = report["frames"]
    print(f"      genes kept        : {fr['genes_mapped_to_ensembl']}")
    print(f"      genes dropped     : {fr['genes_dropped_no_ensembl']}")
    print(f"      cells             : {fr['cells']}")
    print(f"      median n_counts   : {fr['median_n_counts']:.0f}")
    if fr["cells_with_zero_counts"]:
        print(f"      WARNING: {fr['cells_with_zero_counts']} cells have zero "
              f"total counts and will not tokenise")

    # ---- write ----------------------------------------------------------
    print("\n[4/4] Writing .h5ad")
    h5ad_path = out / "geneformer_input.h5ad"
    ok = write_h5ad(X, var_df, obs_df, h5ad_path)
    if ok:
        print(f"      wrote {h5ad_path.name}")

    report["output"] = {
        "h5ad_written": ok,
        "h5ad_path": str(h5ad_path) if ok else None,
        "model_series_note": (
            "Tokenise with model_input_size=4096 and special_token=True for the "
            "104M (V2) model; use model_version='V2' in EmbExtractor."
        ),
    }
    io_utils.save_json(report, out / "geneformer_prep_report.json")

    print(f"\n{'=' * 74}")
    print("DONE")
    print("=" * 74)
    print(f"  count source : {report['counts']['count_source']}")
    print(f"  genes        : {fr['genes_mapped_to_ensembl']} "
          f"(Ensembl-mapped, in canonical order)")
    print(f"  cells        : {fr['cells']}")
    if report["counts"]["count_source"] == "reconstructed_from_log_tpm":
        print("\n  NOTE: pseudo-counts were used. For a reported result, add "
              "DepMap's\n  expected-counts file and re-run. See README.")
    print("\n  Next: tokenise + extract embeddings on Kaggle "
          "(run_geneformer_embeddings.py)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
