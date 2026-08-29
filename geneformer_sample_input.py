"""
geneformer_sample_input.py
==========================
Deterministic Geneformer-input builder for ONE external RNA-seq sample -- the
Phase 2 demo tumour, BG003082. This is the *repo-local half* of the BG003082
embedding pipeline: it turns the committed GCT into a tokeniser-ready AnnData
(or, where ``anndata`` is absent, into the three validated component frames plus
a hard error telling you where to finish). It never calls Geneformer, never
touches a GPU, and never writes into the frozen Phase 1 artifact set.

    py geneformer_sample_input.py                 # build + validate, print a summary
    py geneformer_sample_input.py --json          # ...and dump the provenance record
    py geneformer_sample_input.py --write-h5ad    # also write the .h5ad (needs anndata)
    py geneformer_sample_input.py --self-test     # offline synthetic + real-data checks

Why this module exists
----------------------
`prepare_geneformer_input.py` builds a Geneformer input for the *1,140-line
DepMap matrix* and is hard-wired to `io_utils.load_matrix(PROCESSED_DIR /
"expression")`. It has no path for a single external row. `sample_profile.py`
reconciles one external GCT into the frozen 18,460-gene canonical space but
deliberately stops before any model concern: no Ensembl re-keying, no AnnData,
nothing written to `data/processed/`. This module is the missing bridge between
the two, for BG003082 only.

What it does, precisely
-----------------------
1. Calls ``sample_profile.load_external_sample(log2_transform=False)`` to get the
   canonical-space **linear TPM** vector for BG003082: 18,460 values in canonical
   order, 18,427 finite, 33 ``NaN`` (unresolved in the sample). Linear TPM is
   used **directly** as the Geneformer pseudo-count row. This is an approximation
   -- Geneformer wants raw counts -- and it is the *same* approximation the
   frozen 1,140 training embeddings used: those went
   ``log2(TPM+1) -> 2**x - 1``, which round-trips to linear TPM (float-exact to
   ~1e-9). Using linear TPM here is therefore the consistency-preserving choice,
   not a shortcut. It is disclosed in the report and the provenance record and in
   `capstone/geneformer-bg003082-feasibility.md`.
2. Re-keys canonical Entrez -> Ensembl gene ID through the committed
   ``data/processed/ensembl_map.csv`` (NCBI ``gene2ensembl`` provenance;
   18,459/18,460 Entrez carry an Ensembl ID). The one Entrez with no Ensembl ID
   is 79400 (NOX5) -- and NOX5 is **already** one of `sample_profile`'s 33
   unresolved genes, so re-keying drops **no additional gene**. The builder
   asserts this coincidence: any canonical gene that has a finite sample value
   but no Ensembl ID would be a 34th drop, and that is a hard error.
3. Assembles three frames, keyed on Ensembl ID, exactly as
   `prepare_geneformer_input.build_frames` keys them:
     * ``X``   -- 1 x 18,427, index ``["BG003082"]``, columns unique Ensembl IDs
       in canonical order, ``float64``, every value finite and >= 0.
     * ``var`` -- indexed by Ensembl ID (index name ``ensembl_id``) with an
       ``ensembl_id`` column, matching what the Kaggle tokeniser-side var patch
       expects.
     * ``obs`` -- index ``["BG003082"]`` (name ``ModelID``), columns ``ModelID``
       (the tokeniser's ``custom_attr_name_dict`` key) and ``n_counts`` set to
       the **exact retained-row sum** of ``X`` (~972,339; the 2.77% of TPM mass
       on non-coding GCT rows outside the canonical space is not in it).
4. Runs every assertion the task requires (see ``_validate``), then returns.
5. ``--write-h5ad`` assembles the AnnData with ``anndata`` and writes
   ``geneformer_bg003082_input.h5ad`` next to the other processed files. If
   ``anndata`` is not importable it raises ``GeneformerInputDependencyError``
   with a clear message and writes **nothing** -- no placeholder frames, no stub
   file. The frame-building and all validation still run and still pass; only the
   HDF5 serialisation is unavailable.

Determinism
-----------
The builder iterates canonical order taken from ``gene_columns.json``, never the
order of its inputs. It refuses (rather than silently reorders) a ``linear_tpm``
Series whose index is not exactly the canonical label list. `sample_profile`
itself is row-order-independent (its own self-test proves this). The self-test
here confirms the property end to end: two GCTs with identical rows in different
orders produce byte-identical ``X`` / ``var`` / ``obs`` values. Semantic
determinism is what is asserted; HDF5 containers are not byte-compared (timestamps
and internal layout vary between ``anndata`` / ``h5py`` builds).

Not done here
-------------
* No tokenisation, no embedding, no GPU, no network. That half runs only on
  Kaggle/Colab -- see ``capstone/kaggle_bg003082_embedding.py``.
* No symbol fallback. The 33 unresolved canonical genes stay unresolved, exactly
  as `sample_profile` leaves them (see that module's docstring for why the
  symbol pass is rejected).
* Nothing is written into `baseline_results.json`, `head_results.json`,
  `analysis_results.json`, `geneformer_embeddings.csv`, the split files, or any
  other frozen Phase 1 artifact.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import config
import sample_profile


# --------------------------------------------------------------------------
# Errors
# --------------------------------------------------------------------------

class GeneformerInputDependencyError(RuntimeError):
    """Raised when a required optional dependency (anndata) is unavailable."""


class GeneformerInputError(ValueError):
    """Raised when the assembled Geneformer input fails a validation gate."""


# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

_ENSEMBL_RE = re.compile(r"^ENSG\d+$")
_HASH_CHUNK = 1 << 20  # 1 MiB, matching capstone/data-integrity-hashes.md

# Expected reconciliation counts for the committed BG003082 + ensembl_map.csv.
# These are not magic numbers pulled from the air: they are what
# `sample_profile.load_external_sample()` reports on the committed bytes
# (verified: `py sample_profile.py --json`). If a committed input changes, these
# assertions fail loudly and this module -- not a downstream consumer -- is where
# the change surfaces.
EXPECTED_CANONICAL_INPUTS = 18_460
EXPECTED_MAPPED_VALUES = 18_427
EXPECTED_UNRESOLVED_VALUES = 33

# The Geneformer V2 104M model emits 768-dim CLS embeddings (see the committed
# geneformer_embeddings.csv header and the 2026-08-06 Kaggle run).
GENEFORMER_EMB_DIM = 768

PSEUDOCOUNT_DISCLOSURE = (
    "linear TPM used directly as the Geneformer pseudo-count row. Geneformer "
    "expects raw counts; this is an approximation. It is the SAME approximation "
    "the frozen 1,140 training embeddings used (log2(TPM+1) -> 2**x - 1 "
    "round-trips to linear TPM), so it is consistency-preserving, not a "
    "shortcut. Real RSEM expected_count was available for BG003082 and "
    "deliberately not used, for provenance consistency with the training set."
)


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------

def sha256_file(path: str | Path) -> str:
    """SHA-256 of a file, read in 1 MiB chunks (same recipe as the hash table)."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(_HASH_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def _entrez_to_ensembl(ensembl_map: pd.DataFrame) -> dict[str, str]:
    """
    Build a canonical Entrez -> Ensembl gene ID lookup from an ensembl_map frame.

    Enforces the schema `sample_profile` and `prepare_geneformer_input` both
    assume: exactly two columns `entrez,ensembl_id`, `entrez` unique (one row per
    canonical gene), `ensembl_id` unique (the reverse join must be unambiguous),
    and every `ensembl_id` a well-formed bare ENSG accession.
    """
    if list(ensembl_map.columns) != ["entrez", "ensembl_id"]:
        raise GeneformerInputError(
            f"ensembl_map columns are {list(ensembl_map.columns)}, expected "
            f"['entrez', 'ensembl_id']."
        )
    entrez = ensembl_map["entrez"].astype(str)
    ensembl = ensembl_map["ensembl_id"].astype(str)
    if entrez.duplicated().any():
        n = int(entrez.duplicated().sum())
        raise GeneformerInputError(
            f"ensembl_map has {n} duplicate Entrez ID(s); the Entrez->Ensembl "
            f"lookup would be ambiguous."
        )
    if ensembl.duplicated().any():
        n = int(ensembl.duplicated().sum())
        raise GeneformerInputError(
            f"ensembl_map has {n} duplicate Ensembl ID(s); two canonical genes "
            f"would collapse to one Geneformer token column."
        )
    malformed = sorted({e for e in ensembl if not _ENSEMBL_RE.match(e)})
    if malformed:
        raise GeneformerInputError(
            f"ensembl_map has {len(malformed)} malformed Ensembl ID(s) "
            f"(expected /^ENSG\\d+$/), e.g. {malformed[:5]}."
        )
    return dict(zip(entrez, ensembl))


# --------------------------------------------------------------------------
# The frame builder (the testable core -- no anndata, no I/O)
# --------------------------------------------------------------------------

def build_geneformer_frames(
    linear_tpm: pd.Series,
    ensembl_map: pd.DataFrame,
    gene_columns: dict,
    sample_id: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """
    Assemble and validate the three Geneformer-input frames for one sample.

    Parameters
    ----------
    linear_tpm
        Canonical-space **linear TPM**, indexed by the canonical
        ``SYMBOL (ENTREZ)`` labels in canonical order, with ``NaN`` for canonical
        genes the sample did not resolve. This is exactly what
        ``sample_profile.load_external_sample(log2_transform=False)`` returns.
    ensembl_map
        Frame with columns ``entrez,ensembl_id`` (both unique).
    gene_columns
        The ``gene_columns.json`` payload (keys ``entrez_ids``,
        ``canonical_columns``, ...).
    sample_id
        Label for the single row (``"BG003082"``).

    Returns
    -------
    (X, var_df, obs_df, report)
        See the module docstring. ``report`` is a JSON-serialisable dict of every
        reconciliation count and the pseudo-count disclosure.

    Raises
    ------
    GeneformerInputError
        On any schema or validation failure, including the "34th drop" guard: a
        canonical gene with a finite sample value but no Ensembl ID.
    """
    canonical_labels = list(gene_columns["canonical_columns"])
    entrez_ids = [str(e) for e in gene_columns["entrez_ids"]]
    n_canon = len(canonical_labels)
    if len(entrez_ids) != n_canon:
        raise GeneformerInputError(
            f"gene_columns: entrez_ids ({len(entrez_ids)}) and canonical_columns "
            f"({n_canon}) differ in length."
        )

    # Refuse a mis-ordered / mis-labelled input rather than silently reindexing.
    if list(linear_tpm.index) != canonical_labels:
        raise GeneformerInputError(
            "linear_tpm is not indexed by the canonical label list in canonical "
            "order; refusing to reorder silently. Pass the Series from "
            "sample_profile.load_external_sample() unmodified."
        )

    values = linear_tpm.to_numpy(dtype=np.float64)

    # A NaN here means "unresolved in the sample" and is expected (33 of them).
    # An inf or a negative value would be real corruption -> fail loudly.
    if np.isinf(values).any():
        raise GeneformerInputError(
            f"{int(np.isinf(values).sum())} non-finite (inf) value(s) in the "
            f"input TPM vector."
        )
    if np.any(values < 0.0):  # NaN < 0 is False, so NaNs do not trip this
        raise GeneformerInputError(
            f"{int(np.sum(values < 0.0))} negative value(s) in the input TPM "
            f"vector; TPM cannot be negative."
        )

    finite_mask = np.isfinite(values)
    e2e = _entrez_to_ensembl(ensembl_map)
    has_ensembl = np.array([e in e2e for e in entrez_ids], dtype=bool)

    # ---- the "34th drop" guard --------------------------------------------
    # Every canonical gene lacking an Ensembl ID must ALREADY be unresolved in
    # the sample (NaN). Otherwise re-keying would drop a gene that
    # sample_profile counted as mapped, silently making the unresolved set
    # larger than 33.
    extra_drops = [
        entrez_ids[i] for i in range(n_canon)
        if (not has_ensembl[i]) and finite_mask[i]
    ]
    if extra_drops:
        raise GeneformerInputError(
            f"gene(s) with Entrez {extra_drops} have a finite sample value but "
            f"no Ensembl ID in ensembl_map -- re-keying would drop them as a "
            f"34th+ unresolved gene beyond sample_profile's set of "
            f"{EXPECTED_UNRESOLVED_VALUES}. Refusing to proceed."
        )

    resolved_mask = finite_mask & has_ensembl
    n_resolved = int(resolved_mask.sum())
    n_unresolved = int(n_canon - n_resolved)

    kept_ensembl = [e2e[entrez_ids[i]] for i in range(n_canon) if resolved_mask[i]]
    kept_values = values[resolved_mask]

    if len(set(kept_ensembl)) != len(kept_ensembl):
        raise GeneformerInputError(
            "duplicate Ensembl IDs among resolved canonical genes; a token "
            "column would be doubled."
        )
    malformed = sorted({e for e in kept_ensembl if not _ENSEMBL_RE.match(e)})
    if malformed:
        raise GeneformerInputError(
            f"malformed Ensembl ID(s) among resolved genes: {malformed[:5]}."
        )
    if not np.all(np.isfinite(kept_values)):
        raise GeneformerInputError("non-finite value survived into X.")
    if np.any(kept_values < 0.0):
        raise GeneformerInputError("negative value survived into X.")

    # ---- frames ---------------------------------------------------------
    X = pd.DataFrame(
        kept_values.reshape(1, -1),
        index=pd.Index([sample_id], name=config.MODEL_ID),
        columns=kept_ensembl,
        dtype=np.float64,
    )
    n_counts = float(X.to_numpy().sum())

    var_df = pd.DataFrame(
        {"ensembl_id": kept_ensembl},
        index=pd.Index(kept_ensembl, name="ensembl_id"),
    )

    obs_df = pd.DataFrame(
        {"ModelID": [sample_id], "n_counts": [n_counts]},
        index=pd.Index([sample_id], name=config.MODEL_ID),
    )

    report = {
        "sample_id": sample_id,
        "canonical_inputs": n_canon,
        "mapped_values": n_resolved,
        "unresolved_values": n_unresolved,
        "genes_without_ensembl_id_in_map": int((~has_ensembl).sum()),
        "genes_without_ensembl_id_examples": sorted(
            entrez_ids[i] for i in range(n_canon) if not has_ensembl[i]
        )[:10],
        "extra_drops_beyond_sample_profile": len(extra_drops),
        "x_shape": list(X.shape),
        "x_has_nan": bool(np.isnan(X.to_numpy()).any()),
        "x_has_inf": bool(np.isinf(X.to_numpy()).any()),
        "x_has_negative": bool((X.to_numpy() < 0.0).any()),
        "ensembl_ids_unique": len(set(kept_ensembl)) == len(kept_ensembl),
        "ensembl_ids_well_formed": all(_ENSEMBL_RE.match(e) for e in kept_ensembl),
        "n_counts": n_counts,
        "n_counts_equals_retained_row_sum": True,  # by construction; re-checked in _validate
        "symbol_fallback_used": False,
        "pseudocount_disclosure": PSEUDOCOUNT_DISCLOSURE,
        "geneformer_emb_dim_expected": GENEFORMER_EMB_DIM,
    }

    _validate(X, var_df, obs_df, report)
    return X, var_df, obs_df, report


def _validate(
    X: pd.DataFrame, var_df: pd.DataFrame, obs_df: pd.DataFrame, report: dict
) -> None:
    """
    Hard gates. Every one of these corresponds to a requirement in the Phase 2
    task or an invariant a downstream consumer would otherwise trip over.
    """
    # 1. matrix content
    arr = X.to_numpy()
    if np.isnan(arr).any():
        raise GeneformerInputError("X contains NaN.")
    if np.isinf(arr).any():
        raise GeneformerInputError("X contains inf.")
    if (arr < 0.0).any():
        raise GeneformerInputError("X contains a negative value.")
    if X.shape[0] != 1:
        raise GeneformerInputError(f"X must have exactly one row, has {X.shape[0]}.")

    # 2. identifier uniqueness + format (columns of X, index of var)
    cols = list(X.columns)
    if len(set(cols)) != len(cols):
        raise GeneformerInputError("X has duplicate Ensembl-ID columns.")
    if any(not _ENSEMBL_RE.match(c) for c in cols):
        raise GeneformerInputError("X has a malformed Ensembl-ID column.")
    if list(var_df.index) != cols:
        raise GeneformerInputError("var index does not match X columns / order.")
    if "ensembl_id" not in var_df.columns:
        raise GeneformerInputError("var is missing the 'ensembl_id' column.")
    if list(var_df["ensembl_id"]) != cols:
        raise GeneformerInputError("var['ensembl_id'] does not match X columns.")
    if var_df.index.name != "ensembl_id":
        raise GeneformerInputError("var index name must be 'ensembl_id'.")

    # 3. metadata
    if list(obs_df.index) != list(X.index):
        raise GeneformerInputError("obs index does not match X index.")
    if "ModelID" not in obs_df.columns or "n_counts" not in obs_df.columns:
        raise GeneformerInputError("obs must carry 'ModelID' and 'n_counts'.")
    if obs_df["ModelID"].iloc[0] != X.index[0]:
        raise GeneformerInputError("obs['ModelID'] does not match the X row label.")
    retained_row_sum = float(arr.sum())
    if not np.isclose(obs_df["n_counts"].iloc[0], retained_row_sum, rtol=0, atol=1e-6):
        raise GeneformerInputError(
            f"obs['n_counts'] ({obs_df['n_counts'].iloc[0]}) != retained-row sum "
            f"of X ({retained_row_sum})."
        )
    if retained_row_sum <= 0.0:
        raise GeneformerInputError("retained-row sum is not positive.")

    # 4. mapping counts
    if report["canonical_inputs"] != X.shape[1] + report["unresolved_values"]:
        raise GeneformerInputError(
            "canonical_inputs != mapped + unresolved."
        )
    if report["mapped_values"] != X.shape[1]:
        raise GeneformerInputError(
            f"report mapped_values ({report['mapped_values']}) != X columns "
            f"({X.shape[1]})."
        )
    if report["extra_drops_beyond_sample_profile"] != 0:
        raise GeneformerInputError("re-keying dropped a gene beyond the 33.")
    if report["symbol_fallback_used"]:
        raise GeneformerInputError("symbol fallback must not be used.")


def assert_bg003082_counts(report: dict) -> None:
    """
    Assert the committed-input reconciliation counts exactly. Separated out so
    the real-data path and the self-test share one source of truth, and so a
    committed-input change fails here rather than three modules downstream.
    """
    checks = {
        "canonical_inputs": EXPECTED_CANONICAL_INPUTS,
        "mapped_values": EXPECTED_MAPPED_VALUES,
        "unresolved_values": EXPECTED_UNRESOLVED_VALUES,
    }
    for key, want in checks.items():
        got = report.get(key)
        if got != want:
            raise GeneformerInputError(
                f"BG003082 reconciliation drift: {key} = {got}, expected {want}. "
                f"A committed input (GCT or ensembl_map.csv) has changed."
            )
    if report["mapped_values"] + report["unresolved_values"] != report["canonical_inputs"]:
        raise GeneformerInputError("mapped + unresolved != canonical_inputs.")


# --------------------------------------------------------------------------
# Real-data wiring
# --------------------------------------------------------------------------

def build_bg003082_input(
    gct_path: str | Path = config.DEMO_TUMOR_GCT_FILE,
    ensembl_map_path: str | Path = config.ENSEMBL_MAP_FILE,
    gene_columns_path: str | Path | None = None,
    sample_id: str = config.DEMO_TUMOR_SAMPLE_ID,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """
    Build and validate the BG003082 Geneformer-input frames from committed files.

    Returns ``(X, var_df, obs_df, provenance)`` where ``provenance`` merges
    `sample_profile`'s reconciliation record, this module's frame report, and the
    SHA-256 / byte size of every input file. Deterministic: calling it twice on
    the same committed bytes returns semantically identical frames.
    """
    gene_columns_path = Path(
        gene_columns_path
        if gene_columns_path is not None
        else config.PROCESSED_DIR / "gene_columns.json"
    )

    series, sp_prov = sample_profile.load_external_sample(
        sample_id=sample_id,
        gct_path=gct_path,
        ensembl_map_path=ensembl_map_path,
        gene_columns_path=gene_columns_path,
        log2_transform=False,  # linear TPM as the disclosed pseudo-count row
    )
    ensembl_map = pd.read_csv(
        ensembl_map_path, dtype={"entrez": str, "ensembl_id": str}
    )
    gene_columns = json.loads(Path(gene_columns_path).read_text(encoding="utf-8"))

    X, var_df, obs_df, report = build_geneformer_frames(
        series, ensembl_map, gene_columns, sample_id
    )
    assert_bg003082_counts(report)

    # cross-check against sample_profile's independent count
    sp_mapped = sp_prov["reconciliation"]["canonical_genes_mapped"]
    sp_missing = sp_prov["reconciliation"]["canonical_genes_missing"]
    if (report["mapped_values"], report["unresolved_values"]) != (sp_mapped, sp_missing):
        raise GeneformerInputError(
            f"count mismatch vs sample_profile: this module says "
            f"({report['mapped_values']}, {report['unresolved_values']}), "
            f"sample_profile says ({sp_mapped}, {sp_missing})."
        )
    # n_counts must equal nansum of the linear series (the 33 NaNs drop out)
    expected_n_counts = float(np.nansum(series.to_numpy()))
    if not np.isclose(obs_df["n_counts"].iloc[0], expected_n_counts, rtol=0, atol=1e-3):
        raise GeneformerInputError(
            f"n_counts {obs_df['n_counts'].iloc[0]} != nansum(linear series) "
            f"{expected_n_counts}."
        )

    provenance = {
        "sample_id": sample_id,
        "builder": "geneformer_sample_input.build_bg003082_input",
        "inputs": {
            "gct_file": {
                "name": Path(gct_path).name,
                "sha256": sha256_file(gct_path),
                "bytes": Path(gct_path).stat().st_size,
            },
            "ensembl_map_file": {
                "name": Path(ensembl_map_path).name,
                "sha256": sha256_file(ensembl_map_path),
                "bytes": Path(ensembl_map_path).stat().st_size,
                "provenance": "NCBI gene2ensembl (see capstone/data-integrity-hashes.md)",
            },
            "gene_columns_file": {
                "name": Path(gene_columns_path).name,
                "sha256": sha256_file(gene_columns_path),
            },
        },
        "sample_profile_reconciliation": sp_prov["reconciliation"],
        "frame_report": report,
        "n_counts": float(obs_df["n_counts"].iloc[0]),
        "x_shape": list(X.shape),
        "pseudocount_disclosure": PSEUDOCOUNT_DISCLOSURE,
        "not_produced_here": (
            "tokenisation and the 1x768 CLS embedding require Geneformer + a GPU "
            "and are NOT produced by this module. Run "
            "capstone/kaggle_bg003082_embedding.py on Kaggle/Colab. Until then "
            "the repository is in the documented baseline-only fallback state "
            "(see capstone/geneformer-bg003082-feasibility.md)."
        ),
    }
    return X, var_df, obs_df, provenance


# --------------------------------------------------------------------------
# AnnData write (thin, optional, no placeholder on failure)
# --------------------------------------------------------------------------

def write_h5ad(
    X: pd.DataFrame,
    var_df: pd.DataFrame,
    obs_df: pd.DataFrame,
    path: str | Path,
) -> Path:
    """
    Assemble and write the AnnData for the single BG003082 row.

    Raises ``GeneformerInputDependencyError`` -- and writes nothing at all -- if
    ``anndata`` (or ``scipy.sparse``) is unavailable. No component-frame
    placeholder is written; this is a deliberate departure from
    `prepare_geneformer_input.write_h5ad`, which does fall back to frames.
    """
    try:
        import anndata as ad
        from scipy import sparse
    except ImportError as exc:
        raise GeneformerInputDependencyError(
            "anndata (and scipy.sparse) are required to write the .h5ad and are "
            "not importable in this environment. The frames were built and fully "
            "validated; only HDF5 serialisation is unavailable. Finish this step "
            "on a machine with anndata (`pip install anndata`) or run "
            "capstone/kaggle_bg003082_embedding.py on Kaggle. No placeholder "
            "file has been written."
        ) from exc

    path = Path(path)
    adata = ad.AnnData(
        X=sparse.csr_matrix(X.to_numpy(dtype=np.float32)),
        obs=obs_df.copy(),
        var=var_df.copy(),
    )
    adata.write_h5ad(path)
    return path


# --------------------------------------------------------------------------
# Offline self-test (synthetic + real committed data; no anndata / GPU / network)
# --------------------------------------------------------------------------

def _write_gct(path: Path, rows: list[tuple[str, str, str]],
               header: tuple[str, ...] = ("Name", "Description", "S1")) -> None:
    n_rows, n_cols = len(rows), len(header) - 2
    body = "\n".join("\t".join(r) for r in rows)
    path.write_text(
        f"#1.2\n{n_rows}\t{n_cols}\n" + "\t".join(header) + "\n" + body + "\n",
        encoding="utf-8",
    )


def _write_canonical(path: Path, pairs: list[tuple[str, str]]) -> None:
    path.write_text(
        json.dumps({
            "description": "self-test canonical space",
            "n_genes": len(pairs),
            "entrez_ids": [e for _, e in pairs],
            "symbols": [s for s, _ in pairs],
            "canonical_columns": [f"{s} ({e})" for s, e in pairs],
        }),
        encoding="utf-8",
    )


def _write_ensembl_map(path: Path, pairs: list[tuple[str, str]]) -> None:
    path.write_text(
        "\n".join(["entrez,ensembl_id"] + [f"{e},{g}" for e, g in pairs]) + "\n",
        encoding="utf-8",
    )


def _frames_equal(a: tuple, b: tuple) -> bool:
    """Semantic equality of (X, var, obs) triples -- values + labels, not dtype id."""
    Xa, va, oa = a
    Xb, vb, ob = b
    return (
        list(Xa.index) == list(Xb.index)
        and list(Xa.columns) == list(Xb.columns)
        and np.allclose(Xa.to_numpy(), Xb.to_numpy(), rtol=0, atol=0, equal_nan=True)
        and list(va.index) == list(vb.index)
        and list(va["ensembl_id"]) == list(vb["ensembl_id"])
        and list(oa.index) == list(ob.index)
        and list(oa["ModelID"]) == list(ob["ModelID"])
        and np.allclose(oa["n_counts"].to_numpy(), ob["n_counts"].to_numpy())
    )


def _self_test() -> int:
    import random
    import tempfile

    print("Running geneformer_sample_input self-test...")
    tmp = Path(tempfile.mkdtemp(prefix="gf_sample_input_selftest_"))

    # ---- synthetic canonical space of 5 genes ---------------------------
    # GENE5 (entrez 50) is deliberately absent from the ensembl map AND never
    # supplied in the GCT -> it is the synthetic stand-in for NOX5: an
    # unmapped gene that is ALSO unresolved, so no "extra drop".
    canon = tmp / "gene_columns.json"
    _write_canonical(canon, [("GENE1", "10"), ("GENE2", "20"), ("GENE3", "30"),
                             ("GENE4", "40"), ("GENE5", "50")])
    emap = tmp / "ensembl_map.csv"
    _write_ensembl_map(emap, [("10", "ENSG00000000001"), ("20", "ENSG00000000002"),
                              ("30", "ENSG00000000003"), ("40", "ENSG00000000004")])

    rows = [
        ("ENSG00000000001.3", "GENE1", "3.0"),    # -> 3.0
        ("ENSG00000000002.1", "GENE2", "0.0"),    # measured zero -> kept, 0.0
        ("ENSG00000000003", "GENE3", "7.0"),      # bare -> 7.0
        ("ENSG99999999999.2", "SOMETHING", "5.0"),  # unresolved, not canonical
        # GENE4 (entrez 40) has an Ensembl ID in the map but is NOT in the GCT
        #   -> unresolved via NaN, not via missing Ensembl. GENE5 has neither.
    ]
    good = tmp / "good.gct"
    _write_gct(good, rows)

    gene_columns = json.loads(canon.read_text())
    ensembl_map = pd.read_csv(emap, dtype={"entrez": str, "ensembl_id": str})

    series, _ = sample_profile.load_external_sample(
        sample_id="S1", gct_path=good, ensembl_map_path=emap,
        gene_columns_path=canon, log2_transform=False,
    )
    X, var_df, obs_df, report = build_geneformer_frames(
        series, ensembl_map, gene_columns, "S1"
    )
    assert list(X.columns) == ["ENSG00000000001", "ENSG00000000002", "ENSG00000000003"], list(X.columns)
    assert X.shape == (1, 3), X.shape
    assert np.allclose(X.to_numpy().ravel(), [3.0, 0.0, 7.0])
    assert not np.isnan(X.to_numpy()).any()
    assert not np.isinf(X.to_numpy()).any()
    assert (X.to_numpy() >= 0).all()
    assert report["canonical_inputs"] == 5
    assert report["mapped_values"] == 3
    assert report["unresolved_values"] == 2
    assert report["extra_drops_beyond_sample_profile"] == 0
    assert report["symbol_fallback_used"] is False
    assert list(var_df.index) == list(X.columns)
    assert var_df.index.name == "ensembl_id"
    assert list(var_df["ensembl_id"]) == list(X.columns)
    assert obs_df.loc["S1", "ModelID"] == "S1"
    assert np.isclose(obs_df.loc["S1", "n_counts"], 10.0)  # 3 + 0 + 7
    assert np.isclose(obs_df.loc["S1", "n_counts"], float(X.to_numpy().sum()))
    print("  [ok] matrix, metadata, mapping counts, n_counts, identifier checks")

    # ---- the "34th drop" guard fires ----------------------------------
    # An ensembl map that is MISSING entrez 30, whose GENE3 the GCT DID measure.
    emap_missing30 = tmp / "emap_missing30.csv"
    _write_ensembl_map(emap_missing30, [("10", "ENSG00000000001"),
                                        ("20", "ENSG00000000002"),
                                        ("40", "ENSG00000000004")])
    try:
        build_geneformer_frames(
            series, pd.read_csv(emap_missing30, dtype=str), gene_columns, "S1"
        )
        raise AssertionError("missing-Ensembl-for-a-measured-gene should have raised")
    except GeneformerInputError as exc:
        assert "34th" in str(exc) or "beyond" in str(exc), str(exc)
    print("  [ok] extra-drop guard: a measured gene with no Ensembl ID is a hard error")

    # ---- mis-ordered input is rejected, not silently reindexed --------
    shuffled_series = series.iloc[::-1]
    try:
        build_geneformer_frames(shuffled_series, ensembl_map, gene_columns, "S1")
        raise AssertionError("a non-canonical-order Series should have raised")
    except GeneformerInputError as exc:
        assert "canonical" in str(exc)
    print("  [ok] non-canonical-order input rejected")

    # ---- duplicate Ensembl ID in the map is a hard error -------------
    emap_dup = tmp / "emap_dup.csv"
    _write_ensembl_map(emap_dup, [("10", "ENSG00000000001"),
                                  ("20", "ENSG00000000001"),  # collide
                                  ("30", "ENSG00000000003"),
                                  ("40", "ENSG00000000004")])
    try:
        build_geneformer_frames(
            series, pd.read_csv(emap_dup, dtype=str), gene_columns, "S1"
        )
        raise AssertionError("duplicate Ensembl ID should have raised")
    except GeneformerInputError as exc:
        assert "duplicate Ensembl" in str(exc)
    print("  [ok] duplicate Ensembl ID in the map rejected")

    # ---- deterministic preprocessing: GCT row order cannot change output ----
    shuffled_rows = list(rows)
    random.Random(0).shuffle(shuffled_rows)
    shuf = tmp / "shuffled.gct"
    _write_gct(shuf, shuffled_rows)
    series_shuf, _ = sample_profile.load_external_sample(
        sample_id="S1", gct_path=shuf, ensembl_map_path=emap,
        gene_columns_path=canon, log2_transform=False,
    )
    frames_shuf = build_geneformer_frames(series_shuf, ensembl_map, gene_columns, "S1")
    assert _frames_equal((X, var_df, obs_df), frames_shuf[:3]), \
        "GCT row order changed the built frames"
    # also: shuffle the ensembl_map row order
    ensembl_map_shuf = ensembl_map.sample(frac=1.0, random_state=1).reset_index(drop=True)
    frames_emap_shuf = build_geneformer_frames(
        series, ensembl_map_shuf, gene_columns, "S1"
    )
    assert _frames_equal((X, var_df, obs_df), frames_emap_shuf[:3]), \
        "ensembl_map row order changed the built frames"
    print("  [ok] deterministic: neither GCT nor ensembl_map row order affects output")

    # ---- write_h5ad: clear dependency error, no placeholder ----------
    try:
        import anndata  # noqa: F401
        has_anndata = True
    except ImportError:
        has_anndata = False
    target = tmp / "out.h5ad"
    if has_anndata:
        write_h5ad(X, var_df, obs_df, target)
        assert target.is_file()
        print("  [ok] write_h5ad wrote the .h5ad (anndata present)")
    else:
        try:
            write_h5ad(X, var_df, obs_df, target)
            raise AssertionError("write_h5ad should have raised without anndata")
        except GeneformerInputDependencyError as exc:
            assert "placeholder" in str(exc).lower()
        assert not target.exists(), "no file may be written when anndata is absent"
        print("  [ok] write_h5ad raises a clear dependency error and writes nothing")

    # ---- real committed BG003082 data -------------------------------
    frozen_before = _frozen_hashes()
    X1, v1, o1, prov1 = build_bg003082_input()
    X2, v2, o2, prov2 = build_bg003082_input()
    frozen_after = _frozen_hashes()
    assert frozen_before == frozen_after, "a frozen Phase 1 artifact changed during the build"

    assert _frames_equal((X1, v1, o1), (X2, v2, o2)), "real build is not deterministic"
    assert X1.shape == (1, EXPECTED_MAPPED_VALUES), X1.shape
    assert list(X1.index) == ["BG003082"]
    assert prov1["frame_report"]["canonical_inputs"] == EXPECTED_CANONICAL_INPUTS
    assert prov1["frame_report"]["mapped_values"] == EXPECTED_MAPPED_VALUES
    assert prov1["frame_report"]["unresolved_values"] == EXPECTED_UNRESOLVED_VALUES
    assert not np.isnan(X1.to_numpy()).any()
    assert not np.isinf(X1.to_numpy()).any()
    assert (X1.to_numpy() >= 0).all()
    assert len(set(X1.columns)) == X1.shape[1]
    assert all(_ENSEMBL_RE.match(c) for c in X1.columns)
    assert list(v1.index) == list(X1.columns)
    assert np.isclose(o1.loc["BG003082", "n_counts"], float(X1.to_numpy().sum()))
    # NOX5 (79400) is the ONLY canonical Entrez with no Ensembl ID, and it is
    # part of the 33 -- not a 34th drop.
    assert prov1["frame_report"]["genes_without_ensembl_id_in_map"] == 1
    assert prov1["frame_report"]["genes_without_ensembl_id_examples"] == ["79400"]
    assert prov1["frame_report"]["extra_drops_beyond_sample_profile"] == 0
    print(
        f"  [ok] real BG003082 build: X {X1.shape}, "
        f"n_counts {o1.loc['BG003082', 'n_counts']:.1f}, "
        f"{EXPECTED_UNRESOLVED_VALUES} unresolved (NOX5 included, not extra)"
    )

    for child in sorted(tmp.rglob("*"), reverse=True):
        child.unlink() if child.is_file() else child.rmdir()
    tmp.rmdir()
    print("\nSelf-test passed.")
    return 0


def _frozen_hashes() -> dict:
    """SHA-256 of the frozen Phase 1 artifacts this module must never touch."""
    frozen = [
        config.PROCESSED_DIR / "baseline_results.json",
        config.PROCESSED_DIR / "head_results.json",
        config.PROCESSED_DIR / "analysis_results.json",
        config.PROCESSED_DIR / "geneformer_embeddings.csv",
        config.PROCESSED_DIR / "expression.npz",
        config.PROCESSED_DIR / "crispr_effect.npz",
        config.PROCESSED_DIR / "splits.json",
        config.PROCESSED_DIR / "gene_columns.json",
        config.ENSEMBL_MAP_FILE,
        config.DEMO_TUMOR_GCT_FILE,
    ]
    return {p.name: sha256_file(p) for p in frozen if p.is_file()}


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def _print_summary(X, var_df, obs_df, provenance) -> None:
    r = provenance["frame_report"]
    print("=" * 74)
    print(f"GENEFORMER SAMPLE INPUT  --  {provenance['sample_id']}")
    print("=" * 74)
    print(f"  GCT file          : {provenance['inputs']['gct_file']['name']}")
    print(f"    sha256          : {provenance['inputs']['gct_file']['sha256']}")
    print(f"  ensembl_map       : {provenance['inputs']['ensembl_map_file']['name']}")
    print(f"    sha256          : {provenance['inputs']['ensembl_map_file']['sha256']}")
    print(f"    provenance      : {provenance['inputs']['ensembl_map_file']['provenance']}")
    print("-" * 74)
    print(f"  canonical inputs                 : {r['canonical_inputs']}")
    print(f"  mapped values (X columns)        : {r['mapped_values']}")
    print(f"  unresolved values (left out)     : {r['unresolved_values']}")
    print(f"  canonical genes with no Ensembl  : {r['genes_without_ensembl_id_in_map']} "
          f"{r['genes_without_ensembl_id_examples']}  (part of the unresolved set)")
    print(f"  extra drops beyond sample_profile: {r['extra_drops_beyond_sample_profile']}")
    print(f"  symbol fallback used             : {r['symbol_fallback_used']}")
    print("-" * 74)
    print(f"  X shape                          : {tuple(X.shape)}")
    print(f"  X has NaN / inf / negative       : {r['x_has_nan']} / {r['x_has_inf']} / {r['x_has_negative']}")
    print(f"  Ensembl IDs unique / well-formed : {r['ensembl_ids_unique']} / {r['ensembl_ids_well_formed']}")
    print(f"  obs n_counts (retained-row sum)  : {obs_df['n_counts'].iloc[0]:.4f}")
    print(f"  expected CLS embedding dim       : {r['geneformer_emb_dim_expected']}")
    print("-" * 74)
    print(f"  pseudo-count basis : {r['pseudocount_disclosure']}")
    print("=" * 74)
    print("  NOTE: this module does NOT tokenise and does NOT produce an "
          "embedding.")
    print("  Run capstone/kaggle_bg003082_embedding.py on Kaggle/Colab for that "
          "half.")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build + validate the BG003082 Geneformer input (repo-local "
                    "half; no tokenisation, no embedding)."
    )
    parser.add_argument("--self-test", action="store_true",
                        help="Run the offline self-test and exit.")
    parser.add_argument("--json", action="store_true",
                        help="Print the full provenance record as JSON.")
    parser.add_argument("--write-h5ad", action="store_true",
                        help="Also write geneformer_bg003082_input.h5ad (requires "
                             "anndata; errors clearly if absent, writes no "
                             "placeholder).")
    parser.add_argument("--h5ad-path", default=None,
                        help="Override the .h5ad output path.")
    args = parser.parse_args()

    if args.self_test:
        return _self_test()

    X, var_df, obs_df, provenance = build_bg003082_input()
    _print_summary(X, var_df, obs_df, provenance)

    if args.write_h5ad:
        path = Path(args.h5ad_path) if args.h5ad_path else (
            config.PROCESSED_DIR / "geneformer_bg003082_input.h5ad"
        )
        try:
            written = write_h5ad(X, var_df, obs_df, path)
        except GeneformerInputDependencyError as exc:
            print("\n" + "=" * 74)
            print("H5AD NOT WRITTEN")
            print("=" * 74)
            print(f"  {exc}")
            print("\n  Frame-building and every validation gate PASSED. Only the "
                  "HDF5 write is unavailable here.")
            return 1
        print(f"\n  wrote {written}")

    if args.json:
        print(json.dumps(provenance, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
