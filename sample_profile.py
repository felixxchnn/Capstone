"""
sample_profile.py
=================
Load one *external* RNA-seq sample into the frozen canonical gene space.

    py sample_profile.py                # load the committed BG003082 tumour GCT
    py sample_profile.py --json         # ...and dump the full provenance record
    py sample_profile.py --self-test    # offline schema/edge-case checks, no network

Why this module exists
----------------------
Phase 1 froze an 18,460-gene feature space (`data/processed/gene_columns.json`),
keyed on Entrez ID, produced from DepMap's `SYMBOL (ENTREZ)` column labels. Any
expression matrix scored by a Phase 1 model has to be reindexed into *exactly*
that column order or the model silently receives scrambled features and returns
confident nonsense (see `gene_ids.py`'s module docstring).

The Phase 2 demo (`capstone/scope-decisions.md`, 2026-08-25) adds one external
sample the frozen space was never built from: **BG003082**, Sid Sijbrandij's
self-released osteosarcoma primary-tumour RNA-seq (`osteosarc.com`, CC0 1.0),
shipped as a gzip-compressed GCT v1.2 of *linear* gene-level TPM keyed on
*versioned Ensembl* gene IDs. This module is the one place that turns that file
into a canonical-space vector, with every reconciliation count recorded.

What it does, precisely
-----------------------
1. Reads the compressed GCT directly and validates it against the GCT v1.2
   schema: the ``#1.2`` version tag, the ``<n_rows>\\t<n_cols>`` dimension line,
   that the declared dimensions match the actual header width and data-row
   count, and that the required ``Name`` / ``Description`` columns and the
   requested sample column are present. Any violation raises ``GCTFormatError``
   -- it never guesses past a malformed file.
2. Takes the TPM sample column (not any RSEM / expected-count companion file --
   none is committed, and mixing quantifiers is a documented hazard, see
   `prepare_geneformer_input.py`). Non-finite or negative TPM raises
   ``ExternalSampleError``.
3. Strips the Ensembl version suffix (``ENSG00000123456.7`` -> ``ENSG00000123456``)
   where present; anything that is not a recognisable ``ENSG`` accession is left
   untouched and counted as an unexpected identifier.
4. Detects colliding identifiers explicitly and never discards silently:
     * two GCT rows whose *stripped* Ensembl IDs are equal -> summed as linear
       TPM, reported under ``duplicate_external_ids``;
     * two distinct Ensembl IDs that resolve (via ``ensembl_map.csv``) to the
       same canonical Entrez ID -> summed as linear TPM, reported under
       ``canonical_id_collisions``;
     * an ``ensembl_map.csv`` whose ``ensembl_id`` column is not unique ->
       ``ExternalSampleError`` (the Ensembl->Entrez join would be ambiguous).
   Linear TPM is additive across transcripts of one gene; that is the only case
   in which values are ever combined. Rows are never merged merely because they
   share a text symbol.
5. Resolves each stripped Ensembl ID to a canonical Entrez ID through
   ``data/processed/ensembl_map.csv`` (schema ``entrez,ensembl_id``, the cache
   `prepare_geneformer_input.load_ensembl_map` already expects; sourced from
   NCBI ``gene2ensembl`` -- see `capstone/data-integrity-hashes.md`). Ensembl ID
   is the only join key used.
6. Reindexes to the canonical Entrez order read straight from
   ``gene_columns.json`` -- deterministic, independent of GCT row order.
7. Keeps the distinction between *measured* and *absent*: a canonical gene the
   sample measured at TPM 0 becomes ``log2(0 + 1) = 0.0``; a canonical gene the
   sample never mentions (or that does not resolve) stays ``NaN``. Missing
   values are left explicit as ``NaN`` -- the same representation every Phase 1
   evaluation row carries into `baseline.impute_with_train_mean`. This module
   does not impute; that is a model-fit step and belongs downstream.
8. Applies ``log2(TPM + 1)`` (matching `data/processed/expression.npz`, which is
   ``log2(TPM + 1)``; confirmed in `build_dataset.py` and independently from the
   committed value range). Pass ``log2_transform=False`` for the linear vector.
9. Returns ``(pandas.Series, provenance_dict)``. The Series is indexed by the
   canonical ``SYMBOL (ENTREZ)`` labels, in canonical order, named after the
   sample. ``provenance_dict`` records the source file, the mapping file, their
   SHA-256 digests and byte sizes, the GCT header facts, the transformation
   string, and every reconciliation count.

Deliberately *not* done here
----------------------------
* **No symbol fallback.** The plan (`moonlit-dazzling-dream.md`, step 4) allowed
  matching the GCT ``Description`` (symbol) column against `gene_id_map.csv` for
  the residual that Ensembl-ID join misses. For BG003082 that residual is 33
  canonical genes and a symbol pass would rescue only 4 of them (ASPRV1, PAXX,
  NOX5, FAM174C -- see ``symbol_fallback_candidates`` in the provenance record).
  It is not implemented because (a) `capstone/data-integrity-hashes.md` records a
  deliberate decision to keep ``ensembl_map.csv`` single-provenance (NCBI only),
  explicitly *not* patching in NOX5's known Ensembl ID ``ENSG00000255346``; a
  ``Description``-symbol fallback would resurrect exactly that gene through a
  second provenance path, undoing that decision; and (b) the GCT's
  ``Description`` column is not unique -- 231 symbols occur on more than one row
  -- so a symbol pass is materially ambiguous in this file. The 33 unresolved
  canonical genes are left as explicit ``NaN`` and counted. A future
  `sample_profile` extension can revisit this as an additive change.
* **No imputation, no model call, no file written to `data/processed/`.**
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import config


# --------------------------------------------------------------------------
# Errors
# --------------------------------------------------------------------------

class GCTFormatError(ValueError):
    """Raised when a .gct file violates the GCT v1.2 schema."""


class ExternalSampleError(ValueError):
    """Raised when an external expression file fails a sanity or integrity gate."""


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------

_ENSEMBL_VERSIONED_RE = re.compile(r"^(ENSG\d+)\.\d+$")
_ENSEMBL_BARE_RE = re.compile(r"^ENSG\d+$")

_HASH_CHUNK = 1 << 20  # 1 MiB, matching capstone/data-integrity-hashes.md


def sha256_file(path: str | Path) -> str:
    """SHA-256 of a file, read in 1 MiB chunks (same recipe as the hash table)."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(_HASH_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def strip_ensembl_version(identifier: str) -> str:
    """
    ``'ENSG00000123456.7'`` -> ``'ENSG00000123456'``.

    An already-bare ``ENSG`` accession, or any string that is not a recognisable
    Ensembl gene ID, is returned unchanged so the caller can count and report
    identifiers that never matched rather than mangling them.
    """
    match = _ENSEMBL_VERSIONED_RE.match(identifier)
    return match.group(1) if match else identifier


# --------------------------------------------------------------------------
# GCT v1.2 parser
# --------------------------------------------------------------------------

def parse_gct(path: str | Path) -> tuple[pd.DataFrame, dict]:
    """
    Parse a (optionally gzip-compressed) GCT v1.2 file.

    Returns
    -------
    (frame, info)
        ``frame`` has columns ``Name``, ``Description`` then one column per
        sample, all as strings (values are parsed by the caller so that a
        malformed value can be reported against the right sample).
        ``info`` records the version tag, the declared and actual dimensions,
        and the sample column names.

    Raises
    ------
    GCTFormatError
        On any schema violation: wrong version tag, malformed dimension line,
        header width or content that disagrees with the dimension line, a data
        row with the wrong number of fields, or a data-row count that does not
        match the declared row count.
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        text = handle.read()

    # Normalise CRLF / CR line endings, then split on LF only (a bare split
    # avoids str.splitlines() also breaking on U+000B/U+000C/U+0085 etc.).
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    lines = text.split("\n")
    # A trailing newline produces one empty final element; drop trailing blanks.
    while lines and lines[-1] == "":
        lines.pop()

    if len(lines) < 3:
        raise GCTFormatError(
            f"{path.name}: only {len(lines)} non-empty line(s); a GCT needs at "
            f"least a version line, a dimension line and a header."
        )

    version = lines[0].strip()
    if version != "#1.2":
        raise GCTFormatError(
            f"{path.name}: version line is {version!r}, expected '#1.2'. "
            f"Only GCT v1.2 is supported."
        )

    dim_fields = lines[1].split("\t")
    if len(dim_fields) != 2 or not all(f.strip().isdigit() for f in dim_fields):
        raise GCTFormatError(
            f"{path.name}: dimension line is {lines[1]!r}, expected "
            f"'<n_rows>\\t<n_cols>' with two non-negative integers."
        )
    declared_rows, declared_cols = (int(f) for f in dim_fields)

    header = lines[2].split("\t")
    if len(header) != declared_cols + 2:
        raise GCTFormatError(
            f"{path.name}: header has {len(header)} field(s), but the dimension "
            f"line implies {declared_cols + 2} (Name, Description and "
            f"{declared_cols} sample column(s))."
        )
    if header[0] != "Name" or header[1] != "Description":
        raise GCTFormatError(
            f"{path.name}: header begins {header[:2]!r}, expected "
            f"['Name', 'Description']."
        )
    sample_columns = header[2:]

    data_lines = lines[3:]
    if len(data_lines) != declared_rows:
        raise GCTFormatError(
            f"{path.name}: dimension line declares {declared_rows} data row(s) "
            f"but the file has {len(data_lines)}."
        )

    rows: list[list[str]] = []
    width = declared_cols + 2
    for offset, line in enumerate(data_lines):
        fields = line.split("\t")
        if len(fields) != width:
            raise GCTFormatError(
                f"{path.name}: data line {offset + 4} has {len(fields)} "
                f"field(s), expected {width}."
            )
        rows.append(fields)

    frame = pd.DataFrame(rows, columns=header, dtype="object")
    info = {
        "version": version,
        "declared_rows": declared_rows,
        "declared_cols": declared_cols,
        "actual_rows": len(data_lines),
        "sample_columns": list(sample_columns),
    }
    return frame, info


# --------------------------------------------------------------------------
# The loader
# --------------------------------------------------------------------------

def load_external_sample(
    sample_id: str = config.DEMO_TUMOR_SAMPLE_ID,
    gct_path: str | Path = config.DEMO_TUMOR_GCT_FILE,
    ensembl_map_path: str | Path = config.ENSEMBL_MAP_FILE,
    gene_columns_path: str | Path | None = None,
    log2_transform: bool = True,
) -> tuple[pd.Series, dict]:
    """
    Load an external GCT sample into the frozen canonical gene space.

    Parameters
    ----------
    sample_id
        Name of the sample column to read from the GCT. Defaults to the Phase 2
        demo tumour sample.
    gct_path
        Path to the (gzip-compressed) GCT v1.2 file.
    ensembl_map_path
        CSV with columns ``entrez,ensembl_id`` (one row per gene, both unique).
    gene_columns_path
        ``gene_columns.json`` defining the canonical Entrez order. Defaults to
        the committed Phase 1 artifact.
    log2_transform
        If True (default) return ``log2(TPM + 1)``, matching
        ``data/processed/expression.npz``. If False, return linear TPM.

    Returns
    -------
    (series, provenance)
        ``series`` is indexed by the canonical ``SYMBOL (ENTREZ)`` labels in
        canonical order, named ``sample_id``; canonical genes the sample does
        not resolve are ``NaN``. ``provenance`` is a JSON-serialisable dict of
        source identities, digests, GCT facts, the transformation applied and
        every reconciliation count.
    """
    gct_path = Path(gct_path)
    ensembl_map_path = Path(ensembl_map_path)
    gene_columns_path = Path(
        gene_columns_path
        if gene_columns_path is not None
        else config.PROCESSED_DIR / "gene_columns.json"
    )
    for required in (gct_path, ensembl_map_path, gene_columns_path):
        if not required.is_file():
            raise FileNotFoundError(f"required input not found: {required}")

    # ---- 1. parse + schema-validate the GCT --------------------------------
    frame, gct_info = parse_gct(gct_path)
    if sample_id not in gct_info["sample_columns"]:
        raise GCTFormatError(
            f"{gct_path.name}: no sample column named {sample_id!r}; present "
            f"sample column(s): {gct_info['sample_columns']}."
        )

    names = frame["Name"].astype(str)
    descriptions = frame["Description"].astype(str)

    # ---- 2. parse the TPM column, reject non-finite / negative -----------
    tpm_raw = frame[sample_id]
    tpm = pd.to_numeric(tpm_raw, errors="coerce").to_numpy(dtype=np.float64)
    nonfinite = ~np.isfinite(tpm)
    if nonfinite.any():
        bad = np.flatnonzero(nonfinite)[:10]
        raise ExternalSampleError(
            f"{gct_path.name}: {int(nonfinite.sum())} non-finite value(s) in "
            f"sample column {sample_id!r} (0-based data rows "
            f"{[int(b) for b in bad]}, raw {[str(tpm_raw.iloc[b]) for b in bad]}). "
            f"Refusing to guess a replacement."
        )
    negative = tpm < 0.0
    if negative.any():
        bad = np.flatnonzero(negative)[:10]
        raise ExternalSampleError(
            f"{gct_path.name}: {int(negative.sum())} negative TPM value(s) in "
            f"sample column {sample_id!r} (0-based data rows "
            f"{[int(b) for b in bad]}). TPM cannot be negative."
        )

    # ---- 3. normalise Ensembl IDs; classify anything that is not one -----
    versioned_mask = names.str.match(_ENSEMBL_VERSIONED_RE).to_numpy()
    bare_mask = names.str.match(_ENSEMBL_BARE_RE).to_numpy()
    unexpected_mask = ~(versioned_mask | bare_mask)
    unexpected_identifiers = sorted(set(names[unexpected_mask]))

    ens_base = names.map(strip_ensembl_version)

    work = pd.DataFrame(
        {
            "ens": ens_base.to_numpy(),
            "tpm": tpm,
            "recognised": ~unexpected_mask,
        }
    )
    # Identifiers that are not Ensembl accessions cannot join on Ensembl ID.
    work = work[work["recognised"]].drop(columns="recognised")

    # ---- 4a. collapse identical stripped Ensembl IDs (linear TPM sums) ---
    dup_mask = work["ens"].duplicated(keep=False)
    duplicate_external_ids = sorted(set(work.loc[dup_mask, "ens"]))
    collapsed = (
        work.groupby("ens", sort=False, as_index=False)["tpm"].sum()
    )

    # ---- 5. resolve Ensembl ID -> canonical Entrez ID -------------------
    emap = pd.read_csv(
        ensembl_map_path, dtype={"entrez": str, "ensembl_id": str}
    )
    if list(emap.columns) != ["entrez", "ensembl_id"]:
        raise ExternalSampleError(
            f"{ensembl_map_path.name}: columns are {list(emap.columns)}, "
            f"expected ['entrez', 'ensembl_id']."
        )
    if emap["ensembl_id"].duplicated().any():
        offenders = (
            emap.loc[emap["ensembl_id"].duplicated(keep=False)]
            .head(6)
            .to_dict("records")
        )
        raise ExternalSampleError(
            f"{ensembl_map_path.name}: {int(emap['ensembl_id'].duplicated().sum())} "
            f"Ensembl ID(s) map to more than one Entrez ID (e.g. {offenders}); "
            f"the Ensembl->Entrez join would be ambiguous."
        )
    ens_to_entrez = dict(zip(emap["ensembl_id"], emap["entrez"]))

    collapsed["entrez"] = collapsed["ens"].map(ens_to_entrez)
    unresolved_external_rows = int(collapsed["entrez"].isna().sum())
    resolved = collapsed.dropna(subset=["entrez"]).copy()

    # ---- 6. restrict to the canonical space, detect Entrez collisions ---
    gene_columns = json.loads(gene_columns_path.read_text(encoding="utf-8"))
    canonical_entrez = [str(e) for e in gene_columns["entrez_ids"]]
    canonical_labels = list(gene_columns["canonical_columns"])
    if len(canonical_entrez) != len(canonical_labels):
        raise ExternalSampleError(
            f"{gene_columns_path.name}: entrez_ids ({len(canonical_entrez)}) and "
            f"canonical_columns ({len(canonical_labels)}) differ in length."
        )
    canonical_set = set(canonical_entrez)

    in_canon = resolved["entrez"].isin(canonical_set)
    resolved_outside_canonical = sorted(set(resolved.loc[~in_canon, "entrez"]))
    resolved = resolved.loc[in_canon]

    per_entrez = resolved.groupby("entrez", sort=False)
    group_sizes = per_entrez.size()
    canonical_id_collisions = sorted(group_sizes[group_sizes > 1].index)
    linear_by_entrez = per_entrez["tpm"].sum()

    # ---- 7. reindex to canonical order; missing stays explicit ----------
    position = {entrez: i for i, entrez in enumerate(canonical_entrez)}
    linear = np.full(len(canonical_entrez), np.nan, dtype=np.float64)
    for entrez, value in linear_by_entrez.items():
        linear[position[entrez]] = value

    mapped_mask = ~np.isnan(linear)
    canonical_genes_mapped = int(mapped_mask.sum())
    canonical_genes_missing = int((~mapped_mask).sum())
    canonical_genes_measured_zero = int(np.sum(linear[mapped_mask] == 0.0))
    canonical_genes_measured_nonzero = int(np.sum(linear[mapped_mask] > 0.0))

    # ---- 8. log2(TPM + 1), NaN preserved -------------------------------
    if log2_transform:
        values = np.where(mapped_mask, np.log2(linear + 1.0), np.nan)
        transformation = (
            "per-canonical-Entrez sum of linear TPM across version-stripped "
            "Ensembl IDs, then log2(TPM + 1); canonical genes not resolved in "
            "the sample left as NaN"
        )
    else:
        values = linear
        transformation = (
            "per-canonical-Entrez sum of linear TPM across version-stripped "
            "Ensembl IDs; canonical genes not resolved in the sample left as NaN"
        )
    series = pd.Series(
        values, index=canonical_labels, name=sample_id, dtype=np.float64
    )

    # ---- symbol-fallback candidates (identified, deliberately NOT used) --
    missing_entrez = {
        canonical_entrez[i] for i in range(len(canonical_entrez))
        if not mapped_mask[i]
    }
    symbol_to_entrez = {
        label.rsplit(" (", 1)[0]: entrez
        for label, entrez in zip(canonical_labels, canonical_entrez)
    }
    missing_symbols = {
        sym for sym, entrez in symbol_to_entrez.items() if entrez in missing_entrez
    }
    fallback_hits = descriptions[
        (~ens_base.isin(ens_to_entrez.keys())) & descriptions.isin(missing_symbols)
    ]
    symbol_fallback_candidates = sorted(set(fallback_hits))

    # ---- provenance ---------------------------------------------------
    finite_tpm = tpm  # already validated finite and non-negative
    provenance = {
        "sample_id": sample_id,
        "gct_file": {
            "name": gct_path.name,
            "sha256": sha256_file(gct_path),
            "bytes": gct_path.stat().st_size,
        },
        "ensembl_map_file": {
            "name": ensembl_map_path.name,
            "sha256": sha256_file(ensembl_map_path),
            "bytes": ensembl_map_path.stat().st_size,
            "rows": int(emap.shape[0]),
        },
        "gene_columns_file": {
            "name": gene_columns_path.name,
            "sha256": sha256_file(gene_columns_path),
        },
        "gct": {
            "version": gct_info["version"],
            "declared_rows": gct_info["declared_rows"],
            "declared_cols": gct_info["declared_cols"],
            "actual_rows": gct_info["actual_rows"],
            "sample_column": sample_id,
        },
        "transformation": transformation,
        "input_tpm": {
            "rows": int(finite_tpm.size),
            "min": float(finite_tpm.min()),
            "max": float(finite_tpm.max()),
            "sum": float(finite_tpm.sum()),
        },
        "reconciliation": {
            "external_rows": int(names.size),
            "unexpected_identifiers": len(unexpected_identifiers),
            "unexpected_identifier_examples": unexpected_identifiers[:10],
            "versioned_ids_stripped": int(versioned_mask.sum()),
            "bare_ensembl_ids": int(bare_mask.sum()),
            "duplicate_external_ids": len(duplicate_external_ids),
            "duplicate_external_id_examples": duplicate_external_ids[:10],
            "external_rows_after_collapse": int(collapsed.shape[0]),
            "resolved_via_ensembl_map": int(collapsed["entrez"].notna().sum()),
            "unresolved_external_rows": unresolved_external_rows,
            "resolved_outside_canonical": len(resolved_outside_canonical),
            "resolved_outside_canonical_examples": resolved_outside_canonical[:10],
            "canonical_id_collisions": len(canonical_id_collisions),
            "canonical_id_collision_examples": list(canonical_id_collisions[:10]),
            "canonical_genes": len(canonical_entrez),
            "canonical_genes_mapped": canonical_genes_mapped,
            "canonical_genes_missing": canonical_genes_missing,
            "canonical_genes_measured_zero": canonical_genes_measured_zero,
            "canonical_genes_measured_nonzero": canonical_genes_measured_nonzero,
            "symbol_fallback": "not attempted (see module docstring)",
            "symbol_fallback_candidates": symbol_fallback_candidates,
        },
    }
    return series, provenance


# --------------------------------------------------------------------------
# Offline self-test (no network, no committed-artifact mutation)
# --------------------------------------------------------------------------

def _write_gct(
    path: Path,
    rows: list[tuple[str, str, str]],
    *,
    version: str = "#1.2",
    declared: tuple[int, int] | None = None,
    header: tuple[str, ...] = ("Name", "Description", "S1"),
) -> None:
    """Write a tiny uncompressed GCT for the self-test."""
    n_rows = len(rows) if declared is None else declared[0]
    n_cols = len(header) - 2 if declared is None else declared[1]
    body = "\n".join("\t".join(r) for r in rows)
    path.write_text(
        f"{version}\n{n_rows}\t{n_cols}\n" + "\t".join(header) + "\n" + body + "\n",
        encoding="utf-8",
    )


def _write_canonical(path: Path, pairs: list[tuple[str, str]]) -> None:
    """Write a minimal gene_columns.json (symbol, entrez pairs, in order)."""
    symbols = [s for s, _ in pairs]
    entrez = [e for _, e in pairs]
    path.write_text(
        json.dumps(
            {
                "description": "self-test canonical space",
                "n_genes": len(pairs),
                "entrez_ids": entrez,
                "symbols": symbols,
                "canonical_columns": [f"{s} ({e})" for s, e in pairs],
            }
        ),
        encoding="utf-8",
    )


def _write_ensembl_map(path: Path, pairs: list[tuple[str, str]]) -> None:
    """Write a minimal ensembl_map.csv (entrez, ensembl_id)."""
    lines = ["entrez,ensembl_id"] + [f"{e},{g}" for e, g in pairs]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _self_test() -> int:
    import tempfile

    print("Running sample_profile self-test...")
    tmp = Path(tempfile.mkdtemp(prefix="sample_profile_selftest_"))

    # A canonical space of five genes; GENE5 is deliberately never supplied.
    canon = tmp / "gene_columns.json"
    _write_canonical(
        canon,
        [("GENE1", "10"), ("GENE2", "20"), ("GENE3", "30"),
         ("GENE4", "40"), ("GENE5", "50")],
    )
    emap = tmp / "ensembl_map.csv"
    _write_ensembl_map(
        emap,
        [("10", "ENSG0000000001"), ("20", "ENSG0000000002"),
         ("30", "ENSG0000000003"), ("40", "ENSG0000000004")],
        # entrez 50 (GENE5) intentionally absent from the map
    )

    # ---- 1. valid load ------------------------------------------------
    good = tmp / "good.gct"
    _write_gct(
        good,
        [
            ("ENSG0000000001.3", "GENE1", "3.0"),   # -> log2(4) = 2.0
            ("ENSG0000000002.1", "GENE2", "0.0"),   # measured zero -> 0.0
            ("ENSG0000000003", "GENE3", "7.0"),     # already bare -> log2(8) = 3.0
            ("ENSG9999999999.2", "SOMETHING", "5.0"),  # unresolved, not canonical
        ],
        header=("Name", "Description", "S1"),
    )
    series, prov = load_external_sample(
        sample_id="S1", gct_path=good, ensembl_map_path=emap,
        gene_columns_path=canon,
    )
    assert list(series.index) == [
        "GENE1 (10)", "GENE2 (20)", "GENE3 (30)", "GENE4 (40)", "GENE5 (50)"
    ], list(series.index)
    assert np.isclose(series["GENE1 (10)"], 2.0)
    assert np.isclose(series["GENE2 (20)"], 0.0)
    assert np.isclose(series["GENE3 (30)"], 3.0)
    assert np.isnan(series["GENE4 (40)"]) and np.isnan(series["GENE5 (50)"])
    r = prov["reconciliation"]
    assert r["canonical_genes"] == 5
    assert r["canonical_genes_mapped"] == 3
    assert r["canonical_genes_missing"] == 2
    assert r["canonical_genes_measured_zero"] == 1
    assert r["canonical_genes_measured_nonzero"] == 2
    assert r["versioned_ids_stripped"] == 3   # .3 .1 .2  (bare ENSG3 not counted)
    assert r["bare_ensembl_ids"] == 1
    assert r["unresolved_external_rows"] == 1
    print("  [ok] valid load: values, canonical order, measured-zero vs missing")

    # ---- 2. GCT version validation ---------------------------------
    bad_ver = tmp / "bad_version.gct"
    _write_gct(bad_ver, [("ENSG0000000001.1", "GENE1", "1.0")], version="#1.3")
    try:
        parse_gct(bad_ver)
        raise AssertionError("bad version tag should have raised")
    except GCTFormatError as exc:
        assert "#1.2" in str(exc)
    print("  [ok] wrong GCT version tag raises")

    # ---- 3. GCT dimension validation -----------------------------
    bad_rows = tmp / "bad_rows.gct"
    _write_gct(
        bad_rows,
        [("ENSG0000000001.1", "GENE1", "1.0")],
        declared=(9, 1),   # claims 9 rows, has 1
    )
    try:
        parse_gct(bad_rows)
        raise AssertionError("row-count mismatch should have raised")
    except GCTFormatError as exc:
        assert "9" in str(exc)

    bad_cols = tmp / "bad_cols.gct"
    _write_gct(
        bad_cols,
        [("ENSG0000000001.1", "GENE1", "1.0", "extra")],
        declared=(1, 1),   # header/data width 4, dimension line implies 3
        header=("Name", "Description", "S1", "S2"),
    )
    try:
        parse_gct(bad_cols)
        raise AssertionError("header width mismatch should have raised")
    except GCTFormatError:
        pass

    ragged = tmp / "ragged.gct"
    ragged.write_text(
        "#1.2\n2\t1\nName\tDescription\tS1\n"
        "ENSG0000000001.1\tGENE1\t1.0\n"
        "ENSG0000000002.1\tGENE2\n",   # short row
        encoding="utf-8",
    )
    try:
        parse_gct(ragged)
        raise AssertionError("ragged data row should have raised")
    except GCTFormatError as exc:
        assert "data line 5" in str(exc)
    print("  [ok] declared vs actual dimensions and ragged rows raise")

    # ---- 4. negative TPM ----------------------------------------
    neg = tmp / "neg.gct"
    _write_gct(neg, [("ENSG0000000001.1", "GENE1", "-0.5")])
    try:
        load_external_sample(sample_id="S1", gct_path=neg,
                             ensembl_map_path=emap, gene_columns_path=canon)
        raise AssertionError("negative TPM should have raised")
    except ExternalSampleError as exc:
        assert "negative" in str(exc)
    print("  [ok] negative TPM raises")

    # ---- 5. non-finite TPM ------------------------------------
    for token in ("nan", "inf", "NA"):
        nf = tmp / f"nf_{token}.gct"
        _write_gct(nf, [("ENSG0000000001.1", "GENE1", token)])
        try:
            load_external_sample(sample_id="S1", gct_path=nf,
                                 ensembl_map_path=emap, gene_columns_path=canon)
            raise AssertionError(f"{token!r} TPM should have raised")
        except ExternalSampleError as exc:
            assert "non-finite" in str(exc)
    print("  [ok] NaN / inf / non-numeric TPM raise")

    # ---- 6. identifier version stripping ----------------------
    assert strip_ensembl_version("ENSG00000123456.7") == "ENSG00000123456"
    assert strip_ensembl_version("ENSG00000123456.17") == "ENSG00000123456"
    assert strip_ensembl_version("ENSG00000123456") == "ENSG00000123456"
    assert strip_ensembl_version("NOX5") == "NOX5"
    assert strip_ensembl_version("ENST00000123456.1") == "ENST00000123456.1"
    print("  [ok] Ensembl version stripping (and pass-through of non-matches)")

    # ---- 7a. duplicate stripped Ensembl IDs -> summed, reported ---
    dup_ext = tmp / "dup_ext.gct"
    _write_gct(
        dup_ext,
        [
            ("ENSG0000000001.1", "GENE1", "1.0"),
            ("ENSG0000000001.5", "GENE1", "3.0"),   # same gene, different version
        ],
    )
    s_dup, p_dup = load_external_sample(
        sample_id="S1", gct_path=dup_ext, ensembl_map_path=emap,
        gene_columns_path=canon,
    )
    assert np.isclose(s_dup["GENE1 (10)"], np.log2(1.0 + 3.0 + 1.0))  # sum then log2
    assert p_dup["reconciliation"]["duplicate_external_ids"] == 1
    assert "ENSG0000000001" in p_dup["reconciliation"]["duplicate_external_id_examples"]

    # ---- 7b. two Ensembl IDs -> one canonical Entrez -> summed, reported ---
    emap_collide = tmp / "emap_collide.csv"
    _write_ensembl_map(
        emap_collide,
        [("10", "ENSG0000000001"), ("10", "ENSG0000000099"),  # both -> entrez 10
         ("20", "ENSG0000000002"), ("30", "ENSG0000000003"),
         ("40", "ENSG0000000004")],
    )
    collide = tmp / "collide.gct"
    _write_gct(
        collide,
        [
            ("ENSG0000000001.1", "GENE1", "1.0"),
            ("ENSG0000000099.1", "GENE1B", "2.0"),
        ],
    )
    s_col, p_col = load_external_sample(
        sample_id="S1", gct_path=collide, ensembl_map_path=emap_collide,
        gene_columns_path=canon,
    )
    assert np.isclose(s_col["GENE1 (10)"], np.log2(1.0 + 2.0 + 1.0))
    assert p_col["reconciliation"]["canonical_id_collisions"] == 1
    assert "10" in p_col["reconciliation"]["canonical_id_collision_examples"]

    # ---- 7c. non-unique ensembl_id column -> hard error ---
    emap_bad = tmp / "emap_bad.csv"
    _write_ensembl_map(
        emap_bad,
        [("10", "ENSG0000000001"), ("20", "ENSG0000000001"),  # same ENSG, 2 entrez
         ("30", "ENSG0000000003")],
    )
    try:
        load_external_sample(sample_id="S1", gct_path=good,
                             ensembl_map_path=emap_bad, gene_columns_path=canon)
        raise AssertionError("non-unique ensembl_id should have raised")
    except ExternalSampleError as exc:
        assert "ambiguous" in str(exc)
    print("  [ok] duplicate / colliding identifiers handled explicitly")

    # ---- 8. deterministic output order, independent of GCT row order ---
    import random
    shuffled_rows = [
        ("ENSG0000000003", "GENE3", "7.0"),
        ("ENSG0000000001.3", "GENE1", "3.0"),
        ("ENSG9999999999.2", "SOMETHING", "5.0"),
        ("ENSG0000000002.1", "GENE2", "0.0"),
    ]
    random.Random(0).shuffle(shuffled_rows)
    shuf = tmp / "shuffled.gct"
    _write_gct(shuf, shuffled_rows)
    s_shuf, _ = load_external_sample(
        sample_id="S1", gct_path=shuf, ensembl_map_path=emap,
        gene_columns_path=canon,
    )
    assert list(s_shuf.index) == list(series.index)
    assert s_shuf.equals(series) or np.allclose(
        s_shuf.to_numpy(), series.to_numpy(), equal_nan=True
    )
    print("  [ok] canonical output order is deterministic regardless of row order")

    # ---- 9. measured-zero vs missing already asserted in (1); re-confirm mask
    assert bool(np.isnan(series.to_numpy())[3]) and bool(np.isnan(series.to_numpy())[4])
    assert series.to_numpy()[1] == 0.0
    print("  [ok] measured TPM==0 is 0.0, genuinely absent genes are NaN")

    # ---- 10. provenance completeness + internal consistency ---
    required_top = {
        "sample_id", "gct_file", "ensembl_map_file", "gene_columns_file",
        "gct", "transformation", "input_tpm", "reconciliation",
    }
    assert required_top <= set(prov), required_top - set(prov)
    for block in ("gct_file", "ensembl_map_file", "gene_columns_file"):
        assert len(prov[block]["sha256"]) == 64
    rr = prov["reconciliation"]
    assert rr["canonical_genes_mapped"] + rr["canonical_genes_missing"] == rr["canonical_genes"]
    assert (rr["canonical_genes_measured_zero"] + rr["canonical_genes_measured_nonzero"]
            == rr["canonical_genes_mapped"])
    json.dumps(prov)  # must be serialisable
    print("  [ok] provenance record is complete, consistent and serialisable")

    # ---- 11. the real committed inputs load, and nothing frozen mutates ---
    frozen = [
        config.PROCESSED_DIR / "baseline_results.json",
        config.PROCESSED_DIR / "head_results.json",
        config.PROCESSED_DIR / "analysis_results.json",
        config.PROCESSED_DIR / "expression.npz",
        config.PROCESSED_DIR / "expression.labels.json",
        config.PROCESSED_DIR / "crispr_effect.npz",
        config.PROCESSED_DIR / "splits.json",
        config.PROCESSED_DIR / "gene_columns.json",
        config.PROCESSED_DIR / "selective_genes.json",
        config.ENSEMBL_MAP_FILE,
        config.DEMO_TUMOR_GCT_FILE,
    ]
    before = {p: sha256_file(p) for p in frozen if p.is_file()}
    real_series, real_prov = load_external_sample()
    after = {p: sha256_file(p) for p in frozen if p.is_file()}
    assert before == after, "a committed input changed during load()"
    assert real_series.name == config.DEMO_TUMOR_SAMPLE_ID
    assert len(real_series) == real_prov["reconciliation"]["canonical_genes"]
    print(
        "  [ok] real BG003082 load: "
        f"{real_prov['reconciliation']['canonical_genes_mapped']} mapped / "
        f"{real_prov['reconciliation']['canonical_genes_missing']} missing; "
        "no committed input mutated"
    )

    # tidy up
    for child in sorted(tmp.iterdir()):
        child.unlink()
    tmp.rmdir()
    print("\nSelf-test passed.")
    return 0


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def _print_summary(series: pd.Series, provenance: dict) -> None:
    rr = provenance["reconciliation"]
    print("=" * 74)
    print(f"EXTERNAL SAMPLE PROFILE  --  {provenance['sample_id']}")
    print("=" * 74)
    print(f"  GCT file        : {provenance['gct_file']['name']}")
    print(f"    sha256        : {provenance['gct_file']['sha256']}")
    print(f"    GCT version   : {provenance['gct']['version']}")
    print(f"    dimensions    : {provenance['gct']['declared_rows']} rows x "
          f"{provenance['gct']['declared_cols']} sample col(s) "
          f"(actual rows {provenance['gct']['actual_rows']})")
    print(f"  ensembl_map     : {provenance['ensembl_map_file']['name']}  "
          f"({provenance['ensembl_map_file']['rows']} rows)")
    print(f"    sha256        : {provenance['ensembl_map_file']['sha256']}")
    print(f"  transformation  : {provenance['transformation']}")
    print("-" * 74)
    print(f"  external rows                     : {rr['external_rows']}")
    print(f"  unexpected identifiers           : {rr['unexpected_identifiers']}")
    print(f"  versioned IDs stripped           : {rr['versioned_ids_stripped']}")
    print(f"  duplicate external IDs (summed)  : {rr['duplicate_external_ids']}")
    print(f"  resolved via ensembl_map         : {rr['resolved_via_ensembl_map']}")
    print(f"  unresolved external rows         : {rr['unresolved_external_rows']}")
    print(f"  resolved outside canonical space : {rr['resolved_outside_canonical']}")
    print(f"  canonical-Entrez collisions      : {rr['canonical_id_collisions']}")
    print("-" * 74)
    print(f"  canonical genes                  : {rr['canonical_genes']}")
    print(f"    mapped                         : {rr['canonical_genes_mapped']}")
    print(f"      measured TPM == 0            : {rr['canonical_genes_measured_zero']}")
    print(f"      measured TPM  > 0            : {rr['canonical_genes_measured_nonzero']}")
    print(f"    missing / unmappable (NaN)     : {rr['canonical_genes_missing']}")
    print(f"  symbol-fallback candidates (not applied): "
          f"{rr['symbol_fallback_candidates']}")
    print("=" * 74)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Load an external RNA-seq GCT sample into the canonical "
                    "gene space."
    )
    parser.add_argument("--self-test", action="store_true",
                        help="Run the offline self-test and exit.")
    parser.add_argument("--json", action="store_true",
                        help="Print the full provenance record as JSON.")
    parser.add_argument("--linear", action="store_true",
                        help="Return linear TPM instead of log2(TPM + 1).")
    parser.add_argument("--sample-id", default=config.DEMO_TUMOR_SAMPLE_ID)
    args = parser.parse_args()

    if args.self_test:
        return _self_test()

    series, provenance = load_external_sample(
        sample_id=args.sample_id, log2_transform=not args.linear
    )
    _print_summary(series, provenance)
    if args.json:
        print(json.dumps(provenance, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
