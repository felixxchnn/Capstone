"""
config.py
=========
Central configuration for the DepMap perturbation-response engine.

Everything that might change between DepMap releases, between machines, or
between experiments lives here. No other module should hard-code a filename,
a path, or a threshold.

Design notes
------------
* Filenames drift between DepMap quarterly releases (the expression matrix was
  `OmicsExpressionProteinCodingGenesTPMLogp1.csv` in some releases and
  `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` in others). Rather than
  hard-code one, `resolve_file()` searches a list of known aliases and, failing
  that, falls back to a glob pattern. If it still cannot find the file it
  raises a *loud, useful* error that lists what is actually present.
* Set `DEBUG_NROWS` to a small integer to make every loader read only the first
  N rows. Invaluable for fast iteration; set back to None for real runs.
"""

from __future__ import annotations

import os
import glob
from pathlib import Path

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------

# Directory containing this file.
PROJECT_ROOT = Path(__file__).resolve().parent

# Where the raw DepMap CSVs live. Override with the DEPMAP_DATA_DIR
# environment variable, e.g. on Kaggle:
#     os.environ["DEPMAP_DATA_DIR"] = "/kaggle/input/depmap"
_env_raw = os.environ.get("DEPMAP_DATA_DIR")
if _env_raw:
    RAW_DIR = Path(_env_raw)
elif (PROJECT_ROOT / "data" / "raw").is_dir():
    RAW_DIR = PROJECT_ROOT / "data" / "raw"
else:
    # Fall back to "the CSVs sit next to the scripts".
    RAW_DIR = PROJECT_ROOT

# Where processed outputs are written.
PROCESSED_DIR = Path(os.environ.get(
    "DEPMAP_PROCESSED_DIR",
    PROJECT_ROOT / "data" / "processed",
))

# --------------------------------------------------------------------------
# Raw file aliases
# --------------------------------------------------------------------------
# Each entry maps a logical name -> (list of known filenames, glob fallback).
# Order matters: the first match wins.

FILE_ALIASES: dict[str, tuple[list[str], str]] = {
    "expression": (
        [
            "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
            "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
        ],
        "OmicsExpression*TPM*ProteinCoding*.csv",
    ),
    "crispr": (
        [
            "CRISPRGeneEffect.csv",
            "CRISPRGeneEffectUncorrected.csv",
        ],
        "CRISPRGeneEffect*.csv",
    ),
    # Optional: DepMap's gene-level expected-counts file. Referenced by
    # prepare_geneformer_input.py, which needs raw counts, not log-TPM.
    # Deliberately empty: capstone/RESUME-HERE-2026-08-23.md records that the
    # only candidate filenames anyone has ever written down for this file
    # ("OmicsExpressionGenesExpectedCountProfile.csv",
    # "OmicsExpressionGenesExpectedCount.csv") were "guessed at the time and
    # have never been checked against a real DepMap release." Registering
    # unverified names here risks a silent wrong-file match; an empty entry
    # cannot match anything, so resolve_file(required=False) always returns
    # None and prepare_geneformer_input.py falls through to the documented,
    # already-used reconstructed-from-log-TPM path. Fill this in only once a
    # real filename has been confirmed against an actual downloaded file.
    "expression_counts": ([], ""),
    "model": (
        [
            "Model.csv",
            "Model (1).csv",
            "Model__1_.csv",
            "sample_info.csv",
        ],
        "Model*.csv",
    ),
    # --- PRISM drug response (optional; pipeline runs without it) ---
    "prism_matrix": (
        [
            "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv",
            "Repurposing_Public_24Q2_LFC_COLLAPSED.csv",
            "Repurposing_Public_23Q2_LFC_COLLAPSED.csv",
            "primary-screen-replicate-collapsed-logfold-change.csv",
        ],
        "*[Pp]rimary*[Dd]ata*[Mm]atrix*.csv",
    ),
    "prism_compounds": (
        [
            "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv",
            "primary-screen-replicate-treatment-info.csv",
        ],
        "*[Cc]ompound*[Ll]ist*.csv",
    ),
    "prism_cell_lines": (
        [
            "primary-screen-cell-line-info.csv",
        ],
        "*cell*line*info*.csv",
    ),
}


class MissingDataFile(FileNotFoundError):
    """Raised when a required raw file cannot be located."""


class AmbiguousDataFile(Exception):
    """Raised when a glob fallback matches more than one file and the caller
    has asked not to guess (see `resolve_file`'s `reject_ambiguous_glob`)."""


def resolve_file(
    logical_name: str,
    required: bool = True,
    reject_ambiguous_glob: bool = False,
) -> Path | None:
    """
    Locate a raw data file by logical name.

    Parameters
    ----------
    logical_name
        One of the keys of FILE_ALIASES.
    required
        If True, raise MissingDataFile when nothing matches. If False,
        return None instead (used for the optional PRISM files).
    reject_ambiguous_glob
        If True, raise AmbiguousDataFile when the glob fallback matches more
        than one file, instead of silently taking the first (sorted) match.
        Default False preserves existing behaviour for every alias already
        relied on by build_dataset.py; opt in per call site for a newly
        added, less-trusted alias (e.g. one with unverified filenames).

    Returns
    -------
    Path or None
    """
    if logical_name not in FILE_ALIASES:
        raise KeyError(
            f"Unknown logical file name {logical_name!r}. "
            f"Known: {sorted(FILE_ALIASES)}"
        )

    known_names, pattern = FILE_ALIASES[logical_name]

    for name in known_names:
        candidate = RAW_DIR / name
        if candidate.is_file():
            return candidate

    # An empty pattern means "no glob fallback configured" -- glob.glob on a
    # bare directory path matches the directory itself, so this must be
    # skipped rather than passed through.
    matches = sorted(glob.glob(str(RAW_DIR / pattern))) if pattern else []

    if reject_ambiguous_glob and len(matches) > 1:
        raise AmbiguousDataFile(
            f"\nMultiple files match the {logical_name!r} glob pattern "
            f"{pattern!r} and none of the known exact names "
            f"{known_names} matched, so the choice would be a guess:\n"
            + "\n".join(f"    - {Path(m).name}" for m in matches)
            + f"\n\nRename the correct file to one of the known names, or "
              f"narrow FILE_ALIASES[{logical_name!r}] in config.py.\n"
        )

    if matches:
        return Path(matches[0])

    if not required:
        return None

    present = sorted(p.name for p in RAW_DIR.iterdir() if p.is_file()) \
        if RAW_DIR.is_dir() else []
    raise MissingDataFile(
        f"\nCould not find the {logical_name!r} file.\n"
        f"  Looked in : {RAW_DIR}\n"
        f"  Tried     : {known_names}\n"
        f"  Glob      : {pattern}\n"
        f"  Files present in that directory:\n"
        + ("\n".join(f"    - {n}" for n in present) if present
           else "    (directory is empty or does not exist)")
        + "\n\nFix: either rename your file to one of the names above, add its "
          "real name to FILE_ALIASES in config.py, or point DEPMAP_DATA_DIR at "
          "the right folder.\n"
    )


# --------------------------------------------------------------------------
# Column names
# --------------------------------------------------------------------------

# Join key used across every DepMap table.
MODEL_ID = "ModelID"

# Metadata columns that prefix the gene block in the expression matrix.
# These are NOT genes and must be stripped before modelling.
EXPRESSION_META_COLS = [
    "SequencingID",
    "ModelConditionID",
    "ModelID",
    "IsDefaultEntryForMC",
    "IsDefaultEntryForModel",
]

# The flag that marks the canonical profile for each model. The expression
# matrix contains multiple sequencing runs per cell line; failing to filter on
# this produces duplicate ModelIDs, which silently inflate the dataset and leak
# near-identical rows across train/test.
EXPRESSION_DEFAULT_FLAG = "IsDefaultEntryForModel"

# Metadata columns kept from Model.csv.
MODEL_KEEP_COLS = [
    "PatientID",
    "CellLineName",
    "StrippedCellLineName",
    "DepmapModelType",
    "OncotreeLineage",
    "OncotreePrimaryDisease",
    "OncotreeSubtype",
    "OncotreeCode",
    "Age",
    "Sex",
    "PrimaryOrMetastasis",
    "SampleCollectionSite",
    "PediatricModelType",
]

# Grouping variable for leakage-safe splits. Several cell lines can derive from
# the same patient; letting them straddle a split is subtle leakage.
GROUP_COL = "PatientID"

# Stratification variable for splits.
STRATIFY_COL = "OncotreeLineage"

# --------------------------------------------------------------------------
# Selective-essentiality filter
# --------------------------------------------------------------------------
# Chronos gene effect: 0 = no effect, -1 = median of pan-essential genes.
#
# Predicting pan-essential genes is trivial (every line is dependent) and
# predicting never-essential genes is noise. Neither tells you anything about
# context-specific biology, but both inflate aggregate metrics. We therefore
# restrict targets to genes that are essential in *some* lines and not others.

DEPENDENCY_THRESHOLD = -0.5      # a line is "dependent" below this score
MIN_DEPENDENT_LINES = 5          # gene must be a dependency in >= this many lines
MAX_DEPENDENT_FRACTION = 0.90    # ...but not in more than this fraction (excludes pan-essential)
MIN_GENE_EFFECT_STD = 0.15       # gene effect must vary across lines
MAX_GENE_NAN_FRACTION = 0.05     # drop genes missing in >5% of lines

# --------------------------------------------------------------------------
# Expression feature filter
# --------------------------------------------------------------------------

MIN_EXPRESSION_STD = 0.0         # 0.0 keeps all genes; raise to drop flat ones
DROP_ALL_ZERO_EXPRESSION = True  # drop genes that are zero in every line

# --------------------------------------------------------------------------
# Frozen legacy build asymmetry -- do not "fix" without going through
# CLAUDE.md section 13 first (diagnosed 2026-08-26)
# --------------------------------------------------------------------------
# build_dataset.py intersects expression and crispr by Entrez ID, then drops
# all-zero-expression columns from `expression` only (DROP_ALL_ZERO_EXPRESSION
# above); the surviving set becomes the canonical space written to
# gene_columns.json / gene_id_map.csv. The equivalent drop was never applied
# back to the saved `crispr_effect` matrix, so it carries 3 columns outside
# the canonical space. join_report.json's "expression_features" block records
# only the *count* of genes dropped (dropped_all_zero: 3) -- it does not
# record which genes. Their identities were established by direct diff
# against gene_columns.json during manual investigation (2026-08-26) and are
# pinned by name below, so that a future, *different* discrepancy still fails
# checks.py loudly instead of being silently absorbed into a general
# allowance.
#
# Confirmed (against data/processed/selective_genes.json, 2026-08-26): all
# three are absent from the 4,297 published CRISPR targets, and every
# modelling path (baseline.py, train_head.py, analysis.py) always subsets
# crispr via selective_genes, never via the full column set -- this asymmetry
# has never touched any committed result. Re-running build_dataset.py would
# rewrite crispr_effect.npz, gene_id_map.csv, gene_columns.json,
# selective_genes.json and join_report.json together; that is a material
# change to input data under CLAUDE.md section 13 and needs sign-off first.
CRISPR_LEGACY_EXTRA_GENES = frozenset({
    "KRTAP23-1 (337963)",
    "DEFB113 (245927)",
    "KRTAP20-3 (337985)",
})

# --------------------------------------------------------------------------
# Splits
# --------------------------------------------------------------------------

RANDOM_SEED = 20260722
TEST_FRACTION = 0.15
VAL_FRACTION = 0.15   # fraction of the *whole* dataset, taken after test

# If True, entire lineages are held out (a much harder generalisation test).
# Default False: this would remove all Bone/osteosarcoma lines from training,
# which is exactly the biology the engine needs before it is applied to a
# real osteosarcoma tumour.
GROUP_SPLIT_BY_LINEAGE = False

# --------------------------------------------------------------------------
# Baseline model
# --------------------------------------------------------------------------

PCA_COMPONENTS = 200
# Half-decade spacing from 1 to 1e6. The original decade grid from 1 to 1e4
# truncated the search: the inner-CV sweep was monotonically increasing and
# still accelerating at the ceiling (0.1658 -> 0.1659 -> 0.1668 -> 0.1743 ->
# 0.2057), so the selected alpha was the largest offered rather than the best
# available. The ceiling is raised to 1e6 because that is where a ridge on
# this data is shrunk to effectively zero degrees of freedom; nothing above it
# can be optimal.
RIDGE_ALPHAS = [
    1.0, 3.16, 10.0, 31.6, 100.0, 316.0,
    1000.0, 3162.0, 10000.0, 31623.0,
    100000.0, 316228.0, 1000000.0,
]

# --------------------------------------------------------------------------
# Runtime
# --------------------------------------------------------------------------

# Set to an int (e.g. 200) to read only that many rows from each big CSV.
# Use for fast iteration; set to None for a real run.
DEBUG_NROWS: int | None = None

# Store float matrices as float32. Halves memory and disk; log-TPM and Chronos
# scores carry nowhere near float64 precision.
FLOAT_DTYPE = "float32"

# Lineage label used for osteosarcoma in Model.csv (for the OS-specific report).
OSTEOSARCOMA_LINEAGE = "Bone"
OSTEOSARCOMA_DISEASE_KEYWORD = "osteosarcoma"

# --------------------------------------------------------------------------
# Phase 2 -- precision-oncology demo (capstone/scope-decisions.md, 2026-08-25)
# --------------------------------------------------------------------------
# Two samples, two roles -- see CLAUDE.md section 9.7. Neither constant changes
# anything about the Phase 1 pipeline; both are read only by the new demo
# modules (sample_profile.py, evidence.py, case_study.py, report.py).

# Root for external (non-DepMap) data the demo depends on. Small enough files
# under here (the GCT, the license-filtered DGIdb snapshot) are committed
# directly, unlike the gitignored raw DepMap CSVs under RAW_DIR.
EXTERNAL_DIR = PROJECT_ROOT / "data" / "external"
SID_OSTEOSARC_DIR = EXTERNAL_DIR / "sid_osteosarc"
DGIDB_DIR = EXTERNAL_DIR / "dgidb"

# Internal verification anchor: a DepMap osteosarcoma line with real CRISPR
# ground truth, chosen for osteosarcoma identity and complete artifacts, not
# for a favourable predicted score (capstone/scope-decisions.md).
DEMO_VERIFICATION_MODEL_ID = "ACH-000364"  # U-2 OS, val split

# Primary demo sample: Sid Sijbrandij's real primary-tumour RNA-seq
# (osteosarc.com, CC0 1.0). Bulk tumour tissue, not a cultured cell line --
# a real domain shift from every DepMap line the model ever saw. No CRISPR
# ground truth exists for it.
DEMO_TUMOR_SAMPLE_ID = "BG003082"
DEMO_TUMOR_GCT_FILE = SID_OSTEOSARC_DIR / "BG003082.gene_tpm.gct.gz"

# Machine-readable status fields attached to every prediction the demo shows.
# Never inferred at report-render time -- fixed here so the two samples cannot
# accidentally be tagged the same way.
PREDICTION_STATUS_HELD_OUT = "held_out_prediction"
PREDICTION_STATUS_EXPLORATORY = "exploratory_external_prediction"
OUTCOME_STATUS_MEASURED = "measured_crispr"
OUTCOME_STATUS_UNAVAILABLE = "unavailable"

# Entrez <-> Ensembl cross-reference. prepare_geneformer_input.py already
# expects to cache one here (it just never has -- see CLAUDE.md section 7);
# sample_profile.py's Ensembl-first gene-ID reconciliation reads the same file.
ENSEMBL_MAP_FILE = PROCESSED_DIR / "ensembl_map.csv"

# Reconstructed fitted state for the two frozen Phase 1 linear models
# (baseline ridge_pca, Geneformer ridge_head). These are NOT the original
# Phase 1 fitted objects -- those were never serialised. They are
# "reconstructed fitted state at the frozen Phase 1 alpha from the unchanged
# frozen training data": scaler/PCA/ridge re-fit on exactly the committed
# train split, at the alpha READ FROM baseline_results.json / head_results.json
# (no hyperparameter selection), then serialised as plain .npy arrays plus a
# JSON manifest. Built by reconstruct_fitted.py; loaded (no fit()) by
# fitted_artifacts.py. See capstone/data-integrity-hashes.md.
RECONSTRUCTED_FITTED_DIR = PROCESSED_DIR / "reconstructed_fitted"

# Ranking / evidence display thresholds.
TOP_N_DEPENDENCIES = 25     # genes shown per sample in the ranked-dependency table
TOP_K_EVIDENCE_PER_GENE = 5  # interaction records shown per gene, per direction tier

# Phase 2 case-study JSON (written by case_study.py) and the self-contained
# offline HTML report rendered from it (written by report.py). The HTML is a
# committed deliverable; `.gitattributes` marks it `-text` so its SHA-256
# survives a clone byte-exact on any platform.
CASE_STUDY_JSON_FILE = PROCESSED_DIR / "case_study.json"
REPORT_HTML_FILE = PROJECT_ROOT / "phase2_report.html"

# Schema identifiers stamped into the two Phase 2 deliverables. Centralised here
# so case_study.py, report.py and checks.py agree on one source of truth.
CASE_STUDY_SCHEMA_VERSION = "case-study/1"
REPORT_SCHEMA_VERSION = "phase2-report/1"

# Approved SHA-256 of the two committed Phase 2 deliverables. Both regenerate
# byte-identically on any platform (case_study.py writes explicit CRLF bytes and
# stamps the frozen reconstruction environment, not the runtime one; report.py
# writes explicit LF). checks.py sections 9 and 12 gate the committed bytes
# against these; a change here is a deliberate, reviewed re-approval, never a
# side effect of an edit.
CASE_STUDY_JSON_SHA256 = (
    "a962c01a5b65a6ef579ea57dced67048bf9016ba0f66aab2355cf1f054796e8c"
)
REPORT_HTML_SHA256 = (
    "f4a093b04bdda3e573056e2d1e2dbdde86e75cee84adf723b7b94a94dc705163"
)

# --------------------------------------------------------------------------
# Phase 2 -- DGIdb drug-gene interaction *evidence retrieval* layer
# (evidence.py). This is cited-interaction retrieval, NOT drug-response
# prediction and NOT a treatment recommendation: no model is attached, no
# efficacy / clinical relevance / approval / indication / osteosarcoma
# relevance is ever inferred. Distinct from `prism_response` above (which is
# drug-response prediction infrastructure) and from any CPIC/PharmGKB safety
# layer (out of scope -- capstone/scope-decisions.md, CLAUDE.md section 10).
# --------------------------------------------------------------------------

# Authoritative upstream release. Pinned to an immutable tag + per-asset
# byte size and SHA-256 so a refresh cannot silently pull a moving "latest".
DGIDB_RELEASE_TAG = "2026-06b"
DGIDB_RELEASE_PAGE = (
    "https://github.com/dgidb/dgidb-data/releases/tag/2026-06b"
)
_DGIDB_ASSET_BASE = (
    "https://github.com/dgidb/dgidb-data/releases/download/2026-06b"
)
DGIDB_ASSETS: dict[str, dict[str, object]] = {
    "interactions.tsv": {
        "url": f"{_DGIDB_ASSET_BASE}/interactions.tsv",
        "bytes": 16_163_629,
        "sha256":
            "5a798ffdaddd775b8409fa22509d77d55c663ca84625808617a78b1fc75dcc9e",
    },
    "genes.tsv": {
        "url": f"{_DGIDB_ASSET_BASE}/genes.tsv",
        "bytes": 4_248_825,
        "sha256":
            "e3cfcff1001fdad9e0079998e97934c6bc666ea54598b35c570d2c92829689e7",
    },
    "drugs.tsv": {
        "url": f"{_DGIDB_ASSET_BASE}/drugs.tsv",
        "bytes": 7_253_304,
        "sha256":
            "f7d465c153d8235a61b9ad43322c98b1f15ace9577e6dc0250e27f19ac463ff5",
    },
}

# Tracked, committed offline artifacts (filtered snapshot + provenance
# manifest). The unfiltered upstream TSVs are NEVER committed here.
DGIDB_SNAPSHOT_FILE = DGIDB_DIR / f"dgidb_{DGIDB_RELEASE_TAG}.interactions.filtered.tsv"
DGIDB_MANIFEST_FILE = DGIDB_DIR / f"dgidb_{DGIDB_RELEASE_TAG}.manifest.json"

# Attached verbatim to every evidence record returned for display.
DGIDB_EVIDENCE_DISCLAIMER = (
    "A recorded drug–gene interaction does not establish efficacy for this "
    "sample or for osteosarcoma."
)

# The only interaction-direction buckets. Normalised strictly from DGIdb's own
# interaction-type -> directionality vocabulary (evidence.DIRECTIONALITY_BY_TYPE),
# never from drug names, free text, approval, or indication.
DGIDB_DIRECTION_TIERS = ("inhibitory", "activating", "unknown")

# Interaction sources whose redistribution terms are *explicitly verified* as
# compatible with committing the filtered records (CC0 / CC BY / CC BY-SA /
# US-government public domain -- no NonCommercial clause, no "unclear" terms).
# DGIdb aggregates ~15 further interaction sources whose terms are NonCommercial,
# unclear, or absent; those are excluded from the committed snapshot. The full
# per-source decision, licence text and URL live in evidence.SOURCE_LICENCES,
# the manifest, and LICENSES.md. Do NOT claim the whole DGIdb dataset is
# redistributable.
DGIDB_INCLUDED_SOURCES = frozenset({
    "CIViC", "ChEMBL", "GuideToPharmacology", "DoCM", "NCI", "FDA",
})

# --- Provenance: three distinct versions, do NOT conflate them ---
# The release *packages* interaction content dated "Dec-2023"; it is NOT
# June-2026 data. `evidence.py` records all three and never labels the
# interaction records as June 2026.
DGIDB_INTERACTION_DATA_VERSION = "Dec-2023"      # TSV/SQL "# Data version:"
DGIDB_APP_VERSION_TSV_COMMENT = "v.5.0.11"       # TSV "# DGIdb version:"
DGIDB_APP_VERSION_GRAPHQL = "v.5.0.12"           # GraphQL serviceInfo at retrieval

# Fixed retrieval provenance -- a build input, NOT a wall-clock value. Every
# committed artifact must be byte-identical across rebuilds, so nothing in a
# tracked output may contain datetime.now(). Bump this only on a real re-pull of
# the pinned assets. Volatile per-run execution times go to DGIDB_RUNLOG_FILE
# (untracked).
DGIDB_RETRIEVED_UTC = "2026-08-29T00:00:00Z"
DGIDB_RUNLOG_FILE = DGIDB_DIR / "build_runlog.jsonl"   # gitignored

# HGNC ID -> Entrez ID crosswalk. Pinned immutable monthly archive; the gene
# identity join is `DGIdb hgnc:<id>` -> HGNC ID -> Entrez ID. Symbols are used
# only as a consistency check, never as the identity join.
DGIDB_HGNC_ASSET: dict[str, object] = {
    "name": "hgnc_complete_set_2026-06-02.txt",
    "url": (
        "https://storage.googleapis.com/public-download-files/hgnc/archive/"
        "archive/monthly/tsv/hgnc_complete_set_2026-06-02.txt"
    ),
    "version": "HGNC monthly release 2026-06-02",
    "retrieved": "2026-08-29",
    "bytes": 16_738_373,
    "sha256":
        "1f8826fb0b519296233dba3f987bc38efcfb3706032ed0c0685a6936f9b11e93",
    "provenance": (
        "genenames.org immutable monthly archive; aligns with DGIdb 2026-06b's "
        "HGNC source version 20260619 (2026-06-19). HGNC IDs are stable, so "
        "month-level alignment is sufficient for an ID->ID crosswalk."
    ),
    "license": "genenames.org: no restrictions on use of HGNC data",
    "license_url": "https://www.genenames.org/about/license/",
}

# Release-aligned DGIdb SQL dump. Used ONLY as a temporary build input to
# recover claim-level publication (PMID) links -- it is NOT committed. PMIDs are
# NEVER parsed from free text.
DGIDB_SQL_ASSET: dict[str, object] = {
    "name": "dgidb_2026_06b.sql.gz",
    "url": (
        "https://github.com/dgidb/dgidb-data/releases/download/2026-06b/"
        "dgidb_2026_06b.sql.gz"
    ),
    "bytes": 84_060_850,
    "sha256":
        "1d1fae3fa9e4b42a8959ef0a84b7b339a236a59403935dcb7e0f2dba20084435",
    "role": "temporary build input -- NOT committed",
}


def osteosarcoma_mask(metadata) -> "pd.Series":  # noqa: F821 (pandas imported lazily)
    """
    Boolean mask selecting osteosarcoma cell lines from a metadata frame.

    Defined in one place because three modules need it and they must agree.
    DepMap records the label in `OncotreePrimaryDisease` for some lines and in
    `OncotreeSubtype` for others, so both are searched.
    """
    import pandas as pd  # local import keeps config dependency-free at import time

    mask = pd.Series(False, index=metadata.index)
    keyword = OSTEOSARCOMA_DISEASE_KEYWORD
    for column in ("OncotreePrimaryDisease", "OncotreeSubtype"):
        if column in metadata.columns:
            values = metadata[column].astype(str).str.lower()
            mask |= values.str.contains(keyword, na=False)
    return mask


def ensure_dirs() -> None:
    """Create the processed-output directory if it does not exist."""
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)


def describe() -> str:
    """Human-readable summary of the resolved configuration."""
    lines = [
        "Configuration",
        "-------------",
        f"  PROJECT_ROOT  : {PROJECT_ROOT}",
        f"  RAW_DIR       : {RAW_DIR}",
        f"  PROCESSED_DIR : {PROCESSED_DIR}",
        f"  DEBUG_NROWS: int | None = {DEBUG_NROWS}",
        f"  RANDOM_SEED   : {RANDOM_SEED}",
        "",
        "Resolved raw files:",
    ]
    for name in FILE_ALIASES:
        required = not name.startswith("prism")
        try:
            path = resolve_file(name, required=required)
        except MissingDataFile:
            path = None
        status = str(path) if path else "(not found)"
        lines.append(f"  {name:<18}: {status}")
    return "\n".join(lines)


if __name__ == "__main__":
    print(describe())
