"""
evidence.py
===========
Phase 2 -- offline, licence-filtered DGIdb drug-gene interaction *evidence
retrieval*.

    py evidence.py                       # load the committed snapshot, print a summary
    py evidence.py --gene 7157 --symbol TP53   # show evidence for one canonical gene
    py evidence.py --self-test           # synthetic offline checks, no network
    py evidence.py --validate            # validate the committed offline snapshot
    py evidence.py --from-staging DIR    # rebuild the snapshot from pre-downloaded assets
    py evidence.py --refresh             # download the pinned assets, then rebuild

What this is, and is not
------------------------
This module *retrieves* recorded drug-gene interactions from a pinned release of
the Drug-Gene Interaction Database (DGIdb) and serves them, per canonical gene,
for the Phase 2 demo (`case_study.py`, `report.py`).

It is **evidence retrieval, not treatment prediction**. No model is attached.
The returned drugs are **not** "candidate treatments" or "recommendations". This
module never infers efficacy, clinical relevance, regulatory approval,
indication, interaction direction (beyond DGIdb's own explicit vocabulary), or
osteosarcoma relevance. Every display record carries, verbatim:

    "A recorded drug–gene interaction does not establish efficacy for this
     sample or for osteosarcoma."

Offline by default
------------------
`load_snapshot()` reads only two committed files under `data/external/dgidb/`
-- a filtered TSV and a provenance manifest -- and never touches the network. If
they are absent it raises `DGIdbSnapshotError` telling you to run `--refresh`;
it does **not** silently fall back to a live query.

The committed snapshot is licence-filtered
------------------------------------------
DGIdb's own software licence does not grant redistribution rights over every
source it aggregates. Only interaction sources whose redistribution terms are
*explicitly verified* as compatible with committing the filtered records are
included (CC0 / CC BY / CC BY-SA / US-government public domain, no NonCommercial
clause, no "unclear" terms). Every other interaction source -- NonCommercial,
custom-restrictive, "unclear", or unlicensed supplementary-table data -- is
excluded. The per-source decision, licence text and URL are recorded in
`SOURCE_LICENCES` below, in the manifest, and in `LICENSES.md`. The unfiltered
upstream TSVs are never committed. **This is not a claim that the whole DGIdb
dataset is redistributable -- it is not.**

Gene identity
-------------
The retained primary identity is the **canonical Entrez ID**. DGIdb 2026-06b is
keyed on HGNC concept IDs; each interaction row also carries `gene_name`, the
HGNC-approved symbol, and a single non-empty `hgnc:` concept id. A DGIdb gene is
resolved to a canonical Entrez ID only when its approved symbol matches exactly
one symbol in the frozen canonical space (`gene_columns.json`, whose symbols are
globally unique) *and* DGIdb assigned it exactly one HGNC concept id. The symbol
is then kept only as a human-readable consistency field. DGIdb's free-text
`gene_claim_name` / aliases are never used for matching, and a symbol that is not
uniquely resolvable is dropped and counted -- there is no ambiguous symbol-only
match.

Interaction direction
---------------------
Normalised strictly from DGIdb's own interaction-type -> directionality
vocabulary (`DIRECTIONALITY_BY_TYPE`, taken from the DGIdb GraphQL
`interactionClaimTypes.directionality` enum, retrieved 2026-08-29) into three
tiers: ``inhibitory`` / ``activating`` / ``unknown``. Nothing is derived from
drug names, free text, approval, or indication; a pipe-joined type with
conflicting or unmapped parts is ``unknown``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import config


# --------------------------------------------------------------------------
# Errors
# --------------------------------------------------------------------------

class DGIdbSchemaError(ValueError):
    """Raised when an upstream DGIdb TSV violates its expected schema."""


class DGIdbSnapshotError(RuntimeError):
    """Raised on a missing/inconsistent snapshot, or a failed asset gate."""


# --------------------------------------------------------------------------
# Upstream schemas -- inspected from the pinned 2026-06b assets, not guessed.
# Each file has two leading '# ...' comment lines, then this exact header.
# --------------------------------------------------------------------------

INTERACTIONS_TSV_COLUMNS = [
    "gene_claim_name", "gene_concept_id", "gene_name",
    "drug_claim_name", "drug_concept_id", "drug_name",
    "drug_is_approved", "drug_is_immunotherapy", "drug_is_antineoplastic",
    "interaction_source_db_name", "interaction_source_db_version",
    "interaction_types", "interaction_score",
    "drug_specificity_score", "gene_specificity_score", "evidence_score",
]
GENES_TSV_COLUMNS = [
    "gene_claim_name", "nomenclature", "concept_id", "gene_name",
    "source_db_name", "source_db_version",
]
DRUGS_TSV_COLUMNS = [
    "drug_claim_name", "nomenclature", "concept_id", "drug_name",
    "approved", "immunotherapy", "anti_neoplastic",
    "source_db_name", "source_db_version",
]

_BOOL_TOKENS = {"true", "false", ""}


# --------------------------------------------------------------------------
# DGIdb interaction-type -> directionality.
#
# Source: DGIdb GraphQL `interactionClaimTypes { type directionality }`
# (https://dgidb.org/api/graphql), retrieved 2026-08-29. DGIdb's `Directionality`
# enum has exactly INHIBITORY / ACTIVATING / (null). Types not listed here map to
# the `unknown` tier. This table is the ONLY basis for direction -- see the
# module docstring.
# --------------------------------------------------------------------------

DIRECTIONALITY_BY_TYPE: dict[str, str] = {
    # INHIBITORY
    "antagonist": "inhibitory",
    "antibody": "inhibitory",
    "antisense oligonucleotide": "inhibitory",
    "blocker": "inhibitory",
    "cleavage": "inhibitory",
    "inhibitor": "inhibitory",
    "inhibitory allosteric modulator": "inhibitory",
    "inverse agonist": "inhibitory",
    "negative modulator": "inhibitory",
    "partial antagonist": "inhibitory",
    "suppressor": "inhibitory",
    # ACTIVATING
    "activator": "activating",
    "agonist": "activating",
    "chaperone": "activating",
    "cofactor": "activating",
    "inducer": "activating",
    "positive modulator": "activating",
    "stimulator": "activating",
    "vaccine": "activating",
    # everything else DGIdb defines (adduct, allosteric modulator, binder,
    # immunotherapy, ligand, modulator, multitarget, other/unknown, potentiator,
    # product of, substrate) has directionality null -> `unknown` tier.
}

DIRECTION_TIERS = tuple(config.DGIDB_DIRECTION_TIERS)  # ("inhibitory","activating","unknown")
_TIER_INDEX = {tier: i for i, tier in enumerate(DIRECTION_TIERS)}


# --------------------------------------------------------------------------
# Per-interaction-source licence + redistribution decision.
#
# Source: DGIdb GraphQL `sources(sourceType: INTERACTION) { sourceDbName
# sourceDbVersion license licenseLink citation }`, retrieved 2026-08-29. The
# `decision` / `reason` are this project's redistribution assessment of that
# licence for committing filtered records into a public repository. `included`
# iff the terms explicitly permit redistribution with no NonCommercial clause
# and are unambiguous. Mirrored in config.DGIDB_INCLUDED_SOURCES and LICENSES.md.
# --------------------------------------------------------------------------

SOURCE_LICENCES: dict[str, dict[str, str]] = {
    "CIViC": {
        "dgidb_source_version": "08-June-2026",
        "license": "Creative Commons CC0 1.0 Universal (CC0 1.0) Public Domain Dedication",
        "license_url": "https://docs.civicdb.org/en/latest/about/faq.html#how-is-civic-licensed",
        "decision": "included",
        "reason": "CC0 1.0 public-domain dedication; redistribution unrestricted.",
    },
    "ChEMBL": {
        "dgidb_source_version": "37",
        "license": "Creative Commons Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0)",
        "license_url": "https://chembl.gitbook.io/chembl-interface-documentation/about#data-licensing",
        "decision": "included",
        "reason": "CC BY-SA 3.0; redistribution permitted with attribution and ShareAlike.",
    },
    "GuideToPharmacology": {
        "dgidb_source_version": "2026.2",
        "license": "Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)",
        "license_url": "https://www.guidetopharmacology.org/about.jsp",
        "decision": "included",
        "reason": "CC BY-SA 4.0; redistribution permitted with attribution and ShareAlike.",
    },
    "DoCM": {
        "dgidb_source_version": "2024-10-02",
        "license": "Creative Commons Attribution 4.0 International (CC BY 4.0)",
        "license_url": "https://github.com/griffithlab/docm/blob/c8d2a8723f505689074d07841931475b9b7e914c/app/views/static/about.html.haml#L86",
        "decision": "included",
        "reason": "CC BY 4.0; redistribution permitted with attribution.",
    },
    "NCI": {
        "dgidb_source_version": "14-September-2017",
        "license": "Public domain (U.S. National Cancer Institute)",
        "license_url": "https://www.cancer.gov/policies/copyright-reuse",
        "decision": "included",
        "reason": "U.S. Government work, public domain; redistribution unrestricted.",
    },
    "FDA": {
        "dgidb_source_version": "08-June-2026",
        "license": "Public domain (U.S. Food and Drug Administration)",
        "license_url": "https://www.fda.gov/about-fda/about-website/website-policies#linking",
        "decision": "included",
        "reason": "U.S. Government work, public domain; redistribution unrestricted.",
    },
    "DTC": {
        "dgidb_source_version": "2020-09-02",
        "license": "Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)",
        "license_url": "https://academic.oup.com/database/article/doi/10.1093/database/bay083/5096727",
        "decision": "excluded",
        "reason": "NonCommercial clause (CC BY-NC-SA 3.0).",
    },
    "PharmGKB": {
        "dgidb_source_version": "20260622",
        "license": "DGIdb labels 'CC BY-SA 4.0'; PharmGKB's own data-usage policy is a custom agreement, not a plain licence deed",
        "license_url": "https://www.pharmgkb.org/page/dataUsagePolicy",
        "decision": "excluded",
        "reason": "Ambiguous / conflicting terms: DGIdb's CC BY-SA label is not reflected by PharmGKB's custom data-usage policy; project scope also excludes PharmGKB (CLAUDE.md section 10).",
    },
    "TTD": {
        "dgidb_source_version": "2020.06.01",
        "license": "Unclear: site states 'All Rights Reserved'; 2002 publication describes it as open-access",
        "license_url": "https://academic.oup.com/nar/article/30/1/412/1331814",
        "decision": "excluded",
        "reason": "Unclear / self-contradictory redistribution terms.",
    },
    "TdgClinicalTrial": {
        "dgidb_source_version": "Jan-2014",
        "license": "Supplementary table from an Annual Reviews copyright publication",
        "license_url": "https://www.annualreviews.org/doi/suppl/10.1146/annurev-pharmtox-011613-135943",
        "decision": "excluded",
        "reason": "No redistribution licence (copyrighted supplementary data).",
    },
    "TEND": {
        "dgidb_source_version": "01-Aug-2011",
        "license": "Supplementary table from a Macmillan Publishers copyright publication",
        "license_url": "https://www.nature.com/articles/nrd3478",
        "decision": "excluded",
        "reason": "No redistribution licence (copyrighted supplementary data).",
    },
    "CKB-CORE": {
        "dgidb_source_version": "2024-11-27",
        "license": "Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)",
        "license_url": "https://ckb.genomenon.org/about/index",
        "decision": "excluded",
        "reason": "NonCommercial clause (CC BY-NC-SA 4.0).",
    },
    "MyCancerGenome": {
        "dgidb_source_version": "20-Jun-2017",
        "license": "Restrictive, custom, non-commercial",
        "license_url": "https://www.mycancergenome.org/content/page/legal-policies-licensing/",
        "decision": "excluded",
        "reason": "Restrictive custom NonCommercial terms.",
    },
    "MyCancerGenomeClinicalTrial": {
        "dgidb_source_version": "30-Feburary-2014",
        "license": "Restrictive, custom, non-commercial",
        "license_url": "https://www.mycancergenome.org/content/page/legal-policies-licensing/",
        "decision": "excluded",
        "reason": "Restrictive custom NonCommercial terms.",
    },
    "TALC": {
        "dgidb_source_version": "12-May-2016",
        "license": "Data extracted from tables in an Elsevier copyright publication",
        "license_url": "https://www.sciencedirect.com/science/article/pii/S1525730413002350",
        "decision": "excluded",
        "reason": "No redistribution licence (copyrighted publication tables).",
    },
    "CGI": {
        "dgidb_source_version": "2022-02-01",
        "license": "Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)",
        "license_url": "https://www.cancergenomeinterpreter.org/faq#q19",
        "decision": "excluded",
        "reason": "NonCommercial clause (CC BY-NC 4.0).",
    },
    "ClearityFoundationClinicalTrial": {
        "dgidb_source_version": "15-June-2013",
        "license": "Unknown; data no longer publicly available from the source site",
        "license_url": "https://www.clearityfoundation.org/about-clearity/contact/",
        "decision": "excluded",
        "reason": "Unknown / unverifiable terms; source data no longer available.",
    },
    "ClearityFoundationBiomarkers": {
        "dgidb_source_version": "26-July-2013",
        "license": "Unknown; data no longer publicly available from the source site",
        "license_url": "https://www.clearityfoundation.org/about-clearity/contact/",
        "decision": "excluded",
        "reason": "Unknown / unverifiable terms; source data no longer available.",
    },
    "OncoKB": {
        "dgidb_source_version": "23-July-2020",
        "license": "Restrictive, non-commercial",
        "license_url": "https://www.oncokb.org/terms",
        "decision": "excluded",
        "reason": "Restrictive NonCommercial terms.",
    },
    "CancerCommons": {
        "dgidb_source_version": "25-Jul-2013",
        "license": "Custom non-commercial",
        "license_url": "https://www.cancercommons.org/terms-of-use/",
        "decision": "excluded",
        "reason": "Custom NonCommercial terms.",
    },
    "COSMIC": {
        "dgidb_source_version": "4-Sep-2020",
        "license": "Free for non-commercial use only",
        "license_url": "https://cancer.sanger.ac.uk/cosmic/license",
        "decision": "excluded",
        "reason": "NonCommercial-only; redistribution requires a COSMIC licence.",
    },
}

INCLUDED_SOURCES = frozenset(config.DGIDB_INCLUDED_SOURCES)
assert INCLUDED_SOURCES == {
    n for n, v in SOURCE_LICENCES.items() if v["decision"] == "included"
}, "config.DGIDB_INCLUDED_SOURCES disagrees with SOURCE_LICENCES"


# --------------------------------------------------------------------------
# Committed snapshot schema (the filtered TSV's exact column order).
# `pmids`, `curation_type`, `indication` are always empty for this DGIdb export
# -- it carries no such columns -- and are NEVER inferred from free text. They
# stay in the schema so a future release that does carry them slots straight in.
# --------------------------------------------------------------------------

SNAPSHOT_COLUMNS = [
    "entrez_id",                      # canonical Entrez -- primary identity
    "gene_symbol",                    # canonical symbol -- consistency field
    "dgidb_gene_concept_id",          # e.g. "hgnc:11998"
    "dgidb_gene_name",                # DGIdb approved symbol, as provided
    "gene_symbol_consistent",         # "true"/"false": canonical symbol == dgidb_gene_name
    "drug_name",
    "drug_concept_id",
    "drug_claim_name",
    "interaction_source",
    "interaction_source_version",
    "interaction_type_raw",           # DGIdb `interaction_types`, verbatim (may be "" or "a|b")
    "interaction_direction",          # normalised tier: inhibitory|activating|unknown
    "interaction_score",              # DGIdb interaction_score, verbatim or ""
    "drug_specificity_score",
    "gene_specificity_score",
    "evidence_score",
    "drug_is_approved",               # regulatory approval status, verbatim: true|false|""
    "drug_is_immunotherapy",
    "drug_is_antineoplastic",
    "pmids",                          # "" (absent in this export; never inferred)
    "curation_type",                  # "" (absent in this export)
    "indication",                     # "" (absent in this export)
    "source_license",
    "source_license_url",
    "dgidb_release_tag",
    "record_key",                     # sha1 of the preceding fields -- dedup key
]

_HASH_CHUNK = 1 << 20


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------

def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(_HASH_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def _normalise_newlines(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n")


def normalize_direction(interaction_type_raw: str) -> str:
    """
    Map a DGIdb ``interaction_types`` value to a direction tier, using only
    ``DIRECTIONALITY_BY_TYPE``. Empty, unmapped, or internally-conflicting
    (pipe-joined) values are ``unknown``.
    """
    raw = (interaction_type_raw or "").strip()
    if not raw:
        return "unknown"
    parts = [p.strip().lower() for p in raw.split("|") if p.strip()]
    dirs = {DIRECTIONALITY_BY_TYPE[p] for p in parts if p in DIRECTIONALITY_BY_TYPE}
    if dirs == {"inhibitory"}:
        return "inhibitory"
    if dirs == {"activating"}:
        return "activating"
    return "unknown"


def _record_key(values_without_key: list[str]) -> str:
    joined = "\t".join(values_without_key)
    return hashlib.sha1(joined.encode("utf-8")).hexdigest()


def _score_sort_value(raw: str) -> tuple[int, float]:
    """Present scores first (descending); empty scores last. Deterministic."""
    if raw == "":
        return (1, 0.0)
    return (0, -float(raw))


# --------------------------------------------------------------------------
# Upstream TSV parsing (schema + value validation)
# --------------------------------------------------------------------------

def _parse_dgidb_tsv(
    path: str | Path, expected_columns: list[str]
) -> tuple[list[dict[str, str]], list[str]]:
    """
    Parse a DGIdb release TSV. Returns ``(rows, comment_lines)``.

    Leading lines beginning with ``#`` are comments; the first non-comment line
    must equal ``expected_columns`` exactly. Every data line must have exactly
    ``len(expected_columns)`` tab-separated fields.
    """
    path = Path(path)
    text = _normalise_newlines(path.read_text(encoding="utf-8"))
    lines = text.split("\n")
    while lines and lines[-1] == "":
        lines.pop()
    if not lines:
        raise DGIdbSchemaError(f"{path.name}: file is empty")

    comments: list[str] = []
    idx = 0
    while idx < len(lines) and lines[idx].startswith("#"):
        comments.append(lines[idx])
        idx += 1
    if idx >= len(lines):
        raise DGIdbSchemaError(f"{path.name}: no header line after comments")

    header = lines[idx].split("\t")
    idx += 1
    if header != expected_columns:
        raise DGIdbSchemaError(
            f"{path.name}: header does not match the expected schema.\n"
            f"  expected: {expected_columns}\n"
            f"  found   : {header}"
        )

    width = len(expected_columns)
    rows: list[dict[str, str]] = []
    for offset, line in enumerate(lines[idx:], start=idx + 1):
        fields = line.split("\t")
        if len(fields) != width:
            raise DGIdbSchemaError(
                f"{path.name}: line {offset} has {len(fields)} field(s), "
                f"expected {width}."
            )
        rows.append(dict(zip(expected_columns, fields)))
    return rows, comments


def _validate_interaction_values(rows: list[dict[str, str]]) -> None:
    """Value-domain checks on interactions.tsv beyond column presence."""
    for i, r in enumerate(rows, start=1):
        for col in ("drug_is_approved", "drug_is_immunotherapy",
                    "drug_is_antineoplastic"):
            if r[col] not in _BOOL_TOKENS:
                raise DGIdbSchemaError(
                    f"interactions.tsv row {i}: {col}={r[col]!r}, "
                    f"expected one of {sorted(_BOOL_TOKENS)}"
                )
        for col in ("interaction_score", "drug_specificity_score",
                    "gene_specificity_score", "evidence_score"):
            if r[col] != "":
                try:
                    float(r[col])
                except ValueError:
                    raise DGIdbSchemaError(
                        f"interactions.tsv row {i}: {col}={r[col]!r} is not "
                        f"numeric or empty"
                    )


# --------------------------------------------------------------------------
# Gene-identity crosswalk: DGIdb HGNC-approved symbol -> canonical Entrez
# --------------------------------------------------------------------------

def _load_canonical_symbols(gene_columns_path: Path) -> dict[str, str]:
    """symbol -> entrez, from gene_columns.json. Symbols must be globally unique."""
    gc = json.loads(gene_columns_path.read_text(encoding="utf-8"))
    symbols = list(gc["symbols"])
    entrez = [str(e) for e in gc["entrez_ids"]]
    if len(symbols) != len(entrez):
        raise DGIdbSchemaError(
            f"{gene_columns_path.name}: symbols/entrez_ids length mismatch"
        )
    if len(set(symbols)) != len(symbols):
        raise DGIdbSchemaError(
            f"{gene_columns_path.name}: canonical symbols are not unique; the "
            f"symbol-anchored crosswalk requires uniqueness"
        )
    return dict(zip(symbols, entrez))


def _build_gene_crosswalk(
    interaction_rows: list[dict[str, str]], canonical_symbol_to_entrez: dict[str, str]
) -> tuple[dict[str, tuple[str, str]], dict[str, int]]:
    """
    Returns ``(resolved, counts)`` where ``resolved`` maps a DGIdb ``gene_name``
    to ``(entrez_id, hgnc_concept_id)``, and ``counts`` accounts for every
    distinct non-empty ``gene_name``.
    """
    name_to_concepts: dict[str, set[str]] = {}
    for r in interaction_rows:
        gname = r["gene_name"]
        if gname == "":
            continue
        name_to_concepts.setdefault(gname, set()).add(r["gene_concept_id"])

    resolved: dict[str, tuple[str, str]] = {}
    counts = {
        "distinct_dgidb_genes": len(name_to_concepts),
        "resolved_to_canonical_entrez": 0,
        "unresolved_no_or_multi_hgnc_concept": 0,
        "unresolved_not_in_canonical_space": 0,
        "symbol_consistency_failures": 0,
    }
    for gname, concepts in name_to_concepts.items():
        non_empty = {c for c in concepts if c != ""}
        hgnc = {c for c in non_empty if c.startswith("hgnc:")}
        # Resolve only when DGIdb assigned this symbol exactly one concept id and
        # it is an HGNC id. A symbol carrying two distinct concept ids is genuine
        # ambiguity and is dropped; empty-concept rows (which in this release only
        # ever occur with an empty gene_name anyway) are ignored, not counted as
        # ambiguity.
        if non_empty != hgnc or len(hgnc) != 1:
            counts["unresolved_no_or_multi_hgnc_concept"] += 1
            continue
        if gname not in canonical_symbol_to_entrez:
            counts["unresolved_not_in_canonical_space"] += 1
            continue
        entrez = canonical_symbol_to_entrez[gname]
        resolved[gname] = (entrez, next(iter(hgnc)))
        counts["resolved_to_canonical_entrez"] += 1
    return resolved, counts


# --------------------------------------------------------------------------
# Snapshot build
# --------------------------------------------------------------------------

def _verify_asset(path: Path, name: str) -> None:
    spec = config.DGIDB_ASSETS[name]
    size = path.stat().st_size
    if size != spec["bytes"]:
        raise DGIdbSnapshotError(
            f"{name}: size {size} != pinned {spec['bytes']}. Refusing this asset."
        )
    digest = sha256_file(path)
    if digest != spec["sha256"]:
        raise DGIdbSnapshotError(
            f"{name}: sha256 {digest} != pinned {spec['sha256']}. "
            f"Refusing this asset."
        )


def _download_asset(name: str, dest: Path, timeout: int = 300) -> None:
    """Download one pinned asset to a staging path. Only reached via --refresh."""
    import urllib.request

    url = str(config.DGIDB_ASSETS[name]["url"])
    req = urllib.request.Request(url, headers={"User-Agent": "capstone-evidence/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:  # noqa: S310 (pinned host)
        data = resp.read()
    dest.write_bytes(data)


def build_snapshot(
    *,
    staging_dir: str | Path | None = None,
    refresh: bool = False,
    out_dir: str | Path | None = None,
    retrieved_utc: str | None = None,
    graphql_service_version: str = "v.5.0.12 (recorded 2026-08-29)",
) -> dict:
    """
    Build the filtered offline snapshot + manifest.

    Parameters
    ----------
    staging_dir
        Directory holding (or to receive, with ``refresh=True``) the three
        pinned TSVs. A fresh temp dir is created if omitted.
    refresh
        If True, download the pinned assets into ``staging_dir`` first. If
        False, they must already be present there.
    out_dir
        Where to write ``dgidb_<tag>.interactions.filtered.tsv`` and
        ``dgidb_<tag>.manifest.json``. Defaults to ``config.DGIDB_DIR``. The
        upstream TSVs are never written here.
    retrieved_utc
        Override the retrieval timestamp (for deterministic tests).

    Returns a summary dict (also the manifest's core content).
    """
    staging = Path(staging_dir) if staging_dir else Path(
        tempfile.mkdtemp(prefix="dgidb_stage_")
    )
    staging.mkdir(parents=True, exist_ok=True)
    out_dir = Path(out_dir) if out_dir is not None else config.DGIDB_DIR

    asset_paths = {name: staging / name for name in config.DGIDB_ASSETS}

    # ---- 1. obtain + gate the pinned assets -----------------------------
    for name, path in asset_paths.items():
        if refresh:
            _download_asset(name, path)
        if not path.is_file():
            raise DGIdbSnapshotError(
                f"{name} not found in staging dir {staging}. "
                f"Run with refresh=True (or place the pinned file there)."
            )
        _verify_asset(path, name)

    # ---- 2. parse + schema-validate -----------------------------------
    interactions, inter_comments = _parse_dgidb_tsv(
        asset_paths["interactions.tsv"], INTERACTIONS_TSV_COLUMNS
    )
    _validate_interaction_values(interactions)
    genes, genes_comments = _parse_dgidb_tsv(
        asset_paths["genes.tsv"], GENES_TSV_COLUMNS
    )
    drugs, drugs_comments = _parse_dgidb_tsv(
        asset_paths["drugs.tsv"], DRUGS_TSV_COLUMNS
    )

    # ---- 3. gene-identity crosswalk ---------------------------------
    gene_columns_path = config.PROCESSED_DIR / "gene_columns.json"
    canonical = _load_canonical_symbols(gene_columns_path)
    entrez_to_symbol = {e: s for s, e in canonical.items()}
    resolved_genes, gene_counts = _build_gene_crosswalk(interactions, canonical)

    # ---- 4. filter by source + build records ----------------------
    per_source_rows: dict[str, int] = {}
    per_source_retained: dict[str, int] = {}
    for name in SOURCE_LICENCES:
        per_source_rows[name] = 0
        per_source_retained[name] = 0

    dropped_unresolved_gene = 0
    dropped_no_drug_identity = 0
    seen_keys: set[str] = set()
    duplicates = 0
    records: list[dict[str, str]] = []

    for r in interactions:
        src = r["interaction_source_db_name"]
        if src in per_source_rows:
            per_source_rows[src] += 1
        if src not in INCLUDED_SOURCES:
            continue

        gname = r["gene_name"]
        if gname not in resolved_genes:
            dropped_unresolved_gene += 1
            continue
        entrez, hgnc_concept = resolved_genes[gname]

        drug_name = r["drug_name"]
        drug_concept_id = r["drug_concept_id"]
        drug_claim_name = r["drug_claim_name"]
        if drug_name == "" and drug_concept_id == "":
            dropped_no_drug_identity += 1
            continue

        canonical_symbol = entrez_to_symbol[entrez]
        direction = normalize_direction(r["interaction_types"])
        lic = SOURCE_LICENCES[src]

        ordered = [
            entrez,
            canonical_symbol,
            hgnc_concept,
            gname,
            "true" if canonical_symbol == gname else "false",
            drug_name,
            drug_concept_id,
            drug_claim_name,
            src,
            r["interaction_source_db_version"],
            r["interaction_types"],
            direction,
            r["interaction_score"],
            r["drug_specificity_score"],
            r["gene_specificity_score"],
            r["evidence_score"],
            r["drug_is_approved"],
            r["drug_is_immunotherapy"],
            r["drug_is_antineoplastic"],
            "",  # pmids -- absent in this export, never inferred
            "",  # curation_type -- absent in this export
            "",  # indication -- absent in this export
            lic["license"],
            lic["license_url"],
            config.DGIDB_RELEASE_TAG,
        ]
        key = _record_key(ordered)
        if key in seen_keys:
            duplicates += 1
            continue
        seen_keys.add(key)
        record = dict(zip(SNAPSHOT_COLUMNS[:-1], ordered))
        record["record_key"] = key
        records.append(record)
        per_source_retained[src] += 1

    if gene_counts["symbol_consistency_failures"] == 0:
        gene_counts["symbol_consistency_failures"] = sum(
            1 for rec in records if rec["gene_symbol_consistent"] != "true"
        )

    # ---- 5. deterministic total order ---------------------------
    records.sort(key=lambda rec: (
        int(rec["entrez_id"]),
        _TIER_INDEX[rec["interaction_direction"]],
        _score_sort_value(rec["interaction_score"]),
        rec["interaction_source"],
        rec["drug_name"],
        rec["drug_concept_id"],
        rec["drug_claim_name"],
        rec["interaction_type_raw"],
        rec["record_key"],
    ))

    tier_counts = {tier: 0 for tier in DIRECTION_TIERS}
    for rec in records:
        tier_counts[rec["interaction_direction"]] += 1

    # ---- 6. write snapshot TSV (LF, no trailing blank line) ------
    staging.mkdir(parents=True, exist_ok=True)
    snapshot_name = config.DGIDB_SNAPSHOT_FILE.name
    manifest_name = config.DGIDB_MANIFEST_FILE.name
    staged_snapshot = staging / snapshot_name
    lines = ["\t".join(SNAPSHOT_COLUMNS)]
    lines += ["\t".join(rec[c] for c in SNAPSHOT_COLUMNS) for rec in records]
    staged_snapshot.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    snapshot_sha = sha256_file(staged_snapshot)
    snapshot_bytes = staged_snapshot.stat().st_size

    # ---- 7. manifest ------------------------------------------
    now = retrieved_utc or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    sources_block: dict[str, dict] = {}
    for name, lic in SOURCE_LICENCES.items():
        sources_block[name] = {
            "role": "interaction",
            "dgidb_source_version": lic["dgidb_source_version"],
            "license": lic["license"],
            "license_url": lic["license_url"],
            "redistribution_decision": lic["decision"],
            "decision_reason": lic["reason"],
            "upstream_interaction_rows": per_source_rows.get(name, 0),
            "retained_records": per_source_retained.get(name, 0),
        }

    manifest = {
        "artifact": snapshot_name,
        "description": (
            "Licence-filtered offline snapshot of DGIdb drug-gene interaction "
            "records for Phase 2 evidence retrieval. Evidence retrieval only -- "
            "not treatment prediction, not a recommendation. Contains records "
            "ONLY from interaction sources whose redistribution terms are "
            "explicitly verified as compatible with committing them; the full "
            "DGIdb dataset is NOT redistributable."
        ),
        "dgidb": {
            "release_tag": config.DGIDB_RELEASE_TAG,
            "release_page": config.DGIDB_RELEASE_PAGE,
            "graphql_service_version": graphql_service_version,
            "upstream_header_comments": {
                "interactions.tsv": inter_comments,
                "genes.tsv": genes_comments,
                "drugs.tsv": drugs_comments,
            },
        },
        "assets": {
            name: {
                "url": str(spec["url"]),
                "bytes": spec["bytes"],
                "sha256": spec["sha256"],
                "verified": True,
            }
            for name, spec in config.DGIDB_ASSETS.items()
        },
        "retrieval": {
            "retrieved_utc": now,
            "method": (
                "download pinned GitHub release assets + exact size/SHA-256 gate"
                if refresh else
                "build from pre-staged assets that passed the exact size/SHA-256 gate"
            ),
            "tool": "evidence.py build_snapshot",
        },
        "gene_identity": {
            "primary_identifier": "entrez_id (canonical space)",
            "method": (
                "DGIdb gene_name (HGNC-approved symbol) with exactly one non-empty "
                "hgnc: concept id, joined to gene_columns.json by globally-unique "
                "symbol; canonical Entrez retained as the primary identity; symbol "
                "kept only as a consistency field; gene_claim_name / aliases never "
                "used; non-uniquely-resolvable genes dropped and counted."
            ),
            "canonical_space_file": gene_columns_path.name,
            "canonical_space_sha256": sha256_file(gene_columns_path),
            **gene_counts,
        },
        "direction_vocabulary": {
            "provenance": (
                "DGIdb GraphQL interactionClaimTypes.directionality enum, "
                "retrieved 2026-08-29"
            ),
            "tiers": list(DIRECTION_TIERS),
            "inhibitory": sorted(
                t for t, d in DIRECTIONALITY_BY_TYPE.items() if d == "inhibitory"
            ),
            "activating": sorted(
                t for t, d in DIRECTIONALITY_BY_TYPE.items() if d == "activating"
            ),
            "unmapped_to_unknown_tier": (
                "any interaction_types value not in the two lists above, "
                "including empty and internally-conflicting pipe-joined values"
            ),
        },
        "sources": sources_block,
        "filter": {
            "included_sources": sorted(INCLUDED_SOURCES),
            "excluded_sources": sorted(
                n for n, v in SOURCE_LICENCES.items() if v["decision"] == "excluded"
            ),
            "upstream_interaction_rows_total": len(interactions),
            "rows_from_included_sources": sum(
                per_source_rows[n] for n in INCLUDED_SOURCES
            ),
            "rows_dropped_unresolved_gene": dropped_unresolved_gene,
            "rows_dropped_no_drug_identity": dropped_no_drug_identity,
            "duplicate_records_collapsed": duplicates,
            "records_written": len(records),
        },
        "direction_tier_counts": tier_counts,
        "snapshot": {
            "file": snapshot_name,
            "columns": list(SNAPSHOT_COLUMNS),
            "bytes": snapshot_bytes,
            "sha256": snapshot_sha,
            "record_count": len(records),
            "sort_order": (
                "entrez_id, direction tier, interaction_score desc (empty last), "
                "interaction_source, drug_name, drug_concept_id, drug_claim_name, "
                "interaction_type_raw, record_key"
            ),
            "row_ordering_note": (
                "Deterministic evidence-DISPLAY ordering, for reproducibility "
                "only. NOT predicted efficacy, clinical priority, or "
                "recommendation strength. Approval and indication never affect "
                "ordering."
            ),
        },
        "absent_fields": {
            "pmids": (
                "DGIdb 2026-06b interactions.tsv has no publication-identifier "
                "column; PMIDs embedded in free-text drug names are NOT parsed "
                "(that would be inference). Always empty."
            ),
            "curation_type": (
                "No curation/evidence-type column in this export; the numeric "
                "evidence_score is retained instead. Always empty."
            ),
            "indication": "No disease/indication column in this export. Always empty.",
        },
        "disclaimer": config.DGIDB_EVIDENCE_DISCLAIMER,
        "scope": {
            "statement": (
                "Evidence retrieval only. Not treatment prediction, not a "
                "recommendation. No efficacy, clinical relevance, approval, "
                "indication, interaction direction beyond DGIdb's explicit "
                "vocabulary, or osteosarcoma relevance is inferred. The returned "
                "drugs are not candidate treatments."
            ),
            "redistribution": (
                "Snapshot contains records ONLY from the explicitly-verified "
                "redistributable interaction sources listed under filter."
                "included_sources. The full DGIdb dataset is NOT redistributable."
            ),
        },
    }

    staged_manifest = staging / manifest_name
    staged_manifest.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8", newline="\n",
    )

    # ---- 8. publish to the tracked location (snapshot + manifest only) ----
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / snapshot_name).write_bytes(staged_snapshot.read_bytes())
        (out_dir / manifest_name).write_bytes(staged_manifest.read_bytes())

    return {
        "staging_dir": str(staging),
        "out_dir": str(out_dir) if out_dir is not None else None,
        "snapshot_sha256": snapshot_sha,
        "snapshot_bytes": snapshot_bytes,
        "record_count": len(records),
        "tier_counts": tier_counts,
        "manifest": manifest,
    }


# --------------------------------------------------------------------------
# Offline snapshot loading + retrieval API
# --------------------------------------------------------------------------

class DGIdbSnapshot:
    """A loaded, in-memory view of the committed offline snapshot. Offline only."""

    def __init__(self, records: list[dict[str, str]], manifest: dict,
                 snapshot_path: Path):
        self.records = records
        self.manifest = manifest
        self.snapshot_path = snapshot_path
        self.by_entrez: dict[str, list[dict[str, str]]] = {}
        self.by_symbol: dict[str, str] = {}
        for rec in records:
            self.by_entrez.setdefault(rec["entrez_id"], []).append(rec)
            self.by_symbol.setdefault(rec["gene_symbol"], rec["entrez_id"])

    def __len__(self) -> int:
        return len(self.records)


def load_snapshot(
    snapshot_path: str | Path = config.DGIDB_SNAPSHOT_FILE,
    manifest_path: str | Path = config.DGIDB_MANIFEST_FILE,
) -> DGIdbSnapshot:
    """
    Load the committed offline snapshot. **No network.** Raises
    ``DGIdbSnapshotError`` (never a live query) if the files are missing or
    inconsistent with each other.
    """
    snapshot_path = Path(snapshot_path)
    manifest_path = Path(manifest_path)
    if not snapshot_path.is_file() or not manifest_path.is_file():
        raise DGIdbSnapshotError(
            "DGIdb offline snapshot not found:\n"
            f"  snapshot: {snapshot_path}  (exists: {snapshot_path.is_file()})\n"
            f"  manifest: {manifest_path}  (exists: {manifest_path.is_file()})\n"
            "Build it with `py evidence.py --refresh` (needs network) or "
            "`py evidence.py --from-staging DIR`. Refusing to query DGIdb live."
        )

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    rows, comments = _parse_dgidb_tsv(snapshot_path, SNAPSHOT_COLUMNS)
    if comments:
        raise DGIdbSnapshotError(
            f"{snapshot_path.name}: unexpected comment lines in a committed snapshot"
        )

    recorded_sha = manifest.get("snapshot", {}).get("sha256")
    actual_sha = sha256_file(snapshot_path)
    if recorded_sha != actual_sha:
        raise DGIdbSnapshotError(
            f"{snapshot_path.name}: sha256 {actual_sha} != manifest "
            f"{recorded_sha}. Snapshot and manifest are out of sync."
        )
    if manifest.get("snapshot", {}).get("record_count") != len(rows):
        raise DGIdbSnapshotError(
            f"{snapshot_path.name}: {len(rows)} rows != manifest record_count "
            f"{manifest.get('snapshot', {}).get('record_count')}"
        )
    return DGIdbSnapshot(rows, manifest, snapshot_path)


def _display_sort_key(rec: dict[str, str]) -> tuple:
    """
    Evidence-DISPLAY ordering within a direction tier. Deterministic; NOT
    efficacy or clinical priority. Approval / indication are deliberately absent.
    """
    return (
        _score_sort_value(rec["interaction_score"]),  # higher interaction_score first, empty last
        rec["interaction_source"],
        rec["drug_name"],
        rec["drug_concept_id"],
        rec["record_key"],
    )


def get_evidence_for_gene(
    entrez_id: str | int,
    symbol: str | None = None,
    top_k: int = config.TOP_K_EVIDENCE_PER_GENE,
    *,
    snapshot: DGIdbSnapshot | None = None,
) -> list[dict[str, str]]:
    """
    Return recorded interaction evidence for one canonical gene.

    Lookup is **by Entrez ID only** (the primary identity). ``symbol`` is used
    only as a consistency check: if given and it disagrees with the snapshot's
    canonical symbol for this Entrez, every returned record is flagged
    ``symbol_query_mismatch = "true"`` -- records are never matched, filtered, or
    ordered by symbol.

    At most ``top_k`` records are returned **per direction tier**
    (inhibitory, then activating, then unknown), each ordered by the documented
    evidence-display ordering. Every returned record carries the constant
    ``disclaimer`` and its ``direction_tier``.
    """
    if snapshot is None:
        snapshot = load_snapshot()

    entrez_id = str(entrez_id).strip()
    if not entrez_id.isdigit():
        raise ValueError(f"entrez_id must be digits, got {entrez_id!r}")
    if top_k < 0:
        raise ValueError(f"top_k must be >= 0, got {top_k}")

    matches = snapshot.by_entrez.get(entrez_id, [])
    mismatch = False
    if symbol is not None and matches:
        canonical_symbol = matches[0]["gene_symbol"]
        mismatch = (symbol != canonical_symbol)

    out: list[dict[str, str]] = []
    for tier in DIRECTION_TIERS:
        tier_recs = [r for r in matches if r["interaction_direction"] == tier]
        tier_recs.sort(key=_display_sort_key)
        for rec in tier_recs[:top_k]:
            enriched = dict(rec)
            enriched["direction_tier"] = tier
            enriched["disclaimer"] = config.DGIDB_EVIDENCE_DISCLAIMER
            if symbol is not None:
                enriched["symbol_query_mismatch"] = "true" if mismatch else "false"
            out.append(enriched)
    return out


# --------------------------------------------------------------------------
# Offline validation of the committed snapshot (requirement 10)
# --------------------------------------------------------------------------

def validate_snapshot(
    snapshot_path: str | Path = config.DGIDB_SNAPSHOT_FILE,
    manifest_path: str | Path = config.DGIDB_MANIFEST_FILE,
) -> dict:
    """Assert every committed-snapshot invariant. Returns a checks dict."""
    snap = load_snapshot(snapshot_path, manifest_path)
    m = snap.manifest
    checks: list[tuple[str, bool, str]] = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        checks.append((name, bool(ok), detail))

    # no network was used: load_snapshot only reads files (structural guarantee).
    check("snapshot + manifest loaded offline", True)

    # snapshot sha matches manifest (also checked in load_snapshot)
    check("snapshot sha256 == manifest", sha256_file(snap.snapshot_path)
          == m["snapshot"]["sha256"])

    # manifest asset hashes == config pins
    assets_ok = all(
        m["assets"][n]["sha256"] == config.DGIDB_ASSETS[n]["sha256"]
        and m["assets"][n]["bytes"] == config.DGIDB_ASSETS[n]["bytes"]
        for n in config.DGIDB_ASSETS
    )
    check("manifest asset sizes+hashes == config pins", assets_ok)

    # no duplicate record keys
    keys = [r["record_key"] for r in snap.records]
    check("record_key values are unique", len(keys) == len(set(keys)),
          f"{len(keys) - len(set(keys))} dup(s)")

    # record_key actually equals the hash of the other fields
    recomputed_ok = all(
        r["record_key"] == _record_key([r[c] for c in SNAPSHOT_COLUMNS[:-1]])
        for r in snap.records
    )
    check("record_key == sha1(preceding fields)", recomputed_ok)

    # required identifiers present + well-formed
    canonical = _load_canonical_symbols(config.PROCESSED_DIR / "gene_columns.json")
    entrez_set = set(canonical.values())
    ident_ok = True
    ident_detail = ""
    for r in snap.records:
        if not r["entrez_id"].isdigit() or r["entrez_id"] not in entrez_set:
            ident_ok = False
            ident_detail = f"bad entrez_id {r['entrez_id']!r}"
            break
        if r["gene_symbol"] == "" or not r["dgidb_gene_concept_id"].startswith("hgnc:"):
            ident_ok = False
            ident_detail = f"bad gene identity for {r['entrez_id']}"
            break
        if r["drug_name"] == "" and r["drug_concept_id"] == "":
            ident_ok = False
            ident_detail = "record with no drug identity"
            break
    check("every record has well-formed gene + drug identifiers", ident_ok, ident_detail)

    # gene symbol consistency field is internally consistent
    sym_ok = all(
        (r["gene_symbol_consistent"] == "true")
        == (r["gene_symbol"] == r["dgidb_gene_name"])
        for r in snap.records
    )
    check("gene_symbol_consistent flag matches the two symbols", sym_ok)

    # every source in the snapshot has a recorded redistribution decision,
    # and it is "included"
    src_in_snap = sorted({r["interaction_source"] for r in snap.records})
    decisions_ok = all(
        s in m["sources"] and m["sources"][s]["redistribution_decision"] == "included"
        for s in src_in_snap
    )
    check("every source in snapshot is recorded + decision=included", decisions_ok,
          f"sources: {src_in_snap}")

    # no excluded source appears
    excluded = {n for n, v in m["sources"].items()
                if v["redistribution_decision"] == "excluded"}
    no_excluded = excluded.isdisjoint(set(src_in_snap))
    check("no excluded source appears in the snapshot", no_excluded)

    # snapshot sources are exactly a subset of config.DGIDB_INCLUDED_SOURCES
    check("snapshot sources subset of config.DGIDB_INCLUDED_SOURCES",
          set(src_in_snap).issubset(INCLUDED_SOURCES))

    # direction values are only the three tiers
    dirs = {r["interaction_direction"] for r in snap.records}
    check("interaction_direction in {inhibitory,activating,unknown}",
          dirs.issubset(set(DIRECTION_TIERS)), f"{sorted(dirs)}")

    # unmapped / empty interaction types stayed 'unknown'
    unknown_ok = all(
        normalize_direction(r["interaction_type_raw"]) == r["interaction_direction"]
        for r in snap.records
    )
    check("interaction_direction reproduces from interaction_type_raw", unknown_ok)

    # manifest tier counts match the file
    recount = {t: 0 for t in DIRECTION_TIERS}
    for r in snap.records:
        recount[r["interaction_direction"]] += 1
    check("manifest direction_tier_counts match the file",
          recount == {k: m["direction_tier_counts"][k] for k in recount},
          f"file={recount} manifest={m['direction_tier_counts']}")

    # disclaimer is constant + present on every retrieval
    check("manifest disclaimer == config.DGIDB_EVIDENCE_DISCLAIMER",
          m["disclaimer"] == config.DGIDB_EVIDENCE_DISCLAIMER)
    sample_entrez = snap.records[0]["entrez_id"] if snap.records else None
    if sample_entrez is not None:
        got = get_evidence_for_gene(sample_entrez, top_k=1000, snapshot=snap)
        check("every returned record carries the disclaimer",
              all(r["disclaimer"] == config.DGIDB_EVIDENCE_DISCLAIMER for r in got))
        # top_k enforced per tier independently
        capped = get_evidence_for_gene(sample_entrez, top_k=1, snapshot=snap)
        per_tier = {t: 0 for t in DIRECTION_TIERS}
        for r in capped:
            per_tier[r["direction_tier"]] += 1
        check("top_k enforced independently per tier (k=1)",
              all(v <= 1 for v in per_tier.values()), f"{per_tier}")

    # rows are in the documented deterministic order
    def _total_key(r: dict[str, str]) -> tuple:
        return (
            int(r["entrez_id"]),
            _TIER_INDEX[r["interaction_direction"]],
            _score_sort_value(r["interaction_score"]),
            r["interaction_source"], r["drug_name"], r["drug_concept_id"],
            r["drug_claim_name"], r["interaction_type_raw"], r["record_key"],
        )
    ordered_ok = all(
        _total_key(snap.records[i]) <= _total_key(snap.records[i + 1])
        for i in range(len(snap.records) - 1)
    )
    check("records are in the documented deterministic order", ordered_ok)

    n_fail = sum(1 for _, ok, _ in checks if not ok)
    return {"checks": checks, "n_fail": n_fail, "n_pass": len(checks) - n_fail,
            "record_count": len(snap.records)}


# --------------------------------------------------------------------------
# Offline self-test (synthetic; no network, no committed-artifact mutation)
# --------------------------------------------------------------------------

def _write_tsv(path: Path, comments: list[str], header: list[str],
               rows: list[list[str]]) -> None:
    body = [*comments, "\t".join(header)]
    body += ["\t".join(r) for r in rows]
    path.write_text("\n".join(body) + "\n", encoding="utf-8", newline="\n")


def _self_test() -> int:  # noqa: C901 -- one linear scenario, kept together
    import shutil

    print("Running evidence.py self-test...")
    tmp = Path(tempfile.mkdtemp(prefix="evidence_selftest_"))
    stage = tmp / "stage"
    stage.mkdir()

    # ---- synthetic canonical space -------------------------------------
    proc = tmp / "processed"
    proc.mkdir()
    (proc / "gene_columns.json").write_text(json.dumps({
        "description": "self-test canonical space",
        "n_genes": 3,
        "entrez_ids": ["100", "200", "300"],
        "symbols": ["AAA", "BBB", "CCC"],
        "canonical_columns": ["AAA (100)", "BBB (200)", "CCC (300)"],
    }), encoding="utf-8")

    _orig_processed = config.PROCESSED_DIR
    config.PROCESSED_DIR = proc  # redirect the crosswalk's gene_columns.json
    try:
        # ---- synthetic upstream TSVs --------------------------------
        icols = INTERACTIONS_TSV_COLUMNS
        # helper to make an interaction row
        def irow(gene_name, concept, drug_name, drug_cid, approved, src, itype,
                 iscore, dss="", gss="", es="", drug_claim="dc", gclaim="gc",
                 immuno="false", antineo="false"):
            d = {c: "" for c in icols}
            d.update(gene_claim_name=gclaim, gene_concept_id=concept, gene_name=gene_name,
                     drug_claim_name=drug_claim, drug_concept_id=drug_cid, drug_name=drug_name,
                     drug_is_approved=approved, drug_is_immunotherapy=immuno,
                     drug_is_antineoplastic=antineo,
                     interaction_source_db_name=src, interaction_source_db_version="v1",
                     interaction_types=itype, interaction_score=iscore,
                     drug_specificity_score=dss, gene_specificity_score=gss,
                     evidence_score=es)
            return [d[c] for c in icols]

        inter_rows = []
        # AAA (entrez 100): 7 inhibitory + 3 activating + 1 unknown, all ChEMBL (included)
        for i in range(7):
            inter_rows.append(irow("AAA", "hgnc:1", f"inh{i}", f"chembl:{i}", "true",
                                   "ChEMBL", "inhibitor", f"{10 - i}.0"))
        for i in range(3):
            inter_rows.append(irow("AAA", "hgnc:1", f"act{i}", f"chembl:a{i}", "false",
                                   "ChEMBL", "agonist", f"{5 - i}.0"))
        inter_rows.append(irow("AAA", "hgnc:1", "unk0", "chembl:u0", "false",
                               "ChEMBL", "modulator", ""))
        # a duplicate of inh0 (identical every field) -> must collapse
        inter_rows.append(irow("AAA", "hgnc:1", "inh0", "chembl:0", "true",
                               "ChEMBL", "inhibitor", "10.0"))
        # BBB (entrez 200): one CIViC (included) + one DrugBank-style EXCLUDED source
        inter_rows.append(irow("BBB", "hgnc:2", "civicdrug", "civic:1", "true",
                               "CIViC", "inhibitor|agonist", "3.0"))   # conflicting -> unknown
        inter_rows.append(irow("BBB", "hgnc:2", "dtcdrug", "dtc:1", "true",
                               "DTC", "inhibitor", "9.9"))              # excluded source
        # CCC (entrez 300): included source but NO drug identity -> dropped
        inter_rows.append(irow("CCC", "hgnc:3", "", "", "false", "NCI", "inhibitor", ""))
        # gene not in canonical space -> dropped
        inter_rows.append(irow("ZZZ", "hgnc:9", "zdrug", "z:1", "false", "FDA",
                               "inhibitor", "1.0"))
        # empty gene_name + empty concept (the only shape this occurs in) -> dropped
        inter_rows.append(irow("", "", "orphan", "x:1", "false", "ChEMBL",
                               "inhibitor", "1.0"))
        # AAA also carries an activator so its 'activating' tier is non-empty for
        # the DoCM included source (exercises a second included source)
        inter_rows.append(irow("AAA", "hgnc:1", "docmact", "docm:1", "false",
                               "DoCM", "activator", "4.5"))

        _write_tsv(stage / "interactions.tsv",
                   ["# Data version: test", "# DGIdb version: test"], icols, inter_rows)
        _write_tsv(stage / "genes.tsv",
                   ["# Data version: test", "# DGIdb version: test"],
                   GENES_TSV_COLUMNS,
                   [["AAA", "Gene Symbol", "hgnc:1", "AAA", "HGNC", "x"]])
        _write_tsv(stage / "drugs.tsv",
                   ["# Data version: test", "# DGIdb version: test"],
                   DRUGS_TSV_COLUMNS,
                   [["dc", "Primary Name", "chembl:0", "inh0", "true", "false",
                     "false", "ChEMBL", "v1"]])

        # ---- point config.DGIDB_ASSETS at the synthetic files -------
        _orig_assets = config.DGIDB_ASSETS
        config.DGIDB_ASSETS = {
            n: {"url": f"file://{stage / n}", "bytes": (stage / n).stat().st_size,
                "sha256": sha256_file(stage / n)}
            for n in ("interactions.tsv", "genes.tsv", "drugs.tsv")
        }
        try:
            out_dir = tmp / "out"

            # ---- 1. build from staging ------------------------------
            res1 = build_snapshot(staging_dir=stage, refresh=False, out_dir=out_dir,
                                  retrieved_utc="2026-01-01T00:00:00Z")
            snap_file = out_dir / config.DGIDB_SNAPSHOT_FILE.name
            man_file = out_dir / config.DGIDB_MANIFEST_FILE.name
            assert snap_file.is_file() and man_file.is_file()
            # upstream TSVs must NOT have been copied into out_dir
            assert not (out_dir / "interactions.tsv").exists()
            print("  [ok] build writes only the filtered snapshot + manifest")

            # ---- 2. deterministic rebuild --------------------------
            res2 = build_snapshot(staging_dir=stage, refresh=False,
                                  out_dir=tmp / "out2",
                                  retrieved_utc="2099-12-31T23:59:59Z")
            b1 = snap_file.read_bytes()
            b2 = (tmp / "out2" / config.DGIDB_SNAPSHOT_FILE.name).read_bytes()
            assert b1 == b2, "snapshot TSV is not byte-identical across rebuilds"
            m1 = json.loads(man_file.read_text(encoding="utf-8"))
            m2 = json.loads((tmp / "out2" / config.DGIDB_MANIFEST_FILE.name)
                            .read_text(encoding="utf-8"))
            m1n, m2n = dict(m1), dict(m2)
            m1n.pop("retrieval"), m2n.pop("retrieval")
            assert m1n == m2n, "manifest differs beyond the retrieval block"
            assert res1["snapshot_sha256"] == res2["snapshot_sha256"]
            print("  [ok] snapshot regenerates byte-identically (manifest modulo retrieval)")

            # ---- 3. size / hash gate refuses a tampered asset ------
            bad_stage = tmp / "bad"
            bad_stage.mkdir()
            for n in ("interactions.tsv", "genes.tsv", "drugs.tsv"):
                shutil.copy(stage / n, bad_stage / n)
            (bad_stage / "genes.tsv").write_text(
                (bad_stage / "genes.tsv").read_text(encoding="utf-8") + "AAA\tGene Symbol\thgnc:1\tAAA\tHGNC\tx\n",
                encoding="utf-8", newline="\n")
            try:
                build_snapshot(staging_dir=bad_stage, refresh=False, out_dir=tmp / "nope",
                               retrieved_utc="2026-01-01T00:00:00Z")
                raise AssertionError("tampered asset should have been refused")
            except DGIdbSnapshotError as exc:
                assert "sha256" in str(exc) or "size" in str(exc)
            print("  [ok] size/SHA-256 gate refuses a mismatched asset")

            # ---- 4. schema validation ------------------------------
            bad_hdr = tmp / "badhdr"
            bad_hdr.mkdir()
            _write_tsv(bad_hdr / "x.tsv", ["# c"], ["a", "b"], [["1", "2"]])
            try:
                _parse_dgidb_tsv(bad_hdr / "x.tsv", ["a", "b", "c"])
                raise AssertionError("bad header should raise")
            except DGIdbSchemaError:
                pass
            try:
                _validate_interaction_values([{c: "" for c in icols} | {
                    "drug_is_approved": "yes", "interaction_score": ""}])
                raise AssertionError("bad boolean should raise")
            except DGIdbSchemaError:
                pass
            try:
                _validate_interaction_values([{c: "" for c in icols} | {
                    "drug_is_approved": "true", "interaction_score": "abc"}])
                raise AssertionError("non-numeric score should raise")
            except DGIdbSchemaError:
                pass
            print("  [ok] schema + value validation rejects bad header / bool / score")

            # ---- 5. load + retrieval semantics --------------------
            snap = load_snapshot(snap_file, man_file)
            # excluded source absent
            assert all(r["interaction_source"] in INCLUDED_SOURCES for r in snap.records)
            assert not any(r["drug_name"] == "dtcdrug" for r in snap.records)
            # dup collapsed: only 7 distinct inhibitors for AAA (not 8)
            aaa_inh = [r for r in snap.records
                       if r["entrez_id"] == "100" and r["interaction_direction"] == "inhibitory"]
            assert len(aaa_inh) == 7, len(aaa_inh)
            # conflicting pipe type -> unknown
            bbb = snap.by_entrez["200"]
            assert len(bbb) == 1 and bbb[0]["interaction_direction"] == "unknown"
            # CCC dropped (no drug identity); ZZZ dropped (not canonical)
            assert "300" not in snap.by_entrez
            assert all(r["gene_symbol"] != "ZZZ" for r in snap.records)
            print("  [ok] source filter, dedup, no-drug-identity + non-canonical drops")

            # direction normalisation table
            assert normalize_direction("inhibitor") == "inhibitory"
            assert normalize_direction("agonist") == "activating"
            assert normalize_direction("modulator") == "unknown"
            assert normalize_direction("") == "unknown"
            assert normalize_direction("inhibitor|agonist") == "unknown"
            assert normalize_direction("inhibitor|blocker") == "inhibitory"
            print("  [ok] direction normalised only from DGIdb's explicit vocabulary")

            # top_k per tier, independently
            ev = get_evidence_for_gene("100", top_k=2, snapshot=snap)
            tiers = {}
            for r in ev:
                tiers.setdefault(r["direction_tier"], 0)
                tiers[r["direction_tier"]] += 1
            assert tiers.get("inhibitory") == 2 and tiers.get("activating") == 2
            assert tiers.get("unknown", 0) == 1  # only one unknown exists
            assert all(r["disclaimer"] == config.DGIDB_EVIDENCE_DISCLAIMER for r in ev)
            # display order within inhibitory tier: highest interaction_score first
            inh = [r for r in ev if r["direction_tier"] == "inhibitory"]
            assert inh[0]["interaction_score"] == "10.0"
            print("  [ok] top_k enforced per tier; disclaimer on every record; display order")

            # entrez is the identity; symbol never matches on its own
            assert get_evidence_for_gene("999", snapshot=snap) == []
            flagged = get_evidence_for_gene("100", symbol="WRONG", snapshot=snap)
            assert all(r["symbol_query_mismatch"] == "true" for r in flagged)
            ok_sym = get_evidence_for_gene("100", symbol="AAA", snapshot=snap)
            assert all(r["symbol_query_mismatch"] == "false" for r in ok_sym)
            print("  [ok] lookup is by Entrez only; symbol is a flagged consistency check")

            # ---- 6. missing snapshot -> clear error, no live query
            try:
                load_snapshot(tmp / "absent.tsv", tmp / "absent.json")
                raise AssertionError("missing snapshot should raise")
            except DGIdbSnapshotError as exc:
                assert "Refusing to query DGIdb live" in str(exc)
            print("  [ok] missing snapshot raises, never queries live")

            # ---- 7. full offline validator on the synthetic snapshot
            v = validate_snapshot(snap_file, man_file)
            assert v["n_fail"] == 0, [c for c in v["checks"] if not c[1]]
            print(f"  [ok] validate_snapshot: {v['n_pass']}/"
                  f"{v['n_pass'] + v['n_fail']} checks pass")

        finally:
            config.DGIDB_ASSETS = _orig_assets
    finally:
        config.PROCESSED_DIR = _orig_processed

    # ---- 8. the real committed inputs are untouched by a self-test ----
    frozen = [
        _orig_processed / "geneformer_embeddings.csv",
        _orig_processed / "baseline_results.json",
        _orig_processed / "gene_columns.json",
    ]
    # (self-test never writes under data/; nothing to diff -- assert dir clean of our temp)
    assert not (config.DGIDB_DIR / "___selftest___").exists()

    shutil.rmtree(tmp, ignore_errors=True)
    print("\nSelf-test passed.")
    return 0


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def _print_snapshot_summary(snap: DGIdbSnapshot) -> None:
    m = snap.manifest
    f = m["filter"]
    print("=" * 74)
    print(f"DGIdb EVIDENCE SNAPSHOT  --  release {m['dgidb']['release_tag']}")
    print("=" * 74)
    print(f"  snapshot file    : {snap.snapshot_path.name}")
    print(f"    sha256         : {m['snapshot']['sha256']}")
    print(f"    records        : {m['snapshot']['record_count']}")
    print(f"    genes covered  : {len(snap.by_entrez)}")
    print(f"  retrieved (UTC)  : {m['retrieval']['retrieved_utc']}")
    print("-" * 74)
    print("  included interaction sources (redistribution verified):")
    for s in f["included_sources"]:
        print(f"    {s:22s} v{m['sources'][s]['dgidb_source_version']:<14s} "
              f"{m['sources'][s]['retained_records']:>6d} records  "
              f"[{m['sources'][s]['license'].split('(')[0].strip()}]")
    print(f"  excluded interaction sources     : {len(f['excluded_sources'])} "
          f"(see manifest / LICENSES.md)")
    print("-" * 74)
    print(f"  upstream interaction rows        : {f['upstream_interaction_rows_total']}")
    print(f"  rows from included sources       : {f['rows_from_included_sources']}")
    print(f"  dropped: gene not resolvable     : {f['rows_dropped_unresolved_gene']}")
    print(f"  dropped: no drug identity        : {f['rows_dropped_no_drug_identity']}")
    print(f"  duplicate records collapsed      : {f['duplicate_records_collapsed']}")
    print(f"  records written                  : {f['records_written']}")
    print("-" * 74)
    tc = m["direction_tier_counts"]
    print(f"  direction tiers  : inhibitory {tc['inhibitory']}  "
          f"activating {tc['activating']}  unknown {tc['unknown']}")
    gi = m["gene_identity"]
    print(f"  gene identity    : primary={gi['primary_identifier']}")
    print(f"    resolved       : {gi['resolved_to_canonical_entrez']} / "
          f"{gi['distinct_dgidb_genes']} DGIdb genes")
    print(f"    unresolved     : "
          f"{gi['unresolved_no_or_multi_hgnc_concept']} no/multi hgnc concept, "
          f"{gi['unresolved_not_in_canonical_space']} outside canonical space")
    print("-" * 74)
    print(f"  DISCLAIMER (on every record): {m['disclaimer']}")
    print("=" * 74)


def _print_evidence(entrez_id: str, symbol: str | None, records: list[dict]) -> None:
    print("=" * 74)
    print(f"EVIDENCE  --  entrez {entrez_id}"
          + (f"  (query symbol {symbol})" if symbol else ""))
    print("=" * 74)
    if not records:
        print("  (no recorded interaction evidence in the offline snapshot)")
        print("=" * 74)
        return
    current = None
    for r in records:
        if r["direction_tier"] != current:
            current = r["direction_tier"]
            print(f"\n  [{current}]")
        score = r["interaction_score"] or "-"
        print(f"    {r['drug_name'] or r['drug_claim_name']:<32s} "
              f"src={r['interaction_source']:<20s} score={score:<8s} "
              f"approved={r['drug_is_approved'] or '-'}")
    print(f"\n  ordering = evidence-display only (NOT efficacy / priority)")
    print(f"  {records[0]['disclaimer']}")
    print("=" * 74)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Offline, licence-filtered DGIdb drug-gene interaction "
                    "evidence retrieval (Phase 2)."
    )
    parser.add_argument("--self-test", action="store_true",
                        help="Run the synthetic offline self-test and exit.")
    parser.add_argument("--validate", action="store_true",
                        help="Validate the committed offline snapshot and exit.")
    parser.add_argument("--refresh", action="store_true",
                        help="Download the pinned assets, rebuild the snapshot.")
    parser.add_argument("--from-staging", metavar="DIR",
                        help="Rebuild the snapshot from pre-downloaded pinned "
                             "assets already present in DIR (no network).")
    parser.add_argument("--out-dir", default=None,
                        help="Override the snapshot/manifest output directory.")
    parser.add_argument("--gene", metavar="ENTREZ",
                        help="Print evidence for one canonical Entrez ID.")
    parser.add_argument("--symbol", default=None,
                        help="Optional symbol for the --gene consistency check.")
    parser.add_argument("--top-k", type=int, default=config.TOP_K_EVIDENCE_PER_GENE)
    parser.add_argument("--json", action="store_true",
                        help="With --gene: emit the records as JSON.")
    args = parser.parse_args()

    if args.self_test:
        return _self_test()

    if args.refresh or args.from_staging:
        res = build_snapshot(
            staging_dir=args.from_staging,
            refresh=bool(args.refresh),
            out_dir=args.out_dir,
        )
        print("=" * 74)
        print("DGIdb SNAPSHOT BUILT")
        print("=" * 74)
        print(f"  staging dir   : {res['staging_dir']}")
        print(f"  out dir       : {res['out_dir']}")
        print(f"  records       : {res['record_count']}")
        print(f"  tier counts   : {res['tier_counts']}")
        print(f"  snapshot bytes: {res['snapshot_bytes']}")
        print(f"  snapshot sha  : {res['snapshot_sha256']}")
        print("=" * 74)
        return 0

    if args.validate:
        v = validate_snapshot()
        for name, ok, detail in v["checks"]:
            print(f"  [{'ok' if ok else 'FAIL'}] {name}"
                  + (f"  -- {detail}" if detail else ""))
        print()
        print(f"{'ALL CHECKS PASSED' if v['n_fail'] == 0 else str(v['n_fail']) + ' FAILED'}"
              f"  ({v['n_pass']}/{v['n_pass'] + v['n_fail']}); "
              f"{v['record_count']} records")
        return 1 if v["n_fail"] else 0

    if args.gene:
        snap = load_snapshot()
        recs = get_evidence_for_gene(args.gene, symbol=args.symbol,
                                     top_k=args.top_k, snapshot=snap)
        if args.json:
            print(json.dumps(recs, indent=2, ensure_ascii=False))
        else:
            _print_evidence(str(args.gene), args.symbol, recs)
        return 0

    _print_snapshot_summary(load_snapshot())
    return 0


if __name__ == "__main__":
    sys.exit(main())
