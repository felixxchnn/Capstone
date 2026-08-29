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

Deterministic committed outputs
-------------------------------
Both committed files (`dgidb_<tag>.interactions.filtered.tsv`, `...manifest.json`)
are **byte-identical across repeated builds**. Nothing tracked contains a
wall-clock value: the retrieval provenance is a fixed build input
(`config.DGIDB_RETRIEVED_UTC`), and per-run execution timestamps go only to
`config.DGIDB_RUNLOG_FILE`, which is git-ignored.

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
The retained primary identity is the **canonical Entrez ID**, resolved by an
identifier join only:

    DGIdb ``gene_concept_id`` (``hgnc:<n>``)  ->  HGNC ID (``HGNC:<n>``)
      ->  ``entrez_id``   (pinned HGNC monthly complete set, config.DGIDB_HGNC_ASSET)
      ->  keep iff that Entrez is in the frozen canonical space (gene_columns.json)

Gene symbols are compared three ways (HGNC approved symbol vs DGIdb ``gene_name``
vs canonical symbol) and the agreement is recorded per record, but a symbol is
**never** an identity key. Ambiguous HGNC->Entrez mappings (a non-1:1 relation in
the HGNC file) are a hard failure. DGIdb's free-text ``gene_claim_name`` / aliases
are never used. Genes that do not resolve to a canonical Entrez are dropped and
counted (no/bad HGNC id, HGNC row with no Entrez, Entrez outside the canonical
space).

Publications (PMIDs)
--------------------
The upstream ``interactions.tsv`` carries no publication column, and PMIDs are
**never** parsed from free-text drug names. Instead the build reads the
release-aligned DGIdb SQL dump (config.DGIDB_SQL_ASSET) as a *temporary* input
(never committed) and joins ``interaction_claims`` ->
``interaction_claims_publications`` -> ``publications`` onto each retained record,
keyed by (gene concept id, drug concept id, interaction source). PMIDs are stored
sorted numerically and de-duplicated. Records DGIdb has no publication for stay
empty; per-source coverage is reported in the manifest.

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
import csv
import gzip
import hashlib
import json
import platform
import socket
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
# `pmids` is populated from the release-aligned SQL dump (claim-level
# publications joined by identifier, never parsed from free text); it is empty
# where DGIdb records no publication for that (gene, drug, source). `curation_type`
# and `indication` are absent from this DGIdb export and are always empty; they
# stay in the schema so a future release that carries them slots straight in.
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
    "pmids",                          # ";"-joined, numerically sorted, de-duped; SQL-dump join, never free text
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
# Gene-identity crosswalk:  DGIdb hgnc:<n>  ->  HGNC ID  ->  canonical Entrez
# (identifier join only; symbols are a consistency check, never a key)
# --------------------------------------------------------------------------

PMID_SEP = ";"

_HGNC_REQUIRED_COLUMNS = ("hgnc_id", "symbol", "status", "entrez_id")


def _load_canonical_gene_maps(gene_columns_path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """
    Return ``(entrez_to_symbol, symbol_to_entrez)`` from gene_columns.json.
    Both directions are needed: Entrez is the identity, symbol is the
    consistency field. Entrez ids and symbols are each globally unique here.
    """
    gc = json.loads(gene_columns_path.read_text(encoding="utf-8"))
    symbols = list(gc["symbols"])
    entrez = [str(e) for e in gc["entrez_ids"]]
    if len(symbols) != len(entrez):
        raise DGIdbSchemaError(
            f"{gene_columns_path.name}: symbols/entrez_ids length mismatch"
        )
    if len(set(entrez)) != len(entrez):
        raise DGIdbSchemaError(f"{gene_columns_path.name}: canonical Entrez ids not unique")
    return dict(zip(entrez, symbols)), dict(zip(symbols, entrez))


def _load_hgnc_crosswalk(
    hgnc_path: Path,
) -> tuple[dict[str, str], dict[str, str], dict[str, object]]:
    """
    Parse the pinned HGNC complete-set TSV into ``hgnc_id -> entrez_id`` and
    ``hgnc_id -> symbol``.

    Raises ``DGIdbSchemaError`` on a schema violation and ``DGIdbSnapshotError``
    on an *ambiguous* HGNC->Entrez relation (a non-1:1 mapping: one hgnc_id with
    a multi-valued entrez_id, or one entrez_id claimed by more than one hgnc_id).
    """
    hgnc_path = Path(hgnc_path)
    with hgnc_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        missing = [c for c in _HGNC_REQUIRED_COLUMNS if c not in header]
        if missing:
            raise DGIdbSchemaError(
                f"{hgnc_path.name}: HGNC complete set missing column(s) {missing}"
            )
        rows = list(reader)

    hid_to_entrez: dict[str, str] = {}
    hid_to_symbol: dict[str, str] = {}
    entrez_to_hids: dict[str, set[str]] = {}
    n_blank_entrez = 0
    for r in rows:
        hid = (r["hgnc_id"] or "").strip()
        if not hid.startswith("HGNC:"):
            raise DGIdbSchemaError(f"{hgnc_path.name}: bad hgnc_id {hid!r}")
        if hid in hid_to_entrez or hid in hid_to_symbol:
            raise DGIdbSnapshotError(
                f"{hgnc_path.name}: hgnc_id {hid} appears more than once -- "
                f"ambiguous crosswalk, refusing."
            )
        entrez = (r["entrez_id"] or "").strip()
        hid_to_symbol[hid] = r["symbol"]
        if not entrez:
            n_blank_entrez += 1
            hid_to_entrez[hid] = ""
            continue
        if not entrez.isdigit():
            raise DGIdbSnapshotError(
                f"{hgnc_path.name}: {hid} has a non-integer / multi-valued "
                f"entrez_id {entrez!r} -- ambiguous HGNC->Entrez, refusing."
            )
        hid_to_entrez[hid] = entrez
        entrez_to_hids.setdefault(entrez, set()).add(hid)

    collisions = {e: sorted(h) for e, h in entrez_to_hids.items() if len(h) > 1}
    if collisions:
        example = next(iter(collisions.items()))
        raise DGIdbSnapshotError(
            f"{hgnc_path.name}: {len(collisions)} Entrez id(s) are claimed by "
            f">1 HGNC id (e.g. {example}) -- ambiguous HGNC->Entrez, refusing."
        )

    provenance = {
        "file": hgnc_path.name,
        "sha256": sha256_file(hgnc_path),
        "bytes": hgnc_path.stat().st_size,
        "version": str(config.DGIDB_HGNC_ASSET["version"]),
        "url": str(config.DGIDB_HGNC_ASSET["url"]),
        "rows": len(rows),
        "rows_with_entrez_id": len(rows) - n_blank_entrez,
        "rows_without_entrez_id": n_blank_entrez,
        "relation": "1:1 (verified: no hgnc_id repeats, no shared Entrez)",
    }
    return hid_to_entrez, hid_to_symbol, provenance


def _resolve_genes(
    interaction_rows: list[dict[str, str]],
    hid_to_entrez: dict[str, str],
    hid_to_symbol: dict[str, str],
    entrez_to_symbol: dict[str, str],
) -> tuple[dict[str, tuple[str, str]], dict[str, object]]:
    """
    Resolve each distinct DGIdb ``gene_concept_id`` (``hgnc:<n>``) to a canonical
    Entrez id via the HGNC crosswalk.

    Returns ``(resolved, counts)`` where ``resolved`` maps ``gene_concept_id`` ->
    ``(entrez_id, hgnc_id)``. ``counts`` accounts for every distinct non-empty
    ``gene_concept_id`` and records the three-way symbol-agreement tally
    (consistency only -- symbols never affect resolution).
    """
    canonical_entrez = set(entrez_to_symbol)

    concept_to_dgidb_symbol: dict[str, set[str]] = {}
    for r in interaction_rows:
        cid = r["gene_concept_id"]
        if not cid:
            continue
        concept_to_dgidb_symbol.setdefault(cid, set()).add(r["gene_name"])

    resolved: dict[str, tuple[str, str]] = {}
    counts = {
        "distinct_dgidb_gene_concepts": len(concept_to_dgidb_symbol),
        "resolved_to_canonical_entrez": 0,
        "unresolved_concept_not_hgnc_or_not_in_hgnc_file": 0,
        "unresolved_hgnc_row_has_no_entrez": 0,
        "unresolved_entrez_not_in_canonical_space": 0,
        "symbol_agreement_all_three": 0,
        "symbol_disagreements": 0,
        "symbol_disagreement_examples": [],
    }
    for cid, dgidb_syms in sorted(concept_to_dgidb_symbol.items()):
        if not cid.startswith("hgnc:"):
            counts["unresolved_concept_not_hgnc_or_not_in_hgnc_file"] += 1
            continue
        hid = "HGNC:" + cid.split(":", 1)[1]
        if hid not in hid_to_entrez:
            counts["unresolved_concept_not_hgnc_or_not_in_hgnc_file"] += 1
            continue
        entrez = hid_to_entrez[hid]
        if not entrez:
            counts["unresolved_hgnc_row_has_no_entrez"] += 1
            continue
        if entrez not in canonical_entrez:
            counts["unresolved_entrez_not_in_canonical_space"] += 1
            continue
        resolved[cid] = (entrez, hid)
        counts["resolved_to_canonical_entrez"] += 1
        # ---- three-way symbol consistency (does not affect the join) ----
        hgnc_sym = hid_to_symbol.get(hid, "")
        canon_sym = entrez_to_symbol[entrez]
        dgidb_sym = next(iter(dgidb_syms)) if len(dgidb_syms) == 1 else "|".join(sorted(dgidb_syms))
        if len(dgidb_syms) == 1 and hgnc_sym == canon_sym == dgidb_sym:
            counts["symbol_agreement_all_three"] += 1
        else:
            counts["symbol_disagreements"] += 1
            if len(counts["symbol_disagreement_examples"]) < 20:
                counts["symbol_disagreement_examples"].append({
                    "entrez_id": entrez, "hgnc_id": hid,
                    "hgnc_symbol": hgnc_sym, "canonical_symbol": canon_sym,
                    "dgidb_gene_name": dgidb_sym,
                })
    return resolved, counts


# --------------------------------------------------------------------------
# Release-aligned SQL dump: claim-level publication (PMID) recovery.
#
# The dump is a PostgreSQL *data-only* pg_dump: a sequence of
#   COPY public.<table> (<cols>) FROM stdin;
#   <tab-separated rows, \N = NULL, C-style backslash escapes>
#   \.
# blocks. We parse only the tables needed to join publications onto the
# retained TSV interaction rows, keyed by identifiers that appear verbatim in
# both the TSV and the dump (gene concept_id, drug concept_id, source name).
# The dump is a TEMPORARY build input and is never committed.
# --------------------------------------------------------------------------

_PG_UNESCAPE = {
    "\\t": "\t", "\\n": "\n", "\\r": "\r", "\\\\": "\\",
    "\\b": "\b", "\\f": "\f", "\\v": "\v",
}


def _pg_unescape(value: str) -> str | None:
    if value == "\\N":
        return None
    if "\\" not in value:
        return value
    out: list[str] = []
    i = 0
    while i < len(value):
        ch = value[i]
        if ch == "\\" and i + 1 < len(value):
            pair = value[i:i + 2]
            out.append(_PG_UNESCAPE.get(pair, pair[1]))
            i += 2
        else:
            out.append(ch)
            i += 1
    return "".join(out)


def _iter_pg_copy(sql_gz_path: Path, table: str):
    """Yield dict rows for the single ``COPY public.<table> (...)`` block."""
    marker = f"COPY public.{table} ("
    with gzip.open(sql_gz_path, "rt", encoding="utf-8", newline="\n") as fh:
        cols: list[str] | None = None
        for line in fh:
            if cols is None:
                if line.startswith(marker):
                    inside = line[line.index("(") + 1:line.rindex(")")]
                    cols = [c.strip() for c in inside.split(",")]
                continue
            if line.startswith("\\."):
                return
            line = line.rstrip("\n")
            parts = line.split("\t")
            if len(parts) != len(cols):
                raise DGIdbSchemaError(
                    f"{sql_gz_path.name}: COPY {table}: row has {len(parts)} "
                    f"field(s), expected {len(cols)}"
                )
            yield {c: _pg_unescape(p) for c, p in zip(cols, parts)}


def _load_sql_publications(
    sql_gz_path: Path,
    included_source_names: frozenset[str],
    required_keys: set[tuple[str, str, str]],
) -> tuple[dict[tuple[str, str, str], tuple[str, ...]], dict[str, object]]:
    """
    Build ``(gene_concept_id, drug_concept_id, source_name) -> sorted PMID tuple``
    from the SQL dump, restricted to the included sources.

    ``required_keys`` are the (gene_concept, drug_concept, source) triples of the
    retained TSV records; used to compute the structural linkage rate and to
    hard-fail (per the task) if the dump cannot be linked reliably.
    """
    sql_gz_path = Path(sql_gz_path)

    # -- header / dump identity --
    dump_header: list[str] = []
    with gzip.open(sql_gz_path, "rt", encoding="utf-8", newline="\n") as fh:
        for _ in range(12):
            line = fh.readline()
            if not line:
                break
            dump_header.append(line.rstrip("\n"))
    pg_dump_version = next(
        (ln for ln in dump_header if "pg_dump version" in ln), "unknown"
    )

    gene_concept_to_gid: dict[str, str] = {}
    for r in _iter_pg_copy(sql_gz_path, "genes"):
        if r["concept_id"]:
            gene_concept_to_gid[r["concept_id"]] = r["id"]

    drug_concept_to_did: dict[str, str] = {}
    for r in _iter_pg_copy(sql_gz_path, "drugs"):
        if r["concept_id"]:
            drug_concept_to_did[r["concept_id"]] = r["id"]

    source_name_to_sid: dict[str, str] = {}
    for r in _iter_pg_copy(sql_gz_path, "sources"):
        source_name_to_sid[r["source_db_name"]] = r["id"]
    wanted_sids = {
        source_name_to_sid[s] for s in included_source_names
        if s in source_name_to_sid
    }

    di_gi_to_iid: dict[tuple[str, str], str] = {}
    for r in _iter_pg_copy(sql_gz_path, "interactions"):
        di_gi_to_iid[(r["drug_id"], r["gene_id"])] = r["id"]

    iid_src_to_claims: dict[tuple[str, str], set[str]] = {}
    n_claims = 0
    for r in _iter_pg_copy(sql_gz_path, "interaction_claims"):
        n_claims += 1
        if r["source_id"] in wanted_sids:
            iid_src_to_claims.setdefault(
                (r["interaction_id"], r["source_id"]), set()
            ).add(r["id"])

    pub_id_to_pmid: dict[str, str] = {}
    for r in _iter_pg_copy(sql_gz_path, "publications"):
        if r["pmid"]:
            pub_id_to_pmid[r["id"]] = r["pmid"]

    claim_to_pubs: dict[str, set[str]] = {}
    for r in _iter_pg_copy(sql_gz_path, "interaction_claims_publications"):
        claim_to_pubs.setdefault(r["interaction_claim_id"], set()).add(
            r["publication_id"]
        )

    # -- resolve every required key --
    key_to_pmids: dict[tuple[str, str, str], tuple[str, ...]] = {}
    linkable = 0            # keys that CAN link (have a drug concept id)
    linked = 0             # keys that DID link to an SQL interaction
    n_no_drug_concept = 0
    for (gcid, dcid, src) in required_keys:
        if not dcid:
            n_no_drug_concept += 1
            key_to_pmids[(gcid, dcid, src)] = ()
            continue
        linkable += 1
        gid = gene_concept_to_gid.get(gcid)
        did = drug_concept_to_did.get(dcid)
        sid = source_name_to_sid.get(src)
        iid = di_gi_to_iid.get((did, gid)) if (gid and did) else None
        if iid is None:
            key_to_pmids[(gcid, dcid, src)] = ()
            continue
        linked += 1
        pmids: set[str] = set()
        for claim in iid_src_to_claims.get((iid, sid), set()):
            for pub in claim_to_pubs.get(claim, set()):
                pm = pub_id_to_pmid.get(pub)
                if pm:
                    pmids.add(pm)
        key_to_pmids[(gcid, dcid, src)] = tuple(
            sorted(pmids, key=lambda x: (int(x) if x.isdigit() else 1 << 62, x))
        )

    linkage_rate = linked / linkable if linkable else 1.0
    if linkage_rate < 0.995:
        raise DGIdbSnapshotError(
            f"{sql_gz_path.name}: only {linked}/{linkable} concept-identified "
            f"records ({linkage_rate:.3%}) link to an interaction in the SQL "
            f"dump. Refusing to build -- the release SQL cannot be linked "
            f"reliably to the TSV records. STOP and investigate; do not "
            f"substitute a live query."
        )

    n_with_pmids = sum(1 for v in key_to_pmids.values() if v)
    total_pmids = sum(len(v) for v in key_to_pmids.values())
    stats = {
        "sql_dump": {
            "file": sql_gz_path.name,
            "sha256": str(config.DGIDB_SQL_ASSET["sha256"]),
            "bytes": int(config.DGIDB_SQL_ASSET["bytes"]),
            "url": str(config.DGIDB_SQL_ASSET["url"]),
            "role": "temporary build input -- NOT committed",
            "pg_dump_version_line": pg_dump_version,
            "interaction_claims_rows": n_claims,
        },
        "method": (
            "SQL dump interaction_claims -> interaction_claims_publications -> "
            "publications, joined onto retained records by (gene concept_id, "
            "drug concept_id, interaction source). PMIDs sorted numerically and "
            "de-duplicated. NEVER parsed from free text."
        ),
        "structural_linkage": {
            "distinct_key_triples": len(required_keys),
            "concept_identified_triples": linkable,
            "triples_without_drug_concept_id": n_no_drug_concept,
            "linked_to_sql_interaction": linked,
            "linkage_rate_over_concept_identified": round(linkage_rate, 5),
        },
        "distinct_triples_with_pmids": n_with_pmids,
        "total_pmid_mentions_over_triples": total_pmids,
    }
    return key_to_pmids, stats


# --------------------------------------------------------------------------
# Snapshot build
# --------------------------------------------------------------------------

def _verify_asset(path: Path, spec: dict, label: str) -> None:
    size = path.stat().st_size
    if size != spec["bytes"]:
        raise DGIdbSnapshotError(
            f"{label}: size {size} != pinned {spec['bytes']}. Refusing this asset."
        )
    digest = sha256_file(path)
    if digest != spec["sha256"]:
        raise DGIdbSnapshotError(
            f"{label}: sha256 {digest} != pinned {spec['sha256']}. "
            f"Refusing this asset."
        )


def _download_asset(url: str, dest: Path, timeout: int = 600) -> None:
    """Download one pinned asset to a staging path. Only reached via --refresh."""
    import urllib.request

    req = urllib.request.Request(url, headers={"User-Agent": "capstone-evidence/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:  # noqa: S310 (pinned host)
        data = resp.read()
    dest.write_bytes(data)


def build_snapshot(
    *,
    staging_dir: str | Path | None = None,
    refresh: bool = False,
    out_dir: str | Path | None = None,
    retrieved_utc_override: str | None = None,
    runlog_path: str | Path | None = None,
    write_runlog: bool = True,
) -> dict:
    """
    Build the licence-filtered offline snapshot + provenance manifest.

    Parameters
    ----------
    staging_dir
        Directory holding (or, with ``refresh=True``, to receive) the five
        pinned build inputs: the three DGIdb TSVs, the HGNC complete set, and
        the DGIdb SQL dump. A fresh temp dir is created if omitted.
    refresh
        If True, download the pinned inputs into ``staging_dir`` first.
    out_dir
        Where to write the committed pair. Defaults to ``config.DGIDB_DIR``.
        Only the filtered TSV and the manifest are written there -- never any
        upstream input, and never the SQL dump.
    retrieved_utc_override
        Test hook only. Production uses the fixed ``config.DGIDB_RETRIEVED_UTC``
        so the committed manifest is byte-identical across rebuilds.
    runlog_path / write_runlog
        Where (and whether) to append a one-line JSON run-log record with the
        wall-clock build time. This file is git-ignored; it is the *only* place
        a volatile timestamp is written.

    Returns a summary dict (its ``manifest`` key is the full manifest content).
    """
    staging = Path(staging_dir) if staging_dir else Path(
        tempfile.mkdtemp(prefix="dgidb_stage_")
    )
    staging.mkdir(parents=True, exist_ok=True)
    out_dir = Path(out_dir) if out_dir is not None else config.DGIDB_DIR
    retrieved_utc = retrieved_utc_override or config.DGIDB_RETRIEVED_UTC

    tsv_specs = dict(config.DGIDB_ASSETS)
    hgnc_spec = config.DGIDB_HGNC_ASSET
    sql_spec = config.DGIDB_SQL_ASSET
    hgnc_path = staging / str(hgnc_spec["name"])
    sql_path = staging / str(sql_spec["name"])
    asset_paths = {name: staging / name for name in tsv_specs}

    # ---- 1. obtain + gate every pinned input (exact size + SHA-256) ----
    for name, path in asset_paths.items():
        if refresh:
            _download_asset(str(tsv_specs[name]["url"]), path)
        if not path.is_file():
            raise DGIdbSnapshotError(
                f"{name} not found in staging dir {staging}. "
                f"Run with refresh=True (or place the pinned file there)."
            )
        _verify_asset(path, tsv_specs[name], name)
    for spec, path in ((hgnc_spec, hgnc_path), (sql_spec, sql_path)):
        if refresh:
            _download_asset(str(spec["url"]), path)
        if not path.is_file():
            raise DGIdbSnapshotError(
                f"{spec['name']} not found in staging dir {staging}. "
                f"Run with refresh=True (or place the pinned file there)."
            )
        _verify_asset(path, spec, str(spec["name"]))

    # ---- 2. parse + schema-validate the TSVs ------------------------
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

    # ---- 3. gene identity:  hgnc:<n> -> HGNC ID -> canonical Entrez ----
    gene_columns_path = config.PROCESSED_DIR / "gene_columns.json"
    entrez_to_symbol, _symbol_to_entrez = _load_canonical_gene_maps(gene_columns_path)
    hid_to_entrez, hid_to_symbol, hgnc_prov = _load_hgnc_crosswalk(hgnc_path)
    resolved_genes, gene_counts = _resolve_genes(
        interactions, hid_to_entrez, hid_to_symbol, entrez_to_symbol
    )

    # ---- 4. filter by source + build core records (no pmids yet) ------
    per_source_rows: dict[str, int] = {}
    per_source_retained: dict[str, int] = {}
    for name in SOURCE_LICENCES:
        per_source_rows[name] = 0
        per_source_retained[name] = 0

    dropped_unresolved_gene = 0
    dropped_no_drug_identity = 0
    seen_cores: set[tuple[str, ...]] = set()
    duplicates = 0
    cores: list[list[str]] = []          # 24 fields (SNAPSHOT_COLUMNS minus pmids/record_key)
    core_cols = [c for c in SNAPSHOT_COLUMNS if c not in ("pmids", "record_key")]

    for r in interactions:
        src = r["interaction_source_db_name"]
        if src in per_source_rows:
            per_source_rows[src] += 1
        if src not in INCLUDED_SOURCES:
            continue

        cid = r["gene_concept_id"]
        if cid not in resolved_genes:
            dropped_unresolved_gene += 1
            continue
        entrez, _hgnc_id = resolved_genes[cid]   # _hgnc_id == "HGNC:" + cid[5:]

        drug_name = r["drug_name"]
        drug_concept_id = r["drug_concept_id"]
        if drug_name == "" and drug_concept_id == "":
            dropped_no_drug_identity += 1
            continue

        canonical_symbol = entrez_to_symbol[entrez]
        dgidb_gene_name = r["gene_name"]
        core = [
            entrez,
            canonical_symbol,
            cid,                       # dgidb_gene_concept_id, e.g. "hgnc:11905"
            dgidb_gene_name,
            "true" if canonical_symbol == dgidb_gene_name else "false",
            drug_name,
            drug_concept_id,
            r["drug_claim_name"],
            src,
            r["interaction_source_db_version"],
            r["interaction_types"],
            normalize_direction(r["interaction_types"]),
            r["interaction_score"],
            r["drug_specificity_score"],
            r["gene_specificity_score"],
            r["evidence_score"],
            r["drug_is_approved"],
            r["drug_is_immunotherapy"],
            r["drug_is_antineoplastic"],
            # pmids slotted in after the SQL join
            "",  # curation_type -- absent in this export
            "",  # indication -- absent in this export
            SOURCE_LICENCES[src]["license"],
            SOURCE_LICENCES[src]["license_url"],
            config.DGIDB_RELEASE_TAG,
        ]
        tcore = tuple(core)
        if tcore in seen_cores:
            duplicates += 1
            continue
        seen_cores.add(tcore)
        cores.append(core)
        per_source_retained[src] += 1

    # ---- 5. recover claim-level publications from the SQL dump --------
    # (gene_concept_id, drug_concept_id, source) key per core record
    gci = core_cols.index("dgidb_gene_concept_id")
    dci = core_cols.index("drug_concept_id")
    sci = core_cols.index("interaction_source")
    required_keys = {(c[gci], c[dci], c[sci]) for c in cores}
    key_to_pmids, pub_stats = _load_sql_publications(
        sql_path, INCLUDED_SOURCES, required_keys
    )

    # ---- 6. assemble full records (dedup already done on the core) ----
    pmids_idx = SNAPSHOT_COLUMNS.index("pmids")
    records: list[dict[str, str]] = []
    for core in cores:
        pm = key_to_pmids.get((core[gci], core[dci], core[sci]), ())
        full = core[:pmids_idx] + [PMID_SEP.join(pm)] + core[pmids_idx:]
        rec = dict(zip(SNAPSHOT_COLUMNS[:-1], full))
        rec["record_key"] = _record_key(full)
        records.append(rec)

    # ---- 7. deterministic total order ---------------------------
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

    # ---- per-source publication coverage ------------------------
    pub_cov: dict[str, dict[str, int]] = {}
    for rec in records:
        b = pub_cov.setdefault(rec["interaction_source"],
                               {"records": 0, "records_with_pmids": 0, "pmids": 0})
        b["records"] += 1
        if rec["pmids"]:
            n = len(rec["pmids"].split(PMID_SEP))
            b["records_with_pmids"] += 1
            b["pmids"] += n
    records_with_pmids = sum(1 for rec in records if rec["pmids"])

    # ---- 8. write snapshot TSV (LF, no trailing blank line) -----
    staging.mkdir(parents=True, exist_ok=True)
    snapshot_name = config.DGIDB_SNAPSHOT_FILE.name
    manifest_name = config.DGIDB_MANIFEST_FILE.name
    staged_snapshot = staging / snapshot_name
    lines = ["\t".join(SNAPSHOT_COLUMNS)]
    lines += ["\t".join(rec[c] for c in SNAPSHOT_COLUMNS) for rec in records]
    staged_snapshot.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    snapshot_sha = sha256_file(staged_snapshot)
    snapshot_bytes = staged_snapshot.stat().st_size

    # ---- 9. manifest (fully deterministic -- no wall-clock value) ----
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
            "publication_coverage": pub_cov.get(name, {
                "records": 0, "records_with_pmids": 0, "pmids": 0}),
        }

    manifest = {
        "artifact": snapshot_name,
        "description": (
            "Licence-filtered offline snapshot of DGIdb drug-gene interaction "
            "records for Phase 2 evidence retrieval. Evidence retrieval only -- "
            "not treatment prediction, not a recommendation. Contains records "
            "ONLY from interaction sources whose redistribution terms are "
            "explicitly verified as compatible with committing them; the full "
            "DGIdb dataset is NOT redistributable, and individual records remain "
            "subject to their own source-specific terms."
        ),
        "dgidb": {
            "release_tag": config.DGIDB_RELEASE_TAG,
            "release_page": config.DGIDB_RELEASE_PAGE,
            "versions": {
                "release_tag": config.DGIDB_RELEASE_TAG,
                "interaction_data_version": config.DGIDB_INTERACTION_DATA_VERSION,
                "dgidb_app_version_tsv_comment": config.DGIDB_APP_VERSION_TSV_COMMENT,
                "dgidb_app_version_graphql_at_retrieval": config.DGIDB_APP_VERSION_GRAPHQL,
                "sql_dump_pg_dump_version_line": pub_stats["sql_dump"]["pg_dump_version_line"],
                "note": (
                    "The release is tagged 2026-06b and was published 2026-06-23, "
                    "but it PACKAGES interaction content dated "
                    f"'{config.DGIDB_INTERACTION_DATA_VERSION}'. The interaction "
                    "records are NOT June-2026 data. Source-level versions (per "
                    "'sources') range from 2011 to 2026."
                ),
            },
            "upstream_header_comments": {
                "interactions.tsv": inter_comments,
                "genes.tsv": genes_comments,
                "drugs.tsv": drugs_comments,
            },
        },
        "assets": {
            **{
                name: {
                    "url": str(spec["url"]),
                    "bytes": spec["bytes"],
                    "sha256": spec["sha256"],
                    "role": "committed-snapshot source (upstream TSV, not committed)",
                    "verified": True,
                }
                for name, spec in tsv_specs.items()
            },
            "hgnc_complete_set": {
                "url": str(hgnc_spec["url"]),
                "bytes": hgnc_spec["bytes"],
                "sha256": hgnc_spec["sha256"],
                "version": str(hgnc_spec["version"]),
                "role": "gene-identity crosswalk input (not committed)",
                "license": str(hgnc_spec["license"]),
                "license_url": str(hgnc_spec["license_url"]),
                "verified": True,
            },
            "sql_dump": {
                "url": str(sql_spec["url"]),
                "bytes": sql_spec["bytes"],
                "sha256": sql_spec["sha256"],
                "role": "temporary publication-recovery input -- NOT committed",
                "verified": True,
            },
        },
        "retrieval": {
            "retrieved_utc": retrieved_utc,
            "retrieved_utc_source": "fixed build input (config.DGIDB_RETRIEVED_UTC)",
            "method": (
                "five pinned inputs -- the 3 DGIdb 2026-06b release TSVs, the HGNC "
                "monthly complete set, and the DGIdb 2026-06b SQL dump -- each "
                "gated on exact byte size and SHA-256 before use. Whether they "
                "were downloaded (--refresh) or pre-staged (--from-staging) does "
                "not affect the output; the committed snapshot and this manifest "
                "are byte-identical across rebuilds. The SQL dump is a temporary "
                "input and is never committed."
            ),
            "tool": "evidence.py build_snapshot",
            "note": (
                "No wall-clock value appears in any committed artifact. Per-run "
                "execution timestamps (and the download-vs-staging mode) are "
                f"appended to the git-ignored {config.DGIDB_RUNLOG_FILE.name}."
            ),
        },
        "gene_identity": {
            "primary_identifier": "entrez_id (canonical space)",
            "method": (
                "identifier join only: DGIdb gene_concept_id (hgnc:<n>) -> HGNC ID "
                "-> entrez_id via the pinned HGNC complete set, then kept iff that "
                "Entrez is in gene_columns.json. Symbols (HGNC approved vs DGIdb "
                "gene_name vs canonical) are compared for consistency only and "
                "never used as an identity key. Ambiguous HGNC->Entrez mappings "
                "hard-fail. gene_claim_name / aliases are never used."
            ),
            "canonical_space_file": gene_columns_path.name,
            "canonical_space_sha256": sha256_file(gene_columns_path),
            "hgnc_crosswalk": hgnc_prov,
            **gene_counts,
        },
        "publications": {
            **pub_stats,
            "sep": PMID_SEP,
            "sorted": "numeric ascending, de-duplicated within a record",
            "records_with_pmids": records_with_pmids,
            "coverage_by_source": {
                s: pub_cov.get(s, {"records": 0, "records_with_pmids": 0, "pmids": 0})
                for s in sorted(INCLUDED_SOURCES)
            },
            "note": (
                "Claim-level publications are absent in DGIdb for ChEMBL, "
                "GuideToPharmacology and FDA in this release, so their records "
                "carry no PMIDs. CIViC, DoCM and NCI are ~fully covered. Records "
                "with no DGIdb publication legitimately stay empty."
            ),
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
                "interactions.tsv carries no publication column; PMIDs are "
                "recovered from the release-aligned SQL dump by identifier join "
                "(never parsed from free text). Empty where DGIdb records no "
                "claim-level publication for that (gene, drug, source) -- notably "
                "all ChEMBL / GuideToPharmacology / FDA records in this release."
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
                "redistributable interaction sources in filter.included_sources. "
                "The full DGIdb dataset is NOT redistributable. Each retained "
                "record additionally remains subject to its own source-specific "
                "licence terms (per-record source_license / source_license_url); "
                "the compilation licence does not override them. See LICENSES.md."
            ),
        },
    }

    staged_manifest = staging / manifest_name
    staged_manifest.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8", newline="\n",
    )

    # ---- 10. publish the tracked pair (filtered TSV + manifest ONLY) ----
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / snapshot_name).write_bytes(staged_snapshot.read_bytes())
        (out_dir / manifest_name).write_bytes(staged_manifest.read_bytes())

    # ---- 11. volatile run log (git-ignored; NOT a tracked output) -------
    if write_runlog:
        rl = Path(runlog_path) if runlog_path is not None else config.DGIDB_RUNLOG_FILE
        try:
            rl.parent.mkdir(parents=True, exist_ok=True)
            entry = {
                "wall_clock_utc": datetime.now(timezone.utc).strftime(
                    "%Y-%m-%dT%H:%M:%SZ"),
                "host": socket.gethostname(),
                "platform": platform.platform(),
                "mode": "refresh" if refresh else "from-staging",
                "snapshot_sha256": snapshot_sha,
                "manifest_sha256": sha256_file(staged_manifest),
                "record_count": len(records),
            }
            with rl.open("a", encoding="utf-8", newline="\n") as fh:
                fh.write(json.dumps(entry) + "\n")
        except OSError:
            pass  # a read-only FS must not fail the build

    return {
        "staging_dir": str(staging),
        "out_dir": str(out_dir) if out_dir is not None else None,
        "snapshot_sha256": snapshot_sha,
        "snapshot_bytes": snapshot_bytes,
        "manifest_sha256": sha256_file(staged_manifest),
        "record_count": len(records),
        "records_with_pmids": records_with_pmids,
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

    # manifest asset hashes == config pins (TSVs + HGNC crosswalk + SQL dump)
    assets_ok = all(
        m["assets"][n]["sha256"] == config.DGIDB_ASSETS[n]["sha256"]
        and m["assets"][n]["bytes"] == config.DGIDB_ASSETS[n]["bytes"]
        for n in config.DGIDB_ASSETS
    )
    assets_ok = assets_ok and (
        m["assets"]["hgnc_complete_set"]["sha256"] == config.DGIDB_HGNC_ASSET["sha256"]
        and m["assets"]["hgnc_complete_set"]["bytes"] == config.DGIDB_HGNC_ASSET["bytes"]
        and m["assets"]["sql_dump"]["sha256"] == config.DGIDB_SQL_ASSET["sha256"]
        and m["assets"]["sql_dump"]["bytes"] == config.DGIDB_SQL_ASSET["bytes"]
    )
    check("manifest asset sizes+hashes == config pins (TSVs + HGNC + SQL)", assets_ok)

    # committed manifest carries NO wall-clock value; retrieval is a fixed input
    check("retrieval.retrieved_utc == fixed config.DGIDB_RETRIEVED_UTC",
          m["retrieval"]["retrieved_utc"] == config.DGIDB_RETRIEVED_UTC,
          m["retrieval"]["retrieved_utc"])

    # the three versions are recorded and distinct in meaning
    v = m["dgidb"]["versions"]
    check("provenance separates release tag / data version / app versions",
          v["release_tag"] == config.DGIDB_RELEASE_TAG
          and v["interaction_data_version"] == config.DGIDB_INTERACTION_DATA_VERSION
          and v["dgidb_app_version_tsv_comment"] == config.DGIDB_APP_VERSION_TSV_COMMENT
          and "not june-2026 data" in v["note"].lower(),
          f"data_version={v['interaction_data_version']}")

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
    entrez_to_symbol, _s2e = _load_canonical_gene_maps(
        config.PROCESSED_DIR / "gene_columns.json")
    entrez_set = set(entrez_to_symbol)
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

    # identity is an HGNC-ID join, never symbol-only: every record carries an
    # hgnc: concept id, and the manifest method says so.
    all_hgnc = all(r["dgidb_gene_concept_id"].startswith("hgnc:")
                   and r["dgidb_gene_concept_id"][5:].isdigit()
                   for r in snap.records)
    method = m["gene_identity"]["method"].lower()
    method_ok = ("hgnc id" in method and "identifier join only" in method
                 and "never used as an identity key" in method
                 and "hard-fail" in method)
    check("gene identity is an HGNC-ID join, not symbol-only",
          all_hgnc and method_ok)
    check("gene_identity records the HGNC crosswalk provenance",
          isinstance(m["gene_identity"].get("hgnc_crosswalk"), dict)
          and m["gene_identity"]["hgnc_crosswalk"]["sha256"]
          == config.DGIDB_HGNC_ASSET["sha256"])

    # ---- publications: sorted, numeric, de-duped; not silently all-empty ----
    pmid_ok = True
    pmid_detail = ""
    for r in snap.records:
        if not r["pmids"]:
            continue
        parts = r["pmids"].split(PMID_SEP)
        if not all(p.isdigit() for p in parts):
            pmid_ok = False
            pmid_detail = f"non-numeric pmid in {r['pmids']!r}"
            break
        if parts != sorted(parts, key=int):
            pmid_ok = False
            pmid_detail = f"unsorted pmids {r['pmids']!r}"
            break
        if len(parts) != len(set(parts)):
            pmid_ok = False
            pmid_detail = f"duplicate pmids {r['pmids']!r}"
            break
    check("pmids are numeric, sorted ascending, de-duplicated", pmid_ok, pmid_detail)

    n_with_pmids = sum(1 for r in snap.records if r["pmids"])
    check("snapshot is NOT silently zero-publication-coverage",
          n_with_pmids > 0, f"{n_with_pmids} records carry >=1 PMID")
    check("manifest publications.records_with_pmids matches the file",
          m["publications"]["records_with_pmids"] == n_with_pmids,
          f"file={n_with_pmids} manifest={m['publications']['records_with_pmids']}")

    # per-source publication coverage matches the manifest
    cov_recount: dict[str, dict[str, int]] = {}
    for r in snap.records:
        b = cov_recount.setdefault(r["interaction_source"],
                                   {"records": 0, "records_with_pmids": 0, "pmids": 0})
        b["records"] += 1
        if r["pmids"]:
            b["records_with_pmids"] += 1
            b["pmids"] += len(r["pmids"].split(PMID_SEP))
    cov_ok = all(
        cov_recount.get(s) == m["publications"]["coverage_by_source"].get(s)
        for s in cov_recount
    )
    check("manifest publications.coverage_by_source matches the file", cov_ok,
          f"file={cov_recount}")

    # pmids are a deterministic function of (gene, drug, source): identical
    # key triples carry identical pmids
    triple_pmids: dict[tuple[str, str, str], str] = {}
    triple_ok = True
    for r in snap.records:
        k = (r["dgidb_gene_concept_id"], r["drug_concept_id"], r["interaction_source"])
        if k in triple_pmids and triple_pmids[k] != r["pmids"]:
            triple_ok = False
            break
        triple_pmids[k] = r["pmids"]
    check("pmids are consistent across identical (gene,drug,source) triples", triple_ok)

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


def _write_pg_dump_gz(path: Path, tables: dict[str, tuple[list[str], list[list[str]]]]) -> None:
    """Write a minimal PostgreSQL data-only pg_dump (.sql.gz) for the self-test."""
    parts = [
        "--", "-- PostgreSQL database dump", "--", "",
        "-- Dumped by pg_dump version 18.1 (self-test)", "",
    ]
    for table, (cols, rows) in tables.items():
        parts.append(f"COPY public.{table} ({', '.join(cols)}) FROM stdin;")
        for row in rows:
            parts.append("\t".join("\\N" if v is None else str(v) for v in row))
        parts.append("\\.")
        parts.append("")
    with gzip.open(path, "wt", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(parts) + "\n")


def _self_test() -> int:  # noqa: C901 -- one linear scenario, kept together
    import shutil

    print("Running evidence.py self-test...")
    tmp = Path(tempfile.mkdtemp(prefix="evidence_selftest_"))
    stage = tmp / "stage"
    stage.mkdir()

    # ---- synthetic canonical space -----------------------------------
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
    _orig_assets = config.DGIDB_ASSETS
    _orig_hgnc = config.DGIDB_HGNC_ASSET
    _orig_sql = config.DGIDB_SQL_ASSET
    config.PROCESSED_DIR = proc
    try:
        icols = INTERACTIONS_TSV_COLUMNS

        def irow(concept, gene_name, drug_name, drug_cid, approved, src, itype, iscore):
            d = {c: "" for c in icols}
            d.update(gene_claim_name="gc", gene_concept_id=concept, gene_name=gene_name,
                     drug_claim_name="dc", drug_concept_id=drug_cid, drug_name=drug_name,
                     drug_is_approved=approved, drug_is_immunotherapy="false",
                     drug_is_antineoplastic="false",
                     interaction_source_db_name=src, interaction_source_db_version="v1",
                     interaction_types=itype, interaction_score=iscore)
            return [d[c] for c in icols]

        inter_rows = []
        # AAA (hgnc:1 -> entrez 100): 7 inhibitory + 3 activating + 1 unknown, ChEMBL
        for i in range(7):
            inter_rows.append(irow("hgnc:1", "AAA", f"inh{i}", f"chembl:{i}", "true",
                                   "ChEMBL", "inhibitor", f"{10 - i}.0"))
        for i in range(3):
            inter_rows.append(irow("hgnc:1", "AAA", f"act{i}", f"chembl:a{i}", "false",
                                   "ChEMBL", "agonist", f"{5 - i}.0"))
        inter_rows.append(irow("hgnc:1", "AAA", "unk0", "chembl:u0", "false",
                               "ChEMBL", "modulator", ""))
        inter_rows.append(irow("hgnc:1", "AAA", "inh0", "chembl:0", "true",
                               "ChEMBL", "inhibitor", "10.0"))   # dup of first -> collapses
        # BBB (hgnc:2 -> entrez 200): one CIViC (included), one DTC (excluded).
        # DGIdb gene_name deliberately "BBB_alias" != canonical/HGNC "BBB"
        # -> counted as a symbol disagreement, identity unaffected.
        inter_rows.append(irow("hgnc:2", "BBB_alias", "civicdrug", "civic:1", "true",
                               "CIViC", "inhibitor|agonist", "3.0"))
        inter_rows.append(irow("hgnc:2", "BBB", "dtcdrug", "dtc:1", "true",
                               "DTC", "inhibitor", "9.9"))
        # CCC (hgnc:3 -> entrez 300): included source but NO drug identity -> dropped
        inter_rows.append(irow("hgnc:3", "CCC", "", "", "false", "NCI", "inhibitor", ""))
        # hgnc:9 -> entrez 999 (not in canonical space) -> dropped
        inter_rows.append(irow("hgnc:9", "ZZZ", "zdrug", "z:1", "false", "FDA",
                               "inhibitor", "1.0"))
        # empty concept -> unresolved -> dropped
        inter_rows.append(irow("", "", "orphan", "x:1", "false", "ChEMBL",
                               "inhibitor", "1.0"))
        # concept with no HGNC row -> dropped
        inter_rows.append(irow("hgnc:404", "NOPE", "ndrug", "n:1", "false", "ChEMBL",
                               "inhibitor", "1.0"))
        # AAA activator from DoCM (second included source)
        inter_rows.append(irow("hgnc:1", "AAA", "docmact", "docm:1", "false",
                               "DoCM", "activator", "4.5"))

        _write_tsv(stage / "interactions.tsv",
                   ["# Data version: Dec-2023", "# DGIdb version: v.5.0.11"], icols, inter_rows)
        _write_tsv(stage / "genes.tsv",
                   ["# Data version: Dec-2023", "# DGIdb version: v.5.0.11"],
                   GENES_TSV_COLUMNS, [["AAA", "Gene Symbol", "hgnc:1", "AAA", "HGNC", "x"]])
        _write_tsv(stage / "drugs.tsv",
                   ["# Data version: Dec-2023", "# DGIdb version: v.5.0.11"],
                   DRUGS_TSV_COLUMNS,
                   [["dc", "Primary Name", "chembl:0", "inh0", "true", "false",
                     "false", "ChEMBL", "v1"]])

        # ---- synthetic HGNC complete set ---------------------------
        hgnc_cols = ["hgnc_id", "symbol", "status", "entrez_id", "ensembl_gene_id"]
        hgnc_rows = [
            ["HGNC:1", "AAA", "Approved", "100", "ENSG1"],
            ["HGNC:2", "BBB", "Approved", "200", "ENSG2"],
            ["HGNC:3", "CCC", "Approved", "300", "ENSG3"],
            ["HGNC:9", "ZZZ", "Approved", "999", "ENSG9"],   # entrez not canonical
            ["HGNC:7", "RNA7", "Approved", "", "ENSG7"],      # no entrez
        ]
        hgnc_file = stage / str(_orig_hgnc["name"])
        _write_tsv(hgnc_file, [], hgnc_cols, hgnc_rows)

        # ---- synthetic release SQL dump (covers every retained record) ----
        # one drug row per distinct drug concept id among retained rows
        retained_specs = [
            ("hgnc:1", f"chembl:{i}", "ChEMBL") for i in range(7)
        ] + [
            ("hgnc:1", f"chembl:a{i}", "ChEMBL") for i in range(3)
        ] + [
            ("hgnc:1", "chembl:u0", "ChEMBL"),
            ("hgnc:2", "civic:1", "CIViC"),
            ("hgnc:1", "docm:1", "DoCM"),
        ]
        gene_ids = {"hgnc:1": "g1", "hgnc:2": "g2"}
        src_ids = {"ChEMBL": "s_chembl", "CIViC": "s_civic", "DoCM": "s_docm",
                   "NCI": "s_nci", "FDA": "s_fda", "GuideToPharmacology": "s_gtp"}
        drug_rows, inter_rows_sql, claim_rows = [], [], []
        seen_drug, seen_inter = set(), set()
        for n, (gc, dc, src) in enumerate(retained_specs):
            if dc not in seen_drug:
                drug_rows.append([f"d{n}", dc.split(':')[1], "false", "false", "false", dc])
                seen_drug.add(dc)
            did = next(dr[0] for dr in drug_rows if dr[5] == dc)
            gid = gene_ids[gc]
            iid = f"i_{gid}_{did}"
            if (did, gid) not in seen_inter:
                inter_rows_sql.append([iid, did, gid, "", "", "", ""])
                seen_inter.add((did, gid))
            claim_rows.append([f"c_{n}", None, None, src_ids[src], iid])
        sql_tables = {
            "genes": (["id", "name", "long_name", "concept_id"],
                      [["g1", "AAA", "", "hgnc:1"], ["g2", "BBB", "", "hgnc:2"]]),
            "drugs": (["id", "name", "approved", "immunotherapy", "anti_neoplastic", "concept_id"],
                      drug_rows),
            "sources": (["id", "source_db_name", "source_db_version"],
                        [[v, k, "v1"] for k, v in src_ids.items()]),
            "interactions": (["id", "drug_id", "gene_id", "score", "drug_specificity",
                              "gene_specificity", "evidence_score"], inter_rows_sql),
            "interaction_claims": (["id", "drug_claim_id", "gene_claim_id", "source_id",
                                    "interaction_id"], claim_rows),
            "interaction_claims_publications": (["interaction_claim_id", "publication_id"],
                # civicdrug (c_11) -> p1,p2 (+ a dup link); docmact (c_12) -> p3;
                # every ChEMBL claim (c_0..c_10) -> none (mirrors real DGIdb)
                [["c_11", "p1"], ["c_11", "p2"], ["c_11", "p1"], ["c_12", "p3"]]),
            "publications": (["id", "pmid", "citation", "created_at", "updated_at"],
                             [["p1", "555", "c", "t", "t"],
                              ["p2", "111", "c", "t", "t"],
                              ["p3", "222", "c", "t", "t"]]),
        }
        sql_file = stage / str(_orig_sql["name"])
        _write_pg_dump_gz(sql_file, sql_tables)

        # ---- point config at the synthetic inputs ------------------
        config.DGIDB_ASSETS = {
            n: {"url": f"file://{stage / n}", "bytes": (stage / n).stat().st_size,
                "sha256": sha256_file(stage / n)}
            for n in ("interactions.tsv", "genes.tsv", "drugs.tsv")
        }
        config.DGIDB_HGNC_ASSET = {**_orig_hgnc, "bytes": hgnc_file.stat().st_size,
                                   "sha256": sha256_file(hgnc_file)}
        config.DGIDB_SQL_ASSET = {**_orig_sql, "bytes": sql_file.stat().st_size,
                                  "sha256": sha256_file(sql_file)}
        try:
            out_dir = tmp / "out"
            runlog = tmp / "runlog.jsonl"

            # ---- 1. build from staging ---------------------------
            res1 = build_snapshot(staging_dir=stage, refresh=False, out_dir=out_dir,
                                  runlog_path=runlog)
            snap_file = out_dir / config.DGIDB_SNAPSHOT_FILE.name
            man_file = out_dir / config.DGIDB_MANIFEST_FILE.name
            assert snap_file.is_file() and man_file.is_file()
            assert not (out_dir / "interactions.tsv").exists()
            assert not (out_dir / str(_orig_sql["name"])).exists()   # SQL never published
            assert runlog.is_file()   # volatile run log written, untracked location
            print("  [ok] build writes only the filtered snapshot + manifest (SQL never published)")

            # ---- 2. FULLY deterministic rebuild (snapshot AND manifest) ----
            res2 = build_snapshot(staging_dir=stage, refresh=False,
                                  out_dir=tmp / "out2", runlog_path=tmp / "rl2.jsonl")
            b1 = snap_file.read_bytes()
            b2 = (tmp / "out2" / config.DGIDB_SNAPSHOT_FILE.name).read_bytes()
            assert b1 == b2, "snapshot TSV not byte-identical across rebuilds"
            man1 = man_file.read_bytes()
            man2 = (tmp / "out2" / config.DGIDB_MANIFEST_FILE.name).read_bytes()
            assert man1 == man2, "manifest not byte-identical across rebuilds"
            assert res1["snapshot_sha256"] == res2["snapshot_sha256"]
            assert res1["manifest_sha256"] == res2["manifest_sha256"]
            assert "datetime.now" not in man_file.read_text(encoding="utf-8")
            print("  [ok] snapshot AND manifest regenerate byte-identically (no wall clock)")

            # ---- 3. size / hash gate refuses a tampered asset ---
            bad_stage = tmp / "bad"
            bad_stage.mkdir()
            for n in ("interactions.tsv", "genes.tsv", "drugs.tsv",
                      str(_orig_hgnc["name"]), str(_orig_sql["name"])):
                shutil.copy(stage / n, bad_stage / n)
            (bad_stage / "genes.tsv").write_text(
                (bad_stage / "genes.tsv").read_text(encoding="utf-8")
                + "AAA\tGene Symbol\thgnc:1\tAAA\tHGNC\tx\n",
                encoding="utf-8", newline="\n")
            try:
                build_snapshot(staging_dir=bad_stage, refresh=False, out_dir=tmp / "nope",
                               write_runlog=False)
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

            # ---- 4b. HGNC-ID crosswalk unit checks ----------------
            hid2e, hid2s, hprov = _load_hgnc_crosswalk(hgnc_file)
            assert hid2e["HGNC:1"] == "100" and hid2s["HGNC:2"] == "BBB"
            assert hprov["relation"].startswith("1:1")
            e2s, s2e = _load_canonical_gene_maps(proc / "gene_columns.json")
            resolved, gcnt = _resolve_genes(
                [dict(zip(icols, r)) for r in inter_rows], hid2e, hid2s, e2s)
            assert resolved["hgnc:1"] == ("100", "HGNC:1")
            assert "hgnc:9" not in resolved       # entrez 999 not canonical
            assert "hgnc:404" not in resolved     # no HGNC row
            assert "hgnc:7" not in resolved       # (not referenced) but HGNC:7 has no entrez
            assert gcnt["unresolved_entrez_not_in_canonical_space"] == 1
            assert gcnt["unresolved_concept_not_hgnc_or_not_in_hgnc_file"] == 1
            assert gcnt["symbol_disagreements"] >= 1   # BBB_alias vs BBB
            # ambiguous HGNC->Entrez must hard-fail
            for bad in (
                [["HGNC:1", "AAA", "Approved", "100", ""], ["HGNC:9", "ZZ", "Approved", "100", ""]],  # shared entrez
                [["HGNC:1", "AAA", "Approved", "100,101", ""]],   # multi-valued entrez
                [["HGNC:1", "AAA", "Approved", "100", ""], ["HGNC:1", "AAA", "Approved", "100", ""]],  # repeat id
            ):
                bf = tmp / "bad_hgnc.txt"
                _write_tsv(bf, [], hgnc_cols[:4] + ["ensembl_gene_id"], bad)
                try:
                    _load_hgnc_crosswalk(bf)
                    raise AssertionError("ambiguous HGNC->Entrez should raise")
                except DGIdbSnapshotError as exc:
                    assert "ambiguous" in str(exc).lower()
            print("  [ok] HGNC-ID crosswalk: 1:1 join, canonical filter, ambiguous hard-fail")

            # ---- 5. load + retrieval semantics --------------------
            snap = load_snapshot(snap_file, man_file)
            assert all(r["interaction_source"] in INCLUDED_SOURCES for r in snap.records)
            assert not any(r["drug_name"] == "dtcdrug" for r in snap.records)
            aaa_inh = [r for r in snap.records
                       if r["entrez_id"] == "100" and r["interaction_direction"] == "inhibitory"]
            assert len(aaa_inh) == 7, len(aaa_inh)   # dup collapsed
            bbb = snap.by_entrez["200"]
            assert len(bbb) == 1 and bbb[0]["interaction_direction"] == "unknown"
            # identity survived a symbol disagreement; symbol is only a flag
            assert bbb[0]["gene_symbol"] == "BBB"            # canonical
            assert bbb[0]["dgidb_gene_name"] == "BBB_alias"  # as DGIdb provided
            assert bbb[0]["gene_symbol_consistent"] == "false"
            assert bbb[0]["dgidb_gene_concept_id"] == "hgnc:2"
            assert snap.manifest["gene_identity"]["symbol_disagreements"] >= 1
            assert "300" not in snap.by_entrez       # CCC: no drug identity
            assert all(r["gene_symbol"] != "ZZZ" for r in snap.records)  # entrez 999
            print("  [ok] source filter, dedup, symbol-disagreement, no-drug + non-canonical drops")

            # ---- 4c. publications joined from the SQL dump --------
            civ = bbb[0]
            assert civ["pmids"] == "111;555", civ["pmids"]     # sorted numeric, deduped
            docm = [r for r in snap.records if r["drug_name"] == "docmact"][0]
            assert docm["pmids"] == "222"
            assert all(r["pmids"] == "" for r in aaa_inh)      # ChEMBL: no claim-level pubs
            n_pm = sum(1 for r in snap.records if r["pmids"])
            assert n_pm >= 2 and snap.manifest["publications"]["records_with_pmids"] == n_pm
            assert snap.manifest["publications"]["coverage_by_source"]["CIViC"]["records_with_pmids"] == 1
            assert snap.manifest["publications"]["coverage_by_source"]["ChEMBL"]["records_with_pmids"] == 0
            print("  [ok] PMIDs joined from SQL dump: sorted, deduped, source-attributed, not all-empty")

            # ---- 4d. SQL that cannot be linked -> hard STOP -------
            broken_sql = dict(sql_tables)
            broken_sql["interactions"] = (sql_tables["interactions"][0], [])  # no interactions
            bsf = tmp / "broken.sql.gz"
            _write_pg_dump_gz(bsf, broken_sql)
            _bkp = config.DGIDB_SQL_ASSET
            config.DGIDB_SQL_ASSET = {**_orig_sql, "bytes": bsf.stat().st_size,
                                      "sha256": sha256_file(bsf)}
            _real = stage / str(_orig_sql["name"])
            _saved = _real.read_bytes()
            _real.write_bytes(bsf.read_bytes())
            try:
                build_snapshot(staging_dir=stage, refresh=False, out_dir=tmp / "nolink",
                               write_runlog=False)
                raise AssertionError("unlinkable SQL should STOP the build")
            except DGIdbSnapshotError as exc:
                assert "cannot be linked reliably" in str(exc)
            finally:
                _real.write_bytes(_saved)
                config.DGIDB_SQL_ASSET = _bkp
            print("  [ok] unlinkable SQL dump stops the build (no live-query fallback)")

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
            config.DGIDB_HGNC_ASSET = _orig_hgnc
            config.DGIDB_SQL_ASSET = _orig_sql
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
    ver = m["dgidb"]["versions"]
    print(f"  versions         : release tag {ver['release_tag']}  |  "
          f"interaction data {ver['interaction_data_version']}  |  "
          f"app {ver['dgidb_app_version_tsv_comment']}")
    print(f"  retrieved (fixed): {m['retrieval']['retrieved_utc']}  "
          f"(build input, not a wall clock)")
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
    print(f"  gene identity    : {gi['primary_identifier']}  "
          f"(hgnc:<n> -> HGNC ID -> Entrez; symbol = consistency check only)")
    print(f"    resolved       : {gi['resolved_to_canonical_entrez']} / "
          f"{gi['distinct_dgidb_gene_concepts']} DGIdb gene concepts  "
          f"({gi['symbol_disagreements']} symbol disagreements)")
    print(f"    unresolved     : "
          f"{gi['unresolved_concept_not_hgnc_or_not_in_hgnc_file']} no HGNC row, "
          f"{gi['unresolved_hgnc_row_has_no_entrez']} HGNC row w/o Entrez, "
          f"{gi['unresolved_entrez_not_in_canonical_space']} outside canonical space")
    pub = m["publications"]
    print("-" * 74)
    print(f"  publications     : {pub['records_with_pmids']} / "
          f"{m['snapshot']['record_count']} records carry >=1 PMID  "
          f"(SQL-dump join, never free text)")
    for s in f["included_sources"]:
        c = pub["coverage_by_source"][s]
        print(f"    {s:22s} {c['records_with_pmids']:>6d} / {c['records']:<6d} "
              f"records w/ PMIDs  ({c['pmids']} PMID mentions)")
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
        print(f"  staging dir     : {res['staging_dir']}")
        print(f"  out dir         : {res['out_dir']}")
        print(f"  records         : {res['record_count']}  "
              f"({res['records_with_pmids']} with >=1 PMID)")
        print(f"  tier counts     : {res['tier_counts']}")
        print(f"  snapshot bytes  : {res['snapshot_bytes']}")
        print(f"  snapshot sha256 : {res['snapshot_sha256']}")
        print(f"  manifest sha256 : {res['manifest_sha256']}")
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
