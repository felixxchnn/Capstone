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

import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import config
import io_utils

# Phase 2 application layer (sections 9-12). Imported at module load, like
# every other dependency here: if one of these cannot import, checks.py must
# not run at all rather than silently skip its own integrity gate.
import case_study
import evidence
import report
import sample_profile


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
        total = len(self.passed) + len(self.failed)
        print()
        print("=" * 74)
        print(f"  {len(self.passed)}/{total} checks passed")
        print(f"  {len(self.warnings)} warnings")
        print(f"  {len(self.failed)} failures")
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


# --------------------------------------------------------------------------
# Phase 2 application-layer integrity  (sections 9-12)
#
# Standalone helpers, deliberately independent of the Phase 1 loaders above.
# Every Phase 2 check is fail-closed: a missing file, a hash mismatch, or an
# exception in the check body is a hard [FAIL] with a readable reason, never a
# skip and never a fabricated pass. Nothing in this block reads, evaluates, or
# reports the held-out test split.
# --------------------------------------------------------------------------

_HASH_CHUNK = 1 << 20


def _sha256_file(path) -> str:
    """SHA-256 of a file, read in 1 MiB chunks (matches the project hash table)."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(_HASH_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def _reject_json_constant(token: str):
    """json.loads parse_constant hook: NaN / Infinity / -Infinity are forbidden."""
    raise ValueError(f"non-standard JSON constant {token!r} is not allowed here")


def _load_strict_json(path):
    """Load JSON, hard-failing on any NaN / Infinity / -Infinity literal."""
    text = Path(path).read_text(encoding="utf-8")
    return json.loads(text, parse_constant=_reject_json_constant)


def _phase2_application_checks(result: "CheckResult") -> None:
    """
    Sections 9-12: close the validation loop across the dataset, the
    reconstructed models, the Phase 2 case study, the drug-gene evidence
    snapshot, and the offline HTML report.
    """
    approved = {config.DEMO_VERIFICATION_MODEL_ID, config.DEMO_TUMOR_SAMPLE_ID}

    def guarded(name: str, fn) -> None:
        """Run fn() -> (ok, detail); any exception becomes a readable [FAIL]."""
        try:
            ok, detail = fn()
        except Exception as exc:  # noqa: BLE001 -- convert to an explicit failure
            ok, detail = False, f"{type(exc).__name__}: {exc}"
        result.check(name, bool(ok), detail)

    # ---- shared, best-effort loads (failure surfaces per-check) -----------
    try:
        cs = _load_strict_json(config.CASE_STUDY_JSON_FILE)
    except Exception as exc:  # noqa: BLE001
        cs = {"_error": f"{type(exc).__name__}: {exc}"}
    try:
        splits_assign = _load_strict_json(
            config.PROCESSED_DIR / "splits.json")["assignment"]
    except Exception as exc:  # noqa: BLE001
        splits_assign = {"_error": f"{type(exc).__name__}: {exc}"}
    try:
        sp_series, sp_prov = sample_profile.load_external_sample()
    except Exception as exc:  # noqa: BLE001
        sp_series, sp_prov = None, {"_error": f"{type(exc).__name__}: {exc}"}

    # ====================================================================
    _section("9. Phase 2 committed artifact identity")

    artifacts = {
        "case_study.json": config.CASE_STUDY_JSON_FILE,
        "phase2_report.html": config.REPORT_HTML_FILE,
        "DGIdb snapshot TSV": config.DGIDB_SNAPSHOT_FILE,
        "DGIdb manifest JSON": config.DGIDB_MANIFEST_FILE,
        "BG003082 GCT": config.DEMO_TUMOR_GCT_FILE,
        "ensembl_map.csv": config.ENSEMBL_MAP_FILE,
    }

    def _c9_present():
        missing = sorted(n for n, p in artifacts.items() if not Path(p).is_file())
        return (not missing), (
            f"missing: {missing} -- restore from git, do not download replacements"
            if missing else "all present")
    guarded("9.1 all six required Phase 2 artifacts are present", _c9_present)

    def _c9_cs_hash():
        got = _sha256_file(config.CASE_STUDY_JSON_FILE)
        return got == config.CASE_STUDY_JSON_SHA256, f"{got} vs pinned {config.CASE_STUDY_JSON_SHA256}"
    guarded("9.2 case_study.json SHA-256 == config.CASE_STUDY_JSON_SHA256", _c9_cs_hash)

    def _c9_html_hash():
        got = _sha256_file(config.REPORT_HTML_FILE)
        return got == config.REPORT_HTML_SHA256, f"{got} vs pinned {config.REPORT_HTML_SHA256}"
    guarded("9.3 phase2_report.html SHA-256 == config.REPORT_HTML_SHA256", _c9_html_hash)

    def _c9_strict():
        _load_strict_json(config.CASE_STUDY_JSON_FILE)
        return True, "parsed; no NaN / Infinity / -Infinity literals"
    guarded("9.4 case_study.json parses under strict JSON rules", _c9_strict)

    def _c9_schema():
        got = cs.get("schema_version")
        return got == config.CASE_STUDY_SCHEMA_VERSION, f"schema_version={got!r}"
    guarded("9.5 case_study.json schema_version == config.CASE_STUDY_SCHEMA_VERSION", _c9_schema)

    # ====================================================================
    _section("10. Phase 2 sample reconciliation and leakage prevention")

    def _c10_samples():
        s = set(cs.get("samples", {}))
        r = set(cs.get("rankings", {}))
        return (s == approved and r == approved), f"samples={sorted(s)} rankings={sorted(r)}"
    guarded("10.1 case-study samples and rankings are exactly {ACH-000364, BG003082}", _c10_samples)

    def _c10_anchor():
        anchor = config.DEMO_VERIFICATION_MODEL_ID
        sa = cs["samples"][anchor]
        conds = {
            "case-study depmap_split == val": sa["depmap_split"] == "val",
            "splits.json assignment == val": splits_assign.get(anchor) == "val",
            "in_training_split is False": sa["in_training_split"] is False,
            "prediction_status == held_out_prediction":
                sa["prediction_status"] == config.PREDICTION_STATUS_HELD_OUT,
            "outcome_status == measured_crispr":
                sa["outcome_status"] == config.OUTCOME_STATUS_MEASURED,
        }
        bad = [k for k, v in conds.items() if not v]
        return (not bad), f"failed: {bad}"
    guarded("10.2 ACH-000364 is a held-out val sample (held_out_prediction / measured_crispr)", _c10_anchor)

    def _c10_tumor():
        tumor = config.DEMO_TUMOR_SAMPLE_ID
        st = cs["samples"][tumor]
        conds = {
            "absent_from_all_depmap_splits is True":
                st["absent_from_all_depmap_splits"] is True,
            "absent from splits.json assignment": tumor not in splits_assign,
            "analysis_role == exploratory_external_prediction":
                st["analysis_role"] == "exploratory_external_prediction",
            "prediction_status == exploratory_external_prediction":
                st["prediction_status"] == config.PREDICTION_STATUS_EXPLORATORY,
            "outcome_status == unavailable":
                st["outcome_status"] == config.OUTCOME_STATUS_UNAVAILABLE,
        }
        bad = [k for k, v in conds.items() if not v]
        return (not bad), f"failed: {bad}"
    guarded("10.3 BG003082 is absent from every split (exploratory_external_prediction / unavailable)", _c10_tumor)

    def _c10_order():
        if sp_series is None:
            return False, sp_prov["_error"]
        gc = _load_strict_json(config.PROCESSED_DIR / "gene_columns.json")
        canonical = list(gc["canonical_columns"])
        return list(sp_series.index) == canonical, (
            f"|series|={len(sp_series)} |canonical|={len(canonical)}; "
            f"order_match={list(sp_series.index) == canonical}")
    guarded("10.4 fresh sample_profile load: canonical gene index & order == gene_columns.json", _c10_order)

    def _c10_recon_equal():
        if sp_series is None:
            return False, sp_prov["_error"]
        committed = cs["samples"][config.DEMO_TUMOR_SAMPLE_ID]["baseline_input"]["reconciliation"]
        fresh = sp_prov["reconciliation"]
        same = (json.dumps(committed, sort_keys=True) == json.dumps(fresh, sort_keys=True))
        return same, ("" if same else "fresh reconciliation != reconciliation in case_study.json")
    guarded("10.5 fresh reconciliation data exactly equals the reconciliation in case_study.json", _c10_recon_equal)

    def _c10_counts():
        if sp_series is None:
            return False, sp_prov["_error"]
        r = sp_prov["reconciliation"]
        conds = {
            "mapped == 18,427": r["canonical_genes_mapped"] == 18427,
            "missing == 33": r["canonical_genes_missing"] == 33,
            "18,427 + 33 == 18,460 == canonical_genes":
                r["canonical_genes_mapped"] + r["canonical_genes_missing"]
                == 18460 == r["canonical_genes"],
            "measured_zero == 1,407": r["canonical_genes_measured_zero"] == 1407,
            "measured_nonzero == 17,020": r["canonical_genes_measured_nonzero"] == 17020,
            "1,407 + 17,020 == 18,427 == mapped":
                r["canonical_genes_measured_zero"] + r["canonical_genes_measured_nonzero"]
                == 18427 == r["canonical_genes_mapped"],
        }
        bad = [k for k, v in conds.items() if not v]
        return (not bad), f"failed: {bad}"
    guarded("10.6 BG003082 counts reconcile: 18,427 + 33 = 18,460 ; 1,407 + 17,020 = 18,427", _c10_counts)

    def _c10_clean():
        if sp_series is None:
            return False, sp_prov["_error"]
        r = sp_prov["reconciliation"]
        conds = {
            "no canonical-Entrez collisions": r["canonical_id_collisions"] == 0,
            "no duplicate external identifiers": r["duplicate_external_ids"] == 0,
            "no gene-symbol fallback applied":
                str(r["symbol_fallback"]).startswith("not attempted"),
        }
        bad = [k for k, v in conds.items() if not v]
        return (not bad), f"failed: {bad}"
    guarded("10.7 BG003082: no identifier collisions, no duplicate ids, no gene-symbol fallback", _c10_clean)

    # ====================================================================
    _section("11. Phase 2 rankings and drug-gene evidence")

    def _c11_cs_validate():
        rep = case_study.validate(verbose=False)
        return rep["n_fail"] == 0, f"{rep['n_fail']} case_study.validate failure(s)"
    guarded("11.1 case_study.validate(verbose=False) reports zero failures", _c11_cs_validate)

    def _c11_rankings():
        problems = []
        for sample in sorted(approved):
            sb = cs["rankings"][sample]
            if set(sb) != {"ridge_pca", "ridge_head"}:
                problems.append(f"{sample}: models {sorted(sb)}")
                continue
            for mk, mb in sb.items():
                genes = mb["genes"]
                if len(genes) != config.TOP_N_DEPENDENCIES:
                    problems.append(f"{sample}/{mk}: {len(genes)} genes")
                if mb["n_targets_ranked"] != 4297:
                    problems.append(f"{sample}/{mk}: n_targets_ranked={mb['n_targets_ranked']}")
                if [g["rank"] for g in genes] != list(range(1, config.TOP_N_DEPENDENCIES + 1)):
                    problems.append(f"{sample}/{mk}: ranks != 1..{config.TOP_N_DEPENDENCIES}")
                vals = [g["predicted_geneeffect"] for g in genes]
                if any(vals[i] > vals[i + 1] for i in range(len(vals) - 1)):
                    problems.append(f"{sample}/{mk}: predicted GeneEffect not ascending")
        return (not problems), f"{problems}"
    guarded("11.2 each sample: 2 model rankings x 25 targets, ranks 1..25, predicted GeneEffect ascending", _c11_rankings)

    def _displayed_entrez():
        disp = set()
        for sb in cs["rankings"].values():
            for mb in sb.values():
                disp |= {g["entrez_id"] for g in mb["genes"]}
        return disp

    def _c11_keys():
        disp = _displayed_entrez()
        ev = cs["drug_gene_interaction_evidence"]
        keys = set(ev["by_entrez"])
        cov_n = ev["coverage"]["n_distinct_genes"]
        ok = (keys == disp and len(disp) == cov_n and len(disp) == 56)
        return ok, f"|displayed|={len(disp)} |by_entrez|={len(keys)} coverage.n_distinct_genes={cov_n}"
    guarded("11.3 evidence by_entrez keys == union of displayed Entrez IDs (derived; 56 distinct)", _c11_keys)

    def _c11_coverage():
        disp = _displayed_entrez()
        ev = cs["drug_gene_interaction_evidence"]
        cov = ev["coverage"]
        recount = {"cited": 0, "source_only": 0, "none_in_filtered_snapshot": 0}
        for b in ev["by_entrez"].values():
            recount[b["evidence_status"]] = recount.get(b["evidence_status"], 0) + 1
        ok = (cov["n_cited"] + cov["n_source_only"] + cov["n_none_in_filtered_snapshot"]
              == cov["n_distinct_genes"] == len(disp)
              and recount["cited"] == cov["n_cited"]
              and recount["source_only"] == cov["n_source_only"]
              and recount["none_in_filtered_snapshot"] == cov["n_none_in_filtered_snapshot"])
        return ok, f"coverage={cov} recount={recount} |displayed|={len(disp)}"
    guarded("11.4 evidence coverage classes sum to the distinct displayed-gene count", _c11_coverage)

    def _c11_records():
        ev = cs["drug_gene_interaction_evidence"]
        tiers = set(config.DGIDB_DIRECTION_TIERS)
        incl = set(config.DGIDB_INCLUDED_SOURCES)
        disc = config.DGIDB_EVIDENCE_DISCLAIMER
        bad = []
        for ent, bucket in ev["by_entrez"].items():
            if bucket["entrez_id"] != ent:
                bad.append(f"{ent}: bucket entrez_id {bucket['entrez_id']}")
            for r in bucket["records"]:
                if r["entrez_id"] != ent:
                    bad.append(f"{ent}: record entrez_id {r['entrez_id']}")
                if r["direction_tier"] not in tiers:
                    bad.append(f"{ent}: direction_tier {r['direction_tier']}")
                if r["interaction_source"] not in incl:
                    bad.append(f"{ent}: source {r['interaction_source']}")
                if r["disclaimer"] != disc:
                    bad.append(f"{ent}: disclaimer mismatch")
        return (not bad), f"{bad[:8]}"
    guarded("11.5 every evidence record: Entrez matches bucket; tier/source/disclaimer in approved vocab", _c11_records)

    def _c11_hashes():
        ev_ret = cs["drug_gene_interaction_evidence"]["retrieval"]
        snap_disk = _sha256_file(config.DGIDB_SNAPSHOT_FILE)
        man_disk = _sha256_file(config.DGIDB_MANIFEST_FILE)
        man = _load_strict_json(config.DGIDB_MANIFEST_FILE)
        inp = cs["input_artifact_sha256"]
        conds = {
            "evidence.snapshot_sha256 == snapshot on disk": ev_ret["snapshot_sha256"] == snap_disk,
            "evidence.manifest_sha256 == manifest on disk": ev_ret["manifest_sha256"] == man_disk,
            "evidence.snapshot_sha256 == manifest.snapshot.sha256":
                ev_ret["snapshot_sha256"] == man["snapshot"]["sha256"],
            "case_study input hash (snapshot) agrees":
                inp["data/external/dgidb/dgidb_2026-06b.interactions.filtered.tsv"] == snap_disk,
            "case_study input hash (manifest) agrees":
                inp["data/external/dgidb/dgidb_2026-06b.manifest.json"] == man_disk,
        }
        bad = [k for k, v in conds.items() if not v]
        return (not bad), f"failed: {bad}"
    guarded("11.6 case_study evidence snapshot/manifest hashes match the committed DGIdb files", _c11_hashes)

    def _c11_snapshot():
        rep = evidence.validate_snapshot()
        return (rep["n_fail"] == 0 and len(rep["checks"]) == 34), (
            f"{rep['n_fail']} failure(s), {len(rep['checks'])} checks (expected 34)")
    guarded("11.7 evidence.validate_snapshot() reports zero failures (34 checks)", _c11_snapshot)

    # ====================================================================
    _section("12. Phase 2 offline report integrity")

    def _c12_validate():
        rep = report.validate(verbose=False)
        return rep["n_fail"] == 0, f"{rep['n_fail']} report.validate failure(s)"
    guarded("12.1 report.validate(verbose=False) reports zero failures", _c12_validate)

    def _c12_embedded():
        doc = Path(config.REPORT_HTML_FILE).read_text(encoding="utf-8")
        m = re.search(
            r'<script type="application/json" id="case-study-data">(.*?)</script>',
            doc, flags=re.S)
        if not m:
            return False, "no #case-study-data script block found"
        raw = m.group(1)
        if "</" in raw or "<script" in raw.lower():
            return False, "embedded JSON block is not fully escaped"
        unesc = (raw.replace("\\u0026", "&")
                    .replace("\\u003c", "<")
                    .replace("\\u003e", ">"))
        embedded = json.loads(unesc, parse_constant=_reject_json_constant)
        committed = _load_strict_json(config.CASE_STUDY_JSON_FILE)
        same = (json.dumps(embedded, sort_keys=True)
                == json.dumps(committed, sort_keys=True))
        return same, ("" if same else "embedded JSON != committed case_study.json")
    guarded("12.2 JSON embedded in #case-study-data equals the strict-parsed committed case_study.json", _c12_embedded)

    def _c12_offline():
        low = Path(config.REPORT_HTML_FILE).read_text(encoding="utf-8").lower()
        bad = []
        if "<script src" in low:
            bad.append("<script src=>")
        if 'rel="stylesheet"' in low or "rel='stylesheet'" in low:
            bad.append("external stylesheet <link>")
        if re.search(r"@import\s+url\(", low):
            bad.append("@import url(...)")
        if "fetch(" in low:
            bad.append("fetch(")
        if "xmlhttprequest" in low:
            bad.append("XMLHttpRequest")
        return (not bad), f"remote/runtime dependency markers: {bad}"
    guarded("12.3 report has no remote runtime dependency (no script src / link / @import / fetch / XHR)", _c12_offline)

    def _c12_disclaimer():
        doc = Path(config.REPORT_HTML_FILE).read_text(encoding="utf-8")
        return config.DGIDB_EVIDENCE_DISCLAIMER in doc, "fixed non-efficacy disclaimer absent from report"
    guarded("12.4 report contains the fixed non-efficacy disclaimer", _c12_disclaimer)


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
    # Frozen legacy build asymmetry, not a general allowance: build_dataset.py
    # never re-sliced the saved crispr matrix after 3 columns were dropped
    # from expression only. See config.CRISPR_LEGACY_EXTRA_GENES for the full
    # history. This asserts full set equality against canonical plus exactly
    # those 3 named columns -- catching both a new extra column (checked
    # above by the old form) AND a canonical gene silently disappearing from
    # crispr, which a one-sided subset/difference check would miss.
    crispr_genes = set(crispr.columns)
    canonical_genes = set(canonical)
    expected_crispr_genes = canonical_genes | config.CRISPR_LEGACY_EXTRA_GENES
    missing_canonical = canonical_genes - crispr_genes
    unexpected_extra = crispr_genes - expected_crispr_genes
    result.check(
        "crispr columns equal the canonical space plus the named frozen "
        "legacy exception, exactly",
        crispr_genes == expected_crispr_genes,
        f"{len(missing_canonical)} canonical gene(s) missing from crispr "
        f"{sorted(missing_canonical)[:10]}; {len(unexpected_extra)} "
        f"unexpected extra column(s) {sorted(unexpected_extra)[:10]}. "
        f"Expected extras are exactly {sorted(config.CRISPR_LEGACY_EXTRA_GENES)} "
        f"(see config.CRISPR_LEGACY_EXTRA_GENES). This is a NEW discrepancy, "
        f"not the known one -- do not wave it through.",
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
    result.check(
        "the 3 frozen legacy extra crispr columns are absent from the "
        "selective (published) target set",
        set(selective_genes).isdisjoint(config.CRISPR_LEGACY_EXTRA_GENES),
        f"one or more of {sorted(config.CRISPR_LEGACY_EXTRA_GENES)} has "
        f"entered selective_genes.json -- this would mean a published "
        f"result now depends on the legacy build asymmetry.",
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

    # ---- Phase 2 application-layer integrity (sections 9-12) -------------
    _phase2_application_checks(result)

    return result.summary()


if __name__ == "__main__":
    sys.exit(main())
