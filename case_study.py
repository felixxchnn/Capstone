"""
case_study.py
=============
Phase 2 proof-of-concept orchestration.

    py case_study.py --build           # write data/processed/case_study.json
    py case_study.py --validate        # regenerate + byte-compare + re-check invariants
    py case_study.py --self-test       # offline unit + structural checks

What this does
--------------
Runs the two **reconstructed** frozen Phase 1 linear models (`ridge_pca`,
`ridge_head`) over two samples and assembles one deterministic JSON artifact:

* **ACH-000364** (U-2 OS, DepMap `val` split) -- internal *verification anchor*.
  `prediction_status = held_out_prediction`, `outcome_status = measured_crispr`.
  Observed CRISPR values are attached to already-ranked genes **after** ranking,
  as a pipeline-verification example -- never used for selection or ordering.
* **BG003082** (Sid Sijbrandij osteosarcoma primary tumour, `osteosarc.com`,
  CC0 1.0) -- *exploratory external prediction*. Absent from every DepMap split.
  `analysis_role = exploratory_external_prediction`,
  `outcome_status = unavailable`. No observed outcome exists, is loaded, or is
  computed for it.

Plus the locked five-line val-split osteosarcoma descriptive aggregate
(`capstone/scope-decisions.md`, 2026-08-29) and, for the displayed genes,
retrieved **drug-gene interaction evidence** from the committed offline DGIdb
snapshot (retrieval only -- never treatment advice).

Guarantees (enforced by --self-test / --validate)
------------------------------------------------
* This module never imports scikit-learn, never calls ``fit`` / ``fit_transform``,
  never imports ``reconstruct_fitted``. Inference is `fitted_artifacts` only.
  (`baseline.per_target_spearman` -- pure numpy/scipy, no fit -- is used for the
  mandated osteosarcoma metric; that is the sole `baseline` import.)
* Reconstructed models are labelled, verbatim, "reconstructed fitted state at the
  frozen Phase 1 alpha from the unchanged frozen training data" -- not the
  unavailable historical fitted objects.
* Exact feature / embedding / target order is asserted before any inference.
* No observed outcome enters prediction or ranking; evidence retrieval happens
  only after the top-N is frozen and never changes it.
* Output JSON: schema-versioned, stably ordered, fixed-precision floats, no
  wall-clock, no absolute paths, byte-identical on regeneration.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import platform
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

import config
import io_utils
import fitted_artifacts
import sample_profile
import evidence
from fitted_artifacts import RECONSTRUCTION_STATUS
from gene_ids import parse_gene_label
from baseline import per_target_spearman  # pure numpy/scipy; no fit, no sklearn state

SCHEMA_VERSION = "case-study/1"
SOURCE_COMMIT = "d6a9b91148c235b1d1215553a3b46b958bc1b212"

OUT_FILE = config.PROCESSED_DIR / "case_study.json"

# fixed-precision rounding (documented in the artifact)
PRED_DP = 10          # predicted / observed GeneEffect
AGG_DP = 6            # osteosarcoma mean per-target Spearman + delta

_HASH_CHUNK = 1 << 20
_REPO_ROOT = config.PROJECT_ROOT


class CaseStudyError(RuntimeError):
    """Raised on any pre-flight / integrity failure -- always a hard stop."""


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(_HASH_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def _rel(path: str | Path) -> str:
    """Repo-relative POSIX path -- never an absolute local path in the artifact."""
    return Path(path).resolve().relative_to(_REPO_ROOT).as_posix()


def _r(x: float, dp: int) -> float:
    v = round(float(x), dp)
    return 0.0 if v == 0.0 else v          # normalise -0.0


def _write_json_deterministic(obj, path: Path) -> None:
    path.write_text(
        json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=True,
                   allow_nan=False) + "\n",
        encoding="utf-8",
    )


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------

# every committed file this module reads, hashed into the artifact
_INPUT_FILES = {
    "expression.npz": config.PROCESSED_DIR / "expression.npz",
    "expression.labels.json": config.PROCESSED_DIR / "expression.labels.json",
    "crispr_effect.npz": config.PROCESSED_DIR / "crispr_effect.npz",
    "crispr_effect.labels.json": config.PROCESSED_DIR / "crispr_effect.labels.json",
    "selective_genes.json": config.PROCESSED_DIR / "selective_genes.json",
    "splits.json": config.PROCESSED_DIR / "splits.json",
    "model_metadata.csv": config.PROCESSED_DIR / "model_metadata.csv",
    "gene_columns.json": config.PROCESSED_DIR / "gene_columns.json",
    "geneformer_embeddings.csv": config.PROCESSED_DIR / "geneformer_embeddings.csv",
    "ensembl_map.csv": config.ENSEMBL_MAP_FILE,
    "BG003082.gene_tpm.gct.gz": config.DEMO_TUMOR_GCT_FILE,
    "geneformer_bg003082_embedding.csv":
        config.PROCESSED_DIR / "geneformer_bg003082_embedding.csv",
    "geneformer_bg003082_embedding.provenance.json":
        config.PROCESSED_DIR / "geneformer_bg003082_embedding.provenance.json",
    "reconstructed_fitted/baseline_ridge_pca/manifest.json":
        config.RECONSTRUCTED_FITTED_DIR / "baseline_ridge_pca" / "manifest.json",
    "reconstructed_fitted/head_ridge_head/manifest.json":
        config.RECONSTRUCTED_FITTED_DIR / "head_ridge_head" / "manifest.json",
    "dgidb_2026-06b.interactions.filtered.tsv": config.DGIDB_SNAPSHOT_FILE,
    "dgidb_2026-06b.manifest.json": config.DGIDB_MANIFEST_FILE,
}


class Inputs:
    def __init__(self):
        self.baseline = fitted_artifacts.load_baseline_ridge_pca()
        self.head = fitted_artifacts.load_head_ridge()

        self.expression = io_utils.load_matrix(config.PROCESSED_DIR / "expression")
        self.crispr = io_utils.load_matrix(config.PROCESSED_DIR / "crispr_effect")
        self.embeddings = _load_embeddings()
        self.metadata = io_utils.load_table(config.PROCESSED_DIR / "model_metadata")
        self.metadata.index.name = config.MODEL_ID

        splits_payload = io_utils.load_json(config.PROCESSED_DIR / "splits.json")
        self.assignment = pd.Series(splits_payload["assignment"], name="split")

        selective = io_utils.load_json(config.PROCESSED_DIR / "selective_genes.json")
        self.target_labels = [g for g in selective["genes"] if g in self.crispr.columns]

        parsed = [parse_gene_label(t) for t in self.target_labels]
        if any(p is None for p in parsed):
            raise CaseStudyError("a selective-gene label did not parse as 'SYMBOL (ENTREZ)'")
        self.target_entrez = [p.entrez for p in parsed]
        self.target_entrez_int = [int(e) for e in self.target_entrez]
        self.target_symbol = [p.symbol for p in parsed]

        self.input_sha256 = {
            _rel(p): _sha256_file(p) for _, p in sorted(_INPUT_FILES.items())
        }
        self._preflight()

    # -- pre-flight (all hard stops) ------------------------------------
    def _preflight(self) -> None:
        b, h = self.baseline, self.head

        if [str(c) for c in self.expression.columns] != b.feature_names:
            raise CaseStudyError("expression columns != baseline artifact feature order")
        if [str(c) for c in self.embeddings.columns] != h.feature_names:
            raise CaseStudyError("embedding columns != head artifact feature order")
        if not (b.target_names == h.target_names == self.target_labels):
            raise CaseStudyError("target order disagreement (artifacts vs selective_genes)")
        b.assert_feature_order([str(c) for c in self.expression.columns])
        b.assert_target_order(self.target_labels)
        h.assert_feature_order([str(c) for c in self.embeddings.columns])
        h.assert_target_order(self.target_labels)

        anchor = config.DEMO_VERIFICATION_MODEL_ID
        self.anchor_split = self.assignment.get(anchor)
        if self.anchor_split == "train":
            raise CaseStudyError(
                f"{anchor} is in the TRAINING split -- refusing to use it as a "
                f"held-out verification anchor")
        if self.anchor_split != "val":
            raise CaseStudyError(
                f"{anchor} split is {self.anchor_split!r}, expected 'val'")
        for name in ("expression", "crispr", "embeddings"):
            if anchor not in getattr(self, name).index:
                raise CaseStudyError(f"{anchor} missing from {name}")

        tumor = config.DEMO_TUMOR_SAMPLE_ID
        if tumor in self.assignment.index:
            raise CaseStudyError(
                f"{tumor} appears in splits.json assignment -- it must be absent "
                f"from every DepMap split")

        if not isinstance(config.TOP_N_DEPENDENCIES, int) or config.TOP_N_DEPENDENCIES <= 0:
            raise CaseStudyError(
                f"config.TOP_N_DEPENDENCIES is {config.TOP_N_DEPENDENCIES!r}; "
                f"expected a positive int -- refusing to invent a top-N")


def _load_embeddings() -> pd.DataFrame:
    csv = config.PROCESSED_DIR / "geneformer_embeddings.csv"
    frame = pd.read_csv(csv, index_col=0)
    frame.index.name = config.MODEL_ID
    return frame.select_dtypes(include=[np.number])


# --------------------------------------------------------------------------
# ranking
# --------------------------------------------------------------------------

RANKING_RULE = (
    "more negative predicted CRISPR GeneEffect = stronger predicted dependency; "
    "sort ascending by RAW finite float64 predicted GeneEffect, most-negative "
    "first; ascending numeric Entrez ID breaks ONLY exact raw-value ties. "
    "predicted_geneeffect is rounded to %d dp for display AFTER ranks and top-N "
    "membership are frozen" % PRED_DP
)
OBSERVED_RANK_RULE = (
    "1-based position of the gene when all 4,297 targets are sorted ascending by "
    "RAW finite float64 observed CRISPR GeneEffect (most-negative first; "
    "ascending numeric Entrez ID breaks ONLY exact raw-value ties); null where "
    "ACH-000364 has no observed value for that target"
)


def _full_ranking(pred_vec: np.ndarray, entrez_int: list[int]) -> list[int]:
    """
    All target indices, strongest predicted dependency first.

    Sort key is the RAW finite float64 predicted GeneEffect (most-negative
    first); ascending numeric Entrez ID breaks ONLY exact raw-value ties.
    Rounding for display happens later, after ranks and top-N membership are
    frozen. A non-finite prediction (should not occur -- predictions are dot
    products of finite inputs) is pushed to the end in Entrez order so it can
    never enter the top-N.
    """
    finite = [j for j in range(len(pred_vec)) if np.isfinite(pred_vec[j])]
    nonfinite = [j for j in range(len(pred_vec)) if not np.isfinite(pred_vec[j])]
    finite.sort(key=lambda j: (float(pred_vec[j]), entrez_int[j]))
    nonfinite.sort(key=lambda j: entrez_int[j])
    return finite + nonfinite


def _observed_rank_lookup(observed_vec: np.ndarray,
                          entrez_int: list[int]) -> dict[int, int]:
    """
    1-based rank of each target with a finite observed value, over the full
    4,297-target universe. Sort key is the RAW finite float64 observed
    GeneEffect (most-negative first); ascending numeric Entrez ID breaks ONLY
    exact raw-value ties.
    """
    finite = [j for j in range(len(observed_vec)) if np.isfinite(observed_vec[j])]
    finite.sort(key=lambda j: (float(observed_vec[j]), entrez_int[j]))
    return {j: rank for rank, j in enumerate(finite, start=1)}


def build_model_ranking(inp: Inputs, model_key: str, pred_vec: np.ndarray,
                        observed_vec: np.ndarray | None) -> dict:
    order = _full_ranking(pred_vec, inp.target_entrez_int)
    top = order[: config.TOP_N_DEPENDENCIES]
    obs_ranks = (_observed_rank_lookup(observed_vec, inp.target_entrez_int)
                 if observed_vec is not None else None)

    genes = []
    for rank, j in enumerate(top, start=1):
        row = {
            "rank": rank,
            "entrez_id": inp.target_entrez[j],
            "symbol": inp.target_symbol[j],
            "predicted_geneeffect": _r(pred_vec[j], PRED_DP),
        }
        if observed_vec is not None:
            ov = observed_vec[j]
            row["observed_geneeffect"] = (None if not np.isfinite(ov)
                                          else _r(ov, PRED_DP))
            row["observed_rank"] = obs_ranks.get(j)
        genes.append(row)

    block = {
        "model": model_key,
        "model_provenance": RECONSTRUCTION_STATUS,
        "n_targets_ranked": len(order),
        "n_displayed": len(genes),
        "ranking_rule": RANKING_RULE,
        "not_a_recommendation": (
            "ranked genes are predicted dependencies, NOT therapeutic targets, "
            "treatments, or recommendations"),
        "genes": genes,
    }
    if observed_vec is not None:
        block["observed_rank_rule"] = OBSERVED_RANK_RULE
        block["observed_values_attached_after_ranking"] = True
        block["n_targets_with_observed_value"] = int(np.isfinite(observed_vec).sum())
    return block


# --------------------------------------------------------------------------
# drug-gene interaction evidence  (retrieval only)
# --------------------------------------------------------------------------

EVIDENCE_LABEL = "Drug\u2013gene interaction evidence"
PMID_SCOPE_NOTE = (
    "PMIDs are joined at the drug-gene / interaction-source group level "
    "(keyed by gene concept + drug concept + interaction source) and may span "
    "multiple interaction claims; they are drug-gene / source-level citations "
    "and may not specifically support the displayed interaction subtype")


def _shape_evidence_record(rec: dict) -> dict:
    pmids = [p for p in rec.get("pmids", "").split(";") if p]
    cited = bool(pmids)
    out = {
        "entrez_id": rec["entrez_id"],
        "gene_symbol": rec["gene_symbol"],
        "dgidb_gene_name": rec["dgidb_gene_name"],
        "gene_symbol_consistent": rec["gene_symbol_consistent"],
        "symbol_query_mismatch": rec.get("symbol_query_mismatch", "false"),
        "drug_name": rec["drug_name"],
        "drug_concept_id": rec["drug_concept_id"],
        "drug_claim_name": rec["drug_claim_name"],
        "interaction_source": rec["interaction_source"],
        "interaction_source_version": rec["interaction_source_version"],
        "source_license": rec["source_license"],
        "source_license_url": rec["source_license_url"],
        "interaction_type_raw": rec["interaction_type_raw"],
        "interaction_direction": rec["interaction_direction"],
        "direction_tier": rec["direction_tier"],
        "interaction_score": rec["interaction_score"],
        "drug_specificity_score": rec["drug_specificity_score"],
        "gene_specificity_score": rec["gene_specificity_score"],
        "evidence_score": rec["evidence_score"],
        "drug_is_approved": rec["drug_is_approved"],
        "drug_is_immunotherapy": rec["drug_is_immunotherapy"],
        "drug_is_antineoplastic": rec["drug_is_antineoplastic"],
        "curation_type": rec.get("curation_type", ""),
        "indication": rec.get("indication", ""),
        "pmids": pmids,
        "pmid_status": "cited" if cited else "source_only",
        "dgidb_release_tag": rec["dgidb_release_tag"],
        "record_key": rec["record_key"],
        "disclaimer": rec["disclaimer"],
    }
    out["pmid_scope_note"] = PMID_SCOPE_NOTE if cited else (
        "source-only interaction evidence: DGIdb records no claim-level "
        "publication for this drug-gene-source group in the filtered snapshot")
    return out


def collect_evidence(inp: Inputs, displayed: dict) -> dict:
    """`displayed` maps entrez_id (str) -> symbol. Called only after rankings freeze."""
    snap = evidence.load_snapshot()
    by_entrez: dict[str, dict] = {}
    n_cited = n_source_only = n_none = 0
    total_records = total_pmids = 0

    for entrez in sorted(displayed, key=int):
        recs = evidence.get_evidence_for_gene(
            entrez, symbol=displayed[entrez],
            top_k=config.TOP_K_EVIDENCE_PER_GENE, snapshot=snap)
        shaped = [_shape_evidence_record(r) for r in recs]
        if not shaped:
            status = "none_in_filtered_snapshot"
            n_none += 1
        elif any(s["pmid_status"] == "cited" for s in shaped):
            status = "cited"
            n_cited += 1
        else:
            status = "source_only"
            n_source_only += 1
        total_records += len(shaped)
        total_pmids += sum(len(s["pmids"]) for s in shaped)
        by_entrez[entrez] = {
            "entrez_id": entrez,
            "symbol": displayed[entrez],
            "evidence_status": status,
            "n_records": len(shaped),
            "records": shaped,
        }

    return {
        "label": EVIDENCE_LABEL,
        "framing": (
            "retrieval of recorded drug-gene interactions from a licence-filtered "
            "offline DGIdb snapshot. NOT treatments, NOT recommendations, NOT "
            "actionable therapy, NO efficacy claim."),
        "disclaimer": config.DGIDB_EVIDENCE_DISCLAIMER,
        "retrieval": {
            "snapshot_file": _rel(config.DGIDB_SNAPSHOT_FILE),
            "snapshot_sha256": _sha256_file(config.DGIDB_SNAPSHOT_FILE),
            "manifest_sha256": _sha256_file(config.DGIDB_MANIFEST_FILE),
            "method": ("evidence.get_evidence_for_gene(entrez_id, symbol) by Entrez "
                       "ID only; symbol is a flagged consistency check, never a key"),
            "top_k_per_direction_tier": config.TOP_K_EVIDENCE_PER_GENE,
            "direction_tiers": list(config.DGIDB_DIRECTION_TIERS),
            "retrieved_after_top_n_frozen": True,
            "evidence_availability_did_not_affect_selection_or_ranking": True,
        },
        "pmid_scope_note": PMID_SCOPE_NOTE,
        "coverage": {
            "n_distinct_genes": len(by_entrez),
            "n_cited": n_cited,
            "n_source_only": n_source_only,
            "n_none_in_filtered_snapshot": n_none,
            "total_records": total_records,
            "total_pmid_citations": total_pmids,
        },
        "by_entrez": by_entrez,
    }


# --------------------------------------------------------------------------
# locked osteosarcoma five-line descriptive aggregate
# --------------------------------------------------------------------------

def osteosarcoma_aggregate(inp: Inputs) -> dict:
    mask = config.osteosarcoma_mask(inp.metadata)
    cohort = sorted(
        mid for mid in inp.metadata.index[mask]
        if inp.assignment.get(mid) == "val"
    )
    if len(cohort) != 5:
        raise CaseStudyError(
            f"locked osteosarcoma cohort must be n=5, got {len(cohort)}: {cohort}")
    for mid in cohort:
        for name in ("expression", "crispr", "embeddings"):
            if mid not in getattr(inp, name).index:
                raise CaseStudyError(f"cohort line {mid} missing from {name}")

    X_expr = inp.expression.loc[cohort, inp.baseline.feature_names].to_numpy(dtype=float)
    X_emb = inp.embeddings.loc[cohort, inp.head.feature_names].to_numpy(dtype=float)
    Y_obs = inp.crispr.loc[cohort, inp.target_labels].to_numpy(dtype=float)  # NaN kept

    rho_pca = per_target_spearman(Y_obs, inp.baseline.predict(X_expr))
    rho_head = per_target_spearman(Y_obs, inp.head.predict(X_emb))

    finite_pca = np.isfinite(rho_pca)
    finite_head = np.isfinite(rho_head)
    common = finite_pca & finite_head
    n_common = int(common.sum())
    if n_common == 0:
        raise CaseStudyError("no target is finite for both models over the cohort")

    mean_pca = _r(rho_pca[common].mean(), AGG_DP)
    mean_head = _r(rho_head[common].mean(), AGG_DP)

    return {
        "status": (
            "DESCRIPTIVE and UNSTABLE because n=5. No confidence interval, no "
            "significance test. This does NOT replace or restate the frozen "
            "Phase 1 primary result (val per-target Spearman 0.2356 vs 0.2047 "
            "over 170 held-out cell lines); it is shown only so ACH-000364's "
            "single-line result is not read as cherry-picked."),
        "definition_source": "capstone/scope-decisions.md (2026-08-29 entry)",
        "cohort": {
            "predicate": ("config.osteosarcoma_mask(model_metadata) "
                          "[OncotreePrimaryDisease / OncotreeSubtype keyword] "
                          "intersected with splits.json assignment == 'val'"),
            "n": len(cohort),
            "model_ids": cohort,
        },
        "models": ["ridge_pca", "ridge_head"],
        "target_universe": len(inp.target_labels),
        "per_target_metric": (
            "baseline.per_target_spearman(observed, predicted) across the 5 "
            "cohort lines (min_samples=5; a constant column -> NaN); rules "
            "unchanged from Phase 1"),
        "common_finite_target_set": {
            "rule": "targets whose per-target Spearman is finite for BOTH models",
            "n_included": n_common,
            "n_excluded": int((~common).sum()),
            "excluded_ridge_pca_nonfinite": int((~finite_pca).sum()),
            "excluded_ridge_head_nonfinite": int((~finite_head).sum()),
            "excluded_reason": (
                "at n=5, a target is excluded when it has <5 finite "
                "(observed, predicted) pairs or a constant column for either model"),
        },
        "mean_per_target_spearman": {
            "ridge_pca": mean_pca,
            "ridge_head": mean_head,
            "rounding_dp": AGG_DP,
        },
        "delta_ridge_head_minus_ridge_pca": _r(mean_head - mean_pca, AGG_DP),
        "used_to_choose_model_or_alter_rankings": False,
    }


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------

def _reconstructed_model_block(loader, key: str) -> dict:
    m = loader.manifest
    return {
        "model": key,
        "provenance_status": RECONSTRUCTION_STATUS,
        "not_original_fitted_objects": (
            "the historical Phase 1 StandardScaler / PCA / Ridge objects were "
            "never serialised; this is a reconstruction re-fit on exactly the "
            "committed train split at the frozen alpha (read, not re-selected)"),
        "artifact_dir": _rel(loader.artifact_dir),
        "manifest_sha256": _sha256_file(loader.artifact_dir / "manifest.json"),
        "frozen_alpha": {
            "value": loader.alpha,
            "source": m["frozen_alpha"]["source_file"] + " :: "
            + m["frozen_alpha"]["source_json_path"],
            "selection": m["frozen_alpha"]["selection"],
        },
        "pipeline": m["pipeline"],
        "n_features": loader.n_features,
        "n_targets": loader.n_targets,
        "feature_order_sha256": loader.feature_order_sha256,
        "target_order_sha256": loader.target_order_sha256,
        "base_commit": m["base_commit"],
    }


def _bg003082_embedding_block(inp: Inputs) -> tuple[np.ndarray, dict]:
    prov = io_utils.load_json(
        config.PROCESSED_DIR / "geneformer_bg003082_embedding.provenance.json")
    sc = pd.read_csv(
        config.PROCESSED_DIR / "geneformer_bg003082_embedding.csv", index_col=0)
    if list(sc.index) != [config.DEMO_TUMOR_SAMPLE_ID]:
        raise CaseStudyError(f"sidecar index {list(sc.index)} != ['{config.DEMO_TUMOR_SAMPLE_ID}']")
    if [str(c) for c in sc.columns] != inp.head.feature_names:
        raise CaseStudyError("sidecar embedding columns != head artifact feature order")
    arr = sc.to_numpy(dtype=float)
    if arr.shape != (1, 768):
        raise CaseStudyError(f"sidecar shape {arr.shape} != (1, 768)")
    if not np.isfinite(arr).all():
        raise CaseStudyError("sidecar embedding has non-finite values")
    disk_sha = _sha256_file(config.PROCESSED_DIR / "geneformer_bg003082_embedding.csv")
    prov_sha = prov["embedding"]["sidecar_csv_sha256"]
    if disk_sha != prov_sha:
        raise CaseStudyError(
            f"sidecar sha256 {disk_sha} != provenance-recorded {prov_sha}")

    block = {
        "sidecar_file": "geneformer_bg003082_embedding.csv",
        "sidecar_sha256": disk_sha,
        "provenance_file": "geneformer_bg003082_embedding.provenance.json",
        "shape": [1, 768],
        "all_finite": True,
        "model": f"{prov['environment']['model_repo']} / {prov['environment']['model_subdir']}",
        "geneformer_revision_pinned": prov["environment"]["geneformer_revision_pinned"],
        "commensurability_caveats": [
            "SEPARATELY GENERATED external embedding: produced on Kaggle by "
            "capstone/kaggle_bg003082_embedding.py, not part of the frozen "
            "1,140-row geneformer_embeddings.csv.",
            "The historical Phase 1 Geneformer code revision was NOT captured; "
            "this sidecar pins revision 04c2b2e8..., a fresh pin for this run.",
            "Commensurability of this embedding with the Phase 1 training "
            "embeddings is NOT proven (bulk primary-tumour TPM input; NCBI "
            "gene2ensembl map, not the vanished mygene map; fresh revision pin).",
            "The BG003082 prediction is EXPLORATORY and cannot be interpreted "
            "as validated model performance.",
        ],
        "provenance_disclosures": list(prov.get("disclosures", [])),
    }
    return arr, block


def assemble() -> dict:
    inp = Inputs()

    # ---- reconstructed models -------------------------------------------
    models = {
        "ridge_pca": _reconstructed_model_block(inp.baseline, "ridge_pca"),
        "ridge_head": _reconstructed_model_block(inp.head, "ridge_head"),
    }

    # ---- ACH-000364 : held-out verification anchor --------------------
    anchor = config.DEMO_VERIFICATION_MODEL_ID
    X_anchor_expr = inp.expression.loc[[anchor], inp.baseline.feature_names].to_numpy(float)
    X_anchor_emb = inp.embeddings.loc[[anchor], inp.head.feature_names].to_numpy(float)
    pred_anchor = {
        "ridge_pca": inp.baseline.predict(X_anchor_expr)[0],
        "ridge_head": inp.head.predict(X_anchor_emb)[0],
    }
    # observed CRISPR -- loaded and attached ONLY now, after prediction
    observed_anchor = inp.crispr.loc[anchor, inp.target_labels].to_numpy(dtype=float)

    # ---- BG003082 : exploratory external prediction -------------------
    tumor = config.DEMO_TUMOR_SAMPLE_ID
    series, sp_prov = sample_profile.load_external_sample()
    if list(series.index) != inp.baseline.feature_names:
        # label sets must match; reorder to the artifact's exact order
        if set(series.index) != set(inp.baseline.feature_names):
            raise CaseStudyError("BG003082 profile label set != baseline feature set")
        series = series.reindex(inp.baseline.feature_names)
    X_tumor_expr = series.to_numpy(dtype=float).reshape(1, -1)
    n_missing_expr = int(np.isnan(X_tumor_expr).sum())
    if n_missing_expr != sp_prov["reconciliation"]["canonical_genes_missing"]:
        raise CaseStudyError(
            f"BG003082 baseline missing-feature count {n_missing_expr} != "
            f"sample_profile canonical_genes_missing "
            f"{sp_prov['reconciliation']['canonical_genes_missing']}")
    X_tumor_emb, tumor_emb_block = _bg003082_embedding_block(inp)
    pred_tumor = {
        "ridge_pca": inp.baseline.predict(X_tumor_expr)[0],
        "ridge_head": inp.head.predict(X_tumor_emb)[0],
    }

    # ---- rankings (predictions only; no outcome in the ranking path) ---
    rankings = {
        anchor: {
            "ridge_pca": build_model_ranking(inp, "ridge_pca",
                                             pred_anchor["ridge_pca"], observed_anchor),
            "ridge_head": build_model_ranking(inp, "ridge_head",
                                              pred_anchor["ridge_head"], observed_anchor),
        },
        tumor: {
            "ridge_pca": build_model_ranking(inp, "ridge_pca",
                                             pred_tumor["ridge_pca"], None),
            "ridge_head": build_model_ranking(inp, "ridge_head",
                                              pred_tumor["ridge_head"], None),
        },
    }

    # ---- evidence : AFTER rankings are frozen --------------------------
    displayed: dict[str, str] = {}
    for sample_block in rankings.values():
        for model_block in sample_block.values():
            for g in model_block["genes"]:
                displayed[g["entrez_id"]] = g["symbol"]
    evidence_block = collect_evidence(inp, displayed)

    # ---- osteosarcoma locked aggregate ------------------------------
    osteo = osteosarcoma_aggregate(inp)

    # ---- samples metadata ------------------------------------------
    samples = {
        anchor: {
            "role": "verification_anchor",
            "cell_line": "U-2 OS",
            "prediction_status": config.PREDICTION_STATUS_HELD_OUT,
            "outcome_status": config.OUTCOME_STATUS_MEASURED,
            "depmap_split": inp.anchor_split,
            "split_assertion": "asserted == 'val'; hard-stop if 'train'",
            "in_training_split": False,
            "baseline_input": {
                "source": "expression.npz row",
                "n_features": inp.baseline.n_features,
                "feature_order_sha256": inp.baseline.feature_order_sha256,
                "missing_features": 0,
                "imputed_features": 0,
            },
            "head_input": {
                "source": "geneformer_embeddings.csv row (frozen 1,140-row matrix)",
                "n_features": inp.head.n_features,
                "feature_order_sha256": inp.head.feature_order_sha256,
                "missing_features": 0,
            },
            "observed_crispr": {
                "source": "crispr_effect.npz row",
                "role": ("pipeline-verification example only; attached to already-"
                         "ranked genes AFTER ranking; never used for selection, "
                         "model choice, or ordering"),
                "n_targets_with_value": int(np.isfinite(observed_anchor).sum()),
                "n_targets_missing": int((~np.isfinite(observed_anchor)).sum()),
            },
        },
        tumor: {
            "role": "exploratory_external_sample",
            "description": ("Sid Sijbrandij osteosarcoma primary-tumour RNA-seq, "
                            "osteosarc.com, CC0 1.0, resected 2022-12-16; bulk "
                            "tumour tissue, a real domain shift from cultured "
                            "DepMap cell lines"),
            "analysis_role": "exploratory_external_prediction",
            "prediction_status": config.PREDICTION_STATUS_EXPLORATORY,
            "outcome_status": config.OUTCOME_STATUS_UNAVAILABLE,
            "absent_from_all_depmap_splits": True,
            "observed_outcome": ("none exists; none is loaded, invented, or "
                                 "computed for this sample"),
            "baseline_input": {
                "source": "sample_profile.load_external_sample() from the committed GCT",
                "gct_file": sp_prov["gct_file"],
                "ensembl_map_file": sp_prov["ensembl_map_file"],
                "gene_columns_file": sp_prov["gene_columns_file"],
                "transformation": sp_prov["transformation"],
                "n_features": inp.baseline.n_features,
                "feature_order_sha256": inp.baseline.feature_order_sha256,
                "reconciliation": sp_prov["reconciliation"],
                "missing_features_represented_as_nan": n_missing_expr,
                "imputation": ("the reconstructed baseline artifact applies its "
                               "STORED training-mean impute vector to those NaN "
                               "features (same impute step every Phase 1 row uses)"),
            },
            "head_input": tumor_emb_block,
        },
    }

    artifact = {
        "schema_version": SCHEMA_VERSION,
        "artifact": OUT_FILE.name,
        "title": "Phase 2 proof-of-concept case study",
        "description": (
            "Ranked predicted CRISPR dependencies for two samples from the two "
            "reconstructed frozen Phase 1 linear models, plus retrieved drug-gene "
            "interaction evidence and a descriptive five-line osteosarcoma "
            "aggregate. Demonstrates a possible future workflow; it does NOT "
            "predict patient treatment response and makes NO efficacy or "
            "treatment-recommendation claim."),
        "source_commit": SOURCE_COMMIT,
        "generated_by": "case_study.py",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": __import__("scipy").__version__,
            "scikit_learn": "not used by case_study.py (inference is fitted_artifacts only)",
        },
        "methodology": {
            "inference": ("fitted_artifacts.py only -- closed-form array "
                          "arithmetic; case_study.py imports no scikit-learn and "
                          "calls no fit()/fit_transform()"),
            "models": ("ridge_pca and ridge_head, loaded as reconstructed fitted "
                       "state at the frozen Phase 1 alpha from the unchanged "
                       "frozen training data (NOT the unavailable historical "
                       "fitted objects)"),
            "no_leakage": ("no observed CRISPR value enters any prediction or "
                           "ranking; evidence retrieval runs only after the "
                           "top-N is frozen and never changes it"),
            "determinism": ("fixed-precision floats (predictions %d dp, aggregate "
                            "%d dp), sorted keys, fixed list order, no wall-clock, "
                            "no absolute paths, byte-identical on regeneration"
                            % (PRED_DP, AGG_DP)),
            "top_n_source": "config.TOP_N_DEPENDENCIES",
            "top_n": config.TOP_N_DEPENDENCIES,
        },
        "input_artifact_sha256": inp.input_sha256,
        "reconstructed_models": models,
        "samples": samples,
        "rankings": rankings,
        "drug_gene_interaction_evidence": evidence_block,
        "osteosarcoma_validation_aggregate": osteo,
        "limitations": [
            "Both models are RECONSTRUCTED fitted state, not the historical "
            "Phase 1 fitted objects (which were never serialised). They "
            "reproduce every committed Phase 1 validation statistic exactly at "
            "the recorded precision, but remain a reconstruction.",
            "ACH-000364 is ONE held-out cell line: a pipeline-verification "
            "example, not a performance measurement. Its observed CRISPR values "
            "are shown for inspection only.",
            "BG003082 is bulk primary-tumour tissue scored by models trained and "
            "validated only on cultured DepMap cell lines -- a real domain "
            "shift. Its prediction is exploratory; no ground truth exists.",
            "BG003082's Geneformer embedding was generated separately; its "
            "commensurability with the Phase 1 training embeddings is not "
            "proven (see samples.BG003082.head_input.commensurability_caveats).",
            "BG003082 baseline input resolves %d / 18,460 canonical genes; the "
            "remaining 33 are represented as missing and receive the artifact's "
            "stored training-mean imputation."
            % sp_prov["reconciliation"]["canonical_genes_mapped"],
            "The osteosarcoma aggregate is descriptive and unstable (n=5); it is "
            "not a confidence-bounded result and does not replace the frozen "
            "Phase 1 primary comparison.",
            "Drug-gene interaction evidence is retrieval from a licence-filtered "
            "offline DGIdb snapshot: not treatments, not recommendations, no "
            "efficacy claim. PMIDs are drug-gene / source-group level citations.",
        ],
        "disclaimers": [
            config.DGIDB_EVIDENCE_DISCLAIMER,
            "Ranked genes are predicted dependencies, not therapeutic targets or "
            "recommended drugs.",
            "This artifact does not predict patient treatment response or "
            "recommend any therapy.",
        ],
    }
    return artifact


# --------------------------------------------------------------------------
# build / validate / self-test
# --------------------------------------------------------------------------

def build(out_path: Path = OUT_FILE, *, verbose: bool = True) -> dict:
    art = assemble()
    _write_json_deterministic(art, Path(out_path))
    if verbose:
        size = Path(out_path).stat().st_size
        print(f"  wrote {_rel(out_path)}  ({size:,} bytes)  "
              f"sha256 {_sha256_file(out_path)}")
    return art


_PROTECTED = {
    "data/processed/baseline_results.json":
        "b49169bd363a596f400b4faff8c21d354275b70404efe08b9109d38f1bdc0ffd",
    "data/processed/head_results.json":
        "1962206fa17646cbd1fec4b642a577cc2586c09c4cabd980541a7e11a8b6f894",
    "data/processed/analysis_results.json":
        "12431dad60d07f0bd2bea9a680367007c9e030e9f17c5c20ef0b0694dcb548f9",
    "data/processed/expression.npz":
        "3d5bfa0c3430584f8943fd2365be0eecf8b994b38bfc7d491d59d7b9ff251a2d",
    "data/processed/crispr_effect.npz":
        "9214efa3ce172079e6ce4ca78853d8bf92fb8f6d4a55d0c6c71e4653b59e8826",
    "data/processed/geneformer_embeddings.csv":
        "af8ee6d734bea11101d07884f1c72d2b4efaff9875506738a037102a712f1e46",
    "data/processed/splits.json":
        "f1419abc7cbd31efc173a5857bab9eb318b53f8e535a17048bfcf0ea2f70aeef",
    "data/processed/selective_genes.json":
        "68c8fe39ae8965ce20b04f50870609cc21734386ceeff859f4d0bddd2e5bab35",
    "data/processed/gene_columns.json":
        "a4b8069cc93af48f01e745bb1a15f4eaf4a7b67c9f92ca44bef3bb9e44c6d0a1",
    "data/processed/geneformer_bg003082_embedding.csv":
        "06a4ab9f85e5ac908975268ed502912317503ed277d28eeab1663d8305835080",
    "data/external/dgidb/dgidb_2026-06b.interactions.filtered.tsv":
        "f7d2089facc17ddac01e422cab8dc89d48aae463573094490f04bc42ef0a0bee",
    "data/external/dgidb/dgidb_2026-06b.manifest.json":
        "9fb585c723cb2102a7cd335dbfac478b206d91cad04951f8ca7f70f495f6f912",
    "data/external/sid_osteosarc/BG003082.gene_tpm.gct.gz":
        "652011a1cdb8ecf42812cc5fcd6c55947a77995ebe47893e3b307165467bb711",
    "data/processed/reconstructed_fitted/baseline_ridge_pca/manifest.json":
        "133af1d12442775a3d16e223380b753650bdfc86445a8f69aaf650f5399b4efe",
    "data/processed/reconstructed_fitted/head_ridge_head/manifest.json":
        "3fceecc9faec1320048a00e972d4b5a38d7cc3ffaed6ff83dd928e98fa182a05",
}

def _assert_no_forbidden_source() -> None:
    """
    AST check (robust to strings / comments): case_study.py must not import
    scikit-learn or reconstruct_fitted, and must not call .fit() / .fit_transform().
    """
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    bad: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                if a.name.split(".")[0] in ("sklearn", "reconstruct_fitted"):
                    bad.append(f"import {a.name}")
        elif isinstance(node, ast.ImportFrom):
            if (node.module or "").split(".")[0] in ("sklearn", "reconstruct_fitted"):
                bad.append(f"from {node.module} import ...")
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            if node.func.attr in ("fit", "fit_transform"):
                bad.append(f".{node.func.attr}() call")
    if bad:
        raise CaseStudyError(f"case_study.py contains forbidden constructs: {bad}")


def _structural_checks(art: dict) -> list[tuple[str, bool, str]]:
    c: list[tuple[str, bool, str]] = []

    def chk(name, ok, detail=""):
        c.append((name, bool(ok), detail))

    chk("schema_version", art.get("schema_version") == SCHEMA_VERSION)
    chk("source_commit recorded", art.get("source_commit") == SOURCE_COMMIT)

    anchor, tumor = config.DEMO_VERIFICATION_MODEL_ID, config.DEMO_TUMOR_SAMPLE_ID
    sa, st = art["samples"][anchor], art["samples"][tumor]
    chk("ACH-000364 split == val", sa["depmap_split"] == "val", sa["depmap_split"])
    chk("ACH-000364 not in training split", sa["in_training_split"] is False)
    chk("ACH-000364 prediction_status", sa["prediction_status"] == config.PREDICTION_STATUS_HELD_OUT)
    chk("ACH-000364 outcome_status", sa["outcome_status"] == config.OUTCOME_STATUS_MEASURED)
    chk("BG003082 analysis_role", st["analysis_role"] == "exploratory_external_prediction")
    chk("BG003082 outcome_status unavailable", st["outcome_status"] == config.OUTCOME_STATUS_UNAVAILABLE)
    chk("BG003082 absent from all splits", st["absent_from_all_depmap_splits"] is True)
    chk("BG003082 has no observed fields anywhere",
        "observed_crispr" not in st and not _contains_observed(art["rankings"][tumor]))
    chk("BG003082 commensurability caveats present",
        len(st["head_input"]["commensurability_caveats"]) >= 4)

    for mk in ("ridge_pca", "ridge_head"):
        mb = art["reconstructed_models"][mk]
        chk(f"{mk} labelled as reconstructed",
            mb["provenance_status"] == RECONSTRUCTION_STATUS)
        chk(f"{mk} frozen alpha read-not-selected",
            mb["frozen_alpha"]["selection"].startswith("read verbatim"))
    chk("ridge_pca alpha == 100000.0",
        art["reconstructed_models"]["ridge_pca"]["frozen_alpha"]["value"] == 100000.0)
    chk("ridge_head alpha == 3162.0",
        art["reconstructed_models"]["ridge_head"]["frozen_alpha"]["value"] == 3162.0)

    # ranking direction + tie-break + top-N
    n = config.TOP_N_DEPENDENCIES
    for sample, sb in art["rankings"].items():
        for mk, mb in sb.items():
            genes = mb["genes"]
            ok_len = len(genes) == n and mb["n_targets_ranked"] == 4297
            vals = [g["predicted_geneeffect"] for g in genes]
            ok_sorted = all(vals[i] <= vals[i + 1] for i in range(len(vals) - 1))
            ok_ties = all(
                not (vals[i] == vals[i + 1]
                     and int(genes[i]["entrez_id"]) > int(genes[i + 1]["entrez_id"]))
                for i in range(len(vals) - 1))
            ok_rank = [g["rank"] for g in genes] == list(range(1, n + 1))
            chk(f"{sample}/{mk}: top-{n}, ascending, Entrez tie-break, ranks 1..N",
                ok_len and ok_sorted and ok_ties and ok_rank)
    # anchor observed attached; tumor not
    for mk, mb in art["rankings"][anchor].items():
        chk(f"{anchor}/{mk}: observed attached after ranking",
            mb.get("observed_values_attached_after_ranking") is True
            and all("observed_geneeffect" in g for g in mb["genes"]))
    for mk, mb in art["rankings"][tumor].items():
        chk(f"{tumor}/{mk}: no observed field on any gene",
            all("observed_geneeffect" not in g and "observed_rank" not in g
                for g in mb["genes"]))

    ev = art["drug_gene_interaction_evidence"]
    chk("evidence label is 'Drug-gene interaction evidence'",
        ev["label"] == EVIDENCE_LABEL)
    chk("evidence retrieved after top-N frozen",
        ev["retrieval"]["retrieved_after_top_n_frozen"] is True
        and ev["retrieval"]["evidence_availability_did_not_affect_selection_or_ranking"] is True)
    disp = set()
    for sb in art["rankings"].values():
        for mb in sb.values():
            disp |= {g["entrez_id"] for g in mb["genes"]}
    chk("every displayed gene has an evidence entry (even if empty)",
        set(ev["by_entrez"]) == disp,
        f"{len(ev['by_entrez'])} entries / {len(disp)} displayed")
    cov = ev["coverage"]
    chk("evidence coverage sums to distinct genes",
        cov["n_cited"] + cov["n_source_only"] + cov["n_none_in_filtered_snapshot"]
        == cov["n_distinct_genes"] == len(disp))
    statuses = {b["evidence_status"] for b in ev["by_entrez"].values()}
    chk("evidence_status vocabulary",
        statuses <= {"cited", "source_only", "none_in_filtered_snapshot"})
    # uncited records labelled source-only; cited records carry the scope note
    lab_ok = True
    for b in ev["by_entrez"].values():
        for r in b["records"]:
            if r["pmids"] and r["pmid_status"] != "cited":
                lab_ok = False
            if not r["pmids"] and r["pmid_status"] != "source_only":
                lab_ok = False
            if r["pmids"] and "group level" not in r["pmid_scope_note"]:
                lab_ok = False
    chk("PMID / source-only labelling + multi-claim scope note", lab_ok)
    forbidden_words = ("treatment", "recommend", "actionable therapy", "efficacious")
    text = json.dumps(ev).lower()
    # 'no efficacy claim' / 'not treatments' etc. are allowed; check we never
    # assert the positive. crude but catches an accidental framing slip.
    chk("evidence section never asserts treatment/recommendation",
        "recommended drug" not in text and "candidate treatment" not in text
        and "actionable therapy\"" not in text)

    o = art["osteosarcoma_validation_aggregate"]
    chk("osteo cohort n==5 with 5 model ids",
        o["cohort"]["n"] == 5 and len(o["cohort"]["model_ids"]) == 5
        and o["cohort"]["model_ids"] == sorted(o["cohort"]["model_ids"]))
    chk("osteo models are exactly ridge_pca + ridge_head",
        o["models"] == ["ridge_pca", "ridge_head"])
    chk("osteo target universe 4297",
        o["target_universe"] == 4297)
    chk("osteo common set accounting",
        o["common_finite_target_set"]["n_included"]
        + o["common_finite_target_set"]["n_excluded"] == 4297
        and o["common_finite_target_set"]["n_included"] > 0)
    chk("osteo delta == head - pca",
        round(o["mean_per_target_spearman"]["ridge_head"]
              - o["mean_per_target_spearman"]["ridge_pca"], AGG_DP)
        == o["delta_ridge_head_minus_ridge_pca"])
    chk("osteo labelled descriptive/unstable, not a replacement",
        "DESCRIPTIVE" in o["status"] and "n=5" in o["status"]
        and "does NOT replace" in o["status"]
        and o["used_to_choose_model_or_alter_rankings"] is False)

    # no absolute path, no volatile timestamp
    blob = json.dumps(art)
    chk("no absolute local path in artifact",
        "C:\\" not in blob and "C:/" not in blob
        and str(config.PROJECT_ROOT) not in blob
        and str(config.PROJECT_ROOT).replace("\\", "/") not in blob)
    chk("no volatile timestamp key",
        not any(k in art for k in ("generated_at", "timestamp", "run_utc", "wall_clock")))
    return c


def _contains_observed(obj) -> bool:
    blob = json.dumps(obj).lower()
    return "observed_geneeffect" in blob or "observed_rank" in blob or "observed_crispr" in blob


def validate(out_path: Path = OUT_FILE, *, verbose: bool = True) -> dict:
    if not Path(out_path).is_file():
        raise CaseStudyError(f"{out_path} not found -- run --build first")

    # protected artifacts unchanged
    prot_ok = True
    for rel, want in sorted(_PROTECTED.items()):
        got = _sha256_file(_REPO_ROOT / rel)
        if got != want:
            prot_ok = False
            if verbose:
                print(f"  [FAIL] protected artifact changed: {rel}")

    # byte-identical regeneration
    tmp = Path(tempfile.mkdtemp(prefix="case_study_val_"))
    try:
        regen = tmp / "case_study.json"
        build(regen, verbose=False)
        committed = Path(out_path).read_bytes()
        fresh = regen.read_bytes()
        byte_ok = committed == fresh
        # second regeneration -> still identical
        regen2 = tmp / "case_study2.json"
        build(regen2, verbose=False)
        byte_ok = byte_ok and regen2.read_bytes() == fresh
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    art = json.loads(Path(out_path).read_text(encoding="utf-8"))
    checks = _structural_checks(art)
    checks.append(("protected artifacts unchanged", prot_ok, ""))
    checks.append(("committed JSON == fresh regeneration (twice)", byte_ok, ""))
    checks.append(("no sklearn / fit / fit_transform / reconstruct_fitted in source",
                   _no_forbidden_ok(), ""))

    n_fail = sum(1 for _, ok, _ in checks if not ok)
    if verbose:
        for name, ok, detail in checks:
            print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""))
        print()
        print(f"  {'ALL CHECKS PASSED' if n_fail == 0 else str(n_fail) + ' FAILED'}"
              f"  ({len(checks) - n_fail}/{len(checks)})")
    return {"checks": checks, "n_fail": n_fail}


def _no_forbidden_ok() -> bool:
    try:
        _assert_no_forbidden_source()
        return True
    except CaseStudyError:
        return False


def _self_test() -> int:
    print("Running case_study.py self-test...")

    _assert_no_forbidden_source()
    print("  [ok] source contains no sklearn / fit / fit_transform / reconstruct_fitted")

    # ---- ranking direction + Entrez tie-break (pure logic) -----------
    pred = np.array([0.5, -1.0, -1.0, -0.2, -3.0])
    entrez = [50, 30, 20, 10, 99]
    order = _full_ranking(pred, entrez)
    assert order == [4, 2, 1, 3, 0], order          # -3.0 ; (-1.0,e20) ; (-1.0,e30) ; -0.2 ; 0.5
    print("  [ok] ranking: most-negative first, ascending numeric Entrez tie-break")

    # ---- observed-rank tie convention -------------------------------
    obs = np.array([np.nan, -2.0, -2.0, 0.1])
    lut = _observed_rank_lookup(obs, [7, 4, 2, 9])
    assert lut == {2: 1, 1: 2, 3: 3}, lut           # nan -> absent; (-2,e2)=1 (-2,e4)=2 (0.1)=3
    print("  [ok] observed rank: NaN -> null, deterministic (value, Entrez) order")

    # ---- evidence status classification ---------------------------
    assert _shape_evidence_record({
        "pmids": "", **_MIN_REC})["pmid_status"] == "source_only"
    assert _shape_evidence_record({
        "pmids": "9;5", **_MIN_REC})["pmid_status"] == "cited"
    assert _shape_evidence_record({"pmids": "9;5", **_MIN_REC})["pmids"] == ["9", "5"]
    print("  [ok] evidence record shaping: cited vs source-only, PMID list split")

    # ---- real end-to-end build, structural checks, byte-identical ----
    tmp = Path(tempfile.mkdtemp(prefix="case_study_selftest_"))
    try:
        a1 = build(tmp / "a.json", verbose=False)
        b1 = build(tmp / "b.json", verbose=False)
        assert (tmp / "a.json").read_bytes() == (tmp / "b.json").read_bytes(), \
            "case_study.json not byte-identical across regenerations"
        print("  [ok] JSON regenerates byte-identically")
        checks = _structural_checks(a1)
        bad = [name for name, ok, _ in checks if not ok]
        assert not bad, f"structural checks failed: {bad}"
        print(f"  [ok] {len(checks)} structural checks pass on a fresh build")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print("\nSelf-test passed.")
    return 0


_MIN_REC = {
    "entrez_id": "1", "gene_symbol": "X", "dgidb_gene_name": "X",
    "gene_symbol_consistent": "true", "drug_name": "D", "drug_concept_id": "d:1",
    "drug_claim_name": "D", "interaction_source": "CIViC",
    "interaction_source_version": "v", "source_license": "CC0",
    "source_license_url": "u", "interaction_type_raw": "inhibitor",
    "interaction_direction": "inhibitory", "direction_tier": "inhibitory",
    "interaction_score": "1", "drug_specificity_score": "", "gene_specificity_score": "",
    "evidence_score": "", "drug_is_approved": "true", "drug_is_immunotherapy": "false",
    "drug_is_antineoplastic": "false", "dgidb_release_tag": "2026-06b",
    "record_key": "k", "disclaimer": config.DGIDB_EVIDENCE_DISCLAIMER,
}


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description="Phase 2 case-study orchestration.")
    ap.add_argument("--build", action="store_true", help="Write case_study.json.")
    ap.add_argument("--validate", action="store_true",
                    help="Regenerate, byte-compare, re-check invariants.")
    ap.add_argument("--self-test", action="store_true", help="Offline unit + build checks.")
    ap.add_argument("--out", default=str(OUT_FILE))
    args = ap.parse_args()

    if args.self_test:
        return _self_test()

    rc = 0
    did = False
    if args.build:
        did = True
        print("=" * 74)
        print("BUILD case_study.json")
        print("=" * 74)
        build(Path(args.out))
    if args.validate:
        did = True
        print("\n" + "=" * 74)
        print("VALIDATE case_study.json")
        print("=" * 74)
        rep = validate(Path(args.out))
        rc = rc or (0 if rep["n_fail"] == 0 else 1)
    if not did:
        ap.print_help()
        return 2
    return rc


if __name__ == "__main__":
    sys.exit(main())
