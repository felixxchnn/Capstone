// Domain types for the committed case-study artifact (schema `case-study/1`).
// These mirror data/processed/case_study.json. The UI never mutates this data.

export type SampleId = 'ACH-000364' | 'BG003082';
export type ModelId = 'ridge_pca' | 'ridge_head';
export type EvidenceStatus = 'cited' | 'source_only' | 'none_in_filtered_snapshot';
export type DirectionTier = 'inhibitory' | 'activating' | 'unknown';

export interface RankedGene {
  rank: number;
  symbol: string;
  entrez_id: string;
  predicted_geneeffect: number;
  /** Present only for ACH-000364 (attached AFTER ranking; never used to rank). */
  observed_geneeffect?: number;
  /** Present only for ACH-000364. 1-based position among all 4,297 targets. */
  observed_rank?: number;
}

export interface ModelRankingBlock {
  model: string;
  model_provenance: string;
  n_displayed: number;
  n_targets_ranked: number;
  ranking_rule: string;
  not_a_recommendation: string;
  genes: RankedGene[];
  /** ACH-000364 only. */
  observed_rank_rule?: string;
  observed_values_attached_after_ranking?: boolean;
  n_targets_with_observed_value?: number;
}

export type Rankings = Record<SampleId, Record<ModelId, ModelRankingBlock>>;

export interface EvidenceRecord {
  entrez_id: string;
  gene_symbol: string;
  dgidb_gene_name: string;
  drug_name: string;
  drug_claim_name: string;
  drug_concept_id: string;
  interaction_source: string;
  interaction_source_version: string;
  interaction_type_raw: string;
  interaction_direction: string;
  direction_tier: DirectionTier;
  interaction_score: string;
  evidence_score: string;
  drug_is_approved: string;
  drug_is_antineoplastic: string;
  drug_is_immunotherapy: string;
  curation_type: string;
  indication: string;
  pmids: string[];
  pmid_status: string;
  pmid_scope_note: string;
  source_license: string;
  source_license_url: string;
  dgidb_release_tag: string;
  record_key: string;
  disclaimer: string;
  gene_symbol_consistent: string;
  symbol_query_mismatch: string;
}

export interface EvidenceBucket {
  entrez_id: string;
  symbol: string;
  evidence_status: EvidenceStatus;
  n_records: number;
  records: EvidenceRecord[];
}

export interface EvidenceCoverage {
  n_distinct_genes: number;
  n_cited: number;
  n_source_only: number;
  n_none_in_filtered_snapshot: number;
  total_records: number;
  total_pmid_citations: number;
}

export interface EvidenceRetrieval {
  snapshot_file: string;
  snapshot_sha256: string;
  manifest_sha256: string;
  method: string;
  top_k_per_direction_tier: number;
  direction_tiers: DirectionTier[];
  retrieved_after_top_n_frozen: boolean;
  evidence_availability_did_not_affect_selection_or_ranking: boolean;
}

export interface DrugGeneInteractionEvidence {
  label: string;
  framing: string;
  disclaimer: string;
  pmid_scope_note: string;
  retrieval: EvidenceRetrieval;
  coverage: EvidenceCoverage;
  by_entrez: Record<string, EvidenceBucket>;
}

export interface FrozenAlpha {
  value: number;
  selection: string;
  source: string;
}

export interface ReconstructedModel {
  model: string;
  provenance_status: string;
  not_original_fitted_objects: string;
  artifact_dir: string;
  manifest_sha256: string;
  base_commit: string;
  feature_order_sha256: string;
  target_order_sha256: string;
  frozen_alpha: FrozenAlpha;
  pipeline: string;
  n_features: number;
  n_targets: number;
}

export interface OsteoAggregate {
  status: string;
  definition_source: string;
  cohort: { n: number; model_ids: string[]; predicate: string };
  models: string[];
  target_universe: number;
  per_target_metric: string;
  common_finite_target_set: {
    n_included: number;
    n_excluded: number;
    excluded_ridge_pca_nonfinite: number;
    excluded_ridge_head_nonfinite: number;
    rule: string;
    excluded_reason: string;
  };
  mean_per_target_spearman: { ridge_pca: number; ridge_head: number; rounding_dp: number };
  delta_ridge_head_minus_ridge_pca: number;
  used_to_choose_model_or_alter_rankings: boolean;
}

export interface BaselineReconciliation {
  canonical_genes: number;
  canonical_genes_mapped: number;
  canonical_genes_missing: number;
  canonical_genes_measured_zero: number;
  canonical_genes_measured_nonzero: number;
  canonical_id_collisions: number;
  duplicate_external_ids: number;
  symbol_fallback: string;
  symbol_fallback_candidates: string[];
  [key: string]: unknown;
}

export interface SampleAnchor {
  role: string;
  cell_line: string;
  prediction_status: string;
  outcome_status: string;
  depmap_split: string;
  in_training_split: boolean;
  split_assertion: string;
  baseline_input: {
    source: string;
    n_features: number;
    feature_order_sha256: string;
    missing_features: number;
    imputed_features: number;
  };
  head_input: {
    source: string;
    n_features: number;
    feature_order_sha256: string;
    missing_features: number;
  };
  observed_crispr: {
    source: string;
    role: string;
    n_targets_with_value: number;
    n_targets_missing: number;
  };
}

export interface SampleExternal {
  role: string;
  description?: string;
  analysis_role: string;
  prediction_status: string;
  outcome_status: string;
  absent_from_all_depmap_splits: boolean;
  observed_outcome?: string;
  baseline_input: {
    source: string;
    transformation: string;
    n_features: number;
    feature_order_sha256: string;
    missing_features_represented_as_nan: number;
    imputation: string;
    reconciliation: BaselineReconciliation;
    gct_file?: { name: string; sha256: string; bytes: number };
    ensembl_map_file?: { name: string; sha256: string; bytes: number; rows: number };
    gene_columns_file?: { name: string; sha256: string };
  };
  head_input: {
    sidecar_file?: string;
    sidecar_sha256?: string;
    shape?: [number, number];
    all_finite?: boolean;
    model?: string;
    geneformer_revision_pinned?: string;
    commensurability_caveats?: string[];
    provenance_disclosures?: string[];
  };
}

export interface CaseStudy {
  schema_version: 'case-study/1';
  title: string;
  source_commit: string;
  description: string;
  generated_by: string;
  artifact: string;
  environment: Record<string, string>;
  rankings: Rankings;
  drug_gene_interaction_evidence: DrugGeneInteractionEvidence;
  osteosarcoma_validation_aggregate: OsteoAggregate;
  reconstructed_models: Record<ModelId, ReconstructedModel>;
  samples: {
    'ACH-000364': SampleAnchor;
    BG003082: SampleExternal;
  };
  input_artifact_sha256: Record<string, string>;
  methodology: Record<string, unknown>;
  limitations: string[];
  disclaimers: string[];
}

export const SAMPLE_IDS: readonly SampleId[] = ['ACH-000364', 'BG003082'];
export const MODEL_IDS: readonly ModelId[] = ['ridge_pca', 'ridge_head'];
