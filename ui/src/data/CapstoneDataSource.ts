// The stable seam between the UI and its scientific data.
//
// Today the only implementation is `StaticCaseStudyDataSource`, which reads the
// committed `case_study.json`. A future `ApiDataSource` (contract documented in
// `adapters/ApiDataSource.contract.md`) can replace it without touching a single
// UI component: every page and component consumes this interface, never the raw
// JSON.
//
// NOTE: there is deliberately no `runInference` / `predict` method. Predictions
// are already frozen in the committed artifact; the UI shows them, it does not
// produce them. A future backend would serve the SAME committed predictions plus
// provenance, not compute new ones in response to the UI.

import type {
  CaseStudy,
  EvidenceBucket,
  ModelId,
  ModelRankingBlock,
  SampleAnchor,
  SampleExternal,
  SampleId,
} from '@/types/caseStudy';

export interface ProjectSummary {
  title: string;
  description: string;
  schemaVersion: string;
  sourceCommit: string;
  caseStudySha256: string;
  /** Frozen Phase 1 headline (values live in the committed artifact / docs). */
  phase1: {
    nValidationLines: number;
    nTargets: number;
    ridgePcaMeanSpearman: number;
    ridgeHeadMeanSpearman: number;
    deltaHeadMinusBaseline: number;
    conclusion: string;
  };
  osteosarcoma: CaseStudy['osteosarcoma_validation_aggregate'];
  evidenceCoverage: CaseStudy['drug_gene_interaction_evidence']['coverage'];
  models: CaseStudy['reconstructed_models'];
  limitations: string[];
  disclaimers: string[];
  environment: Record<string, string>;
  inputArtifactSha256: Record<string, string>;
}

export interface SampleSummary {
  id: SampleId;
  role: string;
  predictionStatus: string;
  outcomeStatus: string;
  hasObservedOutcome: boolean;
  isExternal: boolean;
  splitStatus: string;
  context: string;
}

export type SampleMetadata =
  | { id: 'ACH-000364'; kind: 'anchor'; data: SampleAnchor }
  | { id: 'BG003082'; kind: 'external'; data: SampleExternal };

/** Descriptive-only structure hint derived from the committed data (no network). */
export interface StructureMetadataHint {
  entrezId: string;
  symbol: string;
  humanTaxonomyId: '9606';
  note: string;
}

export interface CapstoneDataSource {
  /** Human/plain-language project framing plus frozen headline numbers. */
  getProjectSummary(): Promise<ProjectSummary>;

  /** Both samples, each with its explicit, non-equivalent role. */
  getSamples(): Promise<SampleSummary[]>;

  /** One model's independent, frozen top-25 ranking for one sample. Never merged. */
  getModelRanking(sample: SampleId, model: ModelId): Promise<ModelRankingBlock>;

  /** Drug–gene interaction evidence for one Entrez ID (retrieved after ranking). */
  getGeneEvidence(entrezId: string): Promise<EvidenceBucket | null>;

  /** Full provenance / mapping / imputation / domain-shift metadata for a sample. */
  getSampleMetadata(sample: SampleId): Promise<SampleMetadata>;

  /** Identifiers the structure providers need. No expression/prediction data. */
  getStructureMetadata(entrezId: string): Promise<StructureMetadataHint>;

  /** Escape hatch for views that need the whole immutable artifact (read-only). */
  getRawCaseStudy(): Promise<Readonly<CaseStudy>>;
}
