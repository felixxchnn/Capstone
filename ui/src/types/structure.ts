// Types for the protein-structure provider chain:
//   Entrez Gene ID -> UniProt (reviewed, human) -> RCSB experimental candidates
//                  -> AlphaFold predicted model (fallback) -> Mol* viewer

export interface UniProtMapping {
  accession: string;
  entryId: string;
  proteinName: string;
  gene: string | null;
  length: number;
  organismId: number;
  reviewed: boolean;
}

export type StructureKind = 'experimental' | 'predicted';

export interface ExperimentalCandidate {
  kind: 'experimental';
  pdbId: string;
  title: string;
  method: string | null;
  /** Ångström; null for methods with no single resolution. */
  resolutionAngstrom: number | null;
  releaseDate: string | null;
  citationPubMedId: string | null;
  citationDoi: string | null;
  /** URL Mol* loads. mmCIF text from RCSB. */
  modelUrl: string;
  modelFormat: 'mmcif';
}

export interface PredictedModel {
  kind: 'predicted';
  source: 'AlphaFold DB';
  alphafoldId: string;
  uniprotAccession: string;
  version: number;
  /** Mean pLDDT (model confidence), 0–100. */
  meanPlddt: number | null;
  fractionVeryHigh: number | null;
  fractionConfident: number | null;
  modelCreatedDate: string | null;
  modelUrl: string;
  modelFormat: 'mmcif' | 'pdb';
  paeImageUrl: string | null;
}

export interface StructureResolution {
  entrezId: string;
  symbol: string;
  humanTaxonomyId: 9606;
  uniprot: UniProtMapping | null;
  experimental: ExperimentalCandidate[];
  predicted: PredictedModel | null;
  /** which provider steps ran + their outcome, for the status panel */
  steps: StructureStep[];
}

export interface StructureStep {
  id: 'uniprot' | 'rcsb' | 'alphafold';
  label: string;
  status: 'ok' | 'empty' | 'error' | 'skipped';
  detail: string;
}

export type StructureSelection =
  | { kind: 'experimental'; candidate: ExperimentalCandidate }
  | { kind: 'predicted'; model: PredictedModel };
