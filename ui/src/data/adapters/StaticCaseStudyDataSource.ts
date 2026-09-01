// The current CapstoneDataSource implementation: reads the committed
// `case_study.json` (synced verbatim to `src/data/case_study.generated.json`).
// Pure, synchronous under the hood, wrapped in Promises so the interface matches
// a future network adapter.

import generated from '@/data/case_study.generated.json';
import expectedShaRaw from '@/data/case_study.expected-sha256.txt?raw';
import { parseCaseStudy } from '@/data/schemas/caseStudy';
import type {
  CapstoneDataSource,
  ProjectSummary,
  SampleMetadata,
  SampleSummary,
  StructureMetadataHint,
} from '@/data/CapstoneDataSource';
import type { CaseStudy, EvidenceBucket, ModelId, ModelRankingBlock, SampleId } from '@/types/caseStudy';

const PHASE1_CONCLUSION =
  'The Geneformer head did not outperform the expression / PCA baseline (delta -0.0308, 95% CI [-0.0365, -0.0255]).';

export interface StaticDataSourceOptions {
  /**
   * The SHA-256 the synced copy is expected to have (mirrors
   * capstone/data-integrity-hashes.md). Defaults to the pinned value bundled
   * with the app. Tests override this to prove the guard bites.
   */
  expectedSha256?: string;
}

export class StaticCaseStudyDataSource implements CapstoneDataSource {
  private readonly caseStudy: CaseStudy;

  constructor(options: StaticDataSourceOptions = {}) {
    this.caseStudy = parseCaseStudy(generated);
    const expected = (options.expectedSha256 ?? expectedShaRaw).trim();
    // The build-time sync writes case_study.sha256.txt; the value is also
    // pinned here. We can't hash the imported object at runtime deterministically
    // (key order), so the byte-level SHA check lives in
    // src/tests/caseStudyLoad.test.ts against the on-disk copy. Here we assert
    // only that a pinned hash exists and is well formed.
    if (!/^[0-9a-f]{64}$/.test(expected)) {
      throw new Error(
        `StaticCaseStudyDataSource: pinned case-study SHA-256 is malformed (${JSON.stringify(expected)})`,
      );
    }
  }

  async getRawCaseStudy(): Promise<Readonly<CaseStudy>> {
    return this.caseStudy;
  }

  async getProjectSummary(): Promise<ProjectSummary> {
    const cs = this.caseStudy;
    return {
      title: cs.title,
      description: cs.description,
      schemaVersion: cs.schema_version,
      sourceCommit: cs.source_commit,
      caseStudySha256: expectedShaRaw.trim(),
      phase1: {
        nValidationLines: 170,
        nTargets: 4297,
        ridgePcaMeanSpearman: 0.2356,
        ridgeHeadMeanSpearman: 0.2047,
        deltaHeadMinusBaseline: -0.0308,
        conclusion: PHASE1_CONCLUSION,
      },
      osteosarcoma: cs.osteosarcoma_validation_aggregate,
      evidenceCoverage: cs.drug_gene_interaction_evidence.coverage,
      models: cs.reconstructed_models,
      limitations: cs.limitations,
      disclaimers: cs.disclaimers,
      environment: cs.environment,
      inputArtifactSha256: cs.input_artifact_sha256,
    };
  }

  async getSamples(): Promise<SampleSummary[]> {
    const anchor = this.caseStudy.samples['ACH-000364'];
    const external = this.caseStudy.samples.BG003082;
    return [
      {
        id: 'ACH-000364',
        role: anchor.role,
        predictionStatus: anchor.prediction_status,
        outcomeStatus: anchor.outcome_status,
        hasObservedOutcome: true,
        isExternal: false,
        splitStatus: `DepMap ${anchor.depmap_split} split (held out; ${anchor.split_assertion})`,
        context: `${anchor.cell_line} — a cultured osteosarcoma cell line used only as a pipeline-verification anchor.`,
      },
      {
        id: 'BG003082',
        role: external.role,
        predictionStatus: external.prediction_status,
        outcomeStatus: external.outcome_status,
        hasObservedOutcome: false,
        isExternal: true,
        splitStatus: 'Absent from every DepMap split (train / val / test).',
        context:
          'A real primary osteosarcoma tumour (bulk RNA-seq, CC0). Bulk tumour tissue is a real domain shift from the cultured cell lines the models were trained and validated on.',
      },
    ];
  }

  async getModelRanking(sample: SampleId, model: ModelId): Promise<ModelRankingBlock> {
    const block = this.caseStudy.rankings[sample]?.[model];
    if (!block) {
      throw new Error(`No committed ranking for sample=${sample} model=${model}`);
    }
    // Defensive copy so a consumer cannot mutate the immutable source.
    return { ...block, genes: block.genes.map((g) => ({ ...g })) };
  }

  async getGeneEvidence(entrezId: string): Promise<EvidenceBucket | null> {
    const bucket = this.caseStudy.drug_gene_interaction_evidence.by_entrez[entrezId];
    if (!bucket) return null;
    return { ...bucket, records: bucket.records.map((r) => ({ ...r })) };
  }

  async getSampleMetadata(sample: SampleId): Promise<SampleMetadata> {
    if (sample === 'ACH-000364') {
      return { id: 'ACH-000364', kind: 'anchor', data: this.caseStudy.samples['ACH-000364'] };
    }
    return { id: 'BG003082', kind: 'external', data: this.caseStudy.samples.BG003082 };
  }

  async getStructureMetadata(entrezId: string): Promise<StructureMetadataHint> {
    // Find the symbol from any ranking / evidence entry that mentions this Entrez.
    let symbol = entrezId;
    for (const sid of ['ACH-000364', 'BG003082'] as SampleId[]) {
      for (const mid of ['ridge_pca', 'ridge_head'] as ModelId[]) {
        const hit = this.caseStudy.rankings[sid][mid].genes.find((g) => g.entrez_id === entrezId);
        if (hit) symbol = hit.symbol;
      }
    }
    const ev = this.caseStudy.drug_gene_interaction_evidence.by_entrez[entrezId];
    if (ev) symbol = ev.symbol;
    return {
      entrezId,
      symbol,
      humanTaxonomyId: '9606',
      note: 'Structure services are queried with the Entrez Gene ID and human taxonomy only. No expression or prediction data is sent.',
    };
  }
}
