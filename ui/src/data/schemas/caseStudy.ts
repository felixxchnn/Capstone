// Runtime validation for the committed case study. Hand-written (no schema-lib
// dependency) so the UI stays dependency-light. Throws `CaseStudyParseError`
// with a path on the first structural problem.

import type {
  CaseStudy,
  ModelRankingBlock,
  RankedGene,
  EvidenceBucket,
} from '@/types/caseStudy';
import { MODEL_IDS, SAMPLE_IDS } from '@/types/caseStudy';

export class CaseStudyParseError extends Error {
  constructor(
    message: string,
    readonly path: string,
  ) {
    super(`case_study.json invalid at ${path}: ${message}`);
    this.name = 'CaseStudyParseError';
  }
}

function req(cond: unknown, path: string, message: string): asserts cond {
  if (!cond) throw new CaseStudyParseError(message, path);
}

function isRecord(v: unknown): v is Record<string, unknown> {
  return typeof v === 'object' && v !== null && !Array.isArray(v);
}

function num(v: unknown, path: string): number {
  req(typeof v === 'number' && Number.isFinite(v), path, `expected finite number, got ${typeof v}`);
  return v as number;
}

function str(v: unknown, path: string): string {
  req(typeof v === 'string', path, `expected string, got ${typeof v}`);
  return v as string;
}

function parseGene(v: unknown, path: string, expectObserved: boolean): RankedGene {
  req(isRecord(v), path, 'expected object');
  const g: RankedGene = {
    rank: num(v.rank, `${path}.rank`),
    symbol: str(v.symbol, `${path}.symbol`),
    entrez_id: str(v.entrez_id, `${path}.entrez_id`),
    predicted_geneeffect: num(v.predicted_geneeffect, `${path}.predicted_geneeffect`),
  };
  if (expectObserved) {
    // observed values may legitimately be null for a target with no measurement;
    // when present they must be finite numbers.
    if (v.observed_geneeffect !== null && v.observed_geneeffect !== undefined) {
      g.observed_geneeffect = num(v.observed_geneeffect, `${path}.observed_geneeffect`);
    }
    if (v.observed_rank !== null && v.observed_rank !== undefined) {
      g.observed_rank = num(v.observed_rank, `${path}.observed_rank`);
    }
  } else {
    req(
      v.observed_geneeffect === undefined && v.observed_rank === undefined,
      path,
      'external-sample ranking must not carry observed values',
    );
  }
  return g;
}

function parseRankingBlock(v: unknown, path: string, expectObserved: boolean): ModelRankingBlock {
  req(isRecord(v), path, 'expected object');
  const genes = v.genes;
  req(Array.isArray(genes), `${path}.genes`, 'expected array');
  req(genes.length === 25, `${path}.genes`, `expected exactly 25 genes, got ${genes.length}`);
  const parsed = genes.map((g, i) => parseGene(g, `${path}.genes[${i}]`, expectObserved));
  // frozen rank order 1..25, strictly increasing, one-to-one
  parsed.forEach((g, i) => {
    req(g.rank === i + 1, `${path}.genes[${i}].rank`, `expected rank ${i + 1}, got ${g.rank}`);
  });
  // predicted GeneEffect must be non-decreasing with rank (more negative = stronger)
  for (let i = 1; i < parsed.length; i += 1) {
    req(
      parsed[i].predicted_geneeffect >= parsed[i - 1].predicted_geneeffect,
      `${path}.genes[${i}].predicted_geneeffect`,
      'predicted GeneEffect must be ascending with rank',
    );
  }
  return {
    model: str(v.model, `${path}.model`),
    model_provenance: str(v.model_provenance, `${path}.model_provenance`),
    n_displayed: num(v.n_displayed, `${path}.n_displayed`),
    n_targets_ranked: num(v.n_targets_ranked, `${path}.n_targets_ranked`),
    ranking_rule: str(v.ranking_rule, `${path}.ranking_rule`),
    not_a_recommendation: str(v.not_a_recommendation, `${path}.not_a_recommendation`),
    genes: parsed,
    ...(typeof v.observed_rank_rule === 'string' ? { observed_rank_rule: v.observed_rank_rule } : {}),
    ...(typeof v.observed_values_attached_after_ranking === 'boolean'
      ? { observed_values_attached_after_ranking: v.observed_values_attached_after_ranking }
      : {}),
    ...(typeof v.n_targets_with_observed_value === 'number'
      ? { n_targets_with_observed_value: v.n_targets_with_observed_value }
      : {}),
  };
}

function parseEvidenceBucket(v: unknown, path: string): EvidenceBucket {
  req(isRecord(v), path, 'expected object');
  const status = v.evidence_status;
  req(
    status === 'cited' || status === 'source_only' || status === 'none_in_filtered_snapshot',
    `${path}.evidence_status`,
    `unknown status ${JSON.stringify(status)}`,
  );
  const records = v.records;
  req(Array.isArray(records), `${path}.records`, 'expected array');
  return {
    entrez_id: str(v.entrez_id, `${path}.entrez_id`),
    symbol: str(v.symbol, `${path}.symbol`),
    evidence_status: status,
    n_records: num(v.n_records, `${path}.n_records`),
    records: records as EvidenceBucket['records'],
  };
}

export function parseCaseStudy(input: unknown): CaseStudy {
  req(isRecord(input), '$', 'expected top-level object');
  req(
    input.schema_version === 'case-study/1',
    '$.schema_version',
    `expected "case-study/1", got ${JSON.stringify(input.schema_version)}`,
  );

  const rankings = input.rankings;
  req(isRecord(rankings), '$.rankings', 'expected object');
  for (const sid of SAMPLE_IDS) {
    const perSample = rankings[sid];
    req(isRecord(perSample), `$.rankings.${sid}`, 'missing sample');
    for (const mid of MODEL_IDS) {
      parseRankingBlock(
        perSample[mid],
        `$.rankings.${sid}.${mid}`,
        sid === 'ACH-000364',
      );
    }
  }

  const dge = input.drug_gene_interaction_evidence;
  req(isRecord(dge), '$.drug_gene_interaction_evidence', 'expected object');
  const byEntrez = dge.by_entrez;
  req(isRecord(byEntrez), '$.drug_gene_interaction_evidence.by_entrez', 'expected object');
  for (const [k, v] of Object.entries(byEntrez)) {
    parseEvidenceBucket(v, `$.drug_gene_interaction_evidence.by_entrez.${k}`);
  }

  req(isRecord(input.samples), '$.samples', 'expected object');
  req(isRecord(input.reconstructed_models), '$.reconstructed_models', 'expected object');
  req(isRecord(input.osteosarcoma_validation_aggregate), '$.osteosarcoma_validation_aggregate', 'expected object');
  req(Array.isArray(input.limitations), '$.limitations', 'expected array');
  req(Array.isArray(input.disclaimers), '$.disclaimers', 'expected array');

  // Structure checked; the compile-time type asserts the rest of the shape.
  return input as unknown as CaseStudy;
}
