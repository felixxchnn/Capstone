// Scientific formatting helpers. Kept in one module so wording stays consistent
// and correct everywhere (sign conventions, rank phrasing, "predicted" vs
// "observed", "not a therapeutic target").

import type { EvidenceStatus } from '@/types/caseStudy';

/** Fixed display precision matching the committed artifact (predictions 10 dp). */
export function geneEffect(value: number | null | undefined): string {
  if (value === null || value === undefined || !Number.isFinite(value)) return 'n/a';
  return value.toFixed(4);
}

export function geneEffectFull(value: number | null | undefined): string {
  if (value === null || value === undefined || !Number.isFinite(value)) return 'n/a';
  // up to 10 dp, trailing zeros trimmed
  return value
    .toFixed(10)
    .replace(/(\.\d*?)0+$/, '$1')
    .replace(/\.$/, '');
}

export function signedDelta(value: number, dp = 4): string {
  if (!Number.isFinite(value)) return 'n/a';
  const s = value.toFixed(dp);
  return value > 0 ? `+${s}` : s;
}

export function integer(value: number | null | undefined): string {
  if (value === null || value === undefined || !Number.isFinite(value)) return 'n/a';
  return Math.round(value).toLocaleString('en-US');
}

export function spearman(value: number, dp = 4): string {
  return value.toFixed(dp);
}

export const DEPENDENCY_DIRECTION =
  'More-negative predicted CRISPR GeneEffect = stronger predicted dependency.';

export const NOT_A_TARGET =
  'This is a predicted dependency — not a therapeutic target and not a drug recommendation.';

export const NON_CLINICAL =
  'Research demonstration — not clinical guidance.';

export const EVIDENCE_NOT_EFFICACY =
  'A recorded drug–gene interaction is retrieval of prior evidence. It does not establish efficacy for this sample or for osteosarcoma.';

export const STRUCTURE_NOT_EVIDENCE =
  'A protein structure (experimental or predicted) is structural evidence about the protein a gene encodes. It does not prove function, drug response, or therapeutic relevance, and it is not drug-response evidence.';

export function evidenceStatusLabel(status: EvidenceStatus): string {
  switch (status) {
    case 'cited':
      return 'Cited evidence';
    case 'source_only':
      return 'Source-only evidence';
    case 'none_in_filtered_snapshot':
      return 'No record in the filtered snapshot';
  }
}

export function evidenceStatusShort(status: EvidenceStatus): string {
  switch (status) {
    case 'cited':
      return 'cited';
    case 'source_only':
      return 'source-only';
    case 'none_in_filtered_snapshot':
      return 'none';
  }
}

export function pmidUrl(pmid: string): string {
  const clean = pmid.replace(/[^0-9]/g, '');
  return `https://pubmed.ncbi.nlm.nih.gov/${clean}/`;
}

/** yes/true/1 -> "yes"; no/false/0 -> "no"; blank -> "not stated" */
export function flag(value: string | undefined): string {
  const v = (value ?? '').trim().toLowerCase();
  if (v === '' ) return 'not stated';
  if (['true', 'yes', '1'].includes(v)) return 'yes';
  if (['false', 'no', '0'].includes(v)) return 'no';
  return value ?? 'not stated';
}

export function shortSha(sha: string, head = 8): string {
  return sha.length > head ? `${sha.slice(0, head)}…` : sha;
}
