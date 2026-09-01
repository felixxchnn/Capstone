// The app's shared selection state (sample / model / gene / search / evidence
// filter), serialised to and from the URL query string so navigation preserves
// it and views are shareable/bookmarkable.

import type { EvidenceStatus, ModelId, SampleId } from '@/types/caseStudy';
import { MODEL_IDS, SAMPLE_IDS } from '@/types/caseStudy';

export type EvidenceFilter = 'all' | EvidenceStatus;

export interface Selection {
  sample: SampleId;
  model: ModelId;
  /** Selected gene Entrez ID for the detail drawer / structure page, or null. */
  gene: string | null;
  search: string;
  evidence: EvidenceFilter;
}

export const DEFAULT_SELECTION: Selection = {
  sample: 'ACH-000364',
  model: 'ridge_pca',
  gene: null,
  search: '',
  evidence: 'all',
};

const EVIDENCE_VALUES: EvidenceFilter[] = [
  'all',
  'cited',
  'source_only',
  'none_in_filtered_snapshot',
];

function isSample(v: string | null): v is SampleId {
  return v !== null && (SAMPLE_IDS as readonly string[]).includes(v);
}
function isModel(v: string | null): v is ModelId {
  return v !== null && (MODEL_IDS as readonly string[]).includes(v);
}
function isEvidence(v: string | null): v is EvidenceFilter {
  return v !== null && (EVIDENCE_VALUES as string[]).includes(v);
}

export function selectionFromParams(params: URLSearchParams): Selection {
  const sample = params.get('sample');
  const model = params.get('model');
  const evidence = params.get('evidence');
  const gene = params.get('gene');
  return {
    sample: isSample(sample) ? sample : DEFAULT_SELECTION.sample,
    model: isModel(model) ? model : DEFAULT_SELECTION.model,
    gene: gene && /^[0-9]+$/.test(gene) ? gene : null,
    search: params.get('q') ?? '',
    evidence: isEvidence(evidence) ? evidence : 'all',
  };
}

export function selectionToParams(sel: Selection, base?: URLSearchParams): URLSearchParams {
  const p = new URLSearchParams(base ?? undefined);
  const set = (key: string, value: string, isDefault: boolean) => {
    if (isDefault) p.delete(key);
    else p.set(key, value);
  };
  set('sample', sel.sample, sel.sample === DEFAULT_SELECTION.sample);
  set('model', sel.model, sel.model === DEFAULT_SELECTION.model);
  set('evidence', sel.evidence, sel.evidence === 'all');
  set('q', sel.search, sel.search.trim() === '');
  if (sel.gene) p.set('gene', sel.gene);
  else p.delete('gene');
  return p;
}

/** Case-insensitive haystack match against a gene row + its evidence. */
export function matchesSearch(
  needle: string,
  parts: { symbol: string; entrez: string; drugs?: string[]; sources?: string[] },
): boolean {
  const q = needle.trim().toLowerCase();
  if (!q) return true;
  const hay = [
    parts.symbol,
    parts.entrez,
    ...(parts.drugs ?? []),
    ...(parts.sources ?? []),
  ]
    .join(' ')
    .toLowerCase();
  return hay.includes(q);
}

export function isFiltering(sel: Selection): boolean {
  return sel.search.trim() !== '' || sel.evidence !== 'all';
}
