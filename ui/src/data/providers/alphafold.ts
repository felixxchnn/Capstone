// AlphaFold DB provider — the predicted-model fallback for a UniProt accession.
//
// API: https://alphafold.ebi.ac.uk/api-docs  ->  /api/prediction/{accession}
// The response carries the current `cifUrl` / `pdbUrl` (versioned) and
// `globalMetricValue` (mean pLDDT). We NEVER guess file URLs — we read them from
// the API response, as the task requires.

import { getJson } from '@/data/providers/http';
import type { PredictedModel } from '@/types/structure';

const API = 'https://alphafold.ebi.ac.uk/api/prediction';

interface AlphaFoldPrediction {
  entryId?: string;
  modelEntityId?: string;
  uniprotAccession?: string;
  latestVersion?: number;
  globalMetricValue?: number;
  fractionPlddtVeryHigh?: number;
  fractionPlddtConfident?: number;
  modelCreatedDate?: string;
  cifUrl?: string;
  pdbUrl?: string;
  bcifUrl?: string;
  paeImageUrl?: string;
}

export interface AlphaFoldLookupResult {
  model: PredictedModel | null;
  detail: string;
}

export async function findAlphaFoldModel(
  uniprotAccession: string,
  signal?: AbortSignal,
): Promise<AlphaFoldLookupResult> {
  const url = `${API}/${encodeURIComponent(uniprotAccession)}`;
  const arr = await getJson<AlphaFoldPrediction[]>(url, { signal, accept: 'application/json' });
  const p = Array.isArray(arr) ? arr[0] : undefined;
  if (!p || (!p.cifUrl && !p.pdbUrl)) {
    return {
      model: null,
      detail: `No AlphaFold DB model for ${uniprotAccession}.`,
    };
  }
  const modelUrl = p.cifUrl ?? (p.pdbUrl as string);
  const model: PredictedModel = {
    kind: 'predicted',
    source: 'AlphaFold DB',
    alphafoldId: p.modelEntityId ?? p.entryId ?? `AF-${uniprotAccession}-F1`,
    uniprotAccession: p.uniprotAccession ?? uniprotAccession,
    version: p.latestVersion ?? 0,
    meanPlddt: typeof p.globalMetricValue === 'number' ? p.globalMetricValue : null,
    fractionVeryHigh: p.fractionPlddtVeryHigh ?? null,
    fractionConfident: p.fractionPlddtConfident ?? null,
    modelCreatedDate: p.modelCreatedDate?.slice(0, 10) ?? null,
    modelUrl,
    modelFormat: p.cifUrl ? 'mmcif' : 'pdb',
    paeImageUrl: p.paeImageUrl ?? null,
  };
  return {
    model,
    detail: `AlphaFold DB predicted model ${model.alphafoldId} (v${model.version}), mean pLDDT ${
      model.meanPlddt !== null ? model.meanPlddt.toFixed(1) : 'n/a'
    }.`,
  };
}
