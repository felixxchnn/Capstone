// Orchestrates the structure provider chain:
//   Entrez Gene ID
//     -> reviewed human UniProt mapping (taxonomy 9606)
//     -> experimental RCSB PDB candidates (listed, user-selectable)
//     -> AlphaFold predicted model (fallback)
//
// Every step is recorded in `steps` for the status panel. Failures in one step
// never abort the others: a missing UniProt entry still lets nothing downstream
// run; a UniProt hit with no PDB still yields the AlphaFold fallback.
//
// Only identifiers (Entrez, UniProt accession, PDB id) + taxonomy 9606 are sent.

import { ProviderError } from '@/data/providers/http';
import { lookupUniProtByEntrez } from '@/data/providers/uniprot';
import { findExperimentalStructures } from '@/data/providers/rcsb';
import { findAlphaFoldModel } from '@/data/providers/alphafold';
import type { StructureResolution, StructureStep } from '@/types/structure';

function stepFromError(id: StructureStep['id'], label: string, err: unknown): StructureStep {
  if (err instanceof ProviderError) {
    if (err.kind === 'offline') {
      return { id, label, status: 'error', detail: 'Upstream service unreachable (offline or blocked).' };
    }
    if (err.kind === 'timeout') {
      return { id, label, status: 'error', detail: 'Upstream service timed out.' };
    }
    return { id, label, status: 'error', detail: err.message };
  }
  return { id, label, status: 'error', detail: err instanceof Error ? err.message : 'Unknown error.' };
}

export async function resolveStructure(
  entrezId: string,
  symbol: string,
  signal?: AbortSignal,
): Promise<StructureResolution> {
  const steps: StructureStep[] = [];
  const result: StructureResolution = {
    entrezId,
    symbol,
    humanTaxonomyId: 9606,
    uniprot: null,
    experimental: [],
    predicted: null,
    steps,
  };

  // --- 1. UniProt -----------------------------------------------------
  try {
    const uni = await lookupUniProtByEntrez(entrezId, signal);
    result.uniprot = uni.mapping;
    steps.push({
      id: 'uniprot',
      label: 'UniProt mapping (reviewed, human)',
      status: uni.mapping ? 'ok' : 'empty',
      detail: uni.detail,
    });
  } catch (err) {
    steps.push(stepFromError('uniprot', 'UniProt mapping (reviewed, human)', err));
  }

  if (signal?.aborted) return result;

  if (!result.uniprot) {
    steps.push({
      id: 'rcsb',
      label: 'RCSB PDB experimental structures',
      status: 'skipped',
      detail: 'Skipped — no UniProt accession to search by.',
    });
    steps.push({
      id: 'alphafold',
      label: 'AlphaFold DB predicted model',
      status: 'skipped',
      detail: 'Skipped — no UniProt accession to look up.',
    });
    return result;
  }

  const acc = result.uniprot.accession;

  // --- 2. RCSB (experimental) + 3. AlphaFold (predicted) in parallel --
  const [rcsbSettled, afSettled] = await Promise.allSettled([
    findExperimentalStructures(acc, signal),
    findAlphaFoldModel(acc, signal),
  ]);

  if (rcsbSettled.status === 'fulfilled') {
    result.experimental = rcsbSettled.value.candidates;
    steps.push({
      id: 'rcsb',
      label: 'RCSB PDB experimental structures',
      status: rcsbSettled.value.candidates.length > 0 ? 'ok' : 'empty',
      detail: rcsbSettled.value.detail,
    });
  } else {
    steps.push(stepFromError('rcsb', 'RCSB PDB experimental structures', rcsbSettled.reason));
  }

  if (afSettled.status === 'fulfilled') {
    result.predicted = afSettled.value.model;
    steps.push({
      id: 'alphafold',
      label: 'AlphaFold DB predicted model',
      status: afSettled.value.model ? 'ok' : 'empty',
      detail: afSettled.value.detail,
    });
  } else {
    steps.push(stepFromError('alphafold', 'AlphaFold DB predicted model', afSettled.reason));
  }

  return result;
}

/** The default selection: best experimental candidate if any, else AlphaFold. */
export function defaultSelection(res: StructureResolution) {
  if (res.experimental.length > 0) {
    return { kind: 'experimental' as const, candidate: res.experimental[0] };
  }
  if (res.predicted) {
    return { kind: 'predicted' as const, model: res.predicted };
  }
  return null;
}
