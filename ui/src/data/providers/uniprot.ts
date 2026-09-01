// UniProt provider — map an Entrez Gene ID to a reviewed human UniProtKB entry.
//
// Endpoint: UniProt REST search (programmatic access documented at
// https://www.uniprot.org/help/api_queries and https://rest.uniprot.org).
// Query: xref:GeneID-<entrez> AND organism_id:9606 AND reviewed:true
// We send the Entrez ID + human taxonomy only. No expression / prediction data.

import { getJson } from '@/data/providers/http';
import type { UniProtMapping } from '@/types/structure';

const BASE = 'https://rest.uniprot.org/uniprotkb/search';
const FIELDS = 'accession,id,protein_name,gene_primary,length,organism_id,reviewed';

interface UniProtSearchResponse {
  results: Array<{
    primaryAccession: string;
    uniProtkbId: string;
    entryType?: string;
    proteinDescription?: { recommendedName?: { fullName?: { value?: string } } };
    genes?: Array<{ geneName?: { value?: string } }>;
    sequence?: { length?: number };
    organism?: { taxonId?: number };
  }>;
}

export interface UniProtLookupResult {
  mapping: UniProtMapping | null;
  /** all reviewed human hits, best (canonical) first — usually 1 */
  candidates: UniProtMapping[];
  detail: string;
}

export async function lookupUniProtByEntrez(
  entrezId: string,
  signal?: AbortSignal,
): Promise<UniProtLookupResult> {
  const query = `xref:GeneID-${encodeURIComponent(entrezId)} AND organism_id:9606 AND reviewed:true`;
  const url = `${BASE}?query=${encodeURIComponent(query)}&fields=${FIELDS}&format=json&size=5`;
  const data = await getJson<UniProtSearchResponse>(url, { signal, accept: 'application/json' });

  const candidates: UniProtMapping[] = (data.results ?? []).map((r) => ({
    accession: r.primaryAccession,
    entryId: r.uniProtkbId,
    proteinName: r.proteinDescription?.recommendedName?.fullName?.value ?? r.uniProtkbId,
    gene: r.genes?.[0]?.geneName?.value ?? null,
    length: r.sequence?.length ?? 0,
    organismId: r.organism?.taxonId ?? 9606,
    reviewed: (r.entryType ?? '').includes('reviewed'),
  }));

  if (candidates.length === 0) {
    return {
      mapping: null,
      candidates: [],
      detail: `No reviewed human UniProtKB entry cross-references Entrez Gene ${entrezId}.`,
    };
  }
  // Prefer the shortest reviewed entry only as a stable tie-break; UniProt already
  // returns the canonical entry first.
  return {
    mapping: candidates[0],
    candidates,
    detail: `Mapped Entrez ${entrezId} → UniProt ${candidates[0].accession} (${candidates[0].proteinName}).`,
  };
}
