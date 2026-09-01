// RCSB PDB provider — experimental structures for a UniProt accession.
//
// Search API: https://search.rcsb.org/#search-api  (rcsbsearch/v2/query)
// Data API:   https://data.rcsb.org/#data-api      (rest/v1/core/entry/{id})
// Model file: https://files.rcsb.org/download/{ID}.cif  (mmCIF text — Mol* loads this)
//
// We query by UniProt accession + experimental content only; no expression data.

import { getJson } from '@/data/providers/http';
import type { ExperimentalCandidate } from '@/types/structure';

const SEARCH_URL = 'https://search.rcsb.org/rcsbsearch/v2/query';
const ENTRY_URL = 'https://data.rcsb.org/rest/v1/core/entry';

interface SearchResponse {
  total_count?: number;
  result_set?: Array<{ identifier: string; score: number }>;
}

interface EntryResponse {
  struct?: { title?: string };
  exptl?: Array<{ method?: string }>;
  rcsb_entry_info?: {
    resolution_combined?: number[];
    experimental_method?: string;
  };
  rcsb_accession_info?: { initial_release_date?: string };
  rcsb_primary_citation?: {
    pdbx_database_id_PubMed?: number;
    pdbx_database_id_DOI?: string;
  };
}

export interface RcsbLookupResult {
  candidates: ExperimentalCandidate[];
  totalCount: number;
  detail: string;
}

export function rcsbModelUrl(pdbId: string): string {
  return `https://files.rcsb.org/download/${pdbId.toUpperCase()}.cif`;
}

export async function findExperimentalStructures(
  uniprotAccession: string,
  signal: AbortSignal | undefined,
  maxCandidates = 8,
): Promise<RcsbLookupResult> {
  const query = {
    query: {
      type: 'terminal',
      service: 'text',
      parameters: {
        attribute:
          'rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession',
        operator: 'exact_match',
        value: uniprotAccession,
      },
    },
    return_type: 'entry',
    request_options: {
      paginate: { start: 0, rows: maxCandidates },
      // ascending resolution = best experimental structures first;
      // entries with no resolution (e.g. NMR) sort last
      sort: [{ sort_by: 'rcsb_entry_info.resolution_combined', direction: 'asc' }],
      results_content_type: ['experimental'],
    },
  };
  const url = `${SEARCH_URL}?json=${encodeURIComponent(JSON.stringify(query))}`;

  let search: SearchResponse;
  try {
    search = await getJson<SearchResponse>(url, { signal, accept: 'application/json' });
  } catch (err) {
    // RCSB search returns 204 (no content) for zero hits on some deployments.
    if (err && typeof err === 'object' && 'status' in err && (err as { status?: number }).status === 204) {
      return { candidates: [], totalCount: 0, detail: `No experimental PDB entry references ${uniprotAccession}.` };
    }
    throw err;
  }

  const ids = (search.result_set ?? []).map((r) => r.identifier);
  if (ids.length === 0) {
    return {
      candidates: [],
      totalCount: search.total_count ?? 0,
      detail: `No experimental PDB entry references ${uniprotAccession}.`,
    };
  }

  const entries = await Promise.all(
    ids.map(async (id): Promise<ExperimentalCandidate | null> => {
      try {
        const e = await getJson<EntryResponse>(`${ENTRY_URL}/${id}`, {
          signal,
          accept: 'application/json',
        });
        const res = e.rcsb_entry_info?.resolution_combined?.[0];
        return {
          kind: 'experimental',
          pdbId: id.toUpperCase(),
          title: e.struct?.title ?? id.toUpperCase(),
          method: e.exptl?.[0]?.method ?? e.rcsb_entry_info?.experimental_method ?? null,
          resolutionAngstrom: typeof res === 'number' ? res : null,
          releaseDate: e.rcsb_accession_info?.initial_release_date?.slice(0, 10) ?? null,
          citationPubMedId: e.rcsb_primary_citation?.pdbx_database_id_PubMed
            ? String(e.rcsb_primary_citation.pdbx_database_id_PubMed)
            : null,
          citationDoi: e.rcsb_primary_citation?.pdbx_database_id_DOI ?? null,
          modelUrl: rcsbModelUrl(id),
          modelFormat: 'mmcif',
        };
      } catch {
        return null;
      }
    }),
  );

  const candidates = entries.filter((c): c is ExperimentalCandidate => c !== null);
  return {
    candidates,
    totalCount: search.total_count ?? candidates.length,
    detail: `${search.total_count ?? candidates.length} experimental PDB ${
      (search.total_count ?? 0) === 1 ? 'entry' : 'entries'
    } reference ${uniprotAccession}; showing the ${candidates.length} best by resolution.`,
  };
}
