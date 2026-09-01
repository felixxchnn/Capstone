import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { resolveStructure, defaultSelection } from '@/data/providers/structureProvider';

// Deterministic fake for the three structure services. No live network.
const UNIPROT_HIT = {
  results: [
    {
      primaryAccession: 'P24941',
      uniProtkbId: 'CDK2_HUMAN',
      entryType: 'UniProtKB reviewed (Swiss-Prot)',
      proteinDescription: { recommendedName: { fullName: { value: 'Cyclin-dependent kinase 2' } } },
      genes: [{ geneName: { value: 'CDK2' } }],
      sequence: { length: 298 },
      organism: { taxonId: 9606 },
    },
  ],
};
const RCSB_SEARCH = {
  total_count: 2,
  result_set: [
    { identifier: '1AQ1', score: 1 },
    { identifier: '6GUE', score: 1 },
  ],
};
const RCSB_ENTRY = (id: string) => ({
  struct: { title: `Structure ${id}` },
  exptl: [{ method: 'X-RAY DIFFRACTION' }],
  rcsb_entry_info: { resolution_combined: [id === '1AQ1' ? 2.0 : 1.6], experimental_method: 'X-ray' },
  rcsb_accession_info: { initial_release_date: '1997-11-12T00:00:00.000+00:00' },
  rcsb_primary_citation: { pdbx_database_id_PubMed: 9334743, pdbx_database_id_DOI: '10.1000/x' },
});
const AF_PREDICTION = [
  {
    entryId: 'AF-P24941-F1',
    modelEntityId: 'AF-P24941-F1',
    uniprotAccession: 'P24941',
    latestVersion: 6,
    globalMetricValue: 88.44,
    fractionPlddtVeryHigh: 0.69,
    fractionPlddtConfident: 0.19,
    modelCreatedDate: '2025-08-01T00:00:00Z',
    cifUrl: 'https://alphafold.ebi.ac.uk/files/AF-P24941-F1-model_v6.cif',
    pdbUrl: 'https://alphafold.ebi.ac.uk/files/AF-P24941-F1-model_v6.pdb',
    bcifUrl: 'https://alphafold.ebi.ac.uk/files/AF-P24941-F1-model_v6.bcif',
    paeImageUrl: 'https://alphafold.ebi.ac.uk/files/AF-P24941-F1-pae_v6.png',
  },
];

function jsonResponse(body: unknown, status = 200): Response {
  return new Response(JSON.stringify(body), {
    status,
    headers: { 'content-type': 'application/json' },
  });
}

let fetchMock: ReturnType<typeof vi.fn>;

beforeEach(() => {
  fetchMock = vi.fn(async (input: RequestInfo | URL) => {
    const url = String(input);
    if (url.includes('rest.uniprot.org')) return jsonResponse(UNIPROT_HIT);
    if (url.includes('search.rcsb.org')) return jsonResponse(RCSB_SEARCH);
    if (url.includes('data.rcsb.org/rest/v1/core/entry/')) {
      const id = url.split('/').pop()!;
      return jsonResponse(RCSB_ENTRY(id));
    }
    if (url.includes('alphafold.ebi.ac.uk/api/prediction/')) return jsonResponse(AF_PREDICTION);
    throw new Error(`unexpected fetch ${url}`);
  });
  vi.stubGlobal('fetch', fetchMock);
});
afterEach(() => vi.unstubAllGlobals());

describe('structure provider chain', () => {
  it('resolves Entrez -> UniProt -> experimental candidates -> AlphaFold', async () => {
    const res = await resolveStructure('1017', 'CDK2');
    expect(res.uniprot?.accession).toBe('P24941');
    expect(res.experimental).toHaveLength(2);
    // sorted by resolution ascending (1AQ1 2.0, 6GUE 1.6 -> RCSB search order kept)
    expect(res.experimental[0].pdbId).toBe('1AQ1');
    expect(res.experimental[0].method).toBe('X-RAY DIFFRACTION');
    expect(res.experimental[0].resolutionAngstrom).toBe(2.0);
    expect(res.predicted?.alphafoldId).toBe('AF-P24941-F1');
    expect(res.predicted?.meanPlddt).toBe(88.44);
    // model URLs come from the API response, not guessed
    expect(res.predicted?.modelUrl).toBe(AF_PREDICTION[0].cifUrl);
    expect(res.steps.map((s) => s.status)).toEqual(['ok', 'ok', 'ok']);
  });

  it('default selection prefers an experimental structure', async () => {
    const res = await resolveStructure('1017', 'CDK2');
    const sel = defaultSelection(res);
    expect(sel?.kind).toBe('experimental');
  });

  it('falls back to AlphaFold when there is no experimental structure', async () => {
    fetchMock.mockImplementation(async (input: RequestInfo | URL) => {
      const url = String(input);
      if (url.includes('rest.uniprot.org')) return jsonResponse(UNIPROT_HIT);
      if (url.includes('search.rcsb.org')) return jsonResponse({ total_count: 0, result_set: [] });
      if (url.includes('alphafold.ebi.ac.uk/api/prediction/')) return jsonResponse(AF_PREDICTION);
      throw new Error(`unexpected ${url}`);
    });
    const res = await resolveStructure('1017', 'CDK2');
    expect(res.experimental).toHaveLength(0);
    const sel = defaultSelection(res);
    expect(sel?.kind).toBe('predicted');
    expect(res.steps.find((s) => s.id === 'rcsb')?.status).toBe('empty');
  });

  it('reports a mapping failure and skips downstream steps', async () => {
    fetchMock.mockImplementation(async (input: RequestInfo | URL) => {
      const url = String(input);
      if (url.includes('rest.uniprot.org')) return jsonResponse({ results: [] });
      throw new Error(`should not be called: ${url}`);
    });
    const res = await resolveStructure('999999', 'FAKE');
    expect(res.uniprot).toBeNull();
    expect(res.steps.find((s) => s.id === 'uniprot')?.status).toBe('empty');
    expect(res.steps.find((s) => s.id === 'rcsb')?.status).toBe('skipped');
    expect(res.steps.find((s) => s.id === 'alphafold')?.status).toBe('skipped');
    expect(defaultSelection(res)).toBeNull();
  });

  it('surfaces an upstream API error without throwing', async () => {
    fetchMock.mockImplementation(async (input: RequestInfo | URL) => {
      const url = String(input);
      if (url.includes('rest.uniprot.org')) return jsonResponse(UNIPROT_HIT);
      if (url.includes('search.rcsb.org')) return new Response('nope', { status: 503 });
      if (url.includes('alphafold.ebi.ac.uk')) return new Response('nope', { status: 500 });
      throw new Error(`unexpected ${url}`);
    });
    const res = await resolveStructure('1017', 'CDK2');
    expect(res.steps.find((s) => s.id === 'rcsb')?.status).toBe('error');
    expect(res.steps.find((s) => s.id === 'alphafold')?.status).toBe('error');
    // still returns; the panel shows the error, the rest of the app is fine
    expect(res.uniprot?.accession).toBe('P24941');
  });

  it('aborts in-flight requests when the signal is already aborted', async () => {
    const ctrl = new AbortController();
    ctrl.abort();
    const res = await resolveStructure('1017', 'CDK2', ctrl.signal);
    // uniprot step ran (or errored) but nothing downstream — aborted early
    expect(res.experimental).toHaveLength(0);
    expect(res.predicted).toBeNull();
  });
});
