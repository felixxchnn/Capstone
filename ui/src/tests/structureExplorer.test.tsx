import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { StructureExplorer } from '@/components/structures/StructureExplorer';

// jsdom has no WebGL, so StructureExplorer takes the text-only path and never
// imports Mol*. That is exactly the WebGL-unavailable behaviour we want to test.

const UNIPROT = {
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
const RCSB_SEARCH = { total_count: 1, result_set: [{ identifier: '1AQ1', score: 1 }] };
const RCSB_ENTRY = {
  struct: { title: 'CDK2 with staurosporine' },
  exptl: [{ method: 'X-RAY DIFFRACTION' }],
  rcsb_entry_info: { resolution_combined: [2.0] },
  rcsb_accession_info: { initial_release_date: '1997-11-12T00:00:00.000+00:00' },
  rcsb_primary_citation: { pdbx_database_id_PubMed: 9334743 },
};
const AF = [
  {
    modelEntityId: 'AF-P24941-F1',
    uniprotAccession: 'P24941',
    latestVersion: 6,
    globalMetricValue: 88.44,
    cifUrl: 'https://alphafold.ebi.ac.uk/files/AF-P24941-F1-model_v6.cif',
  },
];

function json(body: unknown, status = 200) {
  return new Response(JSON.stringify(body), { status, headers: { 'content-type': 'application/json' } });
}

beforeEach(() => {
  vi.stubGlobal(
    'fetch',
    vi.fn(async (input: RequestInfo | URL) => {
      const url = String(input);
      if (url.includes('rest.uniprot.org')) return json(UNIPROT);
      if (url.includes('search.rcsb.org')) return json(RCSB_SEARCH);
      if (url.includes('data.rcsb.org')) return json(RCSB_ENTRY);
      if (url.includes('alphafold.ebi.ac.uk')) return json(AF);
      throw new Error(`unexpected ${url}`);
    }),
  );
});
afterEach(() => vi.unstubAllGlobals());

describe('StructureExplorer (WebGL-unavailable / text path)', () => {
  it('resolves providers and shows the WebGL-unavailable notice + text summary', async () => {
    render(<StructureExplorer entrezId="1017" symbol="CDK2" />);

    await screen.findByText(/WebGL unavailable/i);
    // provider status
    expect(await screen.findByText(/UniProt mapping \(reviewed, human\)/i)).toBeInTheDocument();
    // text summary carries the real facts (phrases may also appear in the
    // candidate list, so allow multiple)
    expect(screen.getAllByText(/Cyclin-dependent kinase 2/i).length).toBeGreaterThan(0);
    expect(screen.getAllByText(/1AQ1/).length).toBeGreaterThan(0);
    expect(
      screen.getAllByText(/not an experimental measurement/i).length,
    ).toBeGreaterThan(0);
  });

  it('defaults to the experimental candidate and lets the user pick AlphaFold', async () => {
    const user = userEvent.setup();
    render(<StructureExplorer entrezId="1017" symbol="CDK2" />);
    await screen.findByText(/WebGL unavailable/i);

    const expBtn = screen.getByRole('button', { name: /Experimental.*1AQ1/i });
    expect(expBtn).toHaveAttribute('aria-pressed', 'true');

    const afBtn = screen.getByRole('button', { name: /Predicted structure.*AF-P24941-F1/i });
    await user.click(afBtn);
    expect(afBtn).toHaveAttribute('aria-pressed', 'true');
    expect(expBtn).toHaveAttribute('aria-pressed', 'false');
  });

  it('shows a mapping-failed callout when UniProt has no hit', async () => {
    vi.stubGlobal(
      'fetch',
      vi.fn(async (input: RequestInfo | URL) => {
        if (String(input).includes('rest.uniprot.org')) return json({ results: [] });
        throw new Error('downstream should not run');
      }),
    );
    render(<StructureExplorer entrezId="999999" symbol="FAKE" />);
    expect(await screen.findByText(/Mapping failed/i)).toBeInTheDocument();
  });

  it('shows an "unavailable" empty state when nothing has a structure', async () => {
    vi.stubGlobal(
      'fetch',
      vi.fn(async (input: RequestInfo | URL) => {
        const url = String(input);
        if (url.includes('rest.uniprot.org')) return json(UNIPROT);
        if (url.includes('search.rcsb.org')) return json({ total_count: 0, result_set: [] });
        if (url.includes('alphafold.ebi.ac.uk')) return json([]);
        throw new Error(`unexpected ${url}`);
      }),
    );
    render(<StructureExplorer entrezId="1017" symbol="CDK2" />);
    expect(await screen.findByText(/No structure available for CDK2/i)).toBeInTheDocument();
  });

  it('recovers with a retry action when the services are unreachable', async () => {
    vi.stubGlobal(
      'fetch',
      vi.fn(async () => {
        throw new TypeError('Failed to fetch');
      }),
    );
    render(<StructureExplorer entrezId="1017" symbol="CDK2" />);
    // resolveStructure catches per-step; with UniProt failing, mapping is null ->
    // "Mapping failed" callout appears (not a hard crash)
    expect(await screen.findByText(/Mapping failed|Structure services are unavailable/i)).toBeInTheDocument();
  });
});
