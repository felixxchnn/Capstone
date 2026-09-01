// Live structure-integration smoke test (NOT a unit test — hits real APIs).
//
//   node scripts/structure-smoke.mjs [ENTREZ_ID ...]
//
// For each gene it runs the real provider chain (UniProt -> RCSB -> AlphaFold)
// and prints: Entrez ID, symbol, UniProt accession, #experimental candidates,
// the chosen structure source, and whether a model file URL resolved (HTTP 200).
// It does NOT bake any of this volatile metadata into the deterministic Vitest
// suite. Requires network access.
//
// Default genes are real protein-coding genes displayed in the committed
// case_study.json rankings.

import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const here = fileURLToPath(new URL('.', import.meta.url));
const cs = JSON.parse(
  readFileSync(resolve(here, '..', 'src', 'data', 'case_study.generated.json'), 'utf8'),
);

// pull a few Entrez IDs actually shown in the ranking tables
const displayed = new Set();
for (const sample of Object.values(cs.rankings)) {
  for (const model of Object.values(sample)) {
    for (const g of model.genes) displayed.add(`${g.entrez_id}\t${g.symbol}`);
  }
}
const argv = process.argv.slice(2);
const targets = argv.length
  ? argv.map((id) => {
      const hit = [...displayed].find((d) => d.startsWith(id + '\t'));
      return hit ?? `${id}\t(not in ranking)`;
    })
  : ['1017\tCDK2', '79693\tYRDC', '1663\tDDX11', '9134\tCCNE2'].filter((t) =>
      displayed.has(t) || argv.length === 0,
    );

const UA = 'capstone-research-interface structure-smoke (+https://github.com/felixxchnn/Capstone)';

async function jget(url) {
  const r = await fetch(url, { headers: { 'User-Agent': UA, Accept: 'application/json' } });
  if (!r.ok) throw new Error(`HTTP ${r.status} ${url}`);
  return r.json();
}
async function headOk(url) {
  try {
    const r = await fetch(url, { method: 'GET', headers: { 'User-Agent': UA, Range: 'bytes=0-64' } });
    return r.status === 200 || r.status === 206;
  } catch {
    return false;
  }
}

console.log(`Live structure integration smoke — ${new Date().toISOString()}`);
console.log('(volatile metadata; never asserted in the deterministic test suite)\n');

let failures = 0;
for (const t of targets) {
  const [entrez, symbol] = t.split('\t');
  const row = { entrez, symbol, uniprot: null, experimental: 0, source: 'none', modelUrl: null, rendered: false };
  try {
    const uni = await jget(
      `https://rest.uniprot.org/uniprotkb/search?query=${encodeURIComponent(
        `xref:GeneID-${entrez} AND organism_id:9606 AND reviewed:true`,
      )}&fields=accession&format=json&size=1`,
    );
    row.uniprot = uni.results?.[0]?.primaryAccession ?? null;

    if (row.uniprot) {
      const q = {
        query: {
          type: 'terminal',
          service: 'text',
          parameters: {
            attribute:
              'rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession',
            operator: 'exact_match',
            value: row.uniprot,
          },
        },
        return_type: 'entry',
        request_options: {
          paginate: { start: 0, rows: 5 },
          sort: [{ sort_by: 'rcsb_entry_info.resolution_combined', direction: 'asc' }],
          results_content_type: ['experimental'],
        },
      };
      let search = { result_set: [], total_count: 0 };
      try {
        search = await jget(`https://search.rcsb.org/rcsbsearch/v2/query?json=${encodeURIComponent(JSON.stringify(q))}`);
      } catch {
        /* 204 / no hits */
      }
      row.experimental = search.total_count ?? search.result_set?.length ?? 0;

      if (row.experimental > 0) {
        const id = search.result_set[0].identifier;
        row.source = `experimental PDB ${id}`;
        row.modelUrl = `https://files.rcsb.org/download/${id}.cif`;
      } else {
        const af = await jget(`https://alphafold.ebi.ac.uk/api/prediction/${row.uniprot}`);
        if (af?.[0]?.cifUrl) {
          row.source = `AlphaFold ${af[0].modelEntityId ?? af[0].entryId} (predicted, mean pLDDT ${af[0].globalMetricValue})`;
          row.modelUrl = af[0].cifUrl;
        }
      }
      row.rendered = row.modelUrl ? await headOk(row.modelUrl) : false;
    }
  } catch (err) {
    row.error = String(err);
  }

  const ok = row.uniprot && row.source !== 'none' && row.rendered;
  if (!ok) failures += 1;
  console.log(`  [${ok ? 'ok' : 'FAIL'}] ${symbol} (Entrez ${entrez})`);
  console.log(`        UniProt        : ${row.uniprot ?? 'no mapping'}`);
  console.log(`        experimental   : ${row.experimental}`);
  console.log(`        structure src  : ${row.source}`);
  console.log(`        model file     : ${row.modelUrl ?? '—'}  (fetch 200: ${row.rendered})`);
  if (row.error) console.log(`        error          : ${row.error}`);
}

console.log(`\n${failures === 0 ? 'All targets resolved a structure.' : `${failures} target(s) failed.`}`);
process.exit(failures === 0 ? 0 : 1);
