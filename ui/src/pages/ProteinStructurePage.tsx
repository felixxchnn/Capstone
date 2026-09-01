import { useMemo } from 'react';
import { useSelection } from '@/hooks/useSelection';
import { useDataSource } from '@/app/dataSourceContext';
import { useAsync } from '@/hooks/useAsync';
import { StructureExplorer } from '@/components/structures/StructureExplorer';
import { SegmentedControl } from '@/components/common/SegmentedControl';
import { Callout, EmptyState, Spinner } from '@/components/common/primitives';
import { STRUCTURE_NOT_EVIDENCE } from '@/lib/format';
import type { ModelId, SampleId } from '@/types/caseStudy';

export function ProteinStructurePage() {
  const { selection, setSelection } = useSelection();
  const ds = useDataSource();

  const ranking = useAsync(
    (_s) => ds.getModelRanking(selection.sample, selection.model),
    [ds, selection.sample, selection.model],
  );

  const genes = useMemo(
    () => (ranking.status === 'success' ? ranking.data.genes : []),
    [ranking.status, ranking.data],
  );
  const activeEntrez = selection.gene ?? genes[0]?.entrez_id ?? null;
  const activeGene = useMemo(
    () => genes.find((g) => g.entrez_id === activeEntrez) ?? null,
    [genes, activeEntrez],
  );

  return (
    <div className="container page">
      <header className="stack">
        <h1>Protein Structure</h1>
        <p className="lede">
          The experimental or predicted structure of the <strong>protein encoded by</strong>{' '}
          the selected gene — retrieved live from UniProt, RCSB PDB and AlphaFold DB. This is
          structural evidence about the protein, not the physical structure of the gene, and
          it is <strong>separate from prediction</strong>: nothing here changes a ranking.
        </p>
        <Callout tone="caution" title="Structure is not drug-response evidence">
          {STRUCTURE_NOT_EVIDENCE}
        </Callout>
      </header>

      <div className="structure-picker no-print">
        <SegmentedControl<SampleId>
          legend="Sample"
          name="structure-sample"
          value={selection.sample}
          onChange={(s) => setSelection({ sample: s, gene: null })}
          options={[
            { value: 'ACH-000364', label: 'ACH-000364' },
            { value: 'BG003082', label: 'BG003082' },
          ]}
        />
        <SegmentedControl<ModelId>
          legend="Model ranking"
          name="structure-model"
          value={selection.model}
          onChange={(m) => setSelection({ model: m, gene: null })}
          options={[
            { value: 'ridge_pca', label: 'ridge_pca' },
            { value: 'ridge_head', label: 'ridge_head' },
          ]}
        />
        <div className="field">
          <label htmlFor="structure-gene">Gene (from the top-25 ranking)</label>
          <select
            id="structure-gene"
            value={activeEntrez ?? ''}
            onChange={(e) => setSelection({ gene: e.target.value })}
          >
            {genes.map((g) => (
              <option key={g.entrez_id} value={g.entrez_id}>
                #{g.rank} · {g.symbol} ({g.entrez_id})
              </option>
            ))}
          </select>
        </div>
      </div>

      {ranking.status === 'loading' ? <Spinner label="Loading gene list" /> : null}
      {ranking.status === 'error' ? (
        <EmptyState title="Could not load the gene list">{ranking.error.message}</EmptyState>
      ) : null}

      {activeGene ? (
        <StructureExplorer entrezId={activeGene.entrez_id} symbol={activeGene.symbol} />
      ) : ranking.status === 'success' ? (
        <EmptyState title="No gene selected" />
      ) : null}
    </div>
  );
}
