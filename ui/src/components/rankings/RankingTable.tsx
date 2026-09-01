import type { GeneRowView } from '@/hooks/useCaseStudyViews';
import type { ModelId, SampleId } from '@/types/caseStudy';
import { EvidenceStatusPill } from '@/components/evidence/EvidenceStatusPill';
import { geneEffect, geneEffectFull, integer } from '@/lib/format';
import { EmptyState } from '@/components/common/primitives';

/**
 * Ordered dependency table. Rows are always in frozen `rank` order 1..25;
 * `row.visible === false` hides a row with the `hidden` attribute but never
 * changes its position or its `rank`.
 */
export function RankingTable({
  sample,
  model,
  rows,
  hasObserved,
  nTargetsRanked,
  onSelectGene,
  onResetFilters,
  filtering,
}: {
  sample: SampleId;
  model: ModelId;
  rows: GeneRowView[];
  hasObserved: boolean;
  nTargetsRanked: number;
  onSelectGene: (entrezId: string) => void;
  onResetFilters: () => void;
  filtering: boolean;
}) {
  const visibleRows = rows.filter((r) => r.visible);

  if (visibleRows.length === 0) {
    return (
      <EmptyState
        title="No genes match the current search and filters"
        action={
          filtering ? (
            <button type="button" className="btn btn--ghost" onClick={onResetFilters}>
              Reset filters
            </button>
          ) : null
        }
      >
        The frozen top-25 ranking is unchanged — every row is simply hidden by the
        active filter.
      </EmptyState>
    );
  }

  return (
    <div className="table-scroll">
      <table className="data ranking-table">
        <caption>
          Predicted CRISPR gene <strong>dependencies</strong> — sample{' '}
          <span className="mono">{sample}</span>, model <span className="mono">{model}</span>.
          Top 25 of {integer(nTargetsRanked)} ranked targets. Predicted dependencies,
          not therapeutic targets.
        </caption>
        <thead>
          <tr>
            <th scope="col">Rank</th>
            <th scope="col">Gene</th>
            <th scope="col">Entrez</th>
            <th scope="col">
              Predicted GeneEffect
              <span className="th-hint"> (more negative = stronger)</span>
            </th>
            {hasObserved ? (
              <>
                <th scope="col">Observed GeneEffect</th>
                <th scope="col">Observed rank / 4,297</th>
              </>
            ) : null}
            <th scope="col">Evidence</th>
            <th scope="col">
              <span className="sr-only">Actions</span>
            </th>
          </tr>
        </thead>
        <tbody>
          {rows.map((row) => {
            const g = row.gene;
            return (
              <tr key={g.entrez_id} hidden={!row.visible} data-rank={g.rank}>
                <td className="num">
                  <span className="rank-chip">{g.rank}</span>
                </td>
                <td>
                  <button
                    type="button"
                    className="linklike"
                    onClick={() => onSelectGene(g.entrez_id)}
                  >
                    {g.symbol}
                  </button>
                </td>
                <td className="num mono">{g.entrez_id}</td>
                <td className="num" title={geneEffectFull(g.predicted_geneeffect)}>
                  {geneEffect(g.predicted_geneeffect)}
                </td>
                {hasObserved ? (
                  <>
                    <td className="num" title={geneEffectFull(g.observed_geneeffect)}>
                      {g.observed_geneeffect === undefined
                        ? 'not measured'
                        : geneEffect(g.observed_geneeffect)}
                    </td>
                    <td className="num">
                      {g.observed_rank === undefined ? '—' : integer(g.observed_rank)}
                    </td>
                  </>
                ) : null}
                <td>
                  <EvidenceStatusPill status={row.evidenceStatus} />{' '}
                  <span className="muted tiny">
                    {row.evidenceCount} record{row.evidenceCount === 1 ? '' : 's'}
                  </span>
                </td>
                <td>
                  <button
                    type="button"
                    className="btn btn--ghost btn--xs"
                    onClick={() => onSelectGene(g.entrez_id)}
                    aria-label={`Open details for ${g.symbol}`}
                  >
                    Details
                  </button>
                </td>
              </tr>
            );
          })}
        </tbody>
      </table>
    </div>
  );
}
