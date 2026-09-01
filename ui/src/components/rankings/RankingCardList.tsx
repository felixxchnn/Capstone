import type { GeneRowView } from '@/hooks/useCaseStudyViews';
import { EvidenceStatusPill } from '@/components/evidence/EvidenceStatusPill';
import { geneEffect, integer } from '@/lib/format';
import { EmptyState } from '@/components/common/primitives';

/** Mobile presentation of the ranking. Cards, but still in frozen rank order. */
export function RankingCardList({
  rows,
  hasObserved,
  onSelectGene,
  onResetFilters,
  filtering,
}: {
  rows: GeneRowView[];
  hasObserved: boolean;
  onSelectGene: (entrezId: string) => void;
  onResetFilters: () => void;
  filtering: boolean;
}) {
  const visible = rows.filter((r) => r.visible);
  if (visible.length === 0) {
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
        The frozen top-25 ranking is unchanged.
      </EmptyState>
    );
  }
  return (
    <ol className="rank-cards" aria-label="Ranked predicted dependencies">
      {rows.map((row) => {
        const g = row.gene;
        return (
          <li key={g.entrez_id} className="rank-card" hidden={!row.visible} value={g.rank}>
            <div className="rank-card__head">
              <span className="rank-chip">{g.rank}</span>
              <button type="button" className="linklike rank-card__sym" onClick={() => onSelectGene(g.entrez_id)}>
                {g.symbol}
              </button>
              <span className="mono muted tiny">{g.entrez_id}</span>
            </div>
            <dl className="rank-card__kv">
              <div>
                <dt>Predicted</dt>
                <dd className="tnum">{geneEffect(g.predicted_geneeffect)}</dd>
              </div>
              {hasObserved ? (
                <>
                  <div>
                    <dt>Observed</dt>
                    <dd className="tnum">
                      {g.observed_geneeffect === undefined ? 'n/m' : geneEffect(g.observed_geneeffect)}
                    </dd>
                  </div>
                  <div>
                    <dt>Obs. rank</dt>
                    <dd className="tnum">{g.observed_rank === undefined ? '—' : integer(g.observed_rank)}</dd>
                  </div>
                </>
              ) : null}
            </dl>
            <div className="rank-card__foot">
              <EvidenceStatusPill status={row.evidenceStatus} />
              <button
                type="button"
                className="btn btn--ghost btn--xs"
                onClick={() => onSelectGene(g.entrez_id)}
              >
                Details
              </button>
            </div>
          </li>
        );
      })}
    </ol>
  );
}
