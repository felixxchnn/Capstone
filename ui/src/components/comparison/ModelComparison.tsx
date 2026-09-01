import { useMemo } from 'react';
import type { ModelRankingBlock, SampleId } from '@/types/caseStudy';
import { geneEffect, signedDelta } from '@/lib/format';
import { Callout } from '@/components/common/primitives';

interface OverlapRow {
  symbol: string;
  entrez_id: string;
  rankPca: number;
  rankHead: number;
  rankDelta: number;
  predPca: number;
  predHead: number;
  predDelta: number;
}

export function ModelComparison({
  sample,
  pca,
  head,
}: {
  sample: SampleId;
  pca: ModelRankingBlock;
  head: ModelRankingBlock;
}) {
  const { overlap, onlyPca, onlyHead } = useMemo(() => {
    const headByEntrez = new Map(head.genes.map((g) => [g.entrez_id, g]));
    const pcaEntrez = new Set(pca.genes.map((g) => g.entrez_id));
    const rows: OverlapRow[] = [];
    for (const g of pca.genes) {
      const h = headByEntrez.get(g.entrez_id);
      if (!h) continue;
      rows.push({
        symbol: g.symbol,
        entrez_id: g.entrez_id,
        rankPca: g.rank,
        rankHead: h.rank,
        rankDelta: h.rank - g.rank,
        predPca: g.predicted_geneeffect,
        predHead: h.predicted_geneeffect,
        predDelta: h.predicted_geneeffect - g.predicted_geneeffect,
      });
    }
    rows.sort((a, b) => a.rankPca - b.rankPca);
    return {
      overlap: rows,
      onlyPca: pca.genes.filter((g) => !headByEntrez.has(g.entrez_id)),
      onlyHead: head.genes.filter((g) => !pcaEntrez.has(g.entrez_id)),
    };
  }, [pca, head]);

  return (
    <div className="model-comparison stack-lg">
      <Callout tone="info" title="Descriptive UI comparison only">
        Each model keeps its own frozen ranking, rank 1 → 25. The overlap and rank-difference
        figures below are descriptive comparisons computed in the browser — they are{' '}
        <strong>not new performance evaluations</strong> and the two lists are never merged
        into a consensus.
      </Callout>

      <div className="model-comparison__lists grid grid--2">
        {[
          { id: 'ridge_pca', block: pca, title: 'ridge_pca — top 25' },
          { id: 'ridge_head', block: head, title: 'ridge_head — top 25' },
        ].map(({ id, block, title }) => (
          <section key={id} aria-label={title} className="table-scroll">
            <table className="data compact">
              <caption>{title}</caption>
              <thead>
                <tr>
                  <th scope="col">#</th>
                  <th scope="col">Gene</th>
                  <th scope="col">Entrez</th>
                  <th scope="col">Predicted</th>
                </tr>
              </thead>
              <tbody>
                {block.genes.map((g) => {
                  const inOther =
                    id === 'ridge_pca'
                      ? head.genes.some((h) => h.entrez_id === g.entrez_id)
                      : pca.genes.some((p) => p.entrez_id === g.entrez_id);
                  return (
                    <tr key={g.entrez_id} className={inOther ? 'is-shared' : undefined}>
                      <td className="num">
                        <span className="rank-chip">{g.rank}</span>
                      </td>
                      <td>
                        {g.symbol}
                        {inOther ? <span className="shared-dot" title="Also in the other model's top 25" /> : null}
                      </td>
                      <td className="num mono">{g.entrez_id}</td>
                      <td className="num">{geneEffect(g.predicted_geneeffect)}</td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </section>
        ))}
      </div>

      <section aria-label="Gene overlap">
        <h3>
          Genes in both top-25 lists: <strong>{overlap.length}</strong> of 25
        </h3>
        {overlap.length === 0 ? (
          <p className="small muted">
            No gene appears in both models' top 25 for {sample}. The two representations
            surface entirely different predicted dependencies here.
          </p>
        ) : (
          <div className="table-scroll">
            <table className="data compact">
              <caption>
                Shared genes — rank &amp; predicted-value differences (ridge_head − ridge_pca).
                Descriptive only.
              </caption>
              <thead>
                <tr>
                  <th scope="col">Gene</th>
                  <th scope="col">Entrez</th>
                  <th scope="col">ridge_pca rank</th>
                  <th scope="col">ridge_head rank</th>
                  <th scope="col">Rank Δ</th>
                  <th scope="col">Predicted Δ</th>
                </tr>
              </thead>
              <tbody>
                {overlap.map((r) => (
                  <tr key={r.entrez_id}>
                    <td>{r.symbol}</td>
                    <td className="num mono">{r.entrez_id}</td>
                    <td className="num">{r.rankPca}</td>
                    <td className="num">{r.rankHead}</td>
                    <td className="num">{signedDelta(r.rankDelta, 0)}</td>
                    <td className="num">{signedDelta(r.predDelta)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        )}
      </section>

      <div className="grid grid--2">
        <section aria-label="Only in ridge_pca">
          <h3>Only in ridge_pca ({onlyPca.length})</h3>
          <ul className="chip-list">
            {onlyPca.map((g) => (
              <li key={g.entrez_id} className="chip">
                {g.symbol} <span className="mono muted tiny">#{g.rank}</span>
              </li>
            ))}
          </ul>
        </section>
        <section aria-label="Only in ridge_head">
          <h3>Only in ridge_head ({onlyHead.length})</h3>
          <ul className="chip-list">
            {onlyHead.map((g) => (
              <li key={g.entrez_id} className="chip">
                {g.symbol} <span className="mono muted tiny">#{g.rank}</span>
              </li>
            ))}
          </ul>
        </section>
      </div>
    </div>
  );
}
