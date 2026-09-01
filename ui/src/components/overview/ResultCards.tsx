import type { ProjectSummary } from '@/data/CapstoneDataSource';
import { integer, signedDelta, spearman } from '@/lib/format';

/** Frozen Phase 1 headline cards. Values come from the verified summary — the
 *  UI does not invent constants. */
export function ResultCards({ summary }: { summary: ProjectSummary }) {
  const p = summary.phase1;
  const cards = [
    {
      k: 'Validation cell lines',
      v: integer(p.nValidationLines),
      sub: 'held out; patient-grouped, lineage-stratified split',
    },
    {
      k: 'Predicted targets',
      v: integer(p.nTargets),
      sub: 'selective CRISPR gene-dependency targets',
    },
    {
      k: 'ridge_pca mean Spearman',
      v: spearman(p.ridgePcaMeanSpearman),
      sub: 'ridge on 200 PCA components of expression',
      strong: true,
    },
    {
      k: 'ridge_head mean Spearman',
      v: spearman(p.ridgeHeadMeanSpearman),
      sub: 'ridge on frozen 768-dim Geneformer embeddings',
      strong: true,
    },
    {
      k: 'Delta (head − baseline)',
      v: signedDelta(p.deltaHeadMinusBaseline),
      sub: '95% CI [−0.0365, −0.0255], paired bootstrap',
      accentNegative: true,
    },
  ];
  return (
    <div className="result-cards" aria-label="Frozen Phase 1 headline">
      {cards.map((c) => (
        <div
          key={c.k}
          className={`result-card${c.strong ? ' result-card--strong' : ''}${
            c.accentNegative ? ' result-card--delta' : ''
          }`}
        >
          <p className="result-card__value tnum">{c.v}</p>
          <p className="result-card__key">{c.k}</p>
          <p className="result-card__sub">{c.sub}</p>
        </div>
      ))}
    </div>
  );
}
