import type { RankingView } from '@/hooks/useCaseStudyViews';
import type { SampleSummary } from '@/data/CapstoneDataSource';
import { integer } from '@/lib/format';

export function ResultsSummary({
  view,
  sample,
}: {
  view: RankingView;
  sample: SampleSummary;
}) {
  const cited = view.rows.filter((r) => r.evidenceStatus === 'cited').length;
  const sourceOnly = view.rows.filter((r) => r.evidenceStatus === 'source_only').length;
  const none = view.rows.filter((r) => r.evidenceStatus === 'none_in_filtered_snapshot').length;

  const cells: { label: string; value: string; hint?: string }[] = [
    { label: 'Sample', value: view.sample, hint: sample.role.replace(/_/g, ' ') },
    { label: 'Model', value: view.model },
    { label: 'Ranked targets', value: integer(view.block.n_targets_ranked) },
    { label: 'Displayed', value: `${view.block.n_displayed} (top 25)` },
    {
      label: 'Evidence coverage',
      value: `${cited} cited · ${sourceOnly} source-only · ${none} none`,
      hint: 'among the 25 displayed genes',
    },
    {
      label: 'Observed outcomes',
      value: sample.hasObservedOutcome ? 'available (attached after ranking)' : 'unavailable',
    },
  ];

  return (
    <dl className="results-summary" aria-label="Results summary">
      {cells.map((c) => (
        <div key={c.label} className="results-summary__cell">
          <dt>{c.label}</dt>
          <dd>
            <span className="results-summary__value">{c.value}</span>
            {c.hint ? <span className="results-summary__hint">{c.hint}</span> : null}
          </dd>
        </div>
      ))}
    </dl>
  );
}
