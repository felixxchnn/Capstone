import { useSelection } from '@/hooks/useSelection';
import { useDataSource } from '@/app/dataSourceContext';
import { useAsync } from '@/hooks/useAsync';
import { SegmentedControl } from '@/components/common/SegmentedControl';
import { SampleRoleLine } from '@/components/samples/SampleRoleBadge';
import { ModelComparison } from '@/components/comparison/ModelComparison';
import { Callout, EmptyState, Spinner } from '@/components/common/primitives';
import type { SampleId } from '@/types/caseStudy';

export function ModelComparisonPage() {
  const { selection, setSelection } = useSelection();
  const ds = useDataSource();

  const state = useAsync(
    async (_s) => {
      const [pca, head, osteo] = await Promise.all([
        ds.getModelRanking(selection.sample, 'ridge_pca'),
        ds.getModelRanking(selection.sample, 'ridge_head'),
        ds.getProjectSummary(),
      ]);
      return { pca, head, osteo: osteo.osteosarcoma };
    },
    [ds, selection.sample],
  );

  return (
    <div className="container page">
      <header className="stack">
        <h1>Model Comparison</h1>
        <p className="lede">
          <code>ridge_pca</code> and <code>ridge_head</code> side by side for one sample.
          Each keeps its own independent, frozen top-25. There is no consensus list.
        </p>
      </header>

      <SegmentedControl<SampleId>
        legend="Sample"
        name="compare-sample"
        value={selection.sample}
        onChange={(s) => setSelection({ sample: s })}
        options={[
          { value: 'ACH-000364', label: 'ACH-000364', hint: 'validation anchor' },
          { value: 'BG003082', label: 'BG003082', hint: 'exploratory external' },
        ]}
      />
      <SampleRoleLine sample={selection.sample} />

      {state.status === 'loading' ? <Spinner label="Loading both rankings" /> : null}
      {state.status === 'error' ? (
        <EmptyState title="Could not load the comparison">{state.error.message}</EmptyState>
      ) : null}

      {state.status === 'success' ? (
        <>
          <ModelComparison sample={selection.sample} pca={state.data.pca} head={state.data.head} />

          <Callout tone="info" title="The frozen Phase 1 verdict still stands">
            Across the full 170-line validation split, ridge_head did not outperform
            ridge_pca (delta −0.0308). The descriptive osteosarcoma-cohort aggregate
            (n = 5, unstable) points the same way: mean per-target Spearman{' '}
            {state.data.osteo.mean_per_target_spearman.ridge_pca.toFixed(6)} for ridge_pca
            vs {state.data.osteo.mean_per_target_spearman.ridge_head.toFixed(6)} for
            ridge_head. Neither figure was used to choose a model or alter any ranking.
          </Callout>
        </>
      ) : null}
    </div>
  );
}
