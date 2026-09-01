import { useMemo, useState } from 'react';
import { useSelection } from '@/hooks/useSelection';
import { useDataSource } from '@/app/dataSourceContext';
import { useAsync } from '@/hooks/useAsync';
import { useRankingView } from '@/hooks/useCaseStudyViews';
import { useIsMobile } from '@/hooks/useMediaQuery';
import { isFiltering } from '@/lib/selection';
import { DEPENDENCY_DIRECTION, NOT_A_TARGET } from '@/lib/format';
import { Callout, EmptyState, Spinner } from '@/components/common/primitives';
import { SampleModelControls } from '@/components/samples/SampleModelControls';
import { SampleRoleLine } from '@/components/samples/SampleRoleBadge';
import { SearchAndFilters } from '@/components/rankings/SearchAndFilters';
import { ResultsSummary } from '@/components/rankings/ResultsSummary';
import { RankingTable } from '@/components/rankings/RankingTable';
import { RankingCardList } from '@/components/rankings/RankingCardList';
import { EvidenceGroup } from '@/components/evidence/EvidenceGroup';
import { GeneDetailDrawer } from '@/components/genes/GeneDetailDrawer';
import { ExportButtons } from '@/components/export/ExportButtons';

export function DependencyExplorerPage() {
  const { selection, setSelection, resetFilters } = useSelection();
  const ds = useDataSource();
  const isMobile = useIsMobile();

  const summary = useAsync((_s) => ds.getProjectSummary(), [ds]);
  const samples = useAsync((_s) => ds.getSamples(), [ds]);
  const { view, status, error, reload } = useRankingView(
    selection.sample,
    selection.model,
    selection,
  );

  const [drawerGene, setDrawerGene] = useState<string | null>(null);
  const openGene = drawerGene ?? selection.gene;

  const activeRow = useMemo(
    () => view?.rows.find((r) => r.gene.entrez_id === openGene) ?? null,
    [view, openGene],
  );

  const filtering = isFiltering(selection);
  const sampleSummary = samples.data?.find((s) => s.id === selection.sample) ?? null;

  return (
    <div className="container page">
      <header className="stack">
        <h1>Dependency Explorer</h1>
        <p className="lede">
          The frozen top-25 predicted CRISPR gene dependencies for one sample and one
          model, with drug–gene interaction evidence retrieved after ranking. Search and
          filters hide rows; they never re-order or renumber them.
        </p>
        <Callout tone="info" title="Reading the ranking">
          {DEPENDENCY_DIRECTION} {NOT_A_TARGET}
        </Callout>
      </header>

      <SampleModelControls
        sample={selection.sample}
        model={selection.model}
        onSample={(s) => setSelection({ sample: s, gene: null })}
        onModel={(m) => setSelection({ model: m, gene: null })}
      />
      <SampleRoleLine sample={selection.sample} />

      <SearchAndFilters
        search={selection.search}
        evidence={selection.evidence}
        onSearch={(v) => setSelection({ search: v })}
        onEvidence={(v) => setSelection({ evidence: v })}
        onReset={resetFilters}
        visibleCount={view?.visibleCount ?? 0}
        totalCount={view?.rows.length ?? 25}
        filtering={filtering}
      />

      {status === 'loading' ? <Spinner label="Loading ranking" /> : null}

      {status === 'error' ? (
        <EmptyState
          title="Could not load this ranking"
          action={
            <button type="button" className="btn" onClick={reload}>
              Retry
            </button>
          }
        >
          {error?.message}
        </EmptyState>
      ) : null}

      {view && sampleSummary ? (
        <>
          <ResultsSummary view={view} sample={sampleSummary} />

          <ExportButtons
            view={view}
            caseStudySha256={summary.data?.caseStudySha256 ?? ''}
            onlyVisible={filtering}
          />

          <section aria-label="Ranked dependency table">
            {isMobile ? (
              <RankingCardList
                rows={view.rows}
                hasObserved={view.hasObserved}
                onSelectGene={(id) => {
                  setDrawerGene(id);
                  setSelection({ gene: id });
                }}
                onResetFilters={resetFilters}
                filtering={filtering}
              />
            ) : (
              <RankingTable
                sample={view.sample}
                model={view.model}
                rows={view.rows}
                hasObserved={view.hasObserved}
                nTargetsRanked={view.block.n_targets_ranked}
                onSelectGene={(id) => {
                  setDrawerGene(id);
                  setSelection({ gene: id });
                }}
                onResetFilters={resetFilters}
                filtering={filtering}
              />
            )}
            <p className="tiny muted">{view.block.ranking_rule}</p>
          </section>

          <section aria-label="Drug–gene interaction evidence" className="stack">
            <h2>Drug–gene interaction evidence</h2>
            <p className="small muted">
              Grouped beneath each displayed gene, in rank order. Retrieved by Entrez ID{' '}
              <em>after</em> the ranking froze — a gene with more evidence did not get a
              better rank.
            </p>
            {view.rows.filter((r) => r.visible).length === 0 ? (
              <EmptyState title="No evidence groups match the current filters" />
            ) : (
              <div className="evi-group-list">
                {view.rows
                  .filter((r) => r.visible)
                  .map((r) => (
                    <EvidenceGroup key={r.gene.entrez_id} row={r} />
                  ))}
              </div>
            )}
          </section>
        </>
      ) : null}

      <GeneDetailDrawer
        open={Boolean(activeRow)}
        onClose={() => {
          setDrawerGene(null);
          setSelection({ gene: null });
        }}
        row={activeRow}
        sample={selection.sample}
        model={selection.model}
        structureHref={`/structure?${new URLSearchParams({
          sample: selection.sample,
          model: selection.model,
          gene: activeRow?.gene.entrez_id ?? '',
        }).toString()}`}
      />
    </div>
  );
}
