import { useState } from 'react';
import type { RankingView } from '@/hooks/useCaseStudyViews';
import {
  buildExportRows,
  exportFilename,
  rowsToCsv,
  rowsToJsonView,
  triggerDownload,
} from '@/lib/export';

/** CSV + JSON export of the currently displayed ranked rows. If a filter hides
 *  rows, only the visible rows are exported and `filtered` is flagged — but the
 *  frozen `rank` of every exported row is preserved (never renumbered). */
export function ExportButtons({
  view,
  caseStudySha256,
  onlyVisible,
}: {
  view: RankingView;
  caseStudySha256: string;
  onlyVisible: boolean;
}) {
  const [error, setError] = useState<string | null>(null);

  const visibleSet = onlyVisible
    ? new Set(view.rows.filter((r) => r.visible).map((r) => r.gene.entrez_id))
    : undefined;

  const evidenceStatusFor = (entrezId: string) => {
    const row = view.rows.find((r) => r.gene.entrez_id === entrezId);
    return { status: row?.evidenceStatus ?? 'none_in_filtered_snapshot', count: row?.evidenceCount ?? 0 };
  };

  const rows = buildExportRows({
    sample: view.sample,
    model: view.model,
    block: view.block,
    evidenceStatusFor,
    visibleEntrezIds: visibleSet,
  });
  const filtered = Boolean(visibleSet && visibleSet.size < view.block.genes.length);

  const doExport = (kind: 'csv' | 'json') => {
    try {
      setError(null);
      if (rows.length === 0) {
        setError('Nothing to export — every row is hidden by the current filter.');
        return;
      }
      const name = exportFilename(view.sample, view.model, kind, filtered);
      if (kind === 'csv') {
        triggerDownload(name, 'text/csv;charset=utf-8', rowsToCsv(rows));
      } else {
        const payload = rowsToJsonView({
          sample: view.sample,
          model: view.model,
          block: view.block,
          rows,
          caseStudySha256,
          filtered,
        });
        triggerDownload(name, 'application/json', JSON.stringify(payload, null, 2) + '\n');
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Export failed.');
    }
  };

  return (
    <div className="export-buttons no-print">
      <span className="field-label">Export {onlyVisible && filtered ? 'visible' : 'displayed'} rows</span>
      <div className="export-buttons__row">
        <button type="button" className="btn btn--subtle btn--xs" onClick={() => doExport('csv')}>
          CSV
        </button>
        <button type="button" className="btn btn--subtle btn--xs" onClick={() => doExport('json')}>
          JSON
        </button>
        <button type="button" className="btn btn--subtle btn--xs" onClick={() => window.print()}>
          Print / Save as PDF
        </button>
      </div>
      {filtered ? (
        <p className="tiny muted">
          A filter is active: {rows.length} visible row{rows.length === 1 ? '' : 's'} exported,
          each keeping its frozen rank.
        </p>
      ) : null}
      {error ? (
        <p className="tiny" role="alert" style={{ color: 'var(--danger-ink)' }}>
          {error}
        </p>
      ) : null}
    </div>
  );
}
