// CSV + JSON export of the *currently displayed* ranked view.
//
// Invariants:
//  - Rows are emitted in frozen `rank` order (1..25), ALWAYS. A visible-row
//    filter is applied only if the caller passes `visibleEntrezIds`; hidden rows
//    are dropped but the remaining rows keep their original `rank` (never
//    renumbered).
//  - Every row records: sample, model, original rank, symbol, entrez, predicted
//    value, observed availability + value, evidence classification.

import type { EvidenceStatus, ModelId, ModelRankingBlock, SampleId } from '@/types/caseStudy';
import { geneEffectFull } from '@/lib/format';

export interface ExportRow {
  sample: SampleId;
  model: ModelId;
  rank: number;
  symbol: string;
  entrez_id: string;
  predicted_geneeffect: number;
  observed_available: boolean;
  observed_geneeffect: number | null;
  observed_rank: number | null;
  evidence_status: EvidenceStatus;
  evidence_record_count: number;
}

export interface BuildExportRowsArgs {
  sample: SampleId;
  model: ModelId;
  block: ModelRankingBlock;
  evidenceStatusFor: (entrezId: string) => { status: EvidenceStatus; count: number };
  /** If given, only these Entrez IDs are exported (order still by rank). */
  visibleEntrezIds?: ReadonlySet<string>;
}

export function buildExportRows(args: BuildExportRowsArgs): ExportRow[] {
  const { sample, model, block, evidenceStatusFor, visibleEntrezIds } = args;
  return block.genes
    .filter((g) => !visibleEntrezIds || visibleEntrezIds.has(g.entrez_id))
    .slice()
    .sort((a, b) => a.rank - b.rank)
    .map((g) => {
      const ev = evidenceStatusFor(g.entrez_id);
      const hasObs = g.observed_geneeffect !== undefined && g.observed_geneeffect !== null;
      return {
        sample,
        model,
        rank: g.rank,
        symbol: g.symbol,
        entrez_id: g.entrez_id,
        predicted_geneeffect: g.predicted_geneeffect,
        observed_available: hasObs,
        observed_geneeffect: hasObs ? (g.observed_geneeffect as number) : null,
        observed_rank: g.observed_rank ?? null,
        evidence_status: ev.status,
        evidence_record_count: ev.count,
      };
    });
}

function csvCell(value: string | number | boolean | null): string {
  if (value === null) return '';
  const s = String(value);
  return /[",\r\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

export function rowsToCsv(rows: ExportRow[]): string {
  const header = [
    'sample',
    'sample_role_note',
    'model',
    'rank',
    'symbol',
    'entrez_id',
    'predicted_geneeffect',
    'observed_available',
    'observed_geneeffect',
    'observed_rank',
    'evidence_status',
    'evidence_record_count',
  ];
  const roleNote = (s: SampleId) =>
    s === 'ACH-000364'
      ? 'held-out validation anchor; observed CRISPR attached AFTER ranking'
      : 'exploratory external bulk-tumour prediction; no measured CRISPR outcome';
  const lines = [header.join(',')];
  for (const r of rows) {
    lines.push(
      [
        csvCell(r.sample),
        csvCell(roleNote(r.sample)),
        csvCell(r.model),
        csvCell(r.rank),
        csvCell(r.symbol),
        csvCell(r.entrez_id),
        csvCell(geneEffectFull(r.predicted_geneeffect)),
        csvCell(r.observed_available),
        csvCell(r.observed_geneeffect === null ? '' : geneEffectFull(r.observed_geneeffect)),
        csvCell(r.observed_rank),
        csvCell(r.evidence_status),
        csvCell(r.evidence_record_count),
      ].join(','),
    );
  }
  return lines.join('\r\n') + '\r\n';
}

export interface ExportJsonView {
  exported_by: 'capstone-research-interface';
  disclaimer: string;
  case_study_sha256: string;
  view: {
    sample: SampleId;
    sample_role: string;
    model: ModelId;
    ranking_rule: string;
    n_targets_ranked: number;
    n_displayed: number;
    n_exported: number;
    filtered: boolean;
  };
  rows: ExportRow[];
}

export function rowsToJsonView(args: {
  sample: SampleId;
  model: ModelId;
  block: ModelRankingBlock;
  rows: ExportRow[];
  caseStudySha256: string;
  filtered: boolean;
}): ExportJsonView {
  const { sample, model, block, rows, caseStudySha256, filtered } = args;
  return {
    exported_by: 'capstone-research-interface',
    disclaimer:
      'Predicted CRISPR gene dependencies from a reconstructed frozen Phase 1 model. Not therapeutic targets, drug recommendations, or clinical guidance. Ranks are the frozen rank order and are not renumbered by filtering.',
    case_study_sha256: caseStudySha256,
    view: {
      sample,
      sample_role:
        sample === 'ACH-000364'
          ? 'held_out_prediction / measured_crispr (verification anchor)'
          : 'exploratory_external_prediction / unavailable',
      model,
      ranking_rule: block.ranking_rule,
      n_targets_ranked: block.n_targets_ranked,
      n_displayed: block.n_displayed,
      n_exported: rows.length,
      filtered,
    },
    rows,
  };
}

/** Browser download. Uses a Blob + object URL; no network. */
export function triggerDownload(filename: string, mime: string, content: string): void {
  const blob = new Blob([content], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  a.rel = 'noopener';
  document.body.appendChild(a);
  a.click();
  a.remove();
  // give the browser a tick to start the download before revoking
  setTimeout(() => URL.revokeObjectURL(url), 4000);
}

export function exportFilename(sample: SampleId, model: ModelId, ext: 'csv' | 'json', filtered: boolean): string {
  return `capstone_${sample}_${model}${filtered ? '_filtered' : ''}_top25.${ext}`;
}
