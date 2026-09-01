import type { GeneRowView } from '@/hooks/useCaseStudyViews';
import { Disclosure } from '@/components/common/Disclosure';
import { EvidenceStatusPill } from '@/components/evidence/EvidenceStatusPill';
import { EvidenceRecordItem } from '@/components/evidence/EvidenceRecordItem';

/** One gene's drug–gene interaction evidence, grouped under the gene and
 *  retrieved AFTER ranking. Evidence volume never affects rank. */
export function EvidenceGroup({ row, defaultOpen = false }: { row: GeneRowView; defaultOpen?: boolean }) {
  const g = row.gene;
  const records = row.evidence?.records ?? [];
  return (
    <Disclosure
      defaultOpen={defaultOpen && records.length > 0}
      summary={
        <span className="evi-group__summary">
          <span className="evi-group__gene">
            {g.symbol} <span className="mono muted tiny">({g.entrez_id})</span>
          </span>
          <EvidenceStatusPill status={row.evidenceStatus} />
          <span className="muted tiny">
            {records.length} record{records.length === 1 ? '' : 's'}
          </span>
        </span>
      }
    >
      {records.length > 0 ? (
        <ul className="evi-record-list">
          {records.map((r) => (
            <EvidenceRecordItem key={r.record_key} record={r} />
          ))}
        </ul>
      ) : (
        <p className="small muted">
          No drug–gene interaction record for {g.symbol} in the licence-filtered offline
          DGIdb snapshot. This is not an error — the snapshot only carries records from
          six redistribution-verified sources.
        </p>
      )}
    </Disclosure>
  );
}
