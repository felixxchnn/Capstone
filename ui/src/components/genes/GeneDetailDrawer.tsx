import { Link } from 'react-router-dom';
import { Drawer } from '@/components/common/Drawer';
import { PredictedObservedStatus } from '@/components/genes/PredictedObservedStatus';
import { EvidenceStatusPill } from '@/components/evidence/EvidenceStatusPill';
import { EvidenceRecordItem } from '@/components/evidence/EvidenceRecordItem';
import { Callout } from '@/components/common/primitives';
import { EVIDENCE_NOT_EFFICACY, NOT_A_TARGET, evidenceStatusLabel } from '@/lib/format';
import type { GeneRowView } from '@/hooks/useCaseStudyViews';
import type { ModelId, SampleId } from '@/types/caseStudy';

export function GeneDetailDrawer({
  open,
  onClose,
  row,
  sample,
  model,
  structureHref,
}: {
  open: boolean;
  onClose: () => void;
  row: GeneRowView | null;
  sample: SampleId;
  model: ModelId;
  structureHref: string;
}) {
  if (!row) return null;
  const g = row.gene;
  const records = row.evidence?.records ?? [];

  return (
    <Drawer
      open={open}
      onClose={onClose}
      labelledById="gene-detail-title"
      title={
        <>
          {g.symbol}{' '}
          <span className="mono muted" style={{ fontWeight: 400 }}>
            Entrez {g.entrez_id}
          </span>
        </>
      }
    >
      <div className="stack">
        <p className="small">
          Rank <strong>{g.rank}</strong> for sample <span className="mono">{sample}</span>,
          model <span className="mono">{model}</span>.
        </p>
        <Callout tone="info">{NOT_A_TARGET}</Callout>

        <section aria-label="Predicted versus observed">
          <h3>Predicted vs observed</h3>
          <PredictedObservedStatus sample={sample} gene={g} />
        </section>

        <section aria-label="Drug–gene interaction evidence">
          <h3>
            Drug–gene interaction evidence <EvidenceStatusPill status={row.evidenceStatus} />
          </h3>
          <p className="small muted">
            {evidenceStatusLabel(row.evidenceStatus)}. Retrieved by Entrez ID from the
            licence-filtered offline DGIdb snapshot <em>after</em> the ranking was frozen.
            Evidence volume did not affect this gene's rank.
          </p>
          {records.length > 0 ? (
            <ul className="evi-record-list">
              {records.map((r) => (
                <EvidenceRecordItem key={r.record_key} record={r} />
              ))}
            </ul>
          ) : (
            <p className="small muted">
              No record for {g.symbol} in the filtered snapshot.
            </p>
          )}
          <Callout tone="caution" title="Evidence is retrieval, not efficacy">
            {EVIDENCE_NOT_EFFICACY}
          </Callout>
        </section>

        <section aria-label="Encoded protein structure">
          <h3>Encoded protein structure</h3>
          <p className="small muted">
            Opens the experimental or predicted structure of the protein encoded by{' '}
            {g.symbol} (Entrez {g.entrez_id}, human taxonomy 9606). Structural evidence
            only — not drug-response evidence.
          </p>
          <Link className="btn" to={structureHref} onClick={onClose}>
            View encoded protein structure →
          </Link>
        </section>
      </div>
    </Drawer>
  );
}
