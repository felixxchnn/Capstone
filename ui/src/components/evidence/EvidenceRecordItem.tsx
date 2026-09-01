import type { EvidenceRecord } from '@/types/caseStudy';
import { flag, pmidUrl } from '@/lib/format';

export function EvidenceRecordItem({ record }: { record: EvidenceRecord }) {
  const drug = record.drug_name || record.drug_claim_name || 'unnamed compound';
  return (
    <li className="evi-record">
      <p className="evi-record__title">
        <strong>{drug}</strong>{' '}
        <span className="muted">
          — {record.interaction_type_raw || 'interaction'} ·{' '}
          direction: {record.interaction_direction || 'unknown'} (tier{' '}
          {record.direction_tier})
        </span>
      </p>
      <dl className="kv kv--tight">
        <dt>Source</dt>
        <dd>
          {record.interaction_source}
          {record.interaction_source_version ? ` v${record.interaction_source_version}` : ''}
        </dd>
        <dt>Source licence</dt>
        <dd>
          {record.source_license_url ? (
            <a href={record.source_license_url} target="_blank" rel="noopener noreferrer">
              {record.source_license}
            </a>
          ) : (
            record.source_license || 'not stated'
          )}
        </dd>
        <dt>DGIdb regulatory-approval flag</dt>
        <dd>{flag(record.drug_is_approved)}</dd>
        <dt>DGIdb antineoplastic flag</dt>
        <dd>{flag(record.drug_is_antineoplastic)}</dd>
        <dt>DGIdb evidence score</dt>
        <dd>{record.evidence_score || 'none'}</dd>
        <dt>Publications</dt>
        <dd>
          {record.pmids.length > 0 ? (
            <ul className="pmid-list">
              {record.pmids.map((p) => (
                <li key={p}>
                  <a href={pmidUrl(p)} target="_blank" rel="noopener noreferrer">
                    PMID {p}
                  </a>
                </li>
              ))}
            </ul>
          ) : (
            <em className="muted">
              Source-only: DGIdb records no claim-level publication for this
              drug–gene–source group in the filtered snapshot.
            </em>
          )}
        </dd>
      </dl>
      {record.pmids.length > 0 ? (
        <p className="tiny muted">
          Group-level citation: PMIDs are attached at the drug–gene / interaction-source
          group level and may span multiple interaction claims.
        </p>
      ) : null}
      <p className="disclaimer tiny">{record.disclaimer}</p>
    </li>
  );
}
