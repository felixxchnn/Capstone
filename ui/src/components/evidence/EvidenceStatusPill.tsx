import type { EvidenceStatus } from '@/types/caseStudy';
import { evidenceStatusLabel } from '@/lib/format';

export function EvidenceStatusPill({ status }: { status: EvidenceStatus }) {
  return (
    <span className={`status-pill status-pill--${status}`}>
      {evidenceStatusLabel(status)}
    </span>
  );
}
