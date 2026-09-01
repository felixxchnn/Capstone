import type { StructureResolution, StructureStep } from '@/types/structure';

const ICON: Record<StructureStep['status'], string> = {
  ok: '●',
  empty: '○',
  error: '▲',
  skipped: '–',
};

export function StructureProviderStatus({ resolution }: { resolution: StructureResolution }) {
  return (
    <div className="provider-status" aria-label="Structure retrieval status">
      <ol className="provider-status__steps">
        {resolution.steps.map((s) => (
          <li key={s.id} className={`provider-status__step is-${s.status}`}>
            <span className="provider-status__icon" aria-hidden="true">
              {ICON[s.status]}
            </span>
            <span className="provider-status__body">
              <span className="provider-status__label">
                {s.label}
                <span className="sr-only"> — {s.status}</span>
              </span>
              <span className="provider-status__detail">{s.detail}</span>
            </span>
          </li>
        ))}
      </ol>
      <p className="tiny muted">
        Chain: Entrez Gene ID → reviewed human UniProt → experimental RCSB PDB candidates →
        AlphaFold DB predicted model (fallback). Only identifiers and human taxonomy 9606
        are sent to these services.
      </p>
    </div>
  );
}
