import { useId } from 'react';
import type { EvidenceFilter } from '@/lib/selection';

const EVIDENCE_OPTIONS: { value: EvidenceFilter; label: string }[] = [
  { value: 'all', label: 'All evidence' },
  { value: 'cited', label: 'Cited evidence' },
  { value: 'source_only', label: 'Source-only evidence' },
  { value: 'none_in_filtered_snapshot', label: 'No evidence in filtered snapshot' },
];

export function SearchAndFilters({
  search,
  evidence,
  onSearch,
  onEvidence,
  onReset,
  visibleCount,
  totalCount,
  filtering,
}: {
  search: string;
  evidence: EvidenceFilter;
  onSearch: (v: string) => void;
  onEvidence: (v: EvidenceFilter) => void;
  onReset: () => void;
  visibleCount: number;
  totalCount: number;
  filtering: boolean;
}) {
  const searchId = useId();
  const evidenceId = useId();
  return (
    <div className="search-filters no-print" role="search">
      <div className="field search-filters__search">
        <label htmlFor={searchId}>Search gene, Entrez ID, drug, or evidence source</label>
        <input
          id={searchId}
          type="search"
          value={search}
          placeholder="e.g. YRDC, 79693, dinaciclib, ChEMBL"
          autoComplete="off"
          spellCheck={false}
          onChange={(e) => onSearch(e.target.value)}
        />
      </div>
      <div className="field search-filters__evidence">
        <label htmlFor={evidenceId}>Evidence filter</label>
        <select
          id={evidenceId}
          value={evidence}
          onChange={(e) => onEvidence(e.target.value as EvidenceFilter)}
        >
          {EVIDENCE_OPTIONS.map((o) => (
            <option key={o.value} value={o.value}>
              {o.label}
            </option>
          ))}
        </select>
      </div>
      <div className="search-filters__meta">
        <p className="small" aria-live="polite">
          <strong>
            {visibleCount} of {totalCount}
          </strong>{' '}
          genes shown
          {filtering ? (
            <span className="filter-flag" aria-label="Filters are active">
              {' '}
              · filters active
            </span>
          ) : null}
        </p>
        <button
          type="button"
          className="btn btn--subtle btn--xs"
          onClick={onReset}
          disabled={!filtering}
        >
          Reset filters
        </button>
      </div>
    </div>
  );
}
