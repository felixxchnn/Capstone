import type {
  ExperimentalCandidate,
  PredictedModel,
  StructureResolution,
  StructureSelection,
} from '@/types/structure';

function sameSelection(a: StructureSelection | null, b: StructureSelection): boolean {
  if (!a || a.kind !== b.kind) return false;
  if (a.kind === 'experimental' && b.kind === 'experimental') {
    return a.candidate.pdbId === b.candidate.pdbId;
  }
  if (a.kind === 'predicted' && b.kind === 'predicted') {
    return a.model.alphafoldId === b.model.alphafoldId;
  }
  return false;
}

export function CandidateList({
  resolution,
  selection,
  onSelect,
}: {
  resolution: StructureResolution;
  selection: StructureSelection | null;
  onSelect: (s: StructureSelection) => void;
}) {
  const { experimental, predicted } = resolution;

  return (
    <div className="candidate-list">
      <fieldset>
        <legend>
          Experimental structures{' '}
          <span className="muted tiny">
            {experimental.length === 0 ? 'none in the PDB' : `${experimental.length} listed, best resolution first`}
          </span>
        </legend>
        {experimental.length === 0 ? (
          <p className="small muted">
            No experimental PDB entry references this protein. The predicted model below is
            the fallback.
          </p>
        ) : (
          <ul className="candidate-list__items">
            {experimental.map((c) => (
              <li key={c.pdbId}>
                <ExperimentalRow
                  candidate={c}
                  selected={sameSelection(selection, { kind: 'experimental', candidate: c })}
                  onSelect={() => onSelect({ kind: 'experimental', candidate: c })}
                />
              </li>
            ))}
          </ul>
        )}
      </fieldset>

      <fieldset>
        <legend>
          Predicted model <span className="muted tiny">AlphaFold DB — fallback</span>
        </legend>
        {predicted ? (
          <PredictedRow
            model={predicted}
            selected={sameSelection(selection, { kind: 'predicted', model: predicted })}
            onSelect={() => onSelect({ kind: 'predicted', model: predicted })}
          />
        ) : (
          <p className="small muted">No AlphaFold DB model is available for this protein.</p>
        )}
      </fieldset>
    </div>
  );
}

function ExperimentalRow({
  candidate,
  selected,
  onSelect,
}: {
  candidate: ExperimentalCandidate;
  selected: boolean;
  onSelect: () => void;
}) {
  return (
    <button
      type="button"
      className={`candidate-card${selected ? ' is-selected' : ''}`}
      aria-pressed={selected}
      onClick={onSelect}
    >
      <span className="candidate-card__id">
        <span className="badge badge--anchor">Experimental</span>
        <span className="mono">{candidate.pdbId}</span>
      </span>
      <span className="candidate-card__meta">
        {candidate.method ?? 'method n/a'}
        {candidate.resolutionAngstrom !== null ? ` · ${candidate.resolutionAngstrom.toFixed(2)} Å` : ''}
        {candidate.releaseDate ? ` · ${candidate.releaseDate}` : ''}
      </span>
      <span className="candidate-card__title">{candidate.title}</span>
      {candidate.citationPubMedId || candidate.citationDoi ? (
        <span className="candidate-card__cite tiny">
          {candidate.citationPubMedId ? (
            <a
              href={`https://pubmed.ncbi.nlm.nih.gov/${candidate.citationPubMedId}/`}
              target="_blank"
              rel="noopener noreferrer"
              onClick={(e) => e.stopPropagation()}
            >
              PMID {candidate.citationPubMedId}
            </a>
          ) : null}
          {candidate.citationDoi ? (
            <a
              href={`https://doi.org/${candidate.citationDoi}`}
              target="_blank"
              rel="noopener noreferrer"
              onClick={(e) => e.stopPropagation()}
            >
              {' '}
              doi:{candidate.citationDoi}
            </a>
          ) : null}
        </span>
      ) : null}
    </button>
  );
}

function PredictedRow({
  model,
  selected,
  onSelect,
}: {
  model: PredictedModel;
  selected: boolean;
  onSelect: () => void;
}) {
  return (
    <button
      type="button"
      className={`candidate-card candidate-card--predicted${selected ? ' is-selected' : ''}`}
      aria-pressed={selected}
      onClick={onSelect}
    >
      <span className="candidate-card__id">
        <span className="badge badge--external">Predicted structure</span>
        <span className="mono">{model.alphafoldId}</span>
      </span>
      <span className="candidate-card__meta">
        AlphaFold DB v{model.version}
        {model.meanPlddt !== null ? ` · mean pLDDT ${model.meanPlddt.toFixed(1)}` : ''}
        {model.modelCreatedDate ? ` · ${model.modelCreatedDate}` : ''}
      </span>
      <span className="candidate-card__title">
        Predicted model — <strong>not an experimental measurement</strong>. pLDDT is a
        per-residue confidence score, not accuracy.
      </span>
    </button>
  );
}
