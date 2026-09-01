import { useMemo, useState } from 'react';
import { useAsync } from '@/hooks/useAsync';
import { useReducedMotion } from '@/hooks/useMediaQuery';
import { resolveStructure, defaultSelection } from '@/data/providers/structureProvider';
import { detectWebgl } from '@/components/structures/webgl';
import { CandidateList } from '@/components/structures/CandidateList';
import { StructureProviderStatus } from '@/components/structures/StructureProviderStatus';
import { StructureTextSummary } from '@/components/structures/StructureTextSummary';
import { MolstarViewerLazy } from '@/components/structures/MolstarViewerLazy';
import { ErrorBoundary } from '@/components/common/ErrorBoundary';
import { Callout, EmptyState, Spinner } from '@/components/common/primitives';
import type { StructureSelection } from '@/types/structure';

export function StructureExplorer({
  entrezId,
  symbol,
}: {
  entrezId: string;
  symbol: string;
}) {
  const reducedMotion = useReducedMotion();
  const webgl = useMemo(() => detectWebgl(), []);
  const [background, setBackground] = useState<'light' | 'dark'>('dark');
  // The user's explicit candidate pick, tagged with the gene it was made for.
  // A pick for a stale gene is ignored — the default selection wins.
  const [userPick, setUserPick] = useState<{ entrezId: string; selection: StructureSelection } | null>(
    null,
  );
  const [viewerState, setViewerState] = useState<{
    status: 'idle' | 'loading' | 'ready' | 'error';
    message?: string;
  }>({ status: 'idle' });

  // resolve providers; aborts + resets when the gene changes
  const resolved = useAsync(
    (signal) => resolveStructure(entrezId, symbol, signal),
    [entrezId, symbol],
  );

  // Selection is derived, never stored in an effect: user pick for THIS gene,
  // otherwise the provider default (best experimental, else AlphaFold).
  const selection: StructureSelection | null = useMemo(() => {
    if (resolved.status !== 'success') return null;
    if (userPick && userPick.entrezId === entrezId) return userPick.selection;
    return defaultSelection(resolved.data);
  }, [resolved.status, resolved.data, userPick, entrezId]);

  const setSelection = (s: StructureSelection) => {
    setUserPick({ entrezId, selection: s });
    setViewerState({ status: 'loading' });
  };

  if (resolved.status === 'loading') {
    return <Spinner label={`Resolving structure providers for ${symbol}`} />;
  }

  if (resolved.status === 'error') {
    return (
      <EmptyState
        title="Structure services are unavailable"
        action={
          <button type="button" className="btn" onClick={resolved.reload}>
            Retry
          </button>
        }
      >
        {resolved.error.message} — UniProt, RCSB and AlphaFold are contacted directly from
        your browser; a network block or outage stops this panel but nothing else.
      </EmptyState>
    );
  }

  const res = resolved.data;
  const nothingAvailable = res.experimental.length === 0 && !res.predicted;

  return (
    <div className="structure-explorer">
      <div className="structure-explorer__grid">
        <aside className="structure-explorer__side">
          <StructureProviderStatus resolution={res} />
          {res.uniprot === null ? (
            <Callout tone="caution" title="Mapping failed">
              No reviewed human UniProt entry cross-references Entrez {entrezId}. Without a
              UniProt accession there is nothing to look up in the PDB or AlphaFold DB.
            </Callout>
          ) : (
            <CandidateList resolution={res} selection={selection} onSelect={setSelection} />
          )}
          <div className="structure-explorer__viewopts">
            <span className="field-label">Viewer background</span>
            <div className="seg-mini" role="group" aria-label="Viewer background">
              <button
                type="button"
                aria-pressed={background === 'dark'}
                className={background === 'dark' ? 'is-selected' : undefined}
                onClick={() => setBackground('dark')}
              >
                Near-black
              </button>
              <button
                type="button"
                aria-pressed={background === 'light'}
                className={background === 'light' ? 'is-selected' : undefined}
                onClick={() => setBackground('light')}
              >
                Off-white
              </button>
            </div>
            <p className="tiny muted">
              Rotate / zoom / pan, reset camera, representation (cartoon · molecular surface
              · ball-and-stick), colour by chain, show/hide ligands and water, fullscreen
              and screenshot are in the Mol* control panels around the canvas. Use the Mol*{' '}
              <em>Components</em> panel to isolate individual chains or ligands.
            </p>
          </div>
        </aside>

        <div className="structure-explorer__main">
          {nothingAvailable ? (
            <EmptyState title={`No structure available for ${symbol}`}>
              Neither the PDB nor AlphaFold DB has a model for this protein. This does not
              affect any prediction or ranking.
            </EmptyState>
          ) : !webgl.ok ? (
            <Callout tone="caution" title="WebGL unavailable — showing the text summary only">
              {webgl.reason} The interactive 3D viewer needs WebGL. The structure facts
              below are complete without it.
            </Callout>
          ) : selection === null ? (
            <EmptyState title="Select a structure">
              Choose an experimental candidate or the predicted model on the left.
            </EmptyState>
          ) : (
            <>
              {selection.kind === 'predicted' ? (
                <Callout tone="caution" title="Predicted structure">
                  This is an <strong>AlphaFold-predicted</strong> model, not an experimental
                  measurement. Colour reflects per-residue pLDDT confidence
                  {res.predicted?.meanPlddt !== null && res.predicted?.meanPlddt !== undefined
                    ? ` (mean ${res.predicted.meanPlddt.toFixed(1)})`
                    : ''}
                  .
                </Callout>
              ) : (
                <Callout tone="info" title={`Experimental structure ${selection.candidate.pdbId}`}>
                  {selection.candidate.method ?? 'Experimental'}{' '}
                  {selection.candidate.resolutionAngstrom !== null
                    ? `· ${selection.candidate.resolutionAngstrom.toFixed(2)} Å resolution`
                    : ''}
                  . Preferred over the predicted model when available.
                </Callout>
              )}

              {viewerState.status === 'error' ? (
                <Callout tone="danger" title="The viewer could not load this structure">
                  {viewerState.message} Try another candidate, or read the text summary
                  below.
                </Callout>
              ) : null}

              <ErrorBoundary
                label="molstar"
                fallback={(err, reset) => (
                  <Callout tone="danger" title="The 3D viewer crashed">
                    {err.message}
                    <div style={{ marginTop: 'var(--sp-2)' }}>
                      <button type="button" className="btn btn--subtle btn--xs" onClick={reset}>
                        Reload viewer
                      </button>
                    </div>
                  </Callout>
                )}
              >
                <MolstarViewerLazy
                  selection={selection}
                  background={background}
                  reducedMotion={reducedMotion}
                  onLoadState={setViewerState}
                />
              </ErrorBoundary>
            </>
          )}

          <StructureTextSummary resolution={res} selection={selection} />
        </div>
      </div>
    </div>
  );
}
