import { useEffect, useRef, useState } from 'react';
import 'molstar/build/viewer/molstar.css';
import { Viewer } from 'molstar/lib/apps/viewer/app';
import { Color } from 'molstar/lib/mol-util/color';
import type { StructureSelection } from '@/types/structure';
import { Spinner } from '@/components/common/primitives';

// This whole module is behind React.lazy() + an ErrorBoundary + a WebGL guard,
// so it only ever runs in a browser that can render it.

type ViewerLike = InstanceType<typeof Viewer>;

const BG_LIGHT = Color(0xf7fbf8);
const BG_DARK = Color(0x07110b);

export interface MolstarViewerProps {
  selection: StructureSelection;
  /** Preferred canvas background, matched to the app theme. */
  background: 'light' | 'dark';
  reducedMotion: boolean;
  onLoadState: (state: { status: 'loading' | 'ready' | 'error'; message?: string }) => void;
}

export function MolstarViewer({
  selection,
  background,
  reducedMotion,
  onLoadState,
}: MolstarViewerProps) {
  const hostRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<ViewerLike | null>(null);
  const [ready, setReady] = useState(false);
  const loadTokenRef = useRef(0);

  // create the viewer once
  useEffect(() => {
    let disposed = false;
    const host = hostRef.current;
    if (!host) return;
    onLoadState({ status: 'loading' });
    Viewer.create(host, {
      layoutIsExpanded: false,
      layoutShowControls: true,
      layoutShowLeftPanel: true,
      collapseLeftPanel: true,
      layoutShowSequence: true,
      layoutShowLog: false,
      layoutShowRemoteState: false,
      viewportShowExpand: true,
      viewportShowScreenshotControls: true,
      viewportShowReset: true,
      viewportShowControls: true,
      viewportShowAnimation: !reducedMotion,
      viewportShowSelectionMode: true,
      viewportBackgroundColor: background === 'dark' ? '#07110b' : '#f7fbf8',
      pdbProvider: 'rcsb',
    })
      .then((v) => {
        if (disposed) {
          v.dispose();
          return;
        }
        viewerRef.current = v;
        setReady(true);
      })
      .catch((err: unknown) => {
        onLoadState({
          status: 'error',
          message: err instanceof Error ? err.message : 'Mol* failed to initialise.',
        });
      });
    return () => {
      disposed = true;
      viewerRef.current?.dispose();
      viewerRef.current = null;
    };
    // create once; background / reducedMotion changes are applied below
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // (re)load the selected structure whenever it changes
  useEffect(() => {
    const v = viewerRef.current;
    if (!v || !ready) return;
    const token = ++loadTokenRef.current;
    onLoadState({ status: 'loading' });
    const run = async () => {
      try {
        await v.plugin.clear();
        if (selection.kind === 'experimental') {
          await v.loadPdb(selection.candidate.pdbId);
        } else {
          await v.loadAlphaFoldDb(selection.model.uniprotAccession);
        }
        if (loadTokenRef.current !== token) return; // superseded
        v.plugin.canvas3d?.requestCameraReset();
        onLoadState({ status: 'ready' });
      } catch (err) {
        if (loadTokenRef.current !== token) return;
        onLoadState({
          status: 'error',
          message:
            err instanceof Error
              ? `Could not load the structure: ${err.message}`
              : 'Could not load the structure.',
        });
      }
    };
    void run();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [ready, selection.kind, selectionKey(selection)]);

  // background toggle
  useEffect(() => {
    const v = viewerRef.current;
    if (!v || !ready || !v.plugin.canvas3d) return;
    v.plugin.canvas3d.setProps({
      renderer: { backgroundColor: background === 'dark' ? BG_DARK : BG_LIGHT },
    });
    v.plugin.canvas3d.requestDraw();
  }, [ready, background]);

  return (
    <div className="molstar-host-wrap">
      {!ready ? (
        <div className="molstar-host-loading">
          <Spinner label="Starting Mol* viewer" />
        </div>
      ) : null}
      <div ref={hostRef} className="molstar-host" data-testid="molstar-host" />
    </div>
  );
}

function selectionKey(s: StructureSelection): string {
  return s.kind === 'experimental'
    ? `x:${s.candidate.pdbId}`
    : `p:${s.model.uniprotAccession}:${s.model.version}`;
}
