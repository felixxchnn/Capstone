import { lazy, Suspense } from 'react';
import type { MolstarViewerProps } from '@/components/structures/MolstarViewer';
import { Spinner } from '@/components/common/primitives';

// The heavy Mol* bundle loads only when this component first renders.
const MolstarViewer = lazy(() =>
  import('@/components/structures/MolstarViewer').then((m) => ({ default: m.MolstarViewer })),
);

export function MolstarViewerLazy(props: MolstarViewerProps) {
  return (
    <Suspense
      fallback={
        <div className="molstar-host-loading">
          <Spinner label="Loading 3D structure viewer" />
        </div>
      }
    >
      <MolstarViewer {...props} />
    </Suspense>
  );
}
