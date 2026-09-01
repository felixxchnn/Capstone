import { Suspense, lazy } from 'react';
import { Navigate, Route, Routes } from 'react-router-dom';
import { TopNav } from '@/components/navigation/TopNav';
import { NonClinicalBanner } from '@/components/safety/NonClinicalBanner';
import { SkipLink, Spinner } from '@/components/common/primitives';
import { HelixMark } from '@/components/navigation/HelixMark';
import { OverviewPage } from '@/pages/OverviewPage';
import { DependencyExplorerPage } from '@/pages/DependencyExplorerPage';
import { ModelComparisonPage } from '@/pages/ModelComparisonPage';
import { MethodsPage } from '@/pages/MethodsPage';

// Protein structure page pulls in the Mol* provider chain and (lazily) Mol*
// itself — keep it out of the initial route bundle.
const ProteinStructurePage = lazy(() =>
  import('@/pages/ProteinStructurePage').then((m) => ({ default: m.ProteinStructurePage })),
);

export function App() {
  return (
    <>
      <SkipLink targetId="main" />
      <TopNav />
      <main id="main" className="app-main" tabIndex={-1}>
        <div className="container">
          <NonClinicalBanner />
        </div>
        <Suspense
          fallback={
            <div className="container page">
              <Spinner label="Loading view" />
            </div>
          }
        >
          <Routes>
            <Route path="/" element={<OverviewPage />} />
            <Route path="/explore" element={<DependencyExplorerPage />} />
            <Route path="/compare" element={<ModelComparisonPage />} />
            <Route path="/structure" element={<ProteinStructurePage />} />
            <Route path="/methods" element={<MethodsPage />} />
            <Route path="*" element={<Navigate to="/" replace />} />
          </Routes>
        </Suspense>
      </main>
      <footer className="app-footer no-print">
        <div className="container app-footer__inner">
          <HelixMark size={24} />
          <p className="small">
            Connected presentation layer over committed model predictions. No model
            inference runs in this browser. The offline{' '}
            <code>phase2_report.html</code> is the fully static counterpart.
          </p>
        </div>
      </footer>
    </>
  );
}
