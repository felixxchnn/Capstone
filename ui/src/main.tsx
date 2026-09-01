import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';
import { HashRouter } from 'react-router-dom';
import { App } from '@/app/App';
import { DataSourceProvider } from '@/app/Providers';
import { ErrorBoundary } from '@/components/common/ErrorBoundary';
import './styles/base.css';
import './styles/components.css';

const rootEl = document.getElementById('root');
if (!rootEl) throw new Error('Missing #root element');

// HashRouter: the app is a static bundle that must also work when opened from a
// plain file server or a sub-path with no server rewrites.
createRoot(rootEl).render(
  <StrictMode>
    <ErrorBoundary
      label="app-root"
      fallback={(error, reset) => (
        <div className="container page">
          <div className="callout callout--danger" role="alert">
            <p className="callout__title">The interface hit an unexpected error.</p>
            <p className="small">{error.message}</p>
            <button type="button" className="btn" onClick={reset}>
              Try again
            </button>
          </div>
        </div>
      )}
    >
      <HashRouter>
        <DataSourceProvider>
          <App />
        </DataSourceProvider>
      </HashRouter>
    </ErrorBoundary>
  </StrictMode>,
);
