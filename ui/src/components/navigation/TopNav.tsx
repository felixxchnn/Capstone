import { useState } from 'react';
import { NavLink, useLocation } from 'react-router-dom';
import { NAV_ROUTES } from '@/app/routes';
import { HelixMark } from './HelixMark';
import { useIsMobile } from '@/hooks/useMediaQuery';

export function TopNav() {
  const isMobile = useIsMobile();
  const [menuOpen, setMenuOpen] = useState(false);
  const location = useLocation();

  // keep the selection query string when moving between destinations
  const withQuery = (path: string) => `${path}${location.search}`;

  return (
    <header className="app-nav no-print">
      <div className="app-nav__inner container">
        <NavLink to={withQuery('/')} className="app-nav__brand" onClick={() => setMenuOpen(false)}>
          <HelixMark size={30} />
          <span className="app-nav__brand-text">
            <strong>Capstone</strong>
            <span className="app-nav__brand-sub">Research Interface</span>
          </span>
        </NavLink>

        {isMobile ? (
          <>
            <button
              type="button"
              className="btn btn--subtle app-nav__toggle"
              aria-expanded={menuOpen}
              aria-controls="app-nav-menu"
              onClick={() => setMenuOpen((o) => !o)}
            >
              {menuOpen ? 'Close' : 'Menu'}
            </button>
            {menuOpen ? (
              <nav id="app-nav-menu" className="app-nav__menu" aria-label="Primary">
                <ul>
                  {NAV_ROUTES.map((r) => (
                    <li key={r.path}>
                      <NavLink
                        to={withQuery(r.path)}
                        end={r.path === '/'}
                        onClick={() => setMenuOpen(false)}
                        className={({ isActive }) => (isActive ? 'is-active' : undefined)}
                      >
                        <span className="app-nav__menu-label">{r.label}</span>
                        <span className="app-nav__menu-desc">{r.description}</span>
                      </NavLink>
                    </li>
                  ))}
                </ul>
              </nav>
            ) : null}
          </>
        ) : (
          <nav className="app-nav__links" aria-label="Primary">
            {NAV_ROUTES.map((r) => (
              <NavLink
                key={r.path}
                to={withQuery(r.path)}
                end={r.path === '/'}
                className={({ isActive }) => (isActive ? 'is-active' : undefined)}
              >
                {r.label}
              </NavLink>
            ))}
          </nav>
        )}
      </div>
    </header>
  );
}
