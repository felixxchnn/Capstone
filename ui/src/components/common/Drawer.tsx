import { useEffect, useRef, type ReactNode } from 'react';
import { trapFocus } from '@/lib/a11y';

/** Accessible slide-over dialog. role="dialog" + aria-modal, focus trap,
 *  Escape to close, restores focus to the trigger, closes on backdrop click. */
export function Drawer({
  open,
  onClose,
  title,
  children,
  labelledById,
}: {
  open: boolean;
  onClose: () => void;
  title: ReactNode;
  children: ReactNode;
  labelledById?: string;
}) {
  const panelRef = useRef<HTMLDivElement>(null);
  const lastFocused = useRef<HTMLElement | null>(null);
  const headingId = labelledById ?? 'drawer-title';

  useEffect(() => {
    if (!open) return;
    lastFocused.current = document.activeElement as HTMLElement | null;
    const panel = panelRef.current;
    if (!panel) return;
    const cleanup = trapFocus(panel, onClose);
    // focus the close button / first focusable
    const focusTarget =
      panel.querySelector<HTMLElement>('[data-autofocus]') ??
      panel.querySelector<HTMLElement>('button, a[href], input, [tabindex]');
    focusTarget?.focus();
    const prevOverflow = document.body.style.overflow;
    document.body.style.overflow = 'hidden';
    return () => {
      cleanup();
      document.body.style.overflow = prevOverflow;
      lastFocused.current?.focus?.();
    };
  }, [open, onClose]);

  if (!open) return null;

  return (
    <div className="drawer-overlay" onMouseDown={(e) => e.target === e.currentTarget && onClose()}>
      <div
        className="drawer-panel"
        role="dialog"
        aria-modal="true"
        aria-labelledby={headingId}
        ref={panelRef}
      >
        <div className="drawer-header">
          <h2 id={headingId} className="drawer-title">
            {title}
          </h2>
          <button
            type="button"
            className="btn btn--subtle drawer-close"
            onClick={onClose}
            data-autofocus
          >
            <svg width="14" height="14" viewBox="0 0 16 16" aria-hidden="true">
              <path d="M3 3l10 10M13 3L3 13" stroke="currentColor" strokeWidth="2" strokeLinecap="round" />
            </svg>
            Close
          </button>
        </div>
        <div className="drawer-body">{children}</div>
      </div>
    </div>
  );
}
