import { useId, useState, type ReactNode } from 'react';

/** Accessible expand/collapse. A real button toggles an aria-controlled region
 *  (works without JS-free <details> quirks and is easy to drive in tests). */
export function Disclosure({
  summary,
  children,
  defaultOpen = false,
  id,
}: {
  summary: ReactNode;
  children: ReactNode;
  defaultOpen?: boolean;
  id?: string;
}) {
  const reactId = useId();
  const panelId = `${id ?? reactId}-panel`;
  const [open, setOpen] = useState(defaultOpen);
  return (
    <div className="disclosure">
      <button
        type="button"
        className="disclosure__summary"
        aria-expanded={open}
        aria-controls={panelId}
        onClick={() => setOpen((o) => !o)}
      >
        <svg
          className="disclosure__chevron"
          width="14"
          height="14"
          viewBox="0 0 16 16"
          aria-hidden="true"
        >
          <path d="M5 3l6 5-6 5" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
        </svg>
        <span>{summary}</span>
      </button>
      <div id={panelId} className="disclosure__body" hidden={!open}>
        {children}
      </div>
    </div>
  );
}
