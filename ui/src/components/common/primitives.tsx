import type { ReactNode, HTMLAttributes } from 'react';

export function VisuallyHidden({ children }: { children: ReactNode }) {
  return <span className="sr-only">{children}</span>;
}

export function SkipLink({ targetId = 'main' }: { targetId?: string }) {
  return (
    <a className="skip-link" href={`#${targetId}`}>
      Skip to main content
    </a>
  );
}

export function Card({
  title,
  children,
  as: Tag = 'section',
  ...rest
}: {
  title?: ReactNode;
  children: ReactNode;
  as?: 'section' | 'div' | 'article';
} & Omit<HTMLAttributes<HTMLElement>, 'title'>) {
  return (
    <Tag className="card" {...rest}>
      {title ? <h3 className="card__title">{title}</h3> : null}
      {children}
    </Tag>
  );
}

export type CalloutTone = 'info' | 'caution' | 'danger';

export function Callout({
  tone = 'info',
  title,
  children,
  role,
}: {
  tone?: CalloutTone;
  title?: ReactNode;
  children: ReactNode;
  role?: 'note' | 'alert' | 'status';
}) {
  return (
    <div className={`callout callout--${tone}`} role={role ?? (tone === 'danger' ? 'alert' : 'note')}>
      {title ? <p className="callout__title">{title}</p> : null}
      <div>{children}</div>
    </div>
  );
}

export function EmptyState({
  title,
  children,
  action,
}: {
  title: ReactNode;
  children?: ReactNode;
  action?: ReactNode;
}) {
  return (
    <div className="empty-state" role="status">
      <p className="empty-state__title">{title}</p>
      {children ? <p className="small">{children}</p> : null}
      {action ? <div style={{ marginTop: 'var(--sp-3)' }}>{action}</div> : null}
    </div>
  );
}

export function Spinner({ label = 'Loading' }: { label?: string }) {
  return (
    <div className="spinner" role="status" aria-live="polite">
      <svg width="20" height="20" viewBox="0 0 24 24" aria-hidden="true" className="spinner__svg">
        <circle
          cx="12"
          cy="12"
          r="9"
          fill="none"
          stroke="currentColor"
          strokeWidth="3"
          strokeLinecap="round"
          strokeDasharray="42 14"
        />
      </svg>
      <span>{label}…</span>
    </div>
  );
}
