import { Component, type ErrorInfo, type ReactNode } from 'react';

interface Props {
  children: ReactNode;
  /** Rendered instead of children on error. Receives the error + a reset fn. */
  fallback: (error: Error, reset: () => void) => ReactNode;
  /** Optional label for logging context. */
  label?: string;
}

interface State {
  error: Error | null;
}

/** Generic error boundary. Used to isolate the heavy Mol* viewer so a WebGL /
 *  render failure never blanks the whole page. */
export class ErrorBoundary extends Component<Props, State> {
  override state: State = { error: null };

  static getDerivedStateFromError(error: Error): State {
    return { error };
  }

  override componentDidCatch(error: Error, info: ErrorInfo): void {
    console.error(
      `[ErrorBoundary${this.props.label ? ` ${this.props.label}` : ''}]`,
      error,
      info.componentStack,
    );
  }

  reset = (): void => this.setState({ error: null });

  override render(): ReactNode {
    if (this.state.error) {
      return this.props.fallback(this.state.error, this.reset);
    }
    return this.props.children;
  }
}
