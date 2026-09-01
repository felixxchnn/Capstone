import { useEffect, useReducer } from 'react';

export type AsyncState<T> =
  | { status: 'loading'; data: null; error: null }
  | { status: 'success'; data: T; error: null }
  | { status: 'error'; data: null; error: Error };

type Action<T> =
  | { type: 'reset' }
  | { type: 'success'; data: T }
  | { type: 'error'; error: Error };

function reducer<T>(_state: AsyncState<T>, action: Action<T>): AsyncState<T> {
  switch (action.type) {
    case 'reset':
      return { status: 'loading', data: null, error: null };
    case 'success':
      return { status: 'success', data: action.data, error: null };
    case 'error':
      return { status: 'error', data: null, error: action.error };
  }
}

/**
 * Run an async producer, tracking loading / success / error.
 *
 * The producer is given an AbortSignal. When `deps` change (or the component
 * unmounts, or `reload()` is called) the in-flight run is aborted and its result
 * is discarded — no stale writes for a superseded selection.
 *
 * `deps` must include everything the producer closes over (same contract as
 * `useEffect`). The producer is intentionally not stashed in a ref: that would
 * hide a stale-closure bug behind the dependency array.
 */
export function useAsync<T>(
  producer: (signal: AbortSignal) => Promise<T>,
  deps: React.DependencyList,
): AsyncState<T> & { reload: () => void } {
  const [state, dispatch] = useReducer(reducer<T>, {
    status: 'loading',
    data: null,
    error: null,
  } as AsyncState<T>);
  const [nonce, forceReload] = useReducer((n: number) => n + 1, 0);

  useEffect(() => {
    const controller = new AbortController();
    let superseded = false;
    const alive = () => !superseded && !controller.signal.aborted;

    // Reset to "loading" on a microtask so it is never a synchronous setState in
    // the effect body (which cascades renders).
    queueMicrotask(() => {
      if (alive()) dispatch({ type: 'reset' });
    });

    producer(controller.signal)
      .then((data) => {
        if (alive()) dispatch({ type: 'success', data });
      })
      .catch((err: unknown) => {
        if (!alive()) return;
        if (err instanceof DOMException && err.name === 'AbortError') return;
        dispatch({
          type: 'error',
          error: err instanceof Error ? err : new Error(String(err)),
        });
      });

    return () => {
      superseded = true;
      controller.abort();
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [...deps, nonce]);

  return { ...state, reload: () => forceReload() };
}
