import { useCallback, useMemo } from 'react';
import { useSearchParams } from 'react-router-dom';
import {
  selectionFromParams,
  selectionToParams,
  type Selection,
} from '@/lib/selection';

/**
 * The shared sample / model / gene / search / evidence-filter selection,
 * persisted in the URL query string so it survives navigation between pages and
 * is shareable.
 */
export function useSelection(): {
  selection: Selection;
  setSelection: (patch: Partial<Selection>) => void;
  resetFilters: () => void;
} {
  const [params, setParams] = useSearchParams();
  const selection = useMemo(() => selectionFromParams(params), [params]);

  const setSelection = useCallback(
    (patch: Partial<Selection>) => {
      setParams(
        (prev) => {
          const current = selectionFromParams(prev);
          return selectionToParams({ ...current, ...patch }, prev);
        },
        { replace: true },
      );
    },
    [setParams],
  );

  const resetFilters = useCallback(() => {
    setSelection({ search: '', evidence: 'all' });
  }, [setSelection]);

  return { selection, setSelection, resetFilters };
}
