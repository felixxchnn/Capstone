import { useMemo, type ReactNode } from 'react';
import type { CapstoneDataSource } from '@/data/CapstoneDataSource';
import { StaticCaseStudyDataSource } from '@/data/adapters/StaticCaseStudyDataSource';
import { DataSourceContext } from '@/app/dataSourceContext';

// The single seam between UI and scientific data. To move to a Python backend,
// pass an `ApiDataSource` as `value` here (see
// src/data/adapters/ApiDataSource.contract.md). Components call `useDataSource()`
// and never import an adapter directly.
export function DataSourceProvider({
  children,
  value,
}: {
  children: ReactNode;
  value?: CapstoneDataSource;
}) {
  const source = useMemo(() => value ?? new StaticCaseStudyDataSource(), [value]);
  return <DataSourceContext.Provider value={source}>{children}</DataSourceContext.Provider>;
}
