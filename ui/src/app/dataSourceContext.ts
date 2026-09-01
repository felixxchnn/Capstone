import { createContext, useContext } from 'react';
import type { CapstoneDataSource } from '@/data/CapstoneDataSource';

export const DataSourceContext = createContext<CapstoneDataSource | null>(null);

export function useDataSource(): CapstoneDataSource {
  const ctx = useContext(DataSourceContext);
  if (!ctx) {
    throw new Error('useDataSource must be used within <DataSourceProvider>');
  }
  return ctx;
}
