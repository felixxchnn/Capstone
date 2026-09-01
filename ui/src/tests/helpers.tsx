import type { ReactElement } from 'react';
import { render, waitFor, type RenderOptions } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';
import { DataSourceContext } from '@/app/dataSourceContext';
import { StaticCaseStudyDataSource } from '@/data/adapters/StaticCaseStudyDataSource';
import type { CapstoneDataSource } from '@/data/CapstoneDataSource';

export function makeDataSource(): CapstoneDataSource {
  return new StaticCaseStudyDataSource();
}

/** Wait until the ranking table has rendered its 25 rows. */
export async function waitForRankingRows(): Promise<HTMLTableRowElement[]> {
  await waitFor(() => {
    const rows = document.querySelectorAll('table.ranking-table tbody tr');
    if (rows.length !== 25) throw new Error(`ranking table has ${rows.length} rows`);
  });
  return Array.from(document.querySelectorAll('table.ranking-table tbody tr'));
}

export function visibleRankingRows(): HTMLTableRowElement[] {
  return Array.from(
    document.querySelectorAll<HTMLTableRowElement>('table.ranking-table tbody tr'),
  ).filter((r) => !r.hidden);
}

export function allRankingRanks(): number[] {
  return Array.from(document.querySelectorAll('table.ranking-table tbody tr')).map((r) =>
    Number(r.getAttribute('data-rank')),
  );
}

export function renderWithProviders(
  ui: ReactElement,
  opts: { route?: string; dataSource?: CapstoneDataSource } & RenderOptions = {},
) {
  const { route = '/', dataSource = makeDataSource(), ...rest } = opts;
  return render(
    <MemoryRouter initialEntries={[route]}>
      <DataSourceContext.Provider value={dataSource}>{ui}</DataSourceContext.Provider>
    </MemoryRouter>,
    rest,
  );
}
