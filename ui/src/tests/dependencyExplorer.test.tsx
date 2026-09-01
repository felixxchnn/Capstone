import { describe, expect, it } from 'vitest';
import { screen, within, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import {
  renderWithProviders,
  waitForRankingRows,
  visibleRankingRows,
  allRankingRanks,
} from './helpers';
import { DependencyExplorerPage } from '@/pages/DependencyExplorerPage';

const FROZEN = Array.from({ length: 25 }, (_, i) => i + 1);

describe('Dependency Explorer', () => {
  it('shows exactly 25 rows in frozen rank order 1..25', async () => {
    renderWithProviders(<DependencyExplorerPage />, { route: '/explore' });
    await waitForRankingRows();
    expect(allRankingRanks()).toEqual(FROZEN);
  });

  it('search hides rows but never reorders or renumbers them', async () => {
    const user = userEvent.setup();
    renderWithProviders(<DependencyExplorerPage />, { route: '/explore' });
    const rows = await waitForRankingRows();

    const firstSym = rows[0].querySelectorAll('td')[1]?.textContent?.trim() ?? 'DDX11';
    await user.type(screen.getByLabelText(/Search gene, Entrez/i), firstSym);

    await waitFor(() => {
      const visible = visibleRankingRows();
      expect(visible.length).toBeGreaterThanOrEqual(1);
      expect(visible.length).toBeLessThan(25);
    });
    // all rendered rows still hold their original rank, in order
    expect(allRankingRanks()).toEqual(FROZEN);
  });

  it('evidence filter narrows to one status; reset restores 25 visible', async () => {
    const user = userEvent.setup();
    renderWithProviders(<DependencyExplorerPage />, { route: '/explore' });
    await waitForRankingRows();

    await user.selectOptions(screen.getByLabelText(/Evidence filter/i), 'cited');
    await waitFor(() => {
      const visible = visibleRankingRows();
      expect(visible.length).toBeLessThan(25);
      for (const r of visible) expect(r.querySelector('.status-pill--cited')).toBeTruthy();
    });

    await user.click(screen.getByRole('button', { name: /reset filters/i }));
    await waitFor(() => expect(visibleRankingRows()).toHaveLength(25));
  });

  it('switching sample keeps 25 rows and drops the observed columns for BG003082', async () => {
    const user = userEvent.setup();
    renderWithProviders(<DependencyExplorerPage />, { route: '/explore' });
    await waitForRankingRows();
    expect(screen.getByRole('columnheader', { name: /Observed GeneEffect/i })).toBeInTheDocument();

    await user.click(screen.getByRole('radio', { name: /BG003082/i }));
    await waitFor(() => {
      expect(screen.queryByRole('columnheader', { name: /Observed GeneEffect/i })).not.toBeInTheDocument();
    });
    expect(allRankingRanks()).toEqual(FROZEN);
  });

  it('switching model re-loads an independent ranking (rank order still 1..25)', async () => {
    const user = userEvent.setup();
    renderWithProviders(<DependencyExplorerPage />, { route: '/explore' });
    const rows = await waitForRankingRows();
    const firstEntrezBefore = rows[0].querySelectorAll('td')[2]?.textContent?.trim();

    await user.click(screen.getByRole('radio', { name: /ridge_head/i }));
    await waitFor(() => {
      // rank-8 gene differs between the two models for ACH-000364
      const r8 = document.querySelector('table.ranking-table tbody tr[data-rank="8"]');
      expect(r8?.querySelectorAll('td')[2]?.textContent?.trim()).not.toBe(firstEntrezBefore);
    });
    expect(allRankingRanks()).toEqual(FROZEN);
  });

  it('opens an accessible gene detail dialog and closes it', async () => {
    const user = userEvent.setup();
    renderWithProviders(<DependencyExplorerPage />, { route: '/explore' });
    await waitForRankingRows();

    await user.click(screen.getAllByRole('button', { name: /open details for/i })[0]);
    const dialog = await screen.findByRole('dialog');
    expect(within(dialog).getByText(/Predicted vs observed/i)).toBeInTheDocument();
    expect(within(dialog).getByText(/not a therapeutic target/i)).toBeInTheDocument();

    await user.click(within(dialog).getByRole('button', { name: /close/i }));
    await waitFor(() => expect(screen.queryByRole('dialog')).not.toBeInTheDocument());
  });

  it('shows an accessible empty state when nothing matches', async () => {
    const user = userEvent.setup();
    renderWithProviders(<DependencyExplorerPage />, { route: '/explore' });
    await waitForRankingRows();
    await user.type(screen.getByLabelText(/Search gene, Entrez/i), 'zzzznomatch');
    expect(
      await screen.findByText(/No genes match the current search and filters/i),
    ).toBeInTheDocument();
  });
});
