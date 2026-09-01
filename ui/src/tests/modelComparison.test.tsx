import { describe, expect, it } from 'vitest';
import { screen, within } from '@testing-library/react';
import { renderWithProviders } from './helpers';
import { ModelComparisonPage } from '@/pages/ModelComparisonPage';

describe('Model Comparison', () => {
  it('renders two independent top-25 lists, each ranked 1..25, no consensus', async () => {
    renderWithProviders(<ModelComparisonPage />, { route: '/compare' });
    await screen.findByText(/ridge_pca — top 25/i);

    const pcaTable = screen.getByRole('table', { name: /ridge_pca — top 25/i });
    const headTable = screen.getByRole('table', { name: /ridge_head — top 25/i });

    const pcaRanks = within(pcaTable)
      .getAllByRole('row')
      .slice(1)
      .map((r) => r.querySelector('.rank-chip')?.textContent);
    const headRanks = within(headTable)
      .getAllByRole('row')
      .slice(1)
      .map((r) => r.querySelector('.rank-chip')?.textContent);

    expect(pcaRanks).toEqual(Array.from({ length: 25 }, (_, i) => String(i + 1)));
    expect(headRanks).toEqual(Array.from({ length: 25 }, (_, i) => String(i + 1)));

    // exactly the two per-model tables + the shared-genes descriptive table —
    // there is no merged / combined "consensus" ranking table
    const tableCaptions = screen.getAllByRole('table').map((t) => t.querySelector('caption')?.textContent ?? '');
    expect(tableCaptions.some((c) => /ridge_pca — top 25/.test(c))).toBe(true);
    expect(tableCaptions.some((c) => /ridge_head — top 25/.test(c))).toBe(true);
    expect(tableCaptions.some((c) => /combined|consensus|merged/i.test(c))).toBe(false);
    expect(screen.getByText(/not new performance evaluations/i)).toBeInTheDocument();
  });

  it('overlap section is descriptive and labels rank differences', async () => {
    renderWithProviders(<ModelComparisonPage />, { route: '/compare' });
    const overlapHeading = await screen.findByText(/Genes in both top-25 lists/i);
    expect(overlapHeading).toBeInTheDocument();
    expect(screen.getByText(/Descriptive only\./i)).toBeInTheDocument();
  });
});
