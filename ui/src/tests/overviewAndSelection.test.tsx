import { describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { renderWithProviders } from './helpers';
import { OverviewPage } from '@/pages/OverviewPage';
import { selectionFromParams, selectionToParams, isFiltering } from '@/lib/selection';

describe('Overview page', () => {
  it('states the frozen Phase 1 result with real numbers and the honest conclusion', async () => {
    renderWithProviders(<OverviewPage />, { route: '/' });
    expect(await screen.findByText('0.2356')).toBeInTheDocument();
    expect(screen.getByText('0.2047')).toBeInTheDocument();
    expect(screen.getByText('-0.0308')).toBeInTheDocument();
    expect(screen.getByText('170')).toBeInTheDocument();
    expect(screen.getByText('4,297')).toBeInTheDocument();
    expect(
      screen.getByText(/did not outperform the expression baseline/i),
    ).toBeInTheDocument();
  });

  it('does not expose lineage_mean / mlp_head / E1 as selectable models', async () => {
    renderWithProviders(<OverviewPage />, { route: '/' });
    await screen.findByText('0.2356');
    expect(screen.queryByRole('radio', { name: /lineage_mean/i })).not.toBeInTheDocument();
    expect(screen.queryByRole('radio', { name: /mlp_head/i })).not.toBeInTheDocument();
    // they may be *mentioned* in prose
    expect(screen.getByText(/evaluation-only controls/i)).toBeInTheDocument();
  });

  it('labels the two samples with different roles', async () => {
    renderWithProviders(<OverviewPage />, { route: '/' });
    await screen.findByText('0.2356');
    expect(screen.getByText(/Validation anchor/i)).toBeInTheDocument();
    expect(screen.getByText(/Exploratory external/i)).toBeInTheDocument();
  });

  it('carries a non-clinical banner', async () => {
    renderWithProviders(<OverviewPage />, { route: '/' });
    // banner is rendered by App, not the page; test the page's own caution
    expect(
      await screen.findByText(/Research software — not clinical guidance/i),
    ).toBeInTheDocument();
  });
});

describe('selection <-> URL', () => {
  it('round-trips a non-default selection', () => {
    const params = new URLSearchParams('sample=BG003082&model=ridge_head&q=cdk&evidence=cited&gene=1017');
    const sel = selectionFromParams(params);
    expect(sel).toEqual({
      sample: 'BG003082',
      model: 'ridge_head',
      search: 'cdk',
      evidence: 'cited',
      gene: '1017',
    });
    const out = selectionToParams(sel);
    expect(out.get('sample')).toBe('BG003082');
    expect(out.get('model')).toBe('ridge_head');
    expect(out.get('q')).toBe('cdk');
    expect(out.get('evidence')).toBe('cited');
    expect(out.get('gene')).toBe('1017');
  });

  it('omits default values from the query string', () => {
    const out = selectionToParams({
      sample: 'ACH-000364',
      model: 'ridge_pca',
      search: '',
      evidence: 'all',
      gene: null,
    });
    expect(out.toString()).toBe('');
  });

  it('falls back to defaults for bad values', () => {
    const sel = selectionFromParams(new URLSearchParams('sample=bogus&model=consensus&evidence=made-up&gene=abc'));
    expect(sel.sample).toBe('ACH-000364');
    expect(sel.model).toBe('ridge_pca');
    expect(sel.evidence).toBe('all');
    expect(sel.gene).toBeNull();
  });

  it('isFiltering is true only with a search or a non-all evidence filter', () => {
    expect(isFiltering({ sample: 'ACH-000364', model: 'ridge_pca', search: '', evidence: 'all', gene: null })).toBe(false);
    expect(isFiltering({ sample: 'ACH-000364', model: 'ridge_pca', search: 'x', evidence: 'all', gene: null })).toBe(true);
    expect(isFiltering({ sample: 'ACH-000364', model: 'ridge_pca', search: '', evidence: 'cited', gene: null })).toBe(true);
  });
});

// A tiny sanity check that user-event is wired up and the router provider works.
describe('smoke: interactive control', () => {
  it('overview CTA navigates without throwing', async () => {
    const user = userEvent.setup();
    renderWithProviders(<OverviewPage />, { route: '/' });
    await screen.findByText('0.2356');
    await user.click(screen.getByRole('link', { name: /Explore predicted dependencies/i }));
  });
});
