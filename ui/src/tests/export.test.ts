import { describe, expect, it } from 'vitest';
import { StaticCaseStudyDataSource } from '@/data/adapters/StaticCaseStudyDataSource';
import { buildExportRows, rowsToCsv, rowsToJsonView } from '@/lib/export';
import type { EvidenceStatus } from '@/types/caseStudy';

const ds = new StaticCaseStudyDataSource();
const statusFor = () => ({ status: 'none_in_filtered_snapshot' as EvidenceStatus, count: 0 });

describe('export', () => {
  it('CSV export keeps all 25 rows in frozen rank order', async () => {
    const block = await ds.getModelRanking('ACH-000364', 'ridge_pca');
    const rows = buildExportRows({ sample: 'ACH-000364', model: 'ridge_pca', block, evidenceStatusFor: statusFor });
    expect(rows.map((r) => r.rank)).toEqual(Array.from({ length: 25 }, (_, i) => i + 1));
    const csv = rowsToCsv(rows);
    const lines = csv.trim().split('\r\n');
    expect(lines).toHaveLength(26); // header + 25
    expect(lines[0]).toContain('rank');
    expect(lines[0]).toContain('evidence_status');
    // first data line is rank 1
    expect(lines[1].split(',')[3]).toBe('1');
  });

  it('a visible-row filter drops rows but never renumbers ranks', async () => {
    const block = await ds.getModelRanking('ACH-000364', 'ridge_pca');
    // keep only ranks 2, 5, 9
    const keep = new Set([block.genes[1], block.genes[4], block.genes[8]].map((g) => g.entrez_id));
    const rows = buildExportRows({
      sample: 'ACH-000364',
      model: 'ridge_pca',
      block,
      evidenceStatusFor: statusFor,
      visibleEntrezIds: keep,
    });
    expect(rows.map((r) => r.rank)).toEqual([2, 5, 9]); // original ranks preserved, not 1,2,3
  });

  it('CSV includes sample, model, identifiers, predicted, observed availability, evidence class', async () => {
    const block = await ds.getModelRanking('ACH-000364', 'ridge_pca');
    const rows = buildExportRows({ sample: 'ACH-000364', model: 'ridge_pca', block, evidenceStatusFor: statusFor });
    const csv = rowsToCsv(rows);
    expect(csv).toContain('ACH-000364');
    expect(csv).toContain('ridge_pca');
    expect(csv).toContain('observed_available');
    expect(csv).toContain(rows[0].entrez_id);
  });

  it('external sample exports observed_available=false and no observed value', async () => {
    const block = await ds.getModelRanking('BG003082', 'ridge_head');
    const rows = buildExportRows({ sample: 'BG003082', model: 'ridge_head', block, evidenceStatusFor: statusFor });
    expect(rows.every((r) => r.observed_available === false)).toBe(true);
    expect(rows.every((r) => r.observed_geneeffect === null)).toBe(true);
  });

  it('JSON export view records the case-study hash and filter flag', async () => {
    const block = await ds.getModelRanking('ACH-000364', 'ridge_head');
    const rows = buildExportRows({ sample: 'ACH-000364', model: 'ridge_head', block, evidenceStatusFor: statusFor });
    const view = rowsToJsonView({
      sample: 'ACH-000364',
      model: 'ridge_head',
      block,
      rows,
      caseStudySha256: 'abc123',
      filtered: false,
    });
    expect(view.case_study_sha256).toBe('abc123');
    expect(view.view.filtered).toBe(false);
    expect(view.view.n_exported).toBe(25);
    expect(view.disclaimer.toLowerCase()).toContain('not therapeutic targets');
  });
});
