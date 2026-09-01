import { describe, expect, it } from 'vitest';
import { StaticCaseStudyDataSource } from '@/data/adapters/StaticCaseStudyDataSource';

describe('StaticCaseStudyDataSource', () => {
  const ds = new StaticCaseStudyDataSource();

  it('exposes the frozen Phase 1 headline', async () => {
    const s = await ds.getProjectSummary();
    expect(s.phase1.nValidationLines).toBe(170);
    expect(s.phase1.nTargets).toBe(4297);
    expect(s.phase1.ridgePcaMeanSpearman).toBe(0.2356);
    expect(s.phase1.ridgeHeadMeanSpearman).toBe(0.2047);
    expect(s.phase1.deltaHeadMinusBaseline).toBe(-0.0308);
    expect(s.phase1.conclusion.toLowerCase()).toContain('did not outperform');
  });

  it('returns two samples with distinct, non-equivalent roles', async () => {
    const samples = await ds.getSamples();
    expect(samples.map((s) => s.id)).toEqual(['ACH-000364', 'BG003082']);
    const anchor = samples[0];
    const external = samples[1];
    expect(anchor.hasObservedOutcome).toBe(true);
    expect(anchor.isExternal).toBe(false);
    expect(external.hasObservedOutcome).toBe(false);
    expect(external.isExternal).toBe(true);
    expect(anchor.outcomeStatus).not.toBe(external.outcomeStatus);
  });

  it('returns independent 25-row rankings per model, in frozen order', async () => {
    for (const sample of ['ACH-000364', 'BG003082'] as const) {
      for (const model of ['ridge_pca', 'ridge_head'] as const) {
        const block = await ds.getModelRanking(sample, model);
        expect(block.genes).toHaveLength(25);
        expect(block.genes.map((g) => g.rank)).toEqual(
          Array.from({ length: 25 }, (_, i) => i + 1),
        );
      }
    }
  });

  it('never merges the two model rankings (the ordered lists genuinely differ)', async () => {
    const pca = await ds.getModelRanking('ACH-000364', 'ridge_pca');
    const head = await ds.getModelRanking('ACH-000364', 'ridge_head');
    const pcaIds = pca.genes.map((g) => g.entrez_id);
    const headIds = head.genes.map((g) => g.entrez_id);
    expect(pcaIds).not.toEqual(headIds);
    // and the two sets are not identical either
    expect(new Set(pcaIds)).not.toEqual(new Set(headIds));
  });

  it('attaches observed values only for ACH-000364', async () => {
    const anchor = await ds.getModelRanking('ACH-000364', 'ridge_pca');
    const external = await ds.getModelRanking('BG003082', 'ridge_pca');
    expect(anchor.genes[0].observed_geneeffect).toBeTypeOf('number');
    expect(anchor.genes[0].observed_rank).toBeTypeOf('number');
    expect(external.genes[0].observed_geneeffect).toBeUndefined();
    expect(external.genes[0].observed_rank).toBeUndefined();
  });

  it('returns evidence buckets by Entrez, or null', async () => {
    const anchor = await ds.getModelRanking('ACH-000364', 'ridge_pca');
    const withEvidence = await ds.getGeneEvidence(anchor.genes.find((g) => g.symbol === 'DDX11')!.entrez_id);
    expect(withEvidence).not.toBeNull();
    expect(['cited', 'source_only', 'none_in_filtered_snapshot']).toContain(
      withEvidence!.evidence_status,
    );
    expect(await ds.getGeneEvidence('0000-not-real')).toBeNull();
  });

  it('returns a defensive copy (mutation does not leak back)', async () => {
    const a = await ds.getModelRanking('ACH-000364', 'ridge_pca');
    a.genes[0].predicted_geneeffect = 999;
    const b = await ds.getModelRanking('ACH-000364', 'ridge_pca');
    expect(b.genes[0].predicted_geneeffect).not.toBe(999);
  });

  it('rejects a malformed pinned SHA', () => {
    expect(() => new StaticCaseStudyDataSource({ expectedSha256: 'nope' })).toThrow(/malformed/);
  });

  it('structure metadata carries the Entrez id and human taxonomy only', async () => {
    const hint = await ds.getStructureMetadata('1017');
    expect(hint.entrezId).toBe('1017');
    expect(hint.humanTaxonomyId).toBe('9606');
    expect(hint.note.toLowerCase()).toContain('no expression');
  });
});
