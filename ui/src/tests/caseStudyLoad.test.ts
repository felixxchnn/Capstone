import { readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { parseCaseStudy, CaseStudyParseError } from '@/data/schemas/caseStudy';
import generated from '@/data/case_study.generated.json';
import { MODEL_IDS, SAMPLE_IDS } from '@/types/caseStudy';

const uiRoot = resolve(__dirname, '..', '..');
const repoRoot = resolve(uiRoot, '..');

describe('case study data loading', () => {
  it('the synced copy is byte-identical to the committed case_study.json', () => {
    const committed = readFileSync(resolve(repoRoot, 'data', 'processed', 'case_study.json'));
    const synced = readFileSync(resolve(uiRoot, 'src', 'data', 'case_study.generated.json'));
    expect(synced.equals(committed)).toBe(true);
  });

  it('the synced copy matches the pinned SHA-256 (mirrors data-integrity-hashes.md)', () => {
    const synced = readFileSync(resolve(uiRoot, 'src', 'data', 'case_study.generated.json'));
    const sha = createHash('sha256').update(synced).digest('hex');
    const pinned = readFileSync(
      resolve(uiRoot, 'src', 'data', 'case_study.expected-sha256.txt'),
      'utf8',
    ).trim();
    expect(sha).toBe(pinned);
    expect(pinned).toBe('a962c01a5b65a6ef579ea57dced67048bf9016ba0f66aab2355cf1f054796e8c');
  });

  it('parses under the schema guard', () => {
    const cs = parseCaseStudy(generated);
    expect(cs.schema_version).toBe('case-study/1');
    for (const sid of SAMPLE_IDS) {
      for (const mid of MODEL_IDS) {
        expect(cs.rankings[sid][mid].genes).toHaveLength(25);
      }
    }
  });

  it('rejects a malformed case study with a path', () => {
    const bad = structuredClone(generated) as Record<string, unknown>;
    (bad.rankings as Record<string, Record<string, { genes: unknown[] }>>)['ACH-000364'][
      'ridge_pca'
    ].genes.pop();
    expect(() => parseCaseStudy(bad)).toThrow(CaseStudyParseError);
    try {
      parseCaseStudy(bad);
    } catch (err) {
      expect((err as CaseStudyParseError).path).toContain('rankings.ACH-000364.ridge_pca.genes');
    }
  });

  it('rejects an out-of-order rank', () => {
    const bad = structuredClone(generated) as {
      rankings: Record<string, Record<string, { genes: { rank: number }[] }>>;
    };
    bad.rankings['ACH-000364']['ridge_pca'].genes[3].rank = 99;
    expect(() => parseCaseStudy(bad)).toThrow(/expected rank 4/);
  });

  it('rejects observed values on the external sample ranking', () => {
    const bad = structuredClone(generated) as {
      rankings: Record<string, Record<string, { genes: Record<string, unknown>[] }>>;
    };
    bad.rankings['BG003082']['ridge_head'].genes[0].observed_geneeffect = -1;
    expect(() => parseCaseStudy(bad)).toThrow(/must not carry observed values/);
  });
});
