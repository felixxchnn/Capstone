import { useMemo } from 'react';
import { useAsync } from '@/hooks/useAsync';
import { useDataSource } from '@/app/dataSourceContext';
import type {
  EvidenceBucket,
  EvidenceStatus,
  ModelId,
  ModelRankingBlock,
  RankedGene,
  SampleId,
} from '@/types/caseStudy';
import { matchesSearch, type Selection } from '@/lib/selection';

export interface GeneRowView {
  gene: RankedGene;
  evidenceStatus: EvidenceStatus;
  evidenceCount: number;
  evidence: EvidenceBucket | null;
  /** searchable haystack parts */
  drugs: string[];
  sources: string[];
  /** does this row pass the current search + evidence filter? */
  visible: boolean;
}

export interface RankingView {
  sample: SampleId;
  model: ModelId;
  block: ModelRankingBlock;
  rows: GeneRowView[];
  visibleCount: number;
  hasObserved: boolean;
}

/**
 * Load one sample/model ranking plus its per-gene evidence, and compute which
 * rows pass the active search + evidence filter. The frozen `rank` order is
 * NEVER changed here — `rows` is always in rank order 1..25; filtering only
 * flips `visible`.
 */
export function useRankingView(
  sample: SampleId,
  model: ModelId,
  selection: Pick<Selection, 'search' | 'evidence'>,
): ReturnType<typeof useAsync<{ block: ModelRankingBlock; buckets: Record<string, EvidenceBucket | null> }>> & {
  view: RankingView | null;
} {
  const ds = useDataSource();

  const loaded = useAsync(
    async (signal) => {
      const block = await ds.getModelRanking(sample, model);
      const entrezIds = block.genes.map((g) => g.entrez_id);
      const buckets: Record<string, EvidenceBucket | null> = {};
      await Promise.all(
        entrezIds.map(async (id) => {
          if (signal.aborted) return;
          buckets[id] = await ds.getGeneEvidence(id);
        }),
      );
      return { block, buckets };
    },
    [ds, sample, model],
  );

  const view = useMemo<RankingView | null>(() => {
    if (loaded.status !== 'success') return null;
    const { block, buckets } = loaded.data;
    const rows: GeneRowView[] = block.genes.map((gene) => {
      const bucket = buckets[gene.entrez_id] ?? null;
      const status: EvidenceStatus = bucket?.evidence_status ?? 'none_in_filtered_snapshot';
      const drugs = bucket ? bucket.records.map((r) => r.drug_name || r.drug_claim_name) : [];
      const sources = bucket ? bucket.records.map((r) => r.interaction_source) : [];
      const passesSearch = matchesSearch(selection.search, {
        symbol: gene.symbol,
        entrez: gene.entrez_id,
        drugs,
        sources,
      });
      const passesEvidence = selection.evidence === 'all' || selection.evidence === status;
      return {
        gene,
        evidenceStatus: status,
        evidenceCount: bucket?.n_records ?? 0,
        evidence: bucket,
        drugs,
        sources,
        visible: passesSearch && passesEvidence,
      };
    });
    return {
      sample,
      model,
      block,
      rows,
      visibleCount: rows.filter((r) => r.visible).length,
      hasObserved: sample === 'ACH-000364',
    };
  }, [loaded, sample, model, selection.search, selection.evidence]);

  return { ...loaded, view };
}
