import type { RankedGene, SampleId } from '@/types/caseStudy';
import { geneEffect, geneEffectFull, integer } from '@/lib/format';
import { Callout } from '@/components/common/primitives';

/** Predicted-vs-observed panel for one gene.
 *  ACH-000364: committed observed values, with the "attached after ranking" note.
 *  BG003082: outcome unavailable — no invented values, no empty measured chart. */
export function PredictedObservedStatus({
  sample,
  gene,
}: {
  sample: SampleId;
  gene: RankedGene;
}) {
  if (sample === 'BG003082') {
    return (
      <div className="pred-obs pred-obs--unavailable">
        <div className="pred-obs__row">
          <div>
            <span className="pred-obs__label">Predicted GeneEffect</span>
            <span className="pred-obs__value tnum" title={geneEffectFull(gene.predicted_geneeffect)}>
              {geneEffect(gene.predicted_geneeffect)}
            </span>
          </div>
          <div>
            <span className="pred-obs__label">Observed CRISPR outcome</span>
            <span className="pred-obs__value pred-obs__value--muted">Outcome unavailable</span>
          </div>
        </div>
        <Callout tone="caution" title="Exploratory external prediction">
          BG003082 has no CRISPR screen. No observed value is loaded, invented, or
          computed, and this prediction is never compared against a measurement. Bulk
          primary-tumour tissue is a domain shift from the cultured cell lines the models
          saw during training and validation.
        </Callout>
      </div>
    );
  }

  const hasObs = gene.observed_geneeffect !== undefined;
  return (
    <div className="pred-obs">
      <div className="pred-obs__row">
        <div>
          <span className="pred-obs__label">Predicted GeneEffect</span>
          <span className="pred-obs__value tnum" title={geneEffectFull(gene.predicted_geneeffect)}>
            {geneEffect(gene.predicted_geneeffect)}
          </span>
          <span className="pred-obs__sub">predicted rank {gene.rank} of 25 displayed</span>
        </div>
        <div>
          <span className="pred-obs__label">Observed GeneEffect</span>
          <span className="pred-obs__value tnum" title={geneEffectFull(gene.observed_geneeffect)}>
            {hasObs ? geneEffect(gene.observed_geneeffect) : 'not measured'}
          </span>
          <span className="pred-obs__sub">
            {gene.observed_rank === undefined
              ? 'no observed rank'
              : `observed rank ${integer(gene.observed_rank)} of 4,297`}
          </span>
        </div>
      </div>
      <p className="small muted">
        Observed CRISPR GeneEffect values for ACH-000364 are real, and were attached to
        the already-ranked genes <strong>after</strong> prediction and ranking. They never
        influenced selection, model choice, or order. One cell line is a pipeline check,
        not a performance estimate.
      </p>
    </div>
  );
}
