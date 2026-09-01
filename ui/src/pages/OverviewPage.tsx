import { Link } from 'react-router-dom';
import { useDataSource } from '@/app/dataSourceContext';
import { useAsync } from '@/hooks/useAsync';
import { HeroHelix } from '@/components/overview/HeroHelix';
import { ResultCards } from '@/components/overview/ResultCards';
import { Callout, Card, EmptyState, Spinner } from '@/components/common/primitives';
import { SampleRoleBadge } from '@/components/samples/SampleRoleBadge';
import { spearman, signedDelta } from '@/lib/format';

export function OverviewPage() {
  const ds = useDataSource();
  const state = useAsync((_s) => ds.getProjectSummary(), [ds]);

  return (
    <div className="overview">
      <section className="hero" aria-labelledby="hero-title">
        <div className="hero__inner container">
          <div className="hero__copy">
            <p className="hero__eyebrow">Computational biology · research demonstration</p>
            <h1 id="hero-title">
              Do frozen Geneformer embeddings beat a ridge-on-PCA-expression baseline at
              predicting CRISPR gene dependency?
            </h1>
            <p className="hero__lede">
              This interface shows the <strong>committed predictions</strong> from a
              high-school capstone that put that question to a controlled test across
              held-out cancer cell lines. It performs no model inference — every number is
              read from a frozen, hash-pinned artifact.
            </p>
            <p className="hero__answer">
              <span className="hero__answer-tag">Answer</span> No — and the deficit is
              quantified. The Geneformer head did not outperform the expression baseline.
            </p>
            <div className="hero__cta">
              <Link className="btn" to="/explore">
                Explore predicted dependencies →
              </Link>
              <Link className="btn btn--ghost hero__cta-ghost" to="/methods">
                Methods &amp; limitations
              </Link>
            </div>
          </div>
          <div className="hero__art">
            <HeroHelix />
            <p className="sr-only">
              Decorative illustration of a DNA double helix with molecular nodes.
            </p>
          </div>
        </div>
      </section>

      <div className="container page">
        {state.status === 'loading' ? <Spinner label="Loading project summary" /> : null}
        {state.status === 'error' ? (
          <EmptyState title="Could not load the project summary">{state.error.message}</EmptyState>
        ) : null}

        {state.status === 'success' ? (
          <>
            <section aria-labelledby="headline-title" className="stack-lg">
              <h2 id="headline-title">The frozen Phase 1 result</h2>
              <ResultCards summary={state.data} />
              <Callout tone="info" title="What the numbers mean">
                Each model was scored by <strong>per-target Spearman correlation</strong>{' '}
                between predicted and observed CRISPR GeneEffect, computed across the{' '}
                {state.data.phase1.nValidationLines} held-out validation cell lines for each
                of the {state.data.phase1.nTargets.toLocaleString('en-US')} targets, then
                averaged. ridge_pca scored {spearman(state.data.phase1.ridgePcaMeanSpearman)};
                ridge_head scored {spearman(state.data.phase1.ridgeHeadMeanSpearman)}; the
                paired difference is {signedDelta(state.data.phase1.deltaHeadMinusBaseline)}{' '}
                and its 95% confidence interval excludes zero.
              </Callout>
            </section>

            <section aria-labelledby="plain-title" className="stack-lg">
              <h2 id="plain-title">In plain language</h2>
              <div className="grid grid--2">
                <Card title="The research question">
                  <p className="small">
                    Can a large pretrained single-cell transformer (Geneformer), used as a
                    frozen feature extractor, predict which genes a cancer cell line depends
                    on better than a simple linear model on the cell line's gene-expression
                    profile? The comparison is held fixed and pre-specified.
                  </p>
                </Card>
                <Card title="CRISPR GeneEffect">
                  <p className="small">
                    A per-gene number from a genome-wide CRISPR knockout screen. Around 0
                    means losing the gene had little effect; strongly negative means the
                    cell line depended on that gene to survive.{' '}
                    <strong>More negative = stronger dependency.</strong>
                  </p>
                </Card>
                <Card title="ridge_pca and ridge_head">
                  <p className="small">
                    <code>ridge_pca</code> — ridge regression on 200 principal components of
                    the expression matrix. <code>ridge_head</code> — the same ridge
                    regression on the 768-dimensional frozen Geneformer embedding of the
                    same cell line. Same estimator, same split, same metric; only the input
                    representation changes.
                  </p>
                </Card>
                <Card title="Evaluation vs prediction">
                  <p className="small">
                    <strong>Evaluation</strong> compares predictions to measured outcomes on
                    held-out data to estimate performance (Phase 1). A{' '}
                    <strong>prediction</strong> is a model output for a sample where the
                    outcome may be unknown. This interface shows predictions for two
                    samples; only one has measured outcomes, attached after the fact.
                  </p>
                </Card>
              </div>
            </section>

            <section aria-labelledby="samples-title" className="stack-lg">
              <h2 id="samples-title">Two samples, two different roles</h2>
              <div className="grid grid--2">
                <Card title={<>ACH-000364 <SampleRoleBadge sample="ACH-000364" /></>}>
                  <p className="small">
                    A held-out DepMap validation cell line (U-2 OS). Real observed CRISPR
                    outcomes exist and were attached to the ranked genes{' '}
                    <strong>after</strong> prediction and ranking, only to verify the
                    pipeline. One cell line is not a performance estimate.
                  </p>
                </Card>
                <Card title={<>BG003082 <SampleRoleBadge sample="BG003082" /></>}>
                  <p className="small">
                    A real primary osteosarcoma tumour (bulk RNA-seq, CC0), absent from
                    every DepMap split. <strong>No CRISPR outcome exists.</strong> Bulk
                    tumour tissue is a domain shift from the cultured cell lines the models
                    were trained and validated on — every prediction for it is exploratory.
                  </p>
                </Card>
              </div>
            </section>

            <Callout tone="caution" title="Research software — not clinical guidance">
              Nothing here is a treatment recommendation, a treatment-response prediction, a
              diagnosis, or clinical advice. Drug–gene interaction evidence is retrieval of
              prior records and does not establish efficacy. Other evaluation-only controls
              (a lineage-mean baseline, an MLP head, and the E1 random-projection control)
              exist but have no sample-level rankings and are described only in Methods.
            </Callout>
          </>
        ) : null}
      </div>
    </div>
  );
}
