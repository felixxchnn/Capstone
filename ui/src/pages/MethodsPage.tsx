import { useDataSource } from '@/app/dataSourceContext';
import { useAsync } from '@/hooks/useAsync';
import { Disclosure } from '@/components/common/Disclosure';
import { Callout, EmptyState, Spinner } from '@/components/common/primitives';
import { SampleInfoPanel } from '@/components/samples/SampleInfoPanel';
import {
  EVIDENCE_NOT_EFFICACY,
  NON_CLINICAL,
  STRUCTURE_NOT_EVIDENCE,
  shortSha,
} from '@/lib/format';

export function MethodsPage() {
  const ds = useDataSource();
  const state = useAsync(
    async (_s) => {
      const [summary, raw, anchorMeta, externalMeta] = await Promise.all([
        ds.getProjectSummary(),
        ds.getRawCaseStudy(),
        ds.getSampleMetadata('ACH-000364'),
        ds.getSampleMetadata('BG003082'),
      ]);
      return { summary, raw, anchorMeta, externalMeta };
    },
    [ds],
  );

  if (state.status === 'loading') {
    return (
      <div className="container page">
        <Spinner label="Loading methods and provenance" />
      </div>
    );
  }
  if (state.status === 'error') {
    return (
      <div className="container page">
        <EmptyState title="Could not load methods">{state.error.message}</EmptyState>
      </div>
    );
  }

  const { summary, raw, anchorMeta, externalMeta } = state.data;
  const models = summary.models;
  const evRet = raw.drug_gene_interaction_evidence.retrieval;
  const osteo = summary.osteosarcoma;

  return (
    <div className="container page">
      <header className="stack">
        <h1>Methods, Provenance &amp; Limitations</h1>
        <p className="lede">
          Every scientific value in this interface is read from the committed{' '}
          <span className="mono">case_study.json</span> (schema {summary.schemaVersion},
          sha256 <span className="mono">{shortSha(summary.caseStudySha256, 12)}</span>),
          itself built from hash-pinned Phase 1 artifacts. Nothing is recomputed here.
        </p>
      </header>

      <section aria-labelledby="sample-info-title" className="stack-lg">
        <h2 id="sample-info-title">Sample information</h2>
        <div className="grid grid--2">
          <div className="card">
            <h3 className="card__title">ACH-000364</h3>
            <SampleInfoPanel meta={anchorMeta} />
          </div>
          <div className="card">
            <h3 className="card__title">BG003082</h3>
            <SampleInfoPanel meta={externalMeta} />
          </div>
        </div>
      </section>

      <section aria-labelledby="methods-title" className="stack">
        <h2 id="methods-title">Methods &amp; provenance</h2>

        <Disclosure summary="Model pipelines" defaultOpen>
          {(['ridge_pca', 'ridge_head'] as const).map((mk) => (
            <div key={mk} className="stack" style={{ marginBottom: 'var(--sp-3)' }}>
              <h4>{mk}</h4>
              <dl className="kv">
                <dt>Pipeline</dt>
                <dd className="mono">{models[mk].pipeline}</dd>
                <dt>Features / targets</dt>
                <dd>
                  {models[mk].n_features.toLocaleString('en-US')} /{' '}
                  {models[mk].n_targets.toLocaleString('en-US')}
                </dd>
                <dt>Frozen alpha</dt>
                <dd>
                  {models[mk].frozen_alpha.value} — {models[mk].frozen_alpha.selection}; from{' '}
                  <span className="mono">{models[mk].frozen_alpha.source}</span>
                </dd>
                <dt>Artifact manifest sha256</dt>
                <dd className="mono">{models[mk].manifest_sha256}</dd>
              </dl>
            </div>
          ))}
        </Disclosure>

        <Disclosure summary="Reconstructed fitted-state status">
          <p className="small">
            <strong>{models.ridge_pca.provenance_status}.</strong>{' '}
            {models.ridge_pca.not_original_fitted_objects} They reproduce every committed
            Phase 1 validation statistic exactly at the recorded precision (see{' '}
            <span className="mono">reconstruct_fitted.py --validate</span> in the repository).
          </p>
        </Disclosure>

        <Disclosure summary="Identifier mapping &amp; imputation (BG003082)">
          <SampleInfoPanel meta={externalMeta} />
        </Disclosure>

        <Disclosure summary="Evidence snapshot provenance">
          <dl className="kv">
            <dt>Method</dt>
            <dd className="small">{evRet.method}</dd>
            <dt>Retrieved after ranking froze</dt>
            <dd>{String(evRet.retrieved_after_top_n_frozen)}</dd>
            <dt>Availability affected ranking?</dt>
            <dd>
              {String(!evRet.evidence_availability_did_not_affect_selection_or_ranking)} — it
              did not
            </dd>
            <dt>top_k per direction tier</dt>
            <dd>{evRet.top_k_per_direction_tier}</dd>
            <dt>Direction tiers</dt>
            <dd>{evRet.direction_tiers.join(', ')}</dd>
            <dt>Snapshot file</dt>
            <dd className="mono small">{evRet.snapshot_file}</dd>
            <dt>Snapshot sha256</dt>
            <dd className="mono small">{evRet.snapshot_sha256}</dd>
            <dt>Manifest sha256</dt>
            <dd className="mono small">{evRet.manifest_sha256}</dd>
          </dl>
          <p className="small muted">{raw.drug_gene_interaction_evidence.pmid_scope_note}</p>
        </Disclosure>

        <Disclosure summary="Input artifact hashes">
          <dl className="kv kv--tight">
            {Object.entries(summary.inputArtifactSha256).map(([k, v]) => (
              <div key={k} style={{ display: 'contents' }}>
                <dt className="mono small">{k}</dt>
                <dd className="mono small">{v}</dd>
              </div>
            ))}
          </dl>
        </Disclosure>

        <Disclosure summary="Osteosarcoma descriptive aggregate (n = 5)">
          <Callout tone="caution">
            {osteo.status} Not a replacement for the frozen Phase 1 result. Not used to
            choose a model or alter rankings.
          </Callout>
          <dl className="kv">
            <dt>Cohort</dt>
            <dd>
              {osteo.cohort.n} lines — {osteo.cohort.model_ids.join(', ')}
            </dd>
            <dt>Common finite targets</dt>
            <dd>
              {osteo.common_finite_target_set.n_included.toLocaleString('en-US')} included,{' '}
              {osteo.common_finite_target_set.n_excluded} excluded
            </dd>
            <dt>ridge_pca mean per-target Spearman</dt>
            <dd>{osteo.mean_per_target_spearman.ridge_pca.toFixed(6)}</dd>
            <dt>ridge_head mean per-target Spearman</dt>
            <dd>{osteo.mean_per_target_spearman.ridge_head.toFixed(6)}</dd>
            <dt>delta (head − pca)</dt>
            <dd>{osteo.delta_ridge_head_minus_ridge_pca.toFixed(6)}</dd>
          </dl>
        </Disclosure>

        <Disclosure summary="Evidence limitations">
          <ul className="small">
            <li>{EVIDENCE_NOT_EFFICACY}</li>
            <li>
              Only six redistribution-verified DGIdb sources are included; PMID coverage is
              source-skewed (~19% of records carry a PMID). A record with no PMID is{' '}
              <em>source-only</em>, not an error.
            </li>
            <li>
              The interaction data vintage is Dec-2023, not current-year data, despite the
              2026-06b release tag.
            </li>
          </ul>
        </Disclosure>

        <Disclosure summary="Sample limitations">
          <ul className="small">
            {summary.limitations.map((l, i) => (
              <li key={i}>{l}</li>
            ))}
          </ul>
        </Disclosure>

        <Disclosure summary="Reproducibility &amp; environment">
          <dl className="kv">
            {Object.entries(summary.environment).map(([k, v]) => (
              <div key={k} style={{ display: 'contents' }}>
                <dt>{k}</dt>
                <dd className="mono small">{v}</dd>
              </div>
            ))}
            <dt>case-study source commit</dt>
            <dd className="mono small">{summary.sourceCommit}</dd>
          </dl>
          <p className="small muted">
            The offline <span className="mono">phase2_report.html</span> is the fully static,
            deterministically generated counterpart of this interface. The baseline arm
            reproduces end to end from the public repository; the Geneformer arm's
            embeddings do not (Kaggle GPU artifact).
          </p>
        </Disclosure>

        <Disclosure summary="Non-clinical &amp; non-efficacy disclaimers" defaultOpen>
          <ul className="small">
            <li>{NON_CLINICAL}</li>
            <li>{EVIDENCE_NOT_EFFICACY}</li>
            <li>{STRUCTURE_NOT_EVIDENCE}</li>
            {summary.disclaimers.map((d, i) => (
              <li key={i}>{d}</li>
            ))}
          </ul>
        </Disclosure>
      </section>
    </div>
  );
}
