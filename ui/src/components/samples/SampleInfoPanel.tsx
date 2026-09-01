import type { SampleMetadata } from '@/data/CapstoneDataSource';
import { integer } from '@/lib/format';
import { Callout } from '@/components/common/primitives';

export function SampleInfoPanel({ meta }: { meta: SampleMetadata }) {
  if (meta.kind === 'anchor') {
    const d = meta.data;
    return (
      <div className="sample-info">
        <dl className="kv">
          <dt>Sample</dt>
          <dd>
            {meta.id} — {d.cell_line}
          </dd>
          <dt>Sample type</dt>
          <dd>Cultured human cancer cell line (DepMap)</dd>
          <dt>Role</dt>
          <dd>{d.role.replace(/_/g, ' ')} — pipeline-verification anchor</dd>
          <dt>Split status</dt>
          <dd>
            DepMap <strong>{d.depmap_split}</strong> split; held out. {d.split_assertion}
          </dd>
          <dt>Tissue / context</dt>
          <dd>Osteosarcoma; grown in vitro</dd>
          <dt>Expression mapping</dt>
          <dd>
            {integer(d.baseline_input.n_features)} features, {d.baseline_input.missing_features}{' '}
            missing, {d.baseline_input.imputed_features} imputed (direct{' '}
            <span className="mono">expression.npz</span> row)
          </dd>
          <dt>Geneformer input</dt>
          <dd>
            {integer(d.head_input.n_features)}-dim row from the frozen 1,140-line
            embedding matrix
          </dd>
          <dt>Domain shift</dt>
          <dd>None beyond the usual cell-line caveats — same data family as training</dd>
          <dt>Outcome availability</dt>
          <dd>
            Observed CRISPR GeneEffect for {integer(d.observed_crispr.n_targets_with_value)}{' '}
            targets. {d.observed_crispr.role}
          </dd>
        </dl>
      </div>
    );
  }

  const d = meta.data;
  const r = d.baseline_input.reconciliation;
  return (
    <div className="sample-info">
      <dl className="kv">
        <dt>Sample</dt>
        <dd>{meta.id}</dd>
        <dt>Sample type</dt>
        <dd>Real primary tumour — bulk RNA-seq (not a cultured cell line)</dd>
        <dt>Role</dt>
        <dd>{d.analysis_role.replace(/_/g, ' ')}</dd>
        <dt>Split status</dt>
        <dd>Absent from every DepMap split (train / val / test)</dd>
        <dt>Tissue / context</dt>
        <dd>
          Osteosarcoma primary tumour, resected tissue; Sid Sijbrandij self-released
          dataset (CC0)
        </dd>
        <dt>Expression mapping coverage</dt>
        <dd>
          {integer(r.canonical_genes_mapped)} of {integer(r.canonical_genes)} canonical
          genes resolved via Ensembl-ID join ({integer(r.canonical_genes_measured_nonzero)}{' '}
          measured &gt; 0, {integer(r.canonical_genes_measured_zero)} measured zero)
        </dd>
        <dt>Missing genes</dt>
        <dd>
          {r.canonical_genes_missing} left as explicit NaN — never zero-filled. No symbol
          fallback ({r.symbol_fallback}).
        </dd>
        <dt>Imputation behaviour</dt>
        <dd>{d.baseline_input.imputation}</dd>
        <dt>Domain shift</dt>
        <dd>
          Bulk primary-tumour tissue scored with models trained and validated only on
          cultured DepMap cell lines. Real, and never phrased as measured performance.
        </dd>
        <dt>Outcome availability</dt>
        <dd>
          <strong>Unavailable.</strong> No CRISPR screen exists; no observed value is
          loaded, invented, or computed.
        </dd>
      </dl>
      <Callout tone="caution" title="Geneformer sidecar embedding">
        The Geneformer embedding for BG003082 was generated separately (Kaggle GPU run,
        revision pinned <span className="mono">{d.head_input.geneformer_revision_pinned}</span>).
        Commensurability with the historical Phase 1 embeddings is not proven — bulk-tumour
        input, an NCBI rather than the vanished <span className="mono">mygene</span> map, and
        a fresh revision pin.
      </Callout>
    </div>
  );
}
