import type { SampleId } from '@/types/caseStudy';

export function SampleRoleBadge({ sample }: { sample: SampleId }) {
  if (sample === 'ACH-000364') {
    return (
      <span className="badge badge--anchor badge--dot" title="Held-out DepMap validation cell line">
        Validation anchor
      </span>
    );
  }
  return (
    <span className="badge badge--external badge--dot" title="External bulk tumour; no measured outcome">
      Exploratory external
    </span>
  );
}

export function SampleRoleLine({ sample }: { sample: SampleId }) {
  if (sample === 'ACH-000364') {
    return (
      <p className="small muted">
        <strong>ACH-000364 (U-2 OS)</strong> — a held-out DepMap validation-split cell
        line. Its observed CRISPR outcomes are real and were attached to the ranked genes{' '}
        <em>after</em> prediction and ranking, only to check the pipeline.
      </p>
    );
  }
  return (
    <p className="small muted">
      <strong>BG003082</strong> — a real primary osteosarcoma tumour (bulk RNA-seq),
      absent from every DepMap split. <strong>No CRISPR outcome exists</strong> for it,
      and bulk tumour tissue is a domain shift from the cultured cell lines the models
      were trained and validated on. Every number for it is exploratory.
    </p>
  );
}
