import { SegmentedControl } from '@/components/common/SegmentedControl';
import { SampleRoleBadge } from '@/components/samples/SampleRoleBadge';
import type { ModelId, SampleId } from '@/types/caseStudy';

/** Combined sample + model selector. The two model rankings are ALWAYS kept
 *  independent — there is deliberately no "consensus" option. */
export function SampleModelControls({
  sample,
  model,
  onSample,
  onModel,
}: {
  sample: SampleId;
  model: ModelId;
  onSample: (s: SampleId) => void;
  onModel: (m: ModelId) => void;
}) {
  return (
    <div className="sample-model-controls no-print">
      <div className="sample-model-controls__row">
        <SegmentedControl<SampleId>
          legend="Sample"
          name="sample"
          value={sample}
          onChange={onSample}
          options={[
            { value: 'ACH-000364', label: 'ACH-000364', hint: 'validation anchor' },
            { value: 'BG003082', label: 'BG003082', hint: 'exploratory external' },
          ]}
        />
        <SegmentedControl<ModelId>
          legend="Model"
          name="model"
          value={model}
          onChange={onModel}
          describedById="model-independence-note"
          options={[
            { value: 'ridge_pca', label: 'ridge_pca', hint: 'PCA of expression' },
            { value: 'ridge_head', label: 'ridge_head', hint: 'Geneformer embedding' },
          ]}
        />
      </div>
      <div className="sample-model-controls__aside">
        <SampleRoleBadge sample={sample} />
        <p id="model-independence-note" className="tiny muted">
          <code>ridge_pca</code> and <code>ridge_head</code> are shown separately. Their
          rankings are never merged into a consensus list.
        </p>
      </div>
    </div>
  );
}
