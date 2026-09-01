import { NON_CLINICAL } from '@/lib/format';

export function NonClinicalBanner() {
  return (
    <p className="nonclinical-banner" role="note">
      <span className="nonclinical-banner__dot" aria-hidden="true" />
      {NON_CLINICAL} Predicted gene dependencies are not treatment-response
      predictions, drug recommendations, or diagnoses.
    </p>
  );
}
