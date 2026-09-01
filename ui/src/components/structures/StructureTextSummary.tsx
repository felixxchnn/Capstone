import type { StructureResolution, StructureSelection } from '@/types/structure';
import { STRUCTURE_NOT_EVIDENCE } from '@/lib/format';

/** Non-3D text description of the selected structure — the accessible
 *  alternative to the interactive viewer. Always rendered. */
export function StructureTextSummary({
  resolution,
  selection,
}: {
  resolution: StructureResolution;
  selection: StructureSelection | null;
}) {
  const uni = resolution.uniprot;
  return (
    <div className="structure-text-summary">
      <h3>Text summary</h3>
      <dl className="kv">
        <dt>Gene</dt>
        <dd>
          {resolution.symbol} (Entrez {resolution.entrezId}), human — taxonomy{' '}
          {resolution.humanTaxonomyId}
        </dd>
        <dt>Encoded protein</dt>
        <dd>
          {uni
            ? `${uni.proteinName} — UniProt ${uni.accession} (${uni.entryId}), ${uni.length} residues, reviewed`
            : 'No reviewed human UniProt entry was found for this Entrez ID.'}
        </dd>
        <dt>Experimental structures</dt>
        <dd>
          {resolution.experimental.length === 0
            ? 'None found in the PDB for this protein.'
            : `${resolution.experimental.length} candidate experimental PDB ${
                resolution.experimental.length === 1 ? 'entry' : 'entries'
              } (listed by resolution). ` +
              resolution.experimental
                .slice(0, 5)
                .map(
                  (c) =>
                    `${c.pdbId} (${c.method ?? 'method n/a'}${
                      c.resolutionAngstrom !== null ? `, ${c.resolutionAngstrom.toFixed(2)} Å` : ''
                    })`,
                )
                .join('; ')}
        </dd>
        <dt>Predicted model</dt>
        <dd>
          {resolution.predicted
            ? `AlphaFold DB ${resolution.predicted.alphafoldId} (v${resolution.predicted.version}), mean pLDDT ${
                resolution.predicted.meanPlddt !== null
                  ? resolution.predicted.meanPlddt.toFixed(1)
                  : 'n/a'
              } — a predicted structure, not an experimental measurement.`
            : 'No AlphaFold DB model available.'}
        </dd>
        <dt>Currently shown</dt>
        <dd>
          {selection === null
            ? 'Nothing — no structure is available for this protein.'
            : selection.kind === 'experimental'
              ? `Experimental structure ${selection.candidate.pdbId}${
                  selection.candidate.method ? ` (${selection.candidate.method})` : ''
                }${
                  selection.candidate.resolutionAngstrom !== null
                    ? `, ${selection.candidate.resolutionAngstrom.toFixed(2)} Å resolution`
                    : ''
                }.`
              : `AlphaFold DB predicted model ${selection.model.alphafoldId} — predicted, not measured.`}
        </dd>
      </dl>
      <p className="small muted">{STRUCTURE_NOT_EVIDENCE}</p>
    </div>
  );
}
