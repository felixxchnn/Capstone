import { Callout } from '@/components/common/primitives';

// Placeholder — the live protein-structure explorer (UniProt -> RCSB -> AlphaFold
// -> Mol*) lands in the next commit.
export function ProteinStructurePage() {
  return (
    <div className="container page">
      <h1>Protein Structure</h1>
      <Callout tone="info" title="Coming in the next commit">
        The interactive structure explorer for the protein encoded by a selected gene.
      </Callout>
    </div>
  );
}
