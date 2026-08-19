"""
make_fixture.py
===============
Generates a small synthetic dataset that mirrors the exact structure of the
real DepMap files, so the whole pipeline can be tested in seconds.

    python make_fixture.py --out data/fixture

Then point the pipeline at it:

    DEPMAP_DATA_DIR=data/fixture DEPMAP_PROCESSED_DIR=data/fixture_out python build_dataset.py

What it reproduces faithfully
-----------------------------
* Expression CSV with an unnamed integer index, five leading metadata columns,
  and `SYMBOL (ENTREZ)` gene columns.
* **Multiple sequencing profiles per cell line**, only one flagged
  `IsDefaultEntryForModel` -- the duplicate-row trap.
* CRISPR CSV with ModelID as an unnamed first column.
* Model.csv with several cell lines sharing a `PatientID` -- the group-leakage
  trap -- and a Bone/osteosarcoma subset.
* Partially overlapping gene sets between expression and CRISPR, plus a few
  malformed column labels, so the gene reconciliation is exercised.
* Planted signal: a subset of CRISPR targets is a genuine linear function of
  expression, so ridge should beat the mean. If it does not, the pipeline is
  broken rather than the biology being hard.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


LINEAGES = [
    "Bone", "Lung", "Skin", "Bowel", "Breast",
    "Myeloid", "Ovary/Fallopian Tube", "CNS/Brain",
]

DISEASE_BY_LINEAGE = {
    "Bone": "Bone Sarcoma",
    "Lung": "Non-Small Cell Lung Cancer",
    "Skin": "Melanoma",
    "Bowel": "Colorectal Adenocarcinoma",
    "Breast": "Invasive Breast Carcinoma",
    "Myeloid": "Acute Myeloid Leukemia",
    "Ovary/Fallopian Tube": "Ovarian Epithelial Tumor",
    "CNS/Brain": "Diffuse Glioma",
}


def make_gene_labels(n_genes: int, start_entrez: int, rng) -> list[str]:
    """Build `SYMBOL (ENTREZ)` labels with plausible-looking symbols."""
    labels = []
    for i in range(n_genes):
        entrez = start_entrez + i
        symbol = f"GENE{i:04d}"
        labels.append(f"{symbol} ({entrez})")
    return labels


def build(
    out_dir: Path,
    n_patients: int = 90,
    n_genes_shared: int = 300,
    n_genes_expr_only: int = 40,
    n_genes_crispr_only: int = 30,
    n_compounds: int = 60,
    seed: int = 7,
) -> None:
    rng = np.random.default_rng(seed)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------- cell lines
    # Most patients contribute one cell line; some contribute two or three.
    # This is what makes patient-grouped splitting necessary.
    records = []
    model_counter = 1
    for p in range(n_patients):
        patient_id = f"PT-{p:05d}"
        lineage = LINEAGES[p % len(LINEAGES)]
        n_lines = 1 if p % 7 else rng.integers(2, 4)
        for k in range(int(n_lines)):
            model_id = f"ACH-{model_counter:06d}"
            model_counter += 1
            if lineage == "Bone":
                subtype = "Osteosarcoma" if k % 2 == 0 else "Ewing Sarcoma"
                disease = "Bone Sarcoma"
            else:
                subtype = DISEASE_BY_LINEAGE[lineage]
                disease = DISEASE_BY_LINEAGE[lineage]
            records.append(
                {
                    "ModelID": model_id,
                    "PatientID": patient_id,
                    "CellLineName": f"CELL-{model_counter:04d}",
                    "StrippedCellLineName": f"CELL{model_counter:04d}",
                    "DepmapModelType": lineage[:4].upper(),
                    "OncotreeLineage": lineage,
                    "OncotreePrimaryDisease": disease,
                    "OncotreeSubtype": subtype,
                    "OncotreeCode": lineage[:4].upper(),
                    "Age": int(rng.integers(5, 80)),
                    "Sex": "Male" if rng.random() < 0.5 else "Female",
                    "PrimaryOrMetastasis": "Primary" if rng.random() < 0.7 else "Metastatic",
                    "SampleCollectionSite": lineage.lower(),
                    "PediatricModelType": bool(rng.random() < 0.2),
                }
            )

    model_df = pd.DataFrame(records).set_index("ModelID")
    model_ids = model_df.index.tolist()
    n_lines_total = len(model_ids)

    # ------------------------------------------------------------- genes
    shared = make_gene_labels(n_genes_shared, 1000, rng)
    expr_only = make_gene_labels(n_genes_expr_only, 90000, rng)
    crispr_only = make_gene_labels(n_genes_crispr_only, 95000, rng)

    expr_genes = shared + expr_only
    crispr_genes = shared + crispr_only

    # -------------------------------------------------------- expression
    # Latent structure so PCA has something to find, plus a lineage effect.
    n_latent = 8
    latent = rng.normal(size=(n_lines_total, n_latent))
    loadings = rng.normal(size=(n_latent, len(expr_genes))) * 0.8

    lineage_codes = pd.Categorical(model_df["OncotreeLineage"]).codes
    lineage_effect = rng.normal(size=(len(LINEAGES), len(expr_genes))) * 0.6
    lineage_component = lineage_effect[lineage_codes]

    expr_values = (
        3.0
        + latent @ loadings
        + lineage_component
        + rng.normal(scale=0.4, size=(n_lines_total, len(expr_genes)))
    )
    expr_values = np.clip(expr_values, 0.0, None)  # log2(TPM+1) is non-negative

    # A few all-zero genes, to exercise the feature filter.
    expr_values[:, -3:] = 0.0

    expression = pd.DataFrame(expr_values, index=model_ids, columns=expr_genes)

    # ------------------------------------------------------------ crispr
    # Plant real signal: selective dependencies are a linear function of the
    # latent factors, so a model that recovers the latent structure can
    # predict them. Everything else is noise or pan-essential.
    n_selective = 80
    crispr_values = rng.normal(scale=0.05, size=(n_lines_total, len(crispr_genes)))

    signal_loadings = rng.normal(size=(n_latent, n_selective)) * 0.45
    selective_block = latent @ signal_loadings
    selective_block -= selective_block.mean(axis=0)
    selective_block -= 0.45  # push into dependency territory for some lines
    crispr_values[:, :n_selective] = (
        selective_block + rng.normal(scale=0.12, size=(n_lines_total, n_selective))
    )

    # Pan-essential block: every line dependent. Must be excluded by the filter.
    pan_start, pan_end = n_selective, n_selective + 25
    crispr_values[:, pan_start:pan_end] = (
        -1.0 + rng.normal(scale=0.08, size=(n_lines_total, pan_end - pan_start))
    )

    crispr = pd.DataFrame(crispr_values, index=model_ids, columns=crispr_genes)

    # Scattered missingness, plus one gene missing almost everywhere.
    nan_mask = rng.random(crispr.shape) < 0.004
    crispr = crispr.mask(nan_mask)
    crispr.iloc[:, -1] = np.nan

    # ------------------------------------ expression file with duplicates
    # The real matrix carries several sequencing runs per cell line; only one
    # is flagged as the default entry.
    rows = []
    for i, model_id in enumerate(model_ids):
        rows.append(
            {
                "SequencingID": f"CDS-{i:06d}a",
                "ModelConditionID": f"MC-{i:06d}",
                "ModelID": model_id,
                "IsDefaultEntryForMC": True,
                "IsDefaultEntryForModel": True,
                "_values": expr_values[i],
            }
        )
        if i % 5 == 0:  # a second, non-default profile for every fifth line
            jitter = expr_values[i] + rng.normal(scale=0.9, size=expr_values.shape[1])
            rows.append(
                {
                    "SequencingID": f"CDS-{i:06d}b",
                    "ModelConditionID": f"MC-{i:06d}-alt",
                    "ModelID": model_id,
                    "IsDefaultEntryForMC": False,
                    "IsDefaultEntryForModel": False,
                    "_values": np.clip(jitter, 0.0, None),
                }
            )

    meta_frame = pd.DataFrame(
        [{k: v for k, v in r.items() if k != "_values"} for r in rows]
    )
    value_frame = pd.DataFrame(
        np.vstack([r["_values"] for r in rows]), columns=expr_genes
    )
    expression_file = pd.concat(
        [meta_frame.reset_index(drop=True), value_frame.reset_index(drop=True)],
        axis=1,
    )

    # Two malformed gene labels, to exercise the parser's error path.
    renamed = list(expression_file.columns)
    renamed[-1] = "BADGENE_NO_ENTREZ"
    renamed[-2] = "ANOTHER (notanumber)"
    expression_file.columns = renamed

    # Written with an unnamed integer index, exactly like the real file.
    expression_file.to_csv(out_dir / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")

    # ------------------------------------------------------- crispr file
    crispr.to_csv(out_dir / "CRISPRGeneEffect.csv")

    # -------------------------------------------------------- model file
    model_df.to_csv(out_dir / "Model.csv")

    # -------------------------------------------------------- prism file
    # Only a subset of lines is screened, as in reality.
    prism_lines = [m for i, m in enumerate(model_ids) if i % 2 == 0]
    prism_latent = latent[[model_ids.index(m) for m in prism_lines]]
    prism_loadings = rng.normal(size=(n_latent, n_compounds)) * 0.35
    prism_values = (
        prism_latent @ prism_loadings
        + rng.normal(scale=0.25, size=(len(prism_lines), n_compounds))
    )
    prism = pd.DataFrame(
        prism_values,
        index=prism_lines,
        columns=[f"BRD-K{i:06d}-001-01-9" for i in range(n_compounds)],
    )
    prism = prism.mask(rng.random(prism.shape) < 0.01)
    prism.to_csv(out_dir / "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv")

    # ------------------------------------------------------------ report
    n_os = int(model_df["OncotreeSubtype"].str.contains("Osteosarcoma").sum())
    multi = int((model_df["PatientID"].value_counts() > 1).sum())

    print(f"Fixture written to {out_dir}")
    print(f"  cell lines                   : {n_lines_total}")
    print(f"  patients                     : {n_patients}")
    print(f"  patients with multiple lines : {multi}  (group-leakage trap)")
    print(f"  osteosarcoma lines           : {n_os}")
    print(f"  expression rows in file      : {len(rows)}  "
          f"({len(rows) - n_lines_total} non-default duplicates)")
    print(f"  shared genes                 : {n_genes_shared}")
    print(f"  expression-only genes        : {n_genes_expr_only}")
    print(f"  crispr-only genes            : {n_genes_crispr_only}")
    print(f"  planted selective targets    : {n_selective}")
    print(f"  pan-essential decoys         : {pan_end - pan_start}")
    print(f"  prism                        : {len(prism_lines)} lines x "
          f"{n_compounds} compounds")


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate a synthetic DepMap fixture.")
    parser.add_argument("--out", default="data/fixture", help="output directory")
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()
    build(Path(args.out), seed=args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
