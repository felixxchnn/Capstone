# Licences and attribution

## DepMap — data

DepMap, Broad (2026). DepMap Public 26Q1. Dataset. https://depmap.org

Licensed [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Redistribution is
permitted with attribution to the DepMap project.

Release: **DepMap Public 26Q1** (April 2026). Downloaded 2026-07-22.
Source files: `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`,
`CRISPRGeneEffect.csv`, `Model.csv`.

Derived files in `data/processed/` — `expression.npz`, `crispr_effect.npz`,
`model_metadata.csv`, `splits.json`, `gene_columns.json`, `gene_id_map.csv`,
`selective_genes.json` — remain under CC BY 4.0.

Funding acknowledgement (per DepMap): this project is partially funded by the DepMap
Consortium, the Robertson Foundation, the Minderoo Foundation and the Pediatric Cancer
Dependencies Accelerator.

## Geneformer — model

Theodoris, C.V. et al. *Transfer learning enables predictions in network biology.*
Nature (2023). https://huggingface.co/ctheodoris/Geneformer — Apache 2.0.

Checkpoint: `Geneformer-V2-104M_CLcancer` (V2 104M after continual learning on ~14M
cancer transcriptomes). Context 4,096 tokens, `special_token=True`, CLS pooling, layer −1.
Used to generate `data/processed/geneformer_embeddings.csv` (1,140 × 768).

## Osteosarcoma dataset

Sid Sijbrandij, https://osteosarc.com — CC0 1.0. Used as a qualitative inference case
study only; never as validation or test data.

## DGIdb — drug–gene interaction evidence (`evidence.py`)

Cannon M, Stevenson J, Stahl K, Basu R, Coffman A, Kiwala S, McMichael JF, Kuzma K,
Morrissey D, Cotto K, Griffith M, Griffith OL, Wagner AH. *DGIdb 5.0: rebuilding the
drug–gene interaction database for precision medicine and drug discovery platforms.*
Nucleic Acids Research (2024). https://dgidb.org

Upstream release **`2026-06b`** of `https://github.com/dgidb/dgidb-data` (DGIdb service
`v.5.0.12`; TSV export tagged "Data version: Dec-2023 / DGIdb version: v.5.0.11"). The
pinned assets `interactions.tsv`, `genes.tsv`, `drugs.tsv` are downloaded, hash-verified,
and used only to build the committed snapshot; **they are never committed to this repo**.

DGIdb aggregates ~21 interaction sources under **their own separate licences**. DGIdb's
software licence does **not** grant redistribution rights over every aggregated source.
`data/external/dgidb/dgidb_2026-06b.interactions.filtered.tsv` therefore contains records
**only** from the six interaction sources whose redistribution terms are explicitly verified
as compatible with committing filtered records into a public repository. Per-source licence
text, URL, DGIdb version, redistribution decision, reason and record counts are recorded in
`data/external/dgidb/dgidb_2026-06b.manifest.json` (`sources`) and in `evidence.py`
(`SOURCE_LICENCES`). **This is not a claim that the full DGIdb dataset is redistributable —
it is not.**

Included interaction sources and their licences (attribution required):

| Source | Version | Licence | Attribution |
|---|---|---|---|
| **CIViC** | 08-June-2026 | [CC0 1.0](https://docs.civicdb.org/en/latest/about/faq.html#how-is-civic-licensed) — public-domain dedication | Griffith Lab / CIViC — civicdb.org |
| **ChEMBL** | 37 | [CC BY-SA 3.0](https://chembl.gitbook.io/chembl-interface-documentation/about#data-licensing) | ChEMBL, EMBL-EBI |
| **Guide to PHARMACOLOGY** | 2026.2 | [CC BY-SA 4.0](https://www.guidetopharmacology.org/about.jsp) | IUPHAR/BPS Guide to PHARMACOLOGY (GtoPdb) |
| **DoCM** | 2024-10-02 | [CC BY 4.0](https://github.com/griffithlab/docm) | Database of Curated Mutations (DoCM), Griffith Lab |
| **NCI** | 14-September-2017 | Public domain — [NCI copyright/reuse](https://www.cancer.gov/policies/copyright-reuse) | U.S. National Cancer Institute |
| **FDA** | 08-June-2026 | Public domain — [FDA website policy](https://www.fda.gov/about-fda/about-website/website-policies) | U.S. Food and Drug Administration |

**ShareAlike.** Records derived from ChEMBL and Guide to PHARMACOLOGY are redistributed
under CC BY-SA 3.0 and CC BY-SA 4.0 respectively. The committed snapshot as a compilation is
made available under **CC BY-SA 4.0** (the most restrictive applicable ShareAlike term) to
satisfy copyleft; each record additionally names its own `source_license` /
`source_license_url`. CC0, CC BY 4.0 and public-domain records are compatible with inclusion
in a CC BY-SA 4.0 compilation.

Excluded interaction sources (records **not** committed): DTC, TTD, PharmGKB, OncoKB, COSMIC,
CGI, CKB-CORE, CancerCommons, MyCancerGenome, MyCancerGenomeClinicalTrial, TALC, TEND,
TdgClinicalTrial, ClearityFoundationClinicalTrial, ClearityFoundationBiomarkers — each is
NonCommercial, custom-restrictive, "unclear", or copyrighted supplementary-table data. See
the manifest for the individual reason.

The layer is **evidence retrieval, not treatment prediction**: no model, and no inference of
efficacy, clinical relevance, approval, indication, interaction direction (beyond DGIdb's
explicit vocabulary), or osteosarcoma relevance. Every returned display record carries: *"A
recorded drug–gene interaction does not establish efficacy for this sample or for
osteosarcoma."*

## This repository's code

MIT — see `LICENSE`.