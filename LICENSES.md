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

## This repository's code

MIT — see `LICENSE`.