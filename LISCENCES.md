# Licences and attribution

## DepMap (data)

DepMap, Broad (2026). *DepMap Public [RELEASE]*. Figshare+. Dataset.
https://doi.org/[DOI]

Licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Redistribution is permitted with attribution to the DepMap project.

Files in `data/processed/` are derived from this release:
`expression.npz`, `crispr_effect.npz`, `model_metadata.csv`, `splits.json`,
`gene_columns.json`, `gene_id_map.csv`, `selective_genes.json`.

## Geneformer (model)

Theodoris, C.V. et al. *Transfer learning enables predictions in network biology.*
Nature (2023). Checkpoint used: `Geneformer-V2-104M_CLcancer`, 4,096-token context,
CLS pooling from the final layer.
https://huggingface.co/ctheodoris/Geneformer — Apache 2.0.

`data/processed/geneformer_embeddings.csv` (1,140 × 768) was generated with this model.

## Osteosarcoma dataset

Sid Sijbrandij, https://osteosarc.com — CC0 1.0. Used as a qualitative inference
case study only; never as validation or test data.

## This repository's code

See `LICENSE`. [MIT or Apache-2.0 — pick one]
