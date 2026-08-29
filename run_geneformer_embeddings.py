"""
run_geneformer_embeddings.py
============================
Tokenise the prepared AnnData and extract Geneformer cell embeddings.

    RUN THIS ON KAGGLE OR COLAB WITH A GPU. It will not run on the machine used
    to build the dataset -- it needs the `geneformer` package, a GPU, and network
    access, and it downloads a ~400 MB model.

    -----------------------------------------------------------------------
    HONESTY NOTE: unlike the rest of this pipeline, this script could not be
    executed and tested in the environment where it was written (no GPU, no
    geneformer package, no network). It is written against Geneformer's current
    documented API, but that API has changed between versions. If a call
    signature here does not match your installed version, the authority is
    Geneformer's own example notebook:
      https://huggingface.co/ctheodoris/Geneformer/blob/main/examples/
        extract_and_plot_cell_embeddings.ipynb
    Treat any mismatch as this script being out of date, not the notebook.
    -----------------------------------------------------------------------

What it does
------------
1. Downloads the 104M CLcancer model (the cancer-domain V2 model) in plain
   PyTorch form from the `ctheodoris/Geneformer` repo -- NOT the NVIDIA
   TransformerEngine build, which needs an A100/H100. This one runs on a T4.
2. Tokenises `geneformer_input.h5ad` with the V2 settings
   (model_input_size=4096, special_token=True).
3. Extracts one cell embedding per cell line (the CLS-token embedding) with
   EmbExtractor, keyed on ModelID.
4. Saves an embeddings matrix, indexed by ModelID, aligned to the same cell
   lines as crispr_effect and splits.json, ready for train_head.py.

Setup (run once in a Kaggle/Colab cell)
---------------------------------------
    !pip install anndata scanpy mygene
    !pip install git+https://huggingface.co/ctheodoris/Geneformer.git
    # enable a GPU accelerator in the notebook settings

Then, if you have not already produced geneformer_input.h5ad on this machine,
run prepare_geneformer_input.py here first (mygene needs the network that
Kaggle provides).
"""

from __future__ import annotations

import os
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------
# Configuration -- adjust paths for your Kaggle/Colab layout
# --------------------------------------------------------------------------

# The 104M CLcancer model lives in a subdirectory of the ctheodoris/Geneformer
# repo. This is the cancer-continually-pretrained V2 model, in plain PyTorch.
MODEL_REPO = "ctheodoris/Geneformer"
MODEL_SUBDIR = "Geneformer-V2-104M_CLcancer"

# V2 (95M/104M series) tokenisation settings. Do not change these for this
# model: the 30M series would use 2048 / special_token=False instead.
MODEL_INPUT_SIZE = 4096
SPECIAL_TOKEN = True
MODEL_VERSION = "V2"

# Emb extraction. "cls" takes the CLS-token embedding as the cell embedding,
# which is the natural representation for the special-token V2 models.
EMB_MODE = "cls"
EMB_LAYER = -1
FORWARD_BATCH_SIZE = 16   # lower if the T4 runs out of memory
NPROC = 2


def download_model(work_dir: Path) -> Path:
    """Download the CLcancer model subdirectory from the Hub."""
    from huggingface_hub import snapshot_download

    print(f"Downloading {MODEL_REPO}/{MODEL_SUBDIR} ...")
    repo_dir = snapshot_download(
        repo_id=MODEL_REPO,
        allow_patterns=[f"{MODEL_SUBDIR}/*"],
        local_dir=str(work_dir / "geneformer_repo"),
    )
    model_dir = Path(repo_dir) / MODEL_SUBDIR
    if not model_dir.is_dir():
        raise FileNotFoundError(
            f"Expected model directory not found at {model_dir}. "
            f"Check the current subdirectory name in the {MODEL_REPO} repo; "
            f"model directory names occasionally change between releases."
        )
    print(f"  model at: {model_dir}")
    return model_dir


def tokenize(h5ad_path: Path, token_out_dir: Path) -> Path:
    """
    Tokenise the AnnData. Carries ModelID through as a custom attribute so the
    embeddings can be mapped back to cell lines.
    """
    from geneformer import TranscriptomeTokenizer

    token_out_dir.mkdir(parents=True, exist_ok=True)

    # The tokenizer reads a directory of .h5ad/.loom files, so place the input
    # in its own directory.
    input_dir = h5ad_path.parent / "geneformer_tokenize_input"
    input_dir.mkdir(parents=True, exist_ok=True)
    target = input_dir / h5ad_path.name
    if not target.exists():
        target.write_bytes(h5ad_path.read_bytes())

    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={"ModelID": "ModelID"},
        nproc=NPROC,
        model_input_size=MODEL_INPUT_SIZE,
        special_token=SPECIAL_TOKEN,
    )
    tk.tokenize_data(
        str(input_dir),
        str(token_out_dir),
        "depmap",
        file_format="h5ad",
    )
    dataset_path = token_out_dir / "depmap.dataset"
    print(f"  tokenised dataset: {dataset_path}")
    return dataset_path


def extract_embeddings(
    model_dir: Path,
    dataset_path: Path,
    emb_out_dir: Path,
) -> pd.DataFrame:
    """
    Extract one cell embedding per cell line, keyed on ModelID.

    Returns a DataFrame indexed by ModelID; columns are the embedding
    dimensions (768 for the Geneformer V2 104M model; the committed
    geneformer_embeddings.csv and the 2026-08-06 Kaggle run both have 768).
    """
    from geneformer import EmbExtractor

    emb_out_dir.mkdir(parents=True, exist_ok=True)

    embex = EmbExtractor(
        model_type="Pretrained",   # pretrained model, not a classifier
        num_classes=0,
        emb_mode=EMB_MODE,
        max_ncells=None,           # keep every cell line
        emb_layer=EMB_LAYER,
        emb_label=["ModelID"],     # attach ModelID to each embedding row
        forward_batch_size=FORWARD_BATCH_SIZE,
        nproc=NPROC,
        model_version=MODEL_VERSION,
    )
    embs = embex.extract_embs(
        model_directory=str(model_dir),
        input_data_file=str(dataset_path),
        output_directory=str(emb_out_dir),
        output_prefix="depmap_emb",
    )

    # embs is a DataFrame with a ModelID column plus the embedding dimensions.
    if "ModelID" not in embs.columns:
        raise RuntimeError(
            "Extracted embeddings have no ModelID column. Check that the "
            "custom attribute survived tokenisation (custom_attr_name_dict) "
            "and that emb_label=['ModelID'] was set."
        )
    embs = embs.set_index("ModelID")

    # Keep only numeric embedding columns.
    numeric = embs.select_dtypes(include=[np.number])
    numeric.index.name = "ModelID"
    print(f"  embeddings: {numeric.shape[0]} cell lines x "
          f"{numeric.shape[1]} dimensions")
    return numeric


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Extract Geneformer cell embeddings for DepMap cell lines."
    )
    parser.add_argument(
        "--processed-dir",
        default=os.environ.get("DEPMAP_PROCESSED_DIR", "data/processed"),
        help="Directory containing geneformer_input.h5ad and where the "
             "embeddings will be written.",
    )
    parser.add_argument(
        "--work-dir",
        default="geneformer_work",
        help="Scratch directory for the model download and tokenised data.",
    )
    args = parser.parse_args()

    processed_dir = Path(args.processed_dir)
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    h5ad_path = processed_dir / "geneformer_input.h5ad"
    if not h5ad_path.is_file():
        raise SystemExit(
            f"\n{h5ad_path} not found.\n"
            "Run prepare_geneformer_input.py first (on a machine with anndata; "
            "Kaggle works well because it also has the network mygene needs)."
        )

    print("=" * 74)
    print("GENEFORMER EMBEDDING EXTRACTION")
    print("=" * 74)

    print("\n[1/3] Model")
    model_dir = download_model(work_dir)

    print("\n[2/3] Tokenise")
    dataset_path = tokenize(h5ad_path, work_dir / "tokenized")

    print("\n[3/3] Extract embeddings")
    embeddings = extract_embeddings(
        model_dir, dataset_path, work_dir / "embeddings"
    )

    # Save aligned to the rest of the pipeline. train_head.py will load this
    # and pair it with crispr_effect + splits.json, reusing the baseline
    # metric machinery.
    out_csv = processed_dir / "geneformer_embeddings.csv"
    embeddings.to_csv(out_csv)
    try:
        embeddings.to_parquet(processed_dir / "geneformer_embeddings.parquet")
    except Exception:
        pass  # parquet optional; csv is the portable copy

    print(f"\n{'=' * 74}")
    print("DONE")
    print("=" * 74)
    print(f"  wrote {out_csv.name}  "
          f"({embeddings.shape[0]} cell lines x {embeddings.shape[1]} dims)")
    print("\n  Next: train_head.py -- train a head on these embeddings and")
    print("  compare, on the SAME splits, against the ridge baseline. If the")
    print("  embeddings do not beat ridge, that is the finding; report it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
