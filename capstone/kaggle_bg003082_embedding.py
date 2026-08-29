"""
kaggle_bg003082_embedding.py
============================
The GPU half of the BG003082 Geneformer pipeline. RUN THIS ON KAGGLE OR COLAB
WITH A GPU. It will not run on the dataset-build machine -- it needs the
`geneformer` package, a GPU, network access, and a ~400 MB model download.

    -----------------------------------------------------------------------
    HONESTY NOTE (same status as run_geneformer_embeddings.py): this script
    has NOT been executed against a live Geneformer install. It is written
    against Geneformer's documented V2 API and the 2026-08-06 Kaggle notebook
    snapshot (capstone/kaggle_notebook_v2_2026-08-06.ipynb). If a call
    signature here does not match your installed version, the authority is
    Geneformer's own example notebook:
      https://huggingface.co/ctheodoris/Geneformer/blob/main/examples/
        extract_and_plot_cell_embeddings.ipynb
    Nothing in this repository may claim a BG003082 embedding exists, or that
    tokenisation was verified, unless THIS script actually ran and produced
    the sidecar files it describes. Until then the repo is in the documented
    baseline-only fallback state (capstone/geneformer-bg003082-feasibility.md).
    -----------------------------------------------------------------------

What it does
------------
0. Records the exact environment: Python version, package versions, the pinned
   Geneformer/HF commit, GPU name. Writes them into the provenance record.
1. Builds the BG003082 Geneformer input frames with the repo-local, fully
   validated builder (`geneformer_sample_input.build_bg003082_input`), then
   assembles a one-row AnnData. Linear TPM is the disclosed pseudo-count basis
   (identical to how the frozen 1,140 training embeddings were made).
2. Tokenises with the V2 settings: model_input_size=4096, special_token=True.
3. VERIFIES the tokenised sample before extracting anything:
     * ModelID survived tokenisation and equals "BG003082";
     * the token sequence starts with <cls> and ends with <eos>;
     * length <= 4096;
     * every token id is in the model's token vocabulary;
     * a deterministic SHA-256 of the token-id list (recorded, so a re-run on a
       clean environment can be shown to reproduce the exact tokenisation);
     * an independent re-implementation of Geneformer's normalisation
       (X / n_counts * 10_000 / gene_median, then rank-descending) reproduces
       the same ordering of the top genes -- the normalisation is not taken on
       faith from prose.
4. Extracts ONE cell embedding (CLS token, layer -1) and checks it is a finite
   1 x 768 vector.
5. Saves it as a SEPARATE Phase 2 sidecar:
       data/processed/geneformer_bg003082_embedding.csv          (1 x 768)
       data/processed/geneformer_bg003082_embedding.provenance.json
   It never opens data/processed/geneformer_embeddings.csv for writing -- an
   explicit guard asserts the output path is not that file.
6. Records SHA-256 of the sidecar CSV and the full provenance block, ready to
   append to capstone/data-integrity-hashes.md, pinned to the commit that adds
   the sidecar.

Setup (run once in a Kaggle/Colab cell, mirroring the 2026-08-06 notebook)
------------------------------------------------------------------------
    !pip install -q anndata scanpy
    !git lfs install
    !GIT_LFS_SKIP_SMUDGE=1 git clone https://huggingface.co/ctheodoris/Geneformer /kaggle/working/Geneformer
    # PIN the revision -- fill GENEFORMER_REVISION below with this commit:
    !cd /kaggle/working/Geneformer && git rev-parse HEAD
    !cd /kaggle/working/Geneformer && pip install -q .
    !pip install -q transformers==4.49.0
    # then repair the LFS .pkl stubs and apply the tokenizer .iloc patch --
    # see _repair_geneformer_dictionaries() and _patch_tokenizer_iloc() below,
    # or cells 7-8 of capstone/kaggle_notebook_v2_2026-08-06.ipynb.
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Pinning / settings
# --------------------------------------------------------------------------

# FILL THIS IN on Kaggle from `git -C Geneformer rev-parse HEAD` before running.
# Leave as None only for a first exploratory pass; a committed provenance record
# MUST have the real 40-char commit here.
GENEFORMER_REVISION: str | None = None

# The model subdirectory inside the ctheodoris/Geneformer repo. Cancer-domain,
# continually-pretrained V2, plain PyTorch (runs on a T4; NOT the NVIDIA
# TransformerEngine build).
MODEL_REPO = "ctheodoris/Geneformer"
MODEL_SUBDIR = "Geneformer-V2-104M_CLcancer"

# V2 (95M/104M series) tokenisation settings. Do NOT change these for this model.
MODEL_INPUT_SIZE = 4096
SPECIAL_TOKEN = True
MODEL_VERSION = "V2"

# Geneformer's internal per-cell normalisation, stated explicitly so the
# verification step can reproduce it rather than trust prose:
#
#     norm[g] = raw_count[g] / n_counts * TARGET_SUM / gene_median[g]
#
# then genes are ranked by descending norm, non-zero entries kept, the list
# truncated to (MODEL_INPUT_SIZE - 2), <cls> prepended and <eos> appended.
TARGET_SUM = 10_000

# Embedding extraction.
EMB_MODE = "cls"
EMB_LAYER = -1
FORWARD_BATCH_SIZE = 16
NPROC = 2
EXPECTED_EMB_DIM = 768

# Output sidecar (NEVER geneformer_embeddings.csv).
SIDECAR_CSV_NAME = "geneformer_bg003082_embedding.csv"
SIDECAR_PROV_NAME = "geneformer_bg003082_embedding.provenance.json"
FROZEN_EMBEDDING_MATRIX = "geneformer_embeddings.csv"

_HASH_CHUNK = 1 << 20


# --------------------------------------------------------------------------
# Environment capture
# --------------------------------------------------------------------------

def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(_HASH_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def capture_environment() -> dict:
    """Record Python + package versions and the pinned Geneformer revision."""
    versions: dict[str, str] = {}
    for name in ("numpy", "pandas", "torch", "transformers", "anndata",
                 "scanpy", "geneformer", "huggingface_hub", "datasets", "scipy"):
        try:
            mod = __import__(name)
            versions[name] = getattr(mod, "__version__", "unknown")
        except Exception as exc:  # noqa: BLE001 -- report, don't fail
            versions[name] = f"not importable: {exc.__class__.__name__}"

    gpu = "cpu / unknown"
    try:
        import torch
        if torch.cuda.is_available():
            gpu = torch.cuda.get_device_name(0)
    except Exception:  # noqa: BLE001
        pass

    resolved_revision = GENEFORMER_REVISION
    if resolved_revision is None:
        # Best effort: read it from a cloned repo in the working dir.
        for candidate in glob.glob("/kaggle/working/**/Geneformer/.git/HEAD",
                                   recursive=True):
            try:
                head = Path(candidate).read_text().strip()
                if head.startswith("ref:"):
                    ref = head.split()[1]
                    resolved_revision = (
                        Path(candidate).parent / ref
                    ).read_text().strip()
                else:
                    resolved_revision = head
            except Exception:  # noqa: BLE001
                pass

    return {
        "python": sys.version,
        "platform": platform.platform(),
        "packages": versions,
        "gpu": gpu,
        "geneformer_revision_pinned": GENEFORMER_REVISION,
        "geneformer_revision_resolved": resolved_revision,
        "model_repo": MODEL_REPO,
        "model_subdir": MODEL_SUBDIR,
        "model_input_size": MODEL_INPUT_SIZE,
        "special_token": SPECIAL_TOKEN,
        "model_version": MODEL_VERSION,
        "target_sum": TARGET_SUM,
        "emb_mode": EMB_MODE,
        "emb_layer": EMB_LAYER,
    }


# --------------------------------------------------------------------------
# Environment repair (LFS stubs + pandas-2 .iloc patch) -- notebook cells 7-8
# --------------------------------------------------------------------------

def repair_geneformer_dictionaries() -> list[str]:
    """
    Replace Git-LFS pointer stubs shipped in the installed `geneformer` package
    with the real .pkl dictionaries from the Hub. Returns the list repaired.
    """
    import geneformer
    from huggingface_hub import hf_hub_download, list_repo_files

    installed = os.path.dirname(geneformer.__file__)

    def is_stub(path: str) -> bool:
        try:
            with open(path, "rb") as fh:
                return fh.read(20).startswith(b"version https")
        except Exception:  # noqa: BLE001
            return True

    stubs = [p for p in glob.glob(f"{installed}/**/*.pkl", recursive=True)
             if is_stub(p)]
    repo_files = list_repo_files(MODEL_REPO)
    repaired: list[str] = []
    for stub in stubs:
        name = os.path.basename(stub)
        match = [f for f in repo_files if f.endswith(name)]
        if not match:
            print(f"  no repo match for {name}")
            continue
        real = hf_hub_download(MODEL_REPO, match[0])
        Path(stub).write_bytes(Path(real).read_bytes())
        repaired.append(os.path.relpath(stub, installed))
    print(f"  repaired {len(repaired)} dictionary file(s)")
    return repaired


def patch_tokenizer_iloc() -> int:
    """
    pandas >= 2: `Series[int_array]` is label-based; Geneformer's tokenizer
    indexes with a positional array. Force `.iloc`. Returns the substitution
    count (0 means already patched or the pattern moved).
    """
    import geneformer
    tok = os.path.join(os.path.dirname(geneformer.__file__), "tokenizer.py")
    src = Path(tok).read_text()
    before = src
    src = src.replace('["ensembl_id_collapsed"][coding_miRNA_loc]',
                      '["ensembl_id_collapsed"].iloc[coding_miRNA_loc]')
    src = src.replace('["ensembl_id"][coding_miRNA_loc]',
                      '["ensembl_id"].iloc[coding_miRNA_loc]')
    if src != before:
        Path(tok).write_text(src)
    n = src.count('.iloc[coding_miRNA_loc]')
    print(f"  tokenizer.py .iloc occurrences now: {n}")
    return n


# --------------------------------------------------------------------------
# Input assembly
# --------------------------------------------------------------------------

def build_anndata(repo_root: Path):
    """
    Build the one-row BG003082 AnnData via the repo-local validated builder.

    Requires `geneformer_sample_input` and `sample_profile` importable (add the
    repo root to sys.path on Kaggle), plus the committed
    data/external/sid_osteosarc/BG003082.gene_tpm.gct.gz and
    data/processed/ensembl_map.csv + gene_columns.json.
    """
    sys.path.insert(0, str(repo_root))
    import anndata as ad
    from scipy import sparse
    import geneformer_sample_input as gsi

    X, var_df, obs_df, provenance = gsi.build_bg003082_input()
    # var index must be a plain RangeIndex with `ensembl_id` kept as a column
    # (notebook cell 6): Geneformer's tokenizer indexes var positionally.
    var_reset = var_df.reset_index(drop=True)
    var_reset["ensembl_id"] = list(var_df["ensembl_id"])
    var_reset.index = [str(i) for i in range(len(var_reset))]

    adata = ad.AnnData(
        X=sparse.csr_matrix(X.to_numpy(dtype=np.float32)),
        obs=obs_df.copy(),
        var=var_reset,
    )
    return adata, X, var_df, obs_df, provenance


# --------------------------------------------------------------------------
# Verification of the tokenised sample
# --------------------------------------------------------------------------

def _load_pickle(pattern: str):
    import geneformer
    import pickle
    installed = os.path.dirname(geneformer.__file__)
    hits = glob.glob(f"{installed}/**/{pattern}", recursive=True)
    if not hits:
        raise FileNotFoundError(f"no {pattern} under {installed}")
    with open(hits[0], "rb") as fh:
        return pickle.load(fh), hits[0]


def verify_tokenised_sample(dataset_path: Path, X: pd.DataFrame,
                            var_df: pd.DataFrame) -> dict:
    """
    Load the tokenised dataset for the single BG003082 row and check every
    structural property before any embedding is extracted.
    """
    from datasets import load_from_disk

    ds = load_from_disk(str(dataset_path))
    if len(ds) != 1:
        raise RuntimeError(f"expected 1 tokenised cell, got {len(ds)}")
    row = ds[0]

    # ModelID survived tokenisation
    model_id = row.get("ModelID")
    if model_id != "BG003082":
        raise RuntimeError(f"ModelID did not survive tokenisation: {model_id!r}")

    input_ids = list(row["input_ids"])
    token_dict, token_dict_path = _load_pickle("*token_dictionary_gc104M.pkl")
    id_to_tok = {v: k for k, v in token_dict.items()}

    cls_id = token_dict.get("<cls>")
    eos_id = token_dict.get("<eos>")
    if cls_id is None or eos_id is None:
        raise RuntimeError("token dictionary has no <cls>/<eos> entries")

    checks = {
        "n_tokenised_cells": len(ds),
        "model_id": model_id,
        "token_length": len(input_ids),
        "length_within_limit": len(input_ids) <= MODEL_INPUT_SIZE,
        "starts_with_cls": input_ids[0] == cls_id,
        "ends_with_eos": input_ids[-1] == eos_id,
        "all_tokens_in_vocab": all(t in id_to_tok for t in input_ids),
        "n_tokens_out_of_vocab": sum(t not in id_to_tok for t in input_ids),
        "token_id_sha256": hashlib.sha256(
            json.dumps(input_ids).encode("utf-8")
        ).hexdigest(),
        "token_dictionary_path": os.path.basename(token_dict_path),
    }

    # Independent replication of Geneformer's normalisation + ranking, to prove
    # the "X / n_counts * TARGET_SUM / gene_median then rank-descending" recipe
    # is what produced the ordering -- not trusted from documentation.
    gene_median, _ = _load_pickle("*gene_median_dictionary_gc104M.pkl")
    ensembl_ids = list(var_df["ensembl_id"])
    raw = X.to_numpy(dtype=np.float64).ravel()
    n_counts = raw.sum()
    in_vocab = np.array([e in token_dict for e in ensembl_ids])
    med = np.array([gene_median.get(e, np.nan) for e in ensembl_ids], dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        norm = raw / n_counts * TARGET_SUM / med
    keep = in_vocab & np.isfinite(norm) & (norm > 0)
    order = np.argsort(-norm[keep], kind="stable")
    expected_top_tokens = [
        token_dict[ensembl_ids[i]] for i in np.flatnonzero(keep)[order]
    ][: MODEL_INPUT_SIZE - 2]
    # the tokeniser body (between <cls> and <eos>)
    body = input_ids[1:-1] if SPECIAL_TOKEN else input_ids
    n_compare = min(len(body), len(expected_top_tokens), 50)
    checks["independent_norm_top50_matches"] = (
        body[:n_compare] == expected_top_tokens[:n_compare]
    )
    checks["independent_norm_n_compared"] = n_compare

    hard = ["length_within_limit", "starts_with_cls", "ends_with_eos",
            "all_tokens_in_vocab", "independent_norm_top50_matches"]
    failed = [k for k in hard if not checks[k]]
    if model_id != "BG003082":
        failed.append("model_id")
    checks["all_hard_checks_passed"] = not failed
    if failed:
        raise RuntimeError(f"tokenisation verification failed: {failed}\n{checks}")
    return checks


# --------------------------------------------------------------------------
# Model + extraction
# --------------------------------------------------------------------------

def download_model(work_dir: Path) -> Path:
    from huggingface_hub import snapshot_download

    kwargs = dict(
        repo_id=MODEL_REPO,
        allow_patterns=[f"{MODEL_SUBDIR}/*"],
        local_dir=str(work_dir / "geneformer_repo"),
    )
    if GENEFORMER_REVISION:
        kwargs["revision"] = GENEFORMER_REVISION
    repo_dir = snapshot_download(**kwargs)
    model_dir = Path(repo_dir) / MODEL_SUBDIR
    if not model_dir.is_dir():
        raise FileNotFoundError(f"model dir not found at {model_dir}")
    return model_dir


def tokenize(adata, work_dir: Path) -> Path:
    from geneformer import TranscriptomeTokenizer

    in_dir = work_dir / "tok_input"
    in_dir.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(in_dir / "bg003082.h5ad")

    out_dir = work_dir / "tokenized"
    out_dir.mkdir(parents=True, exist_ok=True)
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={"ModelID": "ModelID"},
        nproc=NPROC,
        model_input_size=MODEL_INPUT_SIZE,
        special_token=SPECIAL_TOKEN,
    )
    tk.tokenize_data(str(in_dir), str(out_dir), "bg003082", file_format="h5ad")
    return out_dir / "bg003082.dataset"


def extract_embedding(model_dir: Path, dataset_path: Path,
                      work_dir: Path) -> np.ndarray:
    from geneformer import EmbExtractor

    out_dir = work_dir / "emb"
    out_dir.mkdir(parents=True, exist_ok=True)
    embex = EmbExtractor(
        model_type="Pretrained",
        num_classes=0,
        emb_mode=EMB_MODE,
        max_ncells=None,
        emb_layer=EMB_LAYER,
        emb_label=["ModelID"],
        forward_batch_size=FORWARD_BATCH_SIZE,
        nproc=NPROC,
        model_version=MODEL_VERSION,
    )
    embs = embex.extract_embs(
        model_directory=str(model_dir),
        input_data_file=str(dataset_path),
        output_directory=str(out_dir),
        output_prefix="bg003082_emb",
    )
    if "ModelID" not in embs.columns:
        raise RuntimeError("extracted embeddings have no ModelID column")
    embs = embs.set_index("ModelID")
    numeric = embs.select_dtypes(include=[np.number])
    if numeric.shape != (1, EXPECTED_EMB_DIM):
        raise RuntimeError(
            f"expected a 1 x {EXPECTED_EMB_DIM} embedding, got {numeric.shape}"
        )
    vec = numeric.to_numpy(dtype=np.float64).ravel()
    if not np.all(np.isfinite(vec)):
        raise RuntimeError("embedding contains non-finite values")
    return vec


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Kaggle/Colab GPU step: tokenise + extract the BG003082 "
                    "Geneformer CLS embedding into a Phase 2 sidecar."
    )
    parser.add_argument("--repo-root", default=".",
                        help="Path to the Capstone repo checkout.")
    parser.add_argument("--work-dir", default="geneformer_bg003082_work")
    parser.add_argument("--skip-repair", action="store_true",
                        help="Skip the LFS-stub / tokenizer .iloc repair "
                             "(only if you have already run it this session).")
    args = parser.parse_args()

    repo_root = Path(args.repo_root).resolve()
    processed_dir = repo_root / "data" / "processed"
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 74)
    print("BG003082 GENEFORMER EMBEDDING  (Kaggle/Colab GPU step)")
    print("=" * 74)

    env = capture_environment()
    print(f"  python   : {env['python'].splitlines()[0]}")
    print(f"  gpu      : {env['gpu']}")
    print(f"  geneformer revision (pinned/resolved): "
          f"{env['geneformer_revision_pinned']} / {env['geneformer_revision_resolved']}")
    if env["geneformer_revision_resolved"] is None:
        print("  WARNING: no Geneformer revision recorded. Set GENEFORMER_REVISION "
              "before producing a committed provenance record.")

    if not args.skip_repair:
        print("\n[1/6] Repairing Geneformer LFS dictionaries + tokenizer patch")
        env["dictionaries_repaired"] = repair_geneformer_dictionaries()
        env["tokenizer_iloc_occurrences"] = patch_tokenizer_iloc()

    print("\n[2/6] Building + validating the BG003082 input frames")
    adata, X, var_df, obs_df, build_prov = build_anndata(repo_root)
    print(f"      X {tuple(X.shape)}  n_counts {obs_df['n_counts'].iloc[0]:.4f}")
    print(f"      mapped {build_prov['frame_report']['mapped_values']} / "
          f"unresolved {build_prov['frame_report']['unresolved_values']}")

    print("\n[3/6] Model download")
    model_dir = download_model(work_dir)
    print(f"      {model_dir}")

    print("\n[4/6] Tokenise")
    dataset_path = tokenize(adata, work_dir)
    print(f"      {dataset_path}")

    print("\n[5/6] Verify tokenisation (before extracting anything)")
    tok_checks = verify_tokenised_sample(dataset_path, X, var_df)
    for k, v in tok_checks.items():
        print(f"      {k:34s}: {v}")

    print("\n[6/6] Extract the CLS embedding")
    vec = extract_embedding(model_dir, dataset_path, work_dir)
    print(f"      1 x {vec.size} finite embedding extracted")

    # ---- write the sidecar (NEVER the frozen matrix) --------------------
    out_csv = processed_dir / SIDECAR_CSV_NAME
    assert out_csv.name != FROZEN_EMBEDDING_MATRIX, "refusing to write the frozen matrix"
    assert not str(out_csv).endswith(FROZEN_EMBEDDING_MATRIX)
    frozen_before = None
    frozen_path = processed_dir / FROZEN_EMBEDDING_MATRIX
    if frozen_path.is_file():
        frozen_before = _sha256_file(frozen_path)

    emb_df = pd.DataFrame(
        vec.reshape(1, -1),
        index=pd.Index(["BG003082"], name="ModelID"),
        columns=[str(i) for i in range(vec.size)],
    )
    emb_df.to_csv(out_csv)

    if frozen_before is not None:
        assert _sha256_file(frozen_path) == frozen_before, \
            "the frozen embedding matrix changed -- aborting"

    provenance = {
        "sample_id": "BG003082",
        "produced_by": "capstone/kaggle_bg003082_embedding.py",
        "environment": env,
        "input_build": build_prov,
        "tokenisation_checks": tok_checks,
        "embedding": {
            "shape": [1, int(vec.size)],
            "all_finite": bool(np.all(np.isfinite(vec))),
            "l2_norm": float(np.linalg.norm(vec)),
            "sidecar_csv": SIDECAR_CSV_NAME,
            "sidecar_csv_sha256": _sha256_file(out_csv),
        },
        "frozen_matrix_untouched": {
            "file": FROZEN_EMBEDDING_MATRIX,
            "sha256_before": frozen_before,
            "sha256_after": _sha256_file(frozen_path) if frozen_path.is_file() else None,
        },
        "disclosures": [
            "Bulk primary-tumour TPM, not raw-count scRNA-seq: schema-compatible "
            "but scientifically out of distribution for a model trained on "
            "cultured DepMap cell lines.",
            "Pseudo-count basis is linear TPM, identical to the frozen 1,140 "
            "training embeddings (log2(TPM+1) -> 2**x - 1 round-trips to TPM). "
            "Real RSEM counts were available and deliberately not used.",
            "Entrez->Ensembl map is NCBI gene2ensembl; the training embeddings "
            "used a mygene map that no longer exists in the repo. Provenance "
            "differs even though the join key (Entrez) is the same.",
            "This embedding is NOT reproducible from the public repository "
            "(needs Geneformer + a GPU). Same reproducibility asymmetry as the "
            "Phase 1 Geneformer arm, plus the map-provenance difference above.",
        ],
    }
    prov_path = processed_dir / SIDECAR_PROV_NAME
    prov_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")

    print(f"\n  wrote {out_csv}")
    print(f"  wrote {prov_path}")
    print(f"  sidecar sha256: {provenance['embedding']['sidecar_csv_sha256']}")
    print("\n  Append these to capstone/data-integrity-hashes.md, pinned to the "
          "commit that adds the sidecar:")
    print(f"    {provenance['embedding']['sidecar_csv_sha256']}  "
          f"data/processed/{SIDECAR_CSV_NAME}")
    print("\nDONE. The frozen geneformer_embeddings.csv was not modified.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
