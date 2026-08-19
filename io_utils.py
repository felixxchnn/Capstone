"""
io_utils.py
===========
Matrix persistence with a graceful fallback.

Parquet is the preferred on-disk format (compact, preserves dtypes, fast), but
it requires `pyarrow`, which is not installed everywhere. Rather than making
the pipeline fail on a machine without it, these helpers write Parquet when
possible and fall back to a compressed NumPy `.npz` bundle plus a JSON sidecar
holding the index and column labels.

`save_matrix` / `load_matrix` are symmetric: whatever `save_matrix` wrote,
`load_matrix` will read, without the caller needing to know which format
was used.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

try:  # pragma: no cover - depends on the environment
    import pyarrow  # noqa: F401
    HAVE_PARQUET = True
except Exception:  # pragma: no cover
    HAVE_PARQUET = False


def _npz_paths(stem: Path) -> tuple[Path, Path]:
    return stem.with_suffix(".npz"), stem.with_suffix(".labels.json")


def save_matrix(df: pd.DataFrame, stem: Path) -> Path:
    """
    Persist a numeric DataFrame.

    Parameters
    ----------
    df
        DataFrame with a labelled index and labelled columns. Values should be
        numeric; object columns are not supported by the npz fallback.
    stem
        Path *without* extension, e.g. Path("data/processed/expression").

    Returns
    -------
    Path to the file actually written.
    """
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)

    if HAVE_PARQUET:
        target = stem.with_suffix(".parquet")
        df.to_parquet(target)
        return target

    npz_path, labels_path = _npz_paths(stem)
    np.savez_compressed(npz_path, values=df.to_numpy())
    labels = {
        "index": [str(i) for i in df.index],
        "index_name": df.index.name,
        "columns": [str(c) for c in df.columns],
        "dtype": str(df.to_numpy().dtype),
    }
    labels_path.write_text(json.dumps(labels), encoding="utf-8")
    return npz_path


def load_matrix(stem: Path) -> pd.DataFrame:
    """
    Load a matrix previously written by `save_matrix`.

    Tries Parquet first, then the npz fallback.
    """
    stem = Path(stem)

    parquet_path = stem.with_suffix(".parquet")
    if parquet_path.is_file():
        return pd.read_parquet(parquet_path)

    npz_path, labels_path = _npz_paths(stem)
    if npz_path.is_file() and labels_path.is_file():
        values = np.load(npz_path)["values"]
        labels = json.loads(labels_path.read_text(encoding="utf-8"))
        df = pd.DataFrame(
            values,
            index=pd.Index(labels["index"], name=labels["index_name"]),
            columns=labels["columns"],
        )
        return df

    raise FileNotFoundError(
        f"No saved matrix found for stem {stem}. "
        f"Looked for {parquet_path.name} and {npz_path.name}. "
        f"Did you run build_dataset.py first?"
    )


def save_table(df: pd.DataFrame, stem: Path) -> Path:
    """
    Persist a mixed-dtype table (e.g. cell line metadata).

    Uses Parquet when available, CSV otherwise. CSV is acceptable here because
    metadata tables are small.
    """
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)

    if HAVE_PARQUET:
        target = stem.with_suffix(".parquet")
        df.to_parquet(target)
        return target

    target = stem.with_suffix(".csv")
    df.to_csv(target)
    return target


def load_table(stem: Path, index_col: str | int | None = 0) -> pd.DataFrame:
    """Load a table previously written by `save_table`."""
    stem = Path(stem)

    parquet_path = stem.with_suffix(".parquet")
    if parquet_path.is_file():
        return pd.read_parquet(parquet_path)

    csv_path = stem.with_suffix(".csv")
    if csv_path.is_file():
        return pd.read_csv(csv_path, index_col=index_col, low_memory=False)

    raise FileNotFoundError(
        f"No saved table found for stem {stem}. "
        f"Did you run build_dataset.py first?"
    )


def save_json(obj, path: Path) -> Path:
    """Write an object to JSON with stable formatting."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=False), encoding="utf-8")
    return path


def load_json(path: Path):
    """Read a JSON file."""
    return json.loads(Path(path).read_text(encoding="utf-8"))
