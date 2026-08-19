"""
gene_ids.py
===========
Parsing and reconciliation of DepMap gene column labels.

DepMap labels gene columns as ``SYMBOL (ENTREZ)`` -- for example ``A1BG (1)``
or ``PRRC2B (84726)``. This is convenient because the Entrez ID is a stable,
unambiguous key, whereas HGNC symbols are not: they get renamed over time and
a handful are genuinely ambiguous.

Everything downstream therefore keys on the **Entrez ID**, while keeping the
symbol for human readability.

Why this module exists in its own file
--------------------------------------
The canonical gene ordering produced here is written to disk and reused every
time new expression data is scored by the trained model -- including, later,
a patient tumour profile from an entirely different source. If that data is
not mapped into exactly this column space, in exactly this order, the model
silently receives scrambled features and produces confident nonsense. Getting
this right once, in one place, is worth the separate module.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import pandas as pd

# ``SYMBOL (ENTREZ)`` -- symbol may contain spaces, hyphens, dots, slashes.
_GENE_LABEL_RE = re.compile(r"^\s*(?P<symbol>.+?)\s*\((?P<entrez>\d+)\)\s*$")


@dataclass(frozen=True)
class ParsedGene:
    """A single parsed gene column label."""
    label: str      # original column label, verbatim
    symbol: str     # HGNC symbol
    entrez: str     # Entrez ID as a string (avoids int/float coercion issues)


def parse_gene_label(label: str) -> ParsedGene | None:
    """
    Parse one ``SYMBOL (ENTREZ)`` label.

    Returns None if the label does not match the expected pattern, so callers
    can count and report unparseable columns rather than crashing.
    """
    match = _GENE_LABEL_RE.match(str(label))
    if match is None:
        return None
    return ParsedGene(
        label=str(label),
        symbol=match.group("symbol").strip(),
        entrez=match.group("entrez"),
    )


def parse_gene_labels(labels) -> tuple[list[ParsedGene], list[str]]:
    """
    Parse an iterable of column labels.

    Returns
    -------
    (parsed, unparsed)
        `parsed` is a list of ParsedGene; `unparsed` is a list of the raw
        labels that did not match the pattern.
    """
    parsed: list[ParsedGene] = []
    unparsed: list[str] = []
    for label in labels:
        result = parse_gene_label(label)
        if result is None:
            unparsed.append(str(label))
        else:
            parsed.append(result)
    return parsed, unparsed


def build_gene_frame(labels, source: str) -> pd.DataFrame:
    """
    Build a tidy DataFrame describing the gene columns of one matrix.

    Columns: entrez, symbol, label, source
    Indexed by entrez (duplicates retained so callers can detect them).
    """
    parsed, _unparsed = parse_gene_labels(labels)
    frame = pd.DataFrame(
        {
            "entrez": [g.entrez for g in parsed],
            "symbol": [g.symbol for g in parsed],
            "label": [g.label for g in parsed],
            "source": source,
        }
    )
    return frame


def deduplicate_by_entrez(frame: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """
    Collapse a gene frame to one row per Entrez ID.

    DepMap matrices should already be unique per Entrez, but this guards
    against release quirks. Keeps the first occurrence.

    Returns
    -------
    (deduplicated_frame, duplicate_entrez_ids)
    """
    duplicated_mask = frame["entrez"].duplicated(keep=False)
    duplicate_ids = sorted(frame.loc[duplicated_mask, "entrez"].unique().tolist())
    deduped = frame.drop_duplicates(subset="entrez", keep="first").reset_index(drop=True)
    return deduped, duplicate_ids


@dataclass
class GeneSpace:
    """
    The canonical gene space shared by two matrices.

    Attributes
    ----------
    entrez_ids
        Ordered list of Entrez IDs defining the canonical feature order.
    labels_a / labels_b
        The original column labels in each source matrix, aligned to
        `entrez_ids`. Use these to reindex the raw matrices.
    symbols
        Human-readable symbols aligned to `entrez_ids`.
    report
        Dict of counts describing what was kept and dropped.
    """
    entrez_ids: list[str]
    labels_a: list[str]
    labels_b: list[str]
    symbols: list[str]
    report: dict


def intersect_gene_spaces(
    labels_a,
    labels_b,
    source_a: str = "a",
    source_b: str = "b",
) -> GeneSpace:
    """
    Compute the shared gene space of two DepMap matrices.

    Genes are matched on Entrez ID, not symbol. The resulting order is sorted
    by integer Entrez ID so it is deterministic and reproducible across runs
    and machines.

    Parameters
    ----------
    labels_a, labels_b
        Column labels from each matrix (gene columns only).
    source_a, source_b
        Names used in the report.

    Returns
    -------
    GeneSpace
    """
    parsed_a, unparsed_a = parse_gene_labels(labels_a)
    parsed_b, unparsed_b = parse_gene_labels(labels_b)

    frame_a = build_gene_frame([g.label for g in parsed_a], source_a)
    frame_b = build_gene_frame([g.label for g in parsed_b], source_b)

    frame_a, dup_a = deduplicate_by_entrez(frame_a)
    frame_b, dup_b = deduplicate_by_entrez(frame_b)

    map_a = dict(zip(frame_a["entrez"], frame_a["label"]))
    map_b = dict(zip(frame_b["entrez"], frame_b["label"]))
    symbol_a = dict(zip(frame_a["entrez"], frame_a["symbol"]))
    symbol_b = dict(zip(frame_b["entrez"], frame_b["symbol"]))

    shared = sorted(set(map_a) & set(map_b), key=int)

    # Where the two sources disagree on the symbol for the same Entrez ID,
    # prefer source A but record the disagreement.
    symbol_conflicts = [
        {"entrez": e, source_a: symbol_a[e], source_b: symbol_b[e]}
        for e in shared
        if symbol_a[e] != symbol_b[e]
    ]

    report = {
        f"{source_a}_columns_given": len(list(labels_a)),
        f"{source_a}_parsed": len(parsed_a),
        f"{source_a}_unparsed": len(unparsed_a),
        f"{source_a}_unparsed_examples": unparsed_a[:10],
        f"{source_a}_duplicate_entrez": len(dup_a),
        f"{source_b}_columns_given": len(list(labels_b)),
        f"{source_b}_parsed": len(parsed_b),
        f"{source_b}_unparsed": len(unparsed_b),
        f"{source_b}_unparsed_examples": unparsed_b[:10],
        f"{source_b}_duplicate_entrez": len(dup_b),
        "shared_genes": len(shared),
        f"only_in_{source_a}": len(set(map_a) - set(map_b)),
        f"only_in_{source_b}": len(set(map_b) - set(map_a)),
        "symbol_conflicts": len(symbol_conflicts),
        "symbol_conflict_examples": symbol_conflicts[:10],
    }

    return GeneSpace(
        entrez_ids=shared,
        labels_a=[map_a[e] for e in shared],
        labels_b=[map_b[e] for e in shared],
        symbols=[symbol_a[e] for e in shared],
        report=report,
    )


def canonical_labels(entrez_ids: list[str], symbols: list[str]) -> list[str]:
    """
    Build canonical ``SYMBOL (ENTREZ)`` labels for the shared gene space.

    These are the column names written to disk and expected by the trained
    model at inference time.
    """
    if len(entrez_ids) != len(symbols):
        raise ValueError("entrez_ids and symbols must be the same length")
    return [f"{sym} ({ent})" for sym, ent in zip(symbols, entrez_ids)]


def map_external_matrix(
    external: pd.DataFrame,
    entrez_ids: list[str],
    canonical_cols: list[str],
    fill_value: float = 0.0,
) -> tuple[pd.DataFrame, dict]:
    """
    Reindex an *external* expression matrix into the canonical gene space.

    This is the function that will later accept a patient tumour expression
    profile, quantified by a completely different pipeline, and force it into
    exactly the feature space the model was trained on.

    The external matrix may label its columns as bare symbols, bare Entrez IDs,
    or ``SYMBOL (ENTREZ)``. All three are handled.

    Parameters
    ----------
    external
        Samples (rows) x genes (columns).
    entrez_ids
        Canonical Entrez ID order.
    canonical_cols
        Canonical column labels, aligned to `entrez_ids`.
    fill_value
        Value used for canonical genes absent from the external matrix.

    Returns
    -------
    (mapped_frame, report)
        `mapped_frame` has exactly `canonical_cols` as its columns, in order.
    """
    lookup: dict[str, str] = {}   # entrez -> external column name

    for col in external.columns:
        col_str = str(col)

        parsed = parse_gene_label(col_str)
        if parsed is not None:
            lookup.setdefault(parsed.entrez, col_str)
            continue

        stripped = col_str.strip()
        if stripped.isdigit():                      # bare Entrez ID
            lookup.setdefault(stripped, col_str)

    # Second pass: match remaining canonical genes by symbol.
    symbol_to_entrez = {
        lab.rsplit(" (", 1)[0]: ent
        for lab, ent in zip(canonical_cols, entrez_ids)
    }
    for col in external.columns:
        col_str = str(col).strip()
        if col_str in symbol_to_entrez:
            lookup.setdefault(symbol_to_entrez[col_str], str(col))

    matched = [e for e in entrez_ids if e in lookup]
    missing = [e for e in entrez_ids if e not in lookup]

    mapped = pd.DataFrame(
        fill_value,
        index=external.index,
        columns=canonical_cols,
        dtype="float32",
    )
    if matched:
        source_cols = [lookup[e] for e in matched]
        target_cols = [
            canonical_cols[entrez_ids.index(e)] for e in matched
        ]
        mapped[target_cols] = external[source_cols].to_numpy(dtype="float32")

    report = {
        "canonical_genes": len(entrez_ids),
        "matched": len(matched),
        "missing_filled": len(missing),
        "external_columns": external.shape[1],
        "missing_examples": missing[:10],
    }
    return mapped, report
