"""
splits.py
=========
Creates and freezes the train / validation / test split.

Run after build_dataset.py:

    python splits.py

Output
------
    splits.json   ModelID -> "train" | "val" | "test", plus provenance

Why this is its own module, and why the split is written to disk
---------------------------------------------------------------
The single most common way a first machine-learning project reports a number
that is too good is leakage: information from the test set reaching the model
during training. It produces no error and no warning. It just makes the score
better than the truth, and it is usually found by someone else, in the last
question of a presentation.

Three guards are implemented here.

1. **Split by cell line, never by row.** Obvious once stated, but a random
   row-wise split of a long-format table puts the same cell line on both sides.

2. **Group by patient, not just by cell line.** Several DepMap cell lines can
   derive from the same donor (`PatientID` in Model.csv). Sibling lines from
   one patient are far more similar to each other than to the rest of the
   panel. Letting them straddle the boundary is quiet, real leakage. Grouping
   by patient removes it. This is the guard most projects miss.

3. **Stratify by lineage.** Without it, an unlucky seed can leave a tissue
   type entirely absent from training, or concentrate a rare lineage in test.
   Stratification keeps each split's composition representative.

The assignment is written to `splits.json` and every downstream module reads it
from there. Nothing regenerates it on the fly. This means the split cannot
drift between a training run and an evaluation run, and it means the exact
partition behind any reported number can be inspected months later.

A note on the test set
----------------------
`splits.json` marks the test set, but nothing in this pipeline reads it except
the final evaluation. Use `val` for every decision you make -- feature choices,
hyperparameters, model selection, whether the transformer beat the baseline.
Touch `test` once, at the end, and report whatever it says. Every extra look
at the test set quietly converts it into a validation set.
"""

from __future__ import annotations

import sys
from collections import defaultdict

import numpy as np
import pandas as pd

import config
import io_utils


SPLIT_NAMES = ("train", "val", "test")


def _resolve_groups(metadata: pd.DataFrame) -> tuple[pd.Series, str]:
    """
    Determine the grouping key for each cell line.

    Prefers PatientID. Falls back to ModelID (i.e. no grouping beyond the cell
    line itself) if the column is absent, and treats missing PatientIDs as
    singleton groups so that unknown provenance never merges unrelated lines.
    """
    if config.GROUP_COL in metadata.columns:
        raw = metadata[config.GROUP_COL]
        groups = raw.astype(str).where(raw.notna(), pd.Series(metadata.index, index=metadata.index))
        groups = groups.replace({"nan": None})
        groups = groups.where(groups.notna(), pd.Series(metadata.index, index=metadata.index))
        return groups.astype(str), config.GROUP_COL

    return pd.Series(metadata.index.astype(str), index=metadata.index), config.MODEL_ID


def _resolve_strata(metadata: pd.DataFrame) -> pd.Series:
    """Stratification label per cell line, with a safe default."""
    if config.STRATIFY_COL in metadata.columns:
        return metadata[config.STRATIFY_COL].astype(str).fillna("Unknown")
    return pd.Series("all", index=metadata.index)


def _group_stratum(groups: pd.Series, strata: pd.Series) -> dict[str, str]:
    """
    Assign one stratum label per group.

    If a patient's cell lines somehow span lineages, the most frequent one
    wins; ties are broken alphabetically so the result is deterministic.
    """
    tally: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for model_id in groups.index:
        tally[groups[model_id]][strata[model_id]] += 1

    resolved: dict[str, str] = {}
    for group, counts in tally.items():
        best = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        resolved[group] = best
    return resolved


def assign_splits(
    metadata: pd.DataFrame,
    seed: int = config.RANDOM_SEED,
    test_fraction: float = config.TEST_FRACTION,
    val_fraction: float = config.VAL_FRACTION,
    group_by_lineage: bool = config.GROUP_SPLIT_BY_LINEAGE,
) -> tuple[pd.Series, dict]:
    """
    Assign every cell line to train, val, or test.

    Returns
    -------
    (assignment, report)
        `assignment` is a Series indexed by ModelID with values in SPLIT_NAMES.
    """
    rng = np.random.default_rng(seed)

    groups, group_col = _resolve_groups(metadata)
    strata = _resolve_strata(metadata)

    # ---- Optional: hold out whole lineages -----------------------------
    if group_by_lineage:
        lineages = sorted(strata.unique().tolist())
        shuffled = list(rng.permutation(lineages))
        n_test = max(1, int(round(len(shuffled) * test_fraction)))
        n_val = max(1, int(round(len(shuffled) * val_fraction)))
        test_lineages = set(shuffled[:n_test])
        val_lineages = set(shuffled[n_test:n_test + n_val])

        assignment = pd.Series("train", index=metadata.index, name="split")
        assignment[strata.isin(val_lineages)] = "val"
        assignment[strata.isin(test_lineages)] = "test"

        report = {
            "strategy": "lineage_holdout",
            "seed": seed,
            "group_column": group_col,
            "stratify_column": config.STRATIFY_COL,
            "test_lineages": sorted(test_lineages),
            "val_lineages": sorted(val_lineages),
        }
        return assignment, report

    # ---- Default: group-disjoint, lineage-stratified --------------------
    group_to_stratum = _group_stratum(groups, strata)
    group_sizes = groups.value_counts().to_dict()

    by_stratum: dict[str, list[str]] = defaultdict(list)
    for group, stratum in group_to_stratum.items():
        by_stratum[stratum].append(group)

    split_of_group: dict[str, str] = {}

    for stratum in sorted(by_stratum):
        members = sorted(by_stratum[stratum])
        order = rng.permutation(len(members))
        members = [members[i] for i in order]

        total_lines = sum(group_sizes[g] for g in members)
        target_test = total_lines * test_fraction
        target_val = total_lines * val_fraction

        acc_test = 0.0
        acc_val = 0.0

        for group in members:
            size = group_sizes[group]
            if acc_test + size / 2.0 <= target_test:
                split_of_group[group] = "test"
                acc_test += size
            elif acc_val + size / 2.0 <= target_val:
                split_of_group[group] = "val"
                acc_val += size
            else:
                split_of_group[group] = "train"

    assignment = pd.Series(
        [split_of_group[groups[m]] for m in metadata.index],
        index=metadata.index,
        name="split",
    )

    report = {
        "strategy": "grouped_stratified",
        "seed": seed,
        "group_column": group_col,
        "stratify_column": config.STRATIFY_COL,
        "test_fraction_requested": test_fraction,
        "val_fraction_requested": val_fraction,
        "n_groups": len(group_to_stratum),
        "n_strata": len(by_stratum),
    }
    return assignment, report


def verify_splits(
    assignment: pd.Series,
    metadata: pd.DataFrame,
) -> dict:
    """
    Assert the split is sound. Raises AssertionError on any violation.

    These checks are cheap and they are the difference between a number you
    can defend and a number you cannot.
    """
    checks: dict = {}

    # 1. Every cell line assigned exactly once.
    assert assignment.index.is_unique, "Duplicate ModelIDs in split assignment"
    assert set(assignment.unique()) <= set(SPLIT_NAMES), (
        f"Unexpected split labels: {sorted(set(assignment.unique()))}"
    )
    checks["all_lines_assigned"] = bool(assignment.notna().all())
    assert checks["all_lines_assigned"], "Some cell lines were left unassigned"

    # 2. No cell line in more than one split.
    members = {name: set(assignment.index[assignment == name]) for name in SPLIT_NAMES}
    for a in SPLIT_NAMES:
        for b in SPLIT_NAMES:
            if a < b:
                overlap = members[a] & members[b]
                assert not overlap, (
                    f"{len(overlap)} cell lines appear in both {a} and {b}: "
                    f"{sorted(overlap)[:5]}"
                )
    checks["splits_disjoint"] = True

    # 3. No patient group straddles a split.
    groups, group_col = _resolve_groups(metadata)
    straddling: list[str] = []
    for group, sub in assignment.groupby(groups.reindex(assignment.index)):
        if sub.nunique() > 1:
            straddling.append(str(group))
    checks["group_column"] = group_col
    checks["groups_straddling_splits"] = len(straddling)
    checks["straddling_examples"] = straddling[:5]
    assert not straddling, (
        f"{len(straddling)} groups straddle splits, e.g. {straddling[:5]}. "
        f"This is leakage."
    )

    # 4. No split is empty.
    for name in SPLIT_NAMES:
        assert len(members[name]) > 0, (
            f"Split {name!r} is empty. Adjust TEST_FRACTION / VAL_FRACTION, "
            f"or check that the dataset is large enough."
        )
    checks["split_sizes"] = {name: len(members[name]) for name in SPLIT_NAMES}

    return checks


def summarise(assignment: pd.Series, metadata: pd.DataFrame) -> dict:
    """Per-split composition, including lineage balance and osteosarcoma count."""
    summary: dict = {}
    total = len(assignment)

    for name in SPLIT_NAMES:
        mask = assignment == name
        entry: dict = {
            "n_cell_lines": int(mask.sum()),
            "fraction": round(float(mask.sum()) / total, 4) if total else 0.0,
        }

        if config.STRATIFY_COL in metadata.columns:
            counts = (
                metadata.loc[mask, config.STRATIFY_COL]
                .astype(str)
                .value_counts()
                .to_dict()
            )
            entry["n_lineages"] = len(counts)
            entry["top_lineages"] = dict(list(counts.items())[:5])

        os_mask = config.osteosarcoma_mask(metadata)
        if os_mask.any():
            entry["osteosarcoma_lines"] = int((os_mask & mask).sum())

        summary[name] = entry

    return summary


def main() -> int:
    config.ensure_dirs()
    out = config.PROCESSED_DIR

    print("=" * 74)
    print("BUILDING TRAIN / VAL / TEST SPLIT")
    print("=" * 74)

    try:
        metadata = io_utils.load_table(out / "model_metadata")
    except FileNotFoundError as exc:
        print(f"\n{exc}")
        return 1

    metadata.index.name = config.MODEL_ID
    print(f"\nLoaded metadata for {len(metadata)} cell lines")

    assignment, report = assign_splits(metadata)
    print(f"Strategy      : {report['strategy']}")
    print(f"Grouped by    : {report['group_column']}")
    print(f"Stratified by : {report['stratify_column']}")
    print(f"Seed          : {report['seed']}")

    print("\nVerifying...")
    checks = verify_splits(assignment, metadata)
    print(f"  splits disjoint          : {checks['splits_disjoint']}")
    print(f"  groups straddling splits : {checks['groups_straddling_splits']}")
    print("  All leakage checks passed.")

    summary = summarise(assignment, metadata)
    print("\nComposition:")
    for name in SPLIT_NAMES:
        entry = summary[name]
        line = (f"  {name:<6}: {entry['n_cell_lines']:>5} lines "
                f"({entry['fraction']:.1%})")
        if "osteosarcoma_lines" in entry:
            line += f"   osteosarcoma: {entry['osteosarcoma_lines']}"
        print(line)

    payload = {
        "description": (
            "Frozen train/val/test assignment. Every downstream module reads "
            "this file. Do not regenerate it between a training run and an "
            "evaluation run. Use 'val' for all decisions; touch 'test' once."
        ),
        "provenance": report,
        "verification": checks,
        "summary": summary,
        "assignment": assignment.to_dict(),
    }
    path = io_utils.save_json(payload, out / "splits.json")
    print(f"\nWritten to: {path}")
    print("\nNext: python checks.py")
    return 0


def load_splits() -> pd.Series:
    """Read the frozen assignment. Used by every downstream module."""
    payload = io_utils.load_json(config.PROCESSED_DIR / "splits.json")
    return pd.Series(payload["assignment"], name="split")


if __name__ == "__main__":
    sys.exit(main())
