#!/usr/bin/env python3
"""
Harmonize categorical covariates by collapsing rare levels.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Collapse rare categorical levels in a sample/covariate TSV."
    )
    parser.add_argument("--input", required=True, help="Input TSV path")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument(
        "--categorical",
        default="",
        help="Space-separated categorical covariate column names",
    )
    parser.add_argument(
        "--min-count",
        type=int,
        default=20,
        help="Minimum count required to keep a category level (default: 20)",
    )
    parser.add_argument(
        "--other-suffix",
        default="other",
        help="Suffix for merged rare level label (default: other)",
    )
    args = parser.parse_args()

    cat_cols = [c for c in args.categorical.split() if c]
    in_path = Path(args.input)
    out_path = Path(args.output)
    mapping_path = out_path.with_suffix(".cat_mapping.tsv")
    summary_path = out_path.with_suffix(".cat_summary.tsv")

    try:
        df = pd.read_csv(in_path, sep="\t", engine="c")
    except Exception as exc:
        print(f"ERROR: failed to read input file '{in_path}': {exc}", file=sys.stderr)
        return 1

    if not cat_cols:
        df.to_csv(out_path, sep="\t", index=False)
        pd.DataFrame(columns=["variable", "original_level", "new_level"]).to_csv(
            mapping_path, sep="\t", index=False
        )
        pd.DataFrame(columns=["variable", "level", "count"]).to_csv(
            summary_path, sep="\t", index=False
        )
        return 0

    mapping_rows = []
    summary_rows = []

    for col in cat_cols:
        if col not in df.columns:
            print(f"WARNING: categorical covariate '{col}' not found; skipping", file=sys.stderr)
            continue

        series = df[col]
        non_null = series.dropna().astype(str)
        if non_null.empty:
            continue

        before_counts = non_null.value_counts()
        other_label = f"{col}_{args.other_suffix}"
        rare_levels = set(before_counts[before_counts < args.min_count].index.tolist())

        harmonized = series.astype("string")
        mask = harmonized.isin(list(rare_levels))
        harmonized = harmonized.mask(mask, other_label)
        df[col] = harmonized

        after_counts = df[col].dropna().astype(str).value_counts()

        for level, count in before_counts.items():
            new_level = other_label if level in rare_levels else level
            mapping_rows.append(
                {"variable": col, "original_level": level, "new_level": new_level, "count": int(count)}
            )

        for level, count in after_counts.items():
            summary_rows.append({"variable": col, "level": level, "count": int(count)})

    df.to_csv(out_path, sep="\t", index=False)
    pd.DataFrame(mapping_rows).to_csv(mapping_path, sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
