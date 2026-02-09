#!/usr/bin/env python3
"""
Create a PLINK keep file (FID IID) from a TSV table containing IID/#IID.
Outputs two columns with FID fixed to 0 and IID from input.
"""

import argparse
import csv
import sys


def main() -> int:
    parser = argparse.ArgumentParser(description="Build IID-based PLINK keep file")
    parser.add_argument("--input", required=True, help="Input TSV with IID or #IID column")
    parser.add_argument("--output", required=True, help="Output keep file path")
    args = parser.parse_args()

    with open(args.input, "r", newline="") as fin:
        reader = csv.reader(fin, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            print(f"ERROR: input file is empty: {args.input}", file=sys.stderr)
            return 1

        iid_idx = None
        for i, col in enumerate(header):
            if col in ("IID", "#IID"):
                iid_idx = i
                break

        if iid_idx is None:
            print(f"ERROR: IID/#IID column not found in input file: {args.input}", file=sys.stderr)
            return 1

        with open(args.output, "w", newline="") as fout:
            writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
            n_written = 0
            for row in reader:
                if iid_idx >= len(row):
                    continue
                iid = row[iid_idx].strip()
                if not iid:
                    continue
                writer.writerow(["0", iid])
                n_written += 1

    if n_written == 0:
        print(f"ERROR: no IID values found in input file: {args.input}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
