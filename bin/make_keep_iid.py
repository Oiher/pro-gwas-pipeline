#!/usr/bin/env python3
"""
Create a PLINK2 keep file (FID IID) from a TSV table containing IID/#IID.
Requires --psam to look up the real FID for each IID:
  - psam with #FID + IID columns: uses the real FID value
  - psam with #IID only (no FID): uses '0' as FID (plink2 internal default)
Outputs two-column FID<tab>IID, compatible with all plink2 versions.
"""

import argparse
import csv
import sys


def load_fid_map(psam_path: str) -> dict:
    """Return {IID: FID} from a plink2 .psam file."""
    fid_map = {}
    with open(psam_path, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        h0 = header[0].lstrip("#")
        if h0 == "IID":
            # No FID column — plink2 treats missing FID as '0'
            iid_idx, fid_idx = 0, None
        else:
            # Standard: col0=#FID, col1=IID
            fid_idx, iid_idx = 0, 1
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            iid = row[iid_idx].strip()
            fid = row[fid_idx].strip() if fid_idx is not None else "0"
            fid_map[iid] = fid
    return fid_map


def main() -> int:
    parser = argparse.ArgumentParser(description="Build FID/IID PLINK2 keep file")
    parser.add_argument("--input", required=True, help="Input TSV with IID or #IID column")
    parser.add_argument("--output", required=True, help="Output keep file path (FID<tab>IID)")
    parser.add_argument("--psam", required=True, help="plink2 .psam file for FID lookup")
    args = parser.parse_args()

    fid_map = load_fid_map(args.psam)

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
                fid = fid_map.get(iid, iid)  # fall back to FID=IID if not found in psam
                writer.writerow([fid, iid])
                n_written += 1

    if n_written == 0:
        print(f"ERROR: no IID values found in input file: {args.input}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
