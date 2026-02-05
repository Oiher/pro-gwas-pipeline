#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: normalize_psam_iid_only.sh <file.psam>" >&2
  exit 1
fi

PSAM="$1"
TMP="${PSAM}.iidtmp"

if [ ! -f "${PSAM}" ]; then
  echo "ERROR: psam file not found: ${PSAM}" >&2
  exit 1
fi

if ! awk -F'\t' -v OFS='\t' '
  NR==1 {
    for (i=1; i<=NF; i++) {
      if ($i=="#FID" || $i=="FID") fid=i;
      if ($i=="IID" || $i=="#IID") iid=i;
    }
    if (!iid) {
      print "ERROR: IID column not found in " FILENAME > "/dev/stderr";
      exit 2;
    }

    header="";
    for (i=1; i<=NF; i++) {
      if (i==fid) continue;
      col=$i;
      if (i==iid) col="#IID";
      header = (header=="" ? col : header OFS col);
    }
    print header;
    next;
  }
  {
    if (fid) {
      if (seen[$iid]++) dup=1;
      out="";
      for (i=1; i<=NF; i++) {
        if (i==fid) continue;
        out = (out=="" ? $i : out OFS $i);
      }
      print out;
    } else {
      print $0;
    }
  }
  END {
    if (fid && dup) {
      print "ERROR: duplicate IID values found; cannot drop FID safely." > "/dev/stderr";
      exit 3;
    }
  }
' "${PSAM}" > "${TMP}"; then
  status=$?
  rm -f "${TMP}"
  exit $status
fi

mv "${TMP}" "${PSAM}"
