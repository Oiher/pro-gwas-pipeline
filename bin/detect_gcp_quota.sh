#!/usr/bin/env bash
#
# Look up a GCP project's actual regional compute quota (SSD_TOTAL_GB, CPUS), to plug into
# -profile gcb_scaleable's gcb_ssd_quota_gb/gcb_cpu_quota params instead of assuming
# gcb_final.config's tuned baseline (500GB SSD / 200 CPU in europe-west4).
#
# Usage:
#   bin/detect_gcp_quota.sh <region> [project]
#
# Arguments:
#   region:  GCP region to check, e.g. europe-west4
#   project: GCP project ID (default: $GOOGLE_CLOUD_PROJECT, pre-set on VWB VMs)
#
# Run this on the VWB VM directly (already authenticated) or inside the pipeline container.

set -euo pipefail

REGION="${1:-}"
PROJECT="${2:-${GOOGLE_CLOUD_PROJECT:-}}"

if [ -z "$REGION" ]; then
    echo "Error: region required" >&2
    echo "Usage: $0 <region> [project]" >&2
    echo "Example: $0 europe-west4" >&2
    exit 1
fi

if [ -z "$PROJECT" ]; then
    echo "Error: no project given and \$GOOGLE_CLOUD_PROJECT is not set." >&2
    exit 1
fi

if ! command -v gcloud &> /dev/null; then
    echo "Error: gcloud not found. Run this on the VWB VM or inside the pipeline container." >&2
    exit 1
fi

echo "Looking up quota for project=${PROJECT} region=${REGION}..." >&2

gcloud compute regions describe "$REGION" --project "$PROJECT" --format=json | python3 -c '
import json, sys

data = json.load(sys.stdin)
quotas = {q["metric"]: q for q in data.get("quotas", [])}
ssd = quotas.get("SSD_TOTAL_GB")
cpu = quotas.get("CPUS")

if not ssd or not cpu:
    print("Could not find SSD_TOTAL_GB/CPUS quota metrics in the response.", file=sys.stderr)
    sys.exit(1)

ssd_limit = ssd["limit"]
ssd_usage = ssd["usage"]
cpu_limit = cpu["limit"]
cpu_usage = cpu["usage"]

print(f"SSD_TOTAL_GB: limit={ssd_limit:.0f}  current_usage={ssd_usage:.0f}")
print(f"CPUS:         limit={cpu_limit:.0f}  current_usage={cpu_usage:.0f}")
print()
print("Add to your -params-file (or pass as --param) for -profile gcb_scaleable:")
print(f"  gcb_ssd_quota_gb: {ssd_limit:.0f}")
print(f"  gcb_cpu_quota: {cpu_limit:.0f}")
'
