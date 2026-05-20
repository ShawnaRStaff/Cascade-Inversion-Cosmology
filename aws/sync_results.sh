#!/bin/bash
# sync_results.sh - pull the FSS sweep results back from the EC2 instance.
#
# Uses rsync over SSH. Reads instance address from .aws_state written by
# launch_sweep.sh. Idempotent: can be re-run multiple times to pull
# incremental progress while the sweep is still running.
#
# Usage:
#   ./aws/sync_results.sh

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

if [ ! -f .aws_state ]; then
  echo "ERROR: .aws_state not found. Did you run launch_sweep.sh?"
  exit 1
fi
# shellcheck disable=SC1091
source .aws_state

KEY="${KEY:-cosmology-sim-key.pem}"
if [ ! -f "${KEY}" ]; then
  echo "ERROR: ${KEY} not found in repo root"
  exit 1
fi

echo "Syncing from ${PUBLIC_DNS}..."
echo

# Check sweep status (if available)
STATUS=$(ssh -o StrictHostKeyChecking=accept-new -i "${KEY}" "ubuntu@${PUBLIC_DNS}" \
  "cat /home/ubuntu/SWEEP_STATUS 2>/dev/null || echo unknown" 2>/dev/null || echo "(could not query)")
echo "Remote sweep status: ${STATUS}"
echo

# Sweep outputs (data files)
echo "syncing data/outputs/fss_sweep_*..."
mkdir -p data/outputs
rsync -avz --progress \
  -e "ssh -o StrictHostKeyChecking=accept-new -i ${KEY}" \
  "ubuntu@${PUBLIC_DNS}:/home/ubuntu/repo/data/outputs/fss_sweep_*" \
  data/outputs/ || echo "(no sweep outputs yet)"

# Logs
echo
echo "syncing /home/ubuntu/logs/..."
mkdir -p data/aws_logs
rsync -avz --progress \
  -e "ssh -o StrictHostKeyChecking=accept-new -i ${KEY}" \
  "ubuntu@${PUBLIC_DNS}:/home/ubuntu/logs/" \
  data/aws_logs/ || echo "(no logs yet)"

echo
echo "Sync complete."
if [ "${STATUS}" = "done" ]; then
  echo "Sweep marked DONE remotely. Safe to teardown."
elif [ "${STATUS}" = "starting" ]; then
  echo "Sweep still running. Re-run this script later for more results."
fi
