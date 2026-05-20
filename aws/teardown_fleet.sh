#!/bin/bash
# teardown_fleet.sh - terminate every instance in the fleet.
#
# Reads .aws_fleet_state, terminates all instance IDs in one batch,
# archives the state file.
#
# Usage:
#   ./aws/teardown_fleet.sh           # confirms before terminating
#   ./aws/teardown_fleet.sh --force   # no confirmation

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

if [ ! -f .aws_fleet_state ]; then
  echo "ERROR: .aws_fleet_state not found."
  echo "If you have orphan instances, find them via:"
  echo "  aws ec2 describe-instances --filters Name=tag:Project,Values=cascade-inversion-cosmology Name=instance-state-name,Values=running"
  exit 1
fi

AWS_REGION="${AWS_REGION:-us-west-2}"

FORCE=0
if [ "${1:-}" = "--force" ]; then FORCE=1; fi

# Collect instance IDs from state file.
INSTANCE_IDS=()
while IFS='|' read -r INSTANCE_ID DNS L SEED LAUNCHED; do
  INSTANCE_IDS+=("${INSTANCE_ID}")
done < .aws_fleet_state

if [ "${#INSTANCE_IDS[@]}" -eq 0 ]; then
  echo "No instances in fleet state."
  exit 0
fi

echo "About to terminate ${#INSTANCE_IDS[@]} instances:"
while IFS='|' read -r INSTANCE_ID DNS L SEED LAUNCHED; do
  printf "  %s  L=%-4d  seed=%d  launched=%s\n" "${INSTANCE_ID}" "${L}" "${SEED}" "${LAUNCHED}"
done < .aws_fleet_state
echo

if [ "${FORCE}" -eq 0 ]; then
  read -r -p "Confirm terminate all (y/N)? " ans
  if [ "${ans}" != "y" ] && [ "${ans}" != "Y" ]; then
    echo "aborted."
    exit 0
  fi
fi

# Batch-terminate; AWS handles up to 1000 IDs per call.
echo "Terminating ${#INSTANCE_IDS[@]} instances..."
aws ec2 terminate-instances \
  --region "${AWS_REGION}" \
  --instance-ids "${INSTANCE_IDS[@]}" \
  --output text >/dev/null

echo "Waiting for terminated state..."
aws ec2 wait instance-terminated \
  --region "${AWS_REGION}" \
  --instance-ids "${INSTANCE_IDS[@]}"

# Archive state.
mv .aws_fleet_state ".aws_fleet_state.last_$(date -u +%Y%m%d_%H%M%S)"

echo "All ${#INSTANCE_IDS[@]} instances terminated."
echo "Final cost appears in Billing console within ~24 hours."
