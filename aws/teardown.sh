#!/bin/bash
# teardown.sh - terminate the EC2 instance launched by launch_sweep.sh.
#
# Reads instance ID from .aws_state. After teardown, archives the
# state file to .aws_state.last so we have a record.
#
# Usage:
#   ./aws/teardown.sh           # confirms before terminating
#   ./aws/teardown.sh --force   # no confirmation

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

if [ ! -f .aws_state ]; then
  echo "ERROR: .aws_state not found."
  echo "If you have an orphan instance, find it via:"
  echo "  aws ec2 describe-instances --filters Name=tag:Project,Values=cascade-inversion-cosmology"
  echo "and terminate manually:"
  echo "  aws ec2 terminate-instances --instance-ids <id>"
  exit 1
fi
# shellcheck disable=SC1091
source .aws_state

FORCE=0
if [ "${1:-}" = "--force" ]; then FORCE=1; fi

echo "About to terminate:"
echo "  Instance ID: ${INSTANCE_ID}"
echo "  Public DNS:  ${PUBLIC_DNS}"
echo "  Region:      ${AWS_REGION}"
echo "  Launched at: ${LAUNCHED_AT}"
echo

if [ "${FORCE}" -eq 0 ]; then
  read -r -p "Confirm terminate (y/N)? " ans
  if [ "${ans}" != "y" ] && [ "${ans}" != "Y" ]; then
    echo "aborted."
    exit 0
  fi
fi

echo "terminating ${INSTANCE_ID}..."
aws ec2 terminate-instances \
  --region "${AWS_REGION}" \
  --instance-ids "${INSTANCE_ID}" \
  --output text >/dev/null

aws ec2 wait instance-terminated \
  --region "${AWS_REGION}" \
  --instance-ids "${INSTANCE_ID}"

# Archive state for the record.
mv .aws_state ".aws_state.last_$(date -u +%Y%m%d_%H%M%S)"

echo "done. Instance terminated."
echo "Final cost will appear in the Billing console within ~24 hours."
