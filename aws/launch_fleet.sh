#!/bin/bash
# launch_fleet.sh - launch a fleet of small EC2 instances, one per job.
#
# When the AWS account is restricted to small instance types (Free Plan
# can only use m7i-flex.large / c7i-flex.large), we parallelize by
# launching MANY small instances simultaneously, one per (L, seed) job.
# Each instance runs exactly one job via the --single-job mode of the
# sweep driver.
#
# State is tracked in .aws_fleet_state (one line per instance:
# "instance_id|public_dns|L|seed|launched_at"). This is consumed by
# sync_fleet.sh and teardown_fleet.sh.
#
# Usage:
#   ./aws/launch_fleet.sh           # launches the default job list
#   ./aws/launch_fleet.sh --dry-run # prints job list, no launches
#
# Default job list (matches the FSS sweep plan, minus L=32 which we
# already have from yesterday's shakedown):
#   L=48 x 5 seeds, L=64 x 5 seeds, L=96 x 5 seeds, L=128 x 3 seeds
#   = 18 instances

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

AWS_REGION="${AWS_REGION:-us-west-2}"
KEY_NAME="${KEY_NAME:-cosmology-sim-key}"
SECURITY_GROUP="${SECURITY_GROUP:-cosmology-sim-sg}"
INSTANCE_TYPE="${INSTANCE_TYPE:-c7i-flex.large}"
REPO_URL="${REPO_URL:-https://github.com/ShawnaRStaff/Cascade-Inversion-Cosmology.git}"
BRANCH="${BRANCH:-main}"
SAFETY_HOURS="${SAFETY_HOURS:-120}"  # 5 days, exceeds longest expected job

# Job list: "L:seed_idx" pairs. Seeds are computed by the sweep driver
# as 10000 + L*100 + idx, matching seed_for() in the Python code.
JOBS=(
  # L=48: 5 new seeds
  "48:0"  "48:1"  "48:2"  "48:3"  "48:4"
  # L=64: 5 new seeds
  "64:0"  "64:1"  "64:2"  "64:3"  "64:4"
  # L=96: 5 new seeds
  "96:0"  "96:1"  "96:2"  "96:3"  "96:4"
  # L=128: 3 new seeds
  "128:0" "128:1" "128:2"
)

# Convert L:seed_idx to L:seed (the actual seed number).
declare -A JOB_SPECS  # key: "L_idx", value: "L:seed"
for entry in "${JOBS[@]}"; do
  L="${entry%%:*}"
  idx="${entry##*:}"
  seed=$(( 10000 + L * 100 + idx ))
  JOB_SPECS["${L}_${idx}"]="${L}:${seed}"
done

DRY_RUN=0
if [ "${1:-}" = "--dry-run" ]; then DRY_RUN=1; fi

echo "Fleet launch plan"
echo "  region:         ${AWS_REGION}"
echo "  instance type:  ${INSTANCE_TYPE}"
echo "  safety hours:   ${SAFETY_HOURS}"
echo "  total jobs:     ${#JOBS[@]}"
echo
echo "Jobs (L:seed):"
for entry in "${JOBS[@]}"; do
  L="${entry%%:*}"
  idx="${entry##*:}"
  seed=$(( 10000 + L * 100 + idx ))
  printf "  L=%-4d seed=%d\n" "${L}" "${seed}"
done

if [ "${DRY_RUN}" -eq 1 ]; then
  echo
  echo "Dry-run: not launching."
  exit 0
fi

# Resolve AMI.
AMI_ID=$(aws ec2 describe-images \
  --region "${AWS_REGION}" \
  --owners 099720109477 \
  --filters \
    "Name=name,Values=ubuntu/images/hvm-ssd-gp3/ubuntu-noble-24.04-amd64-server-*" \
    "Name=state,Values=available" \
  --query 'Images | sort_by(@, &CreationDate)[-1].ImageId' \
  --output text)
echo
echo "AMI: ${AMI_ID}"

# Use a single sweep_dir name shared across all instances so all
# results land in the same directory after sync.
STAMP=$(date -u +%Y%m%d_%H%M%S)
SWEEP_SUBDIR="fss_sweep_${STAMP}"

# Truncate fleet state file (and back up if present).
if [ -f .aws_fleet_state ]; then
  mv .aws_fleet_state ".aws_fleet_state.last_$(date -u +%Y%m%d_%H%M%S)"
fi

echo "Launching ${#JOBS[@]} instances..."
LAUNCH_FAILED=0

for entry in "${JOBS[@]}"; do
  L="${entry%%:*}"
  idx="${entry##*:}"
  seed=$(( 10000 + L * 100 + idx ))
  SPEC="${L}:${seed}"

  # Per-job sweep flags for the driver, telling it to do one job and
  # write outputs into the shared sweep subdir.
  JOB_SWEEP_FLAGS="--single-job ${SPEC} --sweep-dir data/outputs/${SWEEP_SUBDIR}"

  # Render bootstrap.sh template for this job.
  USER_DATA=$(mktemp)
  sed \
    -e "s|@REPO_URL@|${REPO_URL}|g" \
    -e "s|@BRANCH@|${BRANCH}|g" \
    -e "s|@SWEEP_FLAGS@|${JOB_SWEEP_FLAGS}|g" \
    -e "s|@SAFETY_HOURS@|${SAFETY_HOURS}|g" \
    aws/bootstrap.sh > "${USER_DATA}"

  printf "  L=%-4d seed=%d ... " "${L}" "${seed}"
  INSTANCE_ID=$(aws ec2 run-instances \
    --region "${AWS_REGION}" \
    --image-id "${AMI_ID}" \
    --instance-type "${INSTANCE_TYPE}" \
    --key-name "${KEY_NAME}" \
    --security-groups "${SECURITY_GROUP}" \
    --user-data "file://${USER_DATA}" \
    --instance-initiated-shutdown-behavior terminate \
    --instance-market-options 'MarketType=spot,SpotOptions={SpotInstanceType=one-time,InstanceInterruptionBehavior=terminate}' \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=cosmology-sim-fleet-L${L}-s${seed}},{Key=Project,Value=cascade-inversion-cosmology},{Key=Fleet,Value=${STAMP}},{Key=L,Value=${L}},{Key=Seed,Value=${seed}}]" \
    --query 'Instances[0].InstanceId' \
    --output text 2>/tmp/launch_err) || {
      echo "FAILED: $(cat /tmp/launch_err)"
      LAUNCH_FAILED=$((LAUNCH_FAILED + 1))
      rm -f "${USER_DATA}"
      continue
    }
  rm -f "${USER_DATA}"
  echo "${INSTANCE_ID}"

  # Append to fleet state. Public DNS comes later, after instances are running.
  echo "${INSTANCE_ID}|pending|${L}|${seed}|$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> .aws_fleet_state

  # Brief sleep to avoid hitting AWS API rate limits.
  sleep 1
done

echo
if [ "${LAUNCH_FAILED}" -gt 0 ]; then
  echo "WARNING: ${LAUNCH_FAILED} launches failed."
fi

# Wait for instances to be running, then fill in PublicDnsName.
echo "Waiting for instances to enter 'running' state and resolving DNS..."
TMP_STATE=$(mktemp)
while IFS='|' read -r INSTANCE_ID DNS L SEED LAUNCHED; do
  if [ "${DNS}" = "pending" ]; then
    aws ec2 wait instance-running --region "${AWS_REGION}" --instance-ids "${INSTANCE_ID}" 2>/dev/null || true
    DNS=$(aws ec2 describe-instances \
      --region "${AWS_REGION}" \
      --instance-ids "${INSTANCE_ID}" \
      --query 'Reservations[0].Instances[0].PublicDnsName' \
      --output text 2>/dev/null || echo "unknown")
  fi
  echo "${INSTANCE_ID}|${DNS}|${L}|${SEED}|${LAUNCHED}" >> "${TMP_STATE}"
done < .aws_fleet_state
mv "${TMP_STATE}" .aws_fleet_state

# Summary
echo
echo "========================================================"
echo "FLEET LAUNCHED"
echo "  Fleet ID:      ${STAMP}"
echo "  Sweep subdir:  data/outputs/${SWEEP_SUBDIR}"
echo "  Instances:     $(wc -l < .aws_fleet_state)"
echo "  State file:    .aws_fleet_state"
echo
echo "Next steps:"
echo "  - Check progress:  ./aws/sync_fleet.sh"
echo "  - Pull results:    ./aws/sync_fleet.sh"
echo "  - Tear down all:   ./aws/teardown_fleet.sh"
echo "========================================================"
