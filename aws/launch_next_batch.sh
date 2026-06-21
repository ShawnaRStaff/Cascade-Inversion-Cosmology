#!/bin/bash
# launch_next_batch.sh - launch the next round of FSS sweep jobs.
#
# New jobs:
#   L=128 seed_idx=3 -> seed=22803  (fresh, reaches 5 total at L=128)
#   L=128 seed_idx=4 -> seed=22804  (fresh, reaches 5 total at L=128)
#   L=192 seed_idx=0 -> seed=29200  (first L=192 point; ~638 h wall time,
#                                     ~5 SAFETY_HOURS cycles to complete)
#
# For the in-progress s22801 resume, run FIRST:
#   ./aws/resume_job.sh 128 22801 fss_sweep_20260521_031056
#
# This script creates its own sweep_dir and overwrites .aws_fleet_state
# (backed up automatically). Run it after s22801 is stable, or accept
# that s22801 runs without monitor auto-recovery (on-demand = no spot risk).
#
# Usage:
#   ./aws/launch_next_batch.sh           # launches jobs
#   ./aws/launch_next_batch.sh --dry-run # prints plan, no launches

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

AWS_REGION="${AWS_REGION:-us-west-2}"
KEY_NAME="${KEY_NAME:-cosmology-sim-key}"
SECURITY_GROUP="${SECURITY_GROUP:-cosmology-sim-sg}"
INSTANCE_TYPE="${INSTANCE_TYPE:-c7i-flex.large}"
REPO_URL="${REPO_URL:-https://github.com/ShawnaRStaff/Cascade-Inversion-Cosmology.git}"
BRANCH="${BRANCH:-main}"
SAFETY_HOURS="${SAFETY_HOURS:-120}"

# L=192 expected wall time ~638h (≈5.3 safety-hour cycles). With the
# checkpoint self-heal, each SAFETY_HOURS kill auto-recovers from the
# latest checkpoint. On-demand for safety (no spot interruption risk).
L192_INSTANCE_TYPE="m7i-flex.large"

# Jobs: seed = 10000 + L*100 + idx (matches seed_for() in Python)
#   L=128 idx=3 -> 10000+12800+3 = 22803
#   L=128 idx=4 -> 10000+12800+4 = 22804
#   L=192 idx=0 -> 10000+19200+0 = 29200
JOBS_128=("128:3" "128:4")
JOBS_192=("192:0")

DRY_RUN=0
if [ "${1:-}" = "--dry-run" ]; then DRY_RUN=1; fi

echo "Next-batch launch plan"
echo "  region:        ${AWS_REGION}"
echo "  L=128 type:    ${INSTANCE_TYPE} (spot)"
echo "  L=192 type:    ${L192_INSTANCE_TYPE} (on-demand, long-running)"
echo "  safety hours:  ${SAFETY_HOURS}"
echo
echo "Jobs:"
for entry in "${JOBS_128[@]}"; do
  L="${entry%%:*}"; idx="${entry##*:}"
  seed=$(( 10000 + L * 100 + idx ))
  printf "  L=%-4d seed=%d  [spot]\n" "${L}" "${seed}"
done
for entry in "${JOBS_192[@]}"; do
  L="${entry%%:*}"; idx="${entry##*:}"
  seed=$(( 10000 + L * 100 + idx ))
  printf "  L=%-4d seed=%d  [on-demand, ~638h wall]\n" "${L}" "${seed}"
done

if [ "${DRY_RUN}" -eq 1 ]; then
  echo
  echo "Dry-run: not launching."
  exit 0
fi

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

STAMP=$(date -u +%Y%m%d_%H%M%S)
SWEEP_SUBDIR="fss_sweep_${STAMP}"

if [ -f .aws_fleet_state ]; then
  mv .aws_fleet_state ".aws_fleet_state.last_$(date -u +%Y%m%d_%H%M%S)"
fi
echo "${SWEEP_SUBDIR}" > .aws_fleet_sweep_dir

launch_one() {
  local L="$1" seed="$2" itype="$3" spot_opts="$4"
  local JOB_FLAGS="--single-job ${L}:${seed} --sweep-dir data/outputs/${SWEEP_SUBDIR}"
  local USER_DATA
  USER_DATA=$(mktemp)
  sed \
    -e "s|@REPO_URL@|${REPO_URL}|g" \
    -e "s|@BRANCH@|${BRANCH}|g" \
    -e "s|@SWEEP_FLAGS@|${JOB_FLAGS}|g" \
    -e "s|@SAFETY_HOURS@|${SAFETY_HOURS}|g" \
    -e "s|@CHECKPOINT_URL_B64@||g" \
    -e "s|@CHECKPOINT_DEST@||g" \
    aws/bootstrap.sh > "${USER_DATA}"

  printf "  L=%-4d seed=%d [%s] ... " "${L}" "${seed}" "${itype}"
  local INSTANCE_ID
  # shellcheck disable=SC2086
  INSTANCE_ID=$(aws ec2 run-instances \
    --region "${AWS_REGION}" \
    --image-id "${AMI_ID}" \
    --instance-type "${itype}" \
    --key-name "${KEY_NAME}" \
    --security-groups "${SECURITY_GROUP}" \
    --user-data "file://${USER_DATA}" \
    --instance-initiated-shutdown-behavior terminate \
    ${spot_opts} \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=cosmology-sim-fleet-L${L}-s${seed}},{Key=Project,Value=cascade-inversion-cosmology},{Key=Fleet,Value=${STAMP}},{Key=L,Value=${L}},{Key=Seed,Value=${seed}}]" \
    --query 'Instances[0].InstanceId' \
    --output text 2>/tmp/launch_err) || {
      echo "FAILED: $(cat /tmp/launch_err)"
      rm -f "${USER_DATA}"
      return 1
    }
  rm -f "${USER_DATA}"
  echo "${INSTANCE_ID}"
  echo "${INSTANCE_ID}|pending|${L}|${seed}|$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> .aws_fleet_state
}

SPOT_OPTS="--instance-market-options MarketType=spot,SpotOptions={SpotInstanceType=one-time,InstanceInterruptionBehavior=terminate}"
ONDEMAND_OPTS=""

echo "Launching..."
for entry in "${JOBS_128[@]}"; do
  L="${entry%%:*}"; idx="${entry##*:}"
  seed=$(( 10000 + L * 100 + idx ))
  launch_one "${L}" "${seed}" "${INSTANCE_TYPE}" "${SPOT_OPTS}" || true
  sleep 1
done
for entry in "${JOBS_192[@]}"; do
  L="${entry%%:*}"; idx="${entry##*:}"
  seed=$(( 10000 + L * 100 + idx ))
  launch_one "${L}" "${seed}" "${L192_INSTANCE_TYPE}" "${ONDEMAND_OPTS}" || true
  sleep 1
done

echo
echo "Resolving DNS for launched instances..."
TMP_STATE=$(mktemp)
while IFS='|' read -r INSTANCE_ID DNS L SEED LAUNCHED; do
  if [ "${DNS}" = "pending" ]; then
    aws ec2 wait instance-running --region "${AWS_REGION}" \
      --instance-ids "${INSTANCE_ID}" 2>/dev/null || true
    DNS=$(aws ec2 describe-instances \
      --region "${AWS_REGION}" \
      --instance-ids "${INSTANCE_ID}" \
      --query 'Reservations[0].Instances[0].PublicDnsName' \
      --output text 2>/dev/null || echo "unknown")
  fi
  echo "${INSTANCE_ID}|${DNS}|${L}|${SEED}|${LAUNCHED}" >> "${TMP_STATE}"
done < .aws_fleet_state
mv "${TMP_STATE}" .aws_fleet_state

echo
echo "========================================================"
echo "BATCH LAUNCHED"
echo "  Fleet stamp:   ${STAMP}"
echo "  Sweep subdir:  data/outputs/${SWEEP_SUBDIR}"
echo "  Instances:     $(wc -l < .aws_fleet_state)"
echo
echo "Next steps:"
echo "  - Monitor:     ./aws/monitor_and_recover.sh"
echo "  - Sync:        ./aws/sync_fleet.sh"
echo "  - s22801 also: ./aws/resume_job.sh 128 22801 fss_sweep_20260521_031056"
echo "========================================================"
