#!/bin/bash
# monitor_and_recover.sh - persistent fleet monitor with auto-recovery.
#
# Runs sync_fleet.sh, then checks for spot-interrupted instances and
# automatically relaunches them on on-demand pricing (no further
# interruption risk). Repeats every 15 minutes.
#
# This script is meant to be run via Monitor (or `nohup ... &`) so it
# survives across syncs without needing user approval on each recovery.
#
# The recovery rule: an instance is "lost" if EITHER
#   (a) AWS state is 'terminated' or 'shutting-down'
#   (b) AWS API returns no record (aged out)
# AND the corresponding final.npz is NOT present locally. In that case,
# the loss was a spot interruption (not normal job completion); we
# relaunch on on-demand.

set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

AWS_REGION="${AWS_REGION:-us-west-2}"
KEY="${KEY:-cosmology-sim-key.pem}"
REPO_URL="https://github.com/ShawnaRStaff/Cascade-Inversion-Cosmology.git"
BRANCH="main"
SAFETY_HOURS=120
AMI_ID="ami-06af99ca40d168b3b"
SLEEP_SECONDS="${SLEEP_SECONDS:-900}"  # 15 min default

relaunch_on_demand() {
  local L="$1"
  local SEED="$2"
  local SWEEP_SUBDIR
  SWEEP_SUBDIR=$(cat .aws_fleet_sweep_dir 2>/dev/null) || return 1
  local JOB_FLAGS="--single-job ${L}:${SEED} --sweep-dir data/outputs/${SWEEP_SUBDIR}"
  local USER_DATA
  USER_DATA=$(mktemp)
  sed -e "s|@REPO_URL@|${REPO_URL}|g" \
      -e "s|@BRANCH@|${BRANCH}|g" \
      -e "s|@SWEEP_FLAGS@|${JOB_FLAGS}|g" \
      -e "s|@SAFETY_HOURS@|${SAFETY_HOURS}|g" \
      aws/bootstrap.sh > "${USER_DATA}"
  local INSTANCE_ID
  INSTANCE_ID=$(aws ec2 run-instances \
    --region "${AWS_REGION}" \
    --image-id "${AMI_ID}" \
    --instance-type m7i-flex.large \
    --key-name cosmology-sim-key \
    --security-groups cosmology-sim-sg \
    --user-data "file://${USER_DATA}" \
    --instance-initiated-shutdown-behavior terminate \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=cosmology-sim-fleet-L${L}-s${SEED}-autoondemand},{Key=Project,Value=cascade-inversion-cosmology},{Key=L,Value=${L}},{Key=Seed,Value=${SEED}},{Key=Pricing,Value=on-demand-autorelaunch}]" \
    --query 'Instances[0].InstanceId' \
    --output text 2>/dev/null)
  rm -f "${USER_DATA}"
  if [[ "${INSTANCE_ID}" == i-* ]]; then
    aws ec2 wait instance-running --region "${AWS_REGION}" --instance-ids "${INSTANCE_ID}" 2>/dev/null
    local DNS
    DNS=$(aws ec2 describe-instances --region "${AWS_REGION}" --instance-ids "${INSTANCE_ID}" --query 'Reservations[0].Instances[0].PublicDnsName' --output text 2>/dev/null)
    echo "${INSTANCE_ID}|${DNS}|${L}|${SEED}|$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> .aws_fleet_state
    echo "${INSTANCE_ID}"
    return 0
  fi
  return 1
}

while true; do
  TS=$(date +"%H:%M:%S")
  # 1. Run sync (preserves data, auto-terminates idle done instances).
  SYNC_OUT=$(./aws/sync_fleet.sh 2>&1 || echo "SYNC_FAILED")
  DONE=$(echo "$SYNC_OUT" | grep -oE "done:\s+[0-9]+" | head -1 | grep -oE "[0-9]+")
  TOTAL=$(echo "$SYNC_OUT" | grep -oE "done:\s+[0-9]+/[0-9]+" | head -1 | grep -oE "[0-9]+$")
  RUNNING=$(echo "$SYNC_OUT" | grep -oE "running:\s+[0-9]+" | head -1 | grep -oE "[0-9]+")
  UNREACH=$(echo "$SYNC_OUT" | grep -oE "unreachable:\s+[0-9]+" | head -1 | grep -oE "[0-9]+")
  TERM=$(echo "$SYNC_OUT" | grep -c "terminating idle")

  # 2. Check for spot-interrupted instances (terminated + no final.npz).
  SWEEP_SUBDIR=$(cat .aws_fleet_sweep_dir 2>/dev/null || echo "")
  RECOVERED=0
  RECOVERY_SUMMARY=""
  if [ -n "${SWEEP_SUBDIR}" ] && [ -f .aws_fleet_state ]; then
    TMP_STATE=$(mktemp)
    while IFS='|' read -r INSTANCE_ID DNS L SEED LAUNCHED; do
      STATE=$(aws ec2 describe-instances --region "${AWS_REGION}" --instance-ids "${INSTANCE_ID}" --query 'Reservations[0].Instances[0].State.Name' --output text 2>/dev/null)
      if [ "${STATE}" = "terminated" ] || [ "${STATE}" = "shutting-down" ] || [ -z "${STATE}" ] || [ "${STATE}" = "None" ]; then
        # Truly gone. Check if final.npz on disk.
        if [ ! -f "data/outputs/${SWEEP_SUBDIR}/L${L}_s${SEED}_final.npz" ]; then
          # No final.npz — this was a spot interruption, not normal completion. Recover.
          NEW_ID=$(relaunch_on_demand "${L}" "${SEED}")
          if [ -n "${NEW_ID}" ]; then
            RECOVERED=$((RECOVERED+1))
            RECOVERY_SUMMARY="${RECOVERY_SUMMARY} L=${L}s=${SEED}->ondemand(${NEW_ID})"
          fi
          # Drop the dead entry; the new one was appended by relaunch_on_demand.
          continue
        fi
        # Final.npz present — normal completion, drop the dead entry.
        continue
      fi
      # Still alive — keep entry.
      echo "${INSTANCE_ID}|${DNS}|${L}|${SEED}|${LAUNCHED}"
    done < .aws_fleet_state > "${TMP_STATE}"
    mv "${TMP_STATE}" .aws_fleet_state
  fi

  # 3. Emit one line summary per cycle.
  RECOVERED_NOTE=""
  if [ "${RECOVERED}" -gt 0 ]; then
    RECOVERED_NOTE=" AUTO-RECOVERED ${RECOVERED}:${RECOVERY_SUMMARY}"
  fi
  echo "[${TS}] sync: ${DONE:-?}/${TOTAL:-?} done, ${RUNNING:-?} running, ${UNREACH:-?} unreachable, ${TERM} self-terminated this cycle${RECOVERED_NOTE}"

  sleep "${SLEEP_SECONDS}"
done
