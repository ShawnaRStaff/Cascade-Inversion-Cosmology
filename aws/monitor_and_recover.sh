#!/bin/bash
# monitor_and_recover.sh - persistent fleet monitor with auto-recovery.
#
# Runs sync_fleet.sh, then checks for lost instances and automatically
# relaunches them on spot pricing (~60% cheaper than on-demand), falling
# back to on-demand only when spot capacity is unavailable. Repeats
# every 15 minutes.
#
# This script is meant to be run via Monitor (or `nohup ... &`) so it
# survives across syncs without needing user approval on each recovery.
#
# The recovery rule: an instance is "lost" if EITHER
#   (a) AWS state is 'terminated' or 'shutting-down'
#   (b) AWS API returns no record (aged out)
# AND the corresponding final.npz is NOT present locally. In that case,
# the loss was an interruption (spot reclaim or the 120 h safety
# shutdown), not normal completion; we relaunch from the latest local
# checkpoint, losing at most the ~10 min since it was written.

set -uo pipefail

# Pure loss-detection decision logic (unit-tested in
# tests/test_monitor_recovery_logic.py). See that file and
# aws/lib_recovery.sh for the 2026-08-04 incident this guards against.
source "$(dirname "${BASH_SOURCE[0]}")/lib_recovery.sh"

# count_running_duplicates L SEED — how many pending/running instances
# already carry this job's tags. Echoes 999 when the API call fails so
# should_relaunch answers "no" (never launch blind).
count_running_duplicates() {
  local l="$1" seed="$2" out rc
  out=$(aws ec2 describe-instances --region "${AWS_REGION}" \
    --filters "Name=tag:Project,Values=cascade-inversion-cosmology" \
              "Name=tag:L,Values=${l}" "Name=tag:Seed,Values=${seed}" \
              "Name=instance-state-name,Values=pending,running" \
    --query 'length(Reservations[].Instances[])' --output text 2>/dev/null)
  rc=$?
  if [ ${rc} -ne 0 ] || ! [[ "${out}" =~ ^[0-9]+$ ]]; then
    echo "999"
  else
    echo "${out}"
  fi
}

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

  # Upload the latest checkpoint to S3 so the relaunched instance can resume
  # rather than re-running from scratch. The presigned URL is base64-encoded
  # to keep it safe for sed substitution (presigned URLs contain '&', '?', '=').
  # The bucket is account-agnostic; if any link in the account-id -> upload ->
  # presign chain fails while a local checkpoint EXISTS, we must abort the
  # relaunch (relaunch_preseed_verdict) — launching without the pre-seed
  # would silently restart the job from scratch.
  local ACCOUNT_ID
  ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text 2>/dev/null)
  local CKPT_LOCAL="data/outputs/${SWEEP_SUBDIR}/L${L}_s${SEED}_ckpt.npz"
  local CKPT_EXISTS="no"
  [ -f "${CKPT_LOCAL}" ] && CKPT_EXISTS="yes"
  local PRESIGNED=""
  local CHECKPOINT_URL_B64=""
  local CHECKPOINT_DEST_VAL=""
  if [ "${CKPT_EXISTS}" = "yes" ] && [ "$(valid_account_id "${ACCOUNT_ID}")" = "yes" ]; then
    local S3_BUCKET="cascade-cosmo-ckpts-${ACCOUNT_ID}"
    local S3_KEY="checkpoints/${SWEEP_SUBDIR}/L${L}_s${SEED}_ckpt.npz"
    aws s3 mb "s3://${S3_BUCKET}" --region "${AWS_REGION}" 2>/dev/null || true
    if aws s3 cp "${CKPT_LOCAL}" "s3://${S3_BUCKET}/${S3_KEY}" \
        --region "${AWS_REGION}" 2>/dev/null; then
      PRESIGNED=$(aws s3 presign "s3://${S3_BUCKET}/${S3_KEY}" \
        --region "${AWS_REGION}" --expires-in 86400 2>/dev/null || true)
    fi
  fi
  case "$(relaunch_preseed_verdict "${CKPT_EXISTS}" "${PRESIGNED}")" in
    abort)
      echo "  [L=${L} s=${SEED}] checkpoint exists but pre-seed failed (account='${ACCOUNT_ID}') — relaunch ABORTED, will retry next cycle" >&2
      return 1
      ;;
    proceed)
      CHECKPOINT_URL_B64=$(echo "${PRESIGNED}" | base64 -w 0)
      CHECKPOINT_DEST_VAL="${CKPT_LOCAL}"
      echo "  checkpoint uploaded: checkpoints/${SWEEP_SUBDIR}/L${L}_s${SEED}_ckpt.npz"
      ;;
    proceed-fresh)
      : # no checkpoint anywhere — a genuine fresh start is correct
      ;;
  esac

  local USER_DATA
  USER_DATA=$(mktemp)
  sed -e "s|@REPO_URL@|${REPO_URL}|g" \
      -e "s|@BRANCH@|${BRANCH}|g" \
      -e "s|@SWEEP_FLAGS@|${JOB_FLAGS}|g" \
      -e "s|@SAFETY_HOURS@|${SAFETY_HOURS}|g" \
      -e "s|@CHECKPOINT_URL_B64@|${CHECKPOINT_URL_B64}|g" \
      -e "s|@CHECKPOINT_DEST@|${CHECKPOINT_DEST_VAL}|g" \
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
      RC=$?
      CLASS=$(classify_instance_state "${RC}" "${STATE}")

      if [ "${CLASS}" = "alive" ]; then
        echo "${INSTANCE_ID}|${DNS}|${L}|${SEED}|${LAUNCHED}"
        continue
      fi

      if [ "${CLASS}" = "indeterminate" ]; then
        # API failure or unexpected state (incl. stopped): we know
        # nothing for certain. Keep the entry, do NOT relaunch, and try
        # again next cycle. Acting here is the 2026-08-04 bug.
        echo "  [L=${L} s=${SEED}] ${INSTANCE_ID} state check indeterminate (rc=${RC} state='${STATE}') — no action" >&2
        echo "${INSTANCE_ID}|${DNS}|${L}|${SEED}|${LAUNCHED}"
        continue
      fi

      # CLASS = lost: API-confirmed terminated/shutting-down/aged-out.
      HAS_FINAL="no"
      [ -f "data/outputs/${SWEEP_SUBDIR}/L${L}_s${SEED}_final.npz" ] && HAS_FINAL="yes"
      DUP_COUNT=$(count_running_duplicates "${L}" "${SEED}")

      if [ "$(should_relaunch "${CLASS}" "${HAS_FINAL}" "${DUP_COUNT}")" = "yes" ]; then
        NEW_ID=$(relaunch_on_demand "${L}" "${SEED}")
        if [ -n "${NEW_ID}" ]; then
          RECOVERED=$((RECOVERED+1))
          RECOVERY_SUMMARY="${RECOVERY_SUMMARY} L=${L}s=${SEED}->ondemand(${NEW_ID})"
        fi
        # Drop the dead entry; the new one was appended by relaunch_on_demand.
        continue
      fi

      if [ "${HAS_FINAL}" = "yes" ]; then
        # Normal completion — drop the dead entry silently.
        continue
      fi

      if [ "${DUP_COUNT}" != "0" ] && [ "${DUP_COUNT}" != "999" ]; then
        # An instance with this job's tags is already running (e.g. a
        # relaunch we lost track of). Adopt the newest one into fleet
        # state instead of launching another duplicate.
        ADOPTED=$(aws ec2 describe-instances --region "${AWS_REGION}" \
          --filters "Name=tag:Project,Values=cascade-inversion-cosmology" \
                    "Name=tag:L,Values=${L}" "Name=tag:Seed,Values=${SEED}" \
                    "Name=instance-state-name,Values=pending,running" \
          --query 'sort_by(Reservations[].Instances[], &LaunchTime)[-1].[InstanceId,PublicDnsName,LaunchTime]' \
          --output text 2>/dev/null)
        if [ -n "${ADOPTED}" ]; then
          A_ID=$(echo "${ADOPTED}" | cut -f1)
          A_DNS=$(echo "${ADOPTED}" | cut -f2)
          A_TS=$(echo "${ADOPTED}" | cut -f3)
          echo "  [L=${L} s=${SEED}] adopting already-running ${A_ID} instead of relaunching" >&2
          echo "${A_ID}|${A_DNS}|${L}|${SEED}|${A_TS}"
          continue
        fi
      fi

      # Could not safely relaunch or adopt (e.g. duplicate count query
      # failed). Keep the dead entry so we retry next cycle rather than
      # forgetting the job.
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
