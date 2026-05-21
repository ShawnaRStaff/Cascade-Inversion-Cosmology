#!/bin/bash
# sync_fleet.sh - pull results from every instance in the fleet.
#
# Idempotent: safe to run multiple times while jobs are still in
# progress; pulls whatever's there each time.
#
# Reads .aws_fleet_state, iterates over each instance, and rsyncs
# /home/ubuntu/repo/data/outputs/fss_sweep_*/ back to local. All
# instances share a sweep subdir name (the timestamp from launch),
# so files merge cleanly into data/outputs/fss_sweep_<stamp>/.

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

if [ ! -f .aws_fleet_state ]; then
  echo "ERROR: .aws_fleet_state not found. Did you run launch_fleet.sh?"
  exit 1
fi

KEY="${KEY:-cosmology-sim-key.pem}"
mkdir -p data/outputs data/aws_logs

DONE_COUNT=0
RUNNING_COUNT=0
UNREACHABLE_COUNT=0
TOTAL_COUNT=$(wc -l < .aws_fleet_state)

while IFS='|' read -r INSTANCE_ID DNS L SEED LAUNCHED; do
  if [ -z "${DNS}" ] || [ "${DNS}" = "unknown" ] || [ "${DNS}" = "pending" ]; then
    printf "  [L=%-4d s=%d] %s no DNS yet\n" "${L}" "${SEED}" "${INSTANCE_ID}"
    UNREACHABLE_COUNT=$((UNREACHABLE_COUNT + 1))
    continue
  fi

  # ssh -n redirects stdin from /dev/null. Without it, ssh consumes the
  # while-loop's stdin (the fleet state file) and the loop exits after
  # the first successful ssh.
  STATUS=$(ssh -n -o StrictHostKeyChecking=accept-new -o ConnectTimeout=10 -i "${KEY}" "ubuntu@${DNS}" \
    'cat /home/ubuntu/SWEEP_STATUS 2>/dev/null || echo unknown' 2>/dev/null \
    || echo "unreachable")

  printf "  [L=%-4d s=%d] %s status=%s ... " "${L}" "${SEED}" "${INSTANCE_ID}" "${STATUS}"

  if [ "${STATUS}" = "unreachable" ]; then
    echo "(unreachable, may have terminated)"
    UNREACHABLE_COUNT=$((UNREACHABLE_COUNT + 1))
    continue
  fi

  if [ "${STATUS}" = "done" ]; then
    DONE_COUNT=$((DONE_COUNT + 1))
  else
    RUNNING_COUNT=$((RUNNING_COUNT + 1))
  fi

  # Pull whatever's there. Single-job mode writes into the
  # data/outputs/fss_sweep_*/ subdir.
  rsync -az --quiet \
    -e "ssh -o StrictHostKeyChecking=accept-new -o ConnectTimeout=10 -i ${KEY}" \
    "ubuntu@${DNS}:/home/ubuntu/repo/data/outputs/fss_sweep_*" \
    data/outputs/ 2>/dev/null || true

  # Also pull job log into aws_logs.
  rsync -az --quiet \
    -e "ssh -o StrictHostKeyChecking=accept-new -o ConnectTimeout=10 -i ${KEY}" \
    "ubuntu@${DNS}:/home/ubuntu/logs/" \
    "data/aws_logs/L${L}_s${SEED}/" 2>/dev/null || true

  # Auto-terminate this instance if its job is done AND the final.npz
  # is on local disk for THIS fleet's sweep dir specifically (not just
  # anywhere — that would let an old fleet's matching filename
  # incorrectly trigger termination of a new fleet instance).
  #
  # The current fleet's sweep subdir is recorded in .aws_fleet_sweep_dir
  # at launch time. If absent (old fleets), fall back to a global find.
  if [ -f .aws_fleet_sweep_dir ]; then
    SWEEP_SUBDIR=$(cat .aws_fleet_sweep_dir)
    FINAL_NPZ="data/outputs/${SWEEP_SUBDIR}/L${L}_s${SEED}_final.npz"
  else
    FINAL_NPZ=$(find data/outputs -name "L${L}_s${SEED}_final.npz" 2>/dev/null | head -1)
  fi
  if [ -n "${FINAL_NPZ}" ] && [ -f "${FINAL_NPZ}" ]; then
    echo "  -> final.npz present locally; terminating idle instance"
    aws ec2 terminate-instances --region "${AWS_REGION:-us-west-2}" --instance-ids "${INSTANCE_ID}" --output text >/dev/null 2>&1 || true
    continue
  fi

  echo "synced"
done < .aws_fleet_state

echo
echo "Fleet status:"
echo "  done:        ${DONE_COUNT}/${TOTAL_COUNT}"
echo "  running:     ${RUNNING_COUNT}/${TOTAL_COUNT}"
echo "  unreachable: ${UNREACHABLE_COUNT}/${TOTAL_COUNT}"
echo
if [ "${DONE_COUNT}" -eq "${TOTAL_COUNT}" ]; then
  echo "All jobs done. Safe to run ./aws/teardown_fleet.sh"
fi
