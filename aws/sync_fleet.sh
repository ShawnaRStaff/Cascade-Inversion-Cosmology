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

  # CRITICAL: always attempt rsync, even if the status check failed.
  # SSH can transiently fail (network blip, AWS API throttling) and we
  # MUST not skip data sync on those. rsync uses its own SSH connection
  # and may succeed even when the status query failed. This is the bug
  # that lost the L=96 NEW-metric data overnight: status check returned
  # "unreachable" intermittently during the 3-hour grace window, sync
  # was skipped for that cycle, and by the time grace expired no sync
  # had successfully pulled the data.
  #
  # We track unreachable only for reporting; we don't gate rsync on it.
  if [ "${STATUS}" = "unreachable" ]; then
    UNREACHABLE_COUNT=$((UNREACHABLE_COUNT + 1))
  elif [ "${STATUS}" = "done" ]; then
    DONE_COUNT=$((DONE_COUNT + 1))
  else
    RUNNING_COUNT=$((RUNNING_COUNT + 1))
  fi

  # ALWAYS try rsync. Data preservation takes priority over status
  # accounting. ConnectTimeout=20 a bit longer than status check.
  RSYNC_OK=0
  if rsync -az --quiet --timeout=60 \
    -e "ssh -o StrictHostKeyChecking=accept-new -o ConnectTimeout=20 -i ${KEY}" \
    "ubuntu@${DNS}:/home/ubuntu/repo/data/outputs/fss_sweep_*" \
    data/outputs/ 2>/dev/null; then
    RSYNC_OK=1
  fi

  # Also pull job log into aws_logs.
  rsync -az --quiet --timeout=30 \
    -e "ssh -o StrictHostKeyChecking=accept-new -o ConnectTimeout=20 -i ${KEY}" \
    "ubuntu@${DNS}:/home/ubuntu/logs/" \
    "data/aws_logs/L${L}_s${SEED}/" 2>/dev/null || true

  # If rsync failed AND status was unreachable, the instance is truly
  # unreachable. If rsync succeeded, the instance was reachable enough
  # to pull data even though status check may have flapped.
  if [ "${STATUS}" = "unreachable" ] && [ "${RSYNC_OK}" -eq 0 ]; then
    echo "(unreachable, rsync also failed)"
    continue
  fi

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
