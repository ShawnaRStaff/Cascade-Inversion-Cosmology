#!/bin/bash
# lib_recovery.sh - pure decision logic for fleet loss detection.
#
# Extracted from monitor_and_recover.sh after the 2026-08-04 incident:
# a transient AWS API failure produced an empty state string, which the
# old inline check read as "record aged out -> relaunch". Four duplicate
# instances were launched in one night; their syncs overwrote good
# checkpoints locally and in S3, and the extra spend drained the
# account's remaining credits days early.
#
# These functions are pure (arguments in, verdict on stdout, no side
# effects, no AWS calls) and are unit-tested by
# tests/test_monitor_recovery_logic.py. Callers do the I/O.

# classify_instance_state API_EXIT_CODE STATE_STRING
#
# Echoes one of: alive | lost | indeterminate
#
#   - Any nonzero API exit code is "indeterminate": we KNOW NOTHING and
#     must not act. This is the 2026-08-04 fix.
#   - "lost" requires an API-confirmed answer: state terminated or
#     shutting-down, or the literal "None" the CLI prints for a record
#     that aged out of the API after termination.
#   - running/pending are "alive".
#   - stopped/stopping (EBS volume still exists) and anything
#     unrecognized — including an empty string on API success — are
#     "indeterminate": log and let a human decide.
classify_instance_state() {
  local rc="$1"
  local state="$2"
  if [ "${rc}" != "0" ]; then
    echo "indeterminate"
    return
  fi
  case "${state}" in
    terminated|shutting-down|None) echo "lost" ;;
    running|pending)               echo "alive" ;;
    *)                             echo "indeterminate" ;;
  esac
}

# should_relaunch CLASSIFICATION HAS_FINAL DUPLICATE_COUNT
#
# Echoes yes | no.
#
# Relaunch only when ALL of:
#   - the instance is confirmed lost (see classify_instance_state);
#   - the job's final.npz is NOT already on local disk (has_final=no) —
#     a finished job needs no recovery;
#   - zero running/pending instances already carry this job's tags.
#     Callers that fail to count duplicates must pass a nonzero
#     sentinel (e.g. 999) so the answer is "no".
should_relaunch() {
  local classification="$1"
  local has_final="$2"
  local dup_count="$3"
  if [ "${classification}" = "lost" ] \
     && [ "${has_final}" = "no" ] \
     && [ "${dup_count}" = "0" ]; then
    echo "yes"
  else
    echo "no"
  fi
}

# relaunch_preseed_verdict CKPT_EXISTS PRESIGNED_URL
#
# Echoes proceed | proceed-fresh | abort.
#
# A relaunch may only start a fresh (from-scratch) instance when there
# is genuinely no local checkpoint to resume from. If a checkpoint
# exists but the upload/presign chain produced no URL (account-id
# lookup failed, bucket name collapsed, S3 error), launching anyway
# would silently discard all progress — abort and retry next cycle.
relaunch_preseed_verdict() {
  local ckpt_exists="$1"
  local presigned_url="$2"
  if [ "${ckpt_exists}" = "no" ]; then
    echo "proceed-fresh"
  elif [ -n "${presigned_url}" ]; then
    echo "proceed"
  else
    echo "abort"
  fi
}

# valid_account_id ACCOUNT_ID — "yes" iff exactly 12 digits (the AWS
# account-id format). Guards bucket names built as
# cascade-cosmo-ckpts-${ACCOUNT_ID} from collapsing when the STS
# lookup fails or prints "None".
valid_account_id() {
  if [[ "$1" =~ ^[0-9]{12}$ ]]; then
    echo "yes"
  else
    echo "no"
  fi
}
