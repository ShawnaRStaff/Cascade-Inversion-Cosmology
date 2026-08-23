"""Tests for the fleet monitor's loss-detection logic (aws/lib_recovery.sh).

The 2026-08-04 incident: a transient AWS API failure left the state
string empty, which the monitor read as "instance record aged out" and
answered with a relaunch — four duplicate instances in one night, each
overwriting checkpoints and burning credits. The fixed logic must:

  1. treat API *failure* (nonzero exit) as indeterminate — never lost;
  2. treat only API-confirmed terminated/shutting-down/aged-out as lost;
  3. refuse to relaunch when a duplicate instance for the same job is
     already running, or when the job's final.npz already exists.

The decision functions are pure bash (data in via args, verdict out via
stdout) so they can be tested here without any AWS access.
"""

import subprocess
from pathlib import Path

LIB = Path(__file__).resolve().parents[1] / "aws" / "lib_recovery.sh"


def call(fn: str, *args: str) -> str:
    out = subprocess.run(
        ["bash", "-c", f'source "{LIB}"; {fn} ' + " ".join(f'"{a}"' for a in args)],
        capture_output=True,
        text=True,
        check=True,
    )
    return out.stdout.strip()


# --- classify_instance_state(api_exit_code, state_string) ---------------

def test_api_failure_is_indeterminate_regardless_of_state_string():
    assert call("classify_instance_state", "1", "") == "indeterminate"
    assert call("classify_instance_state", "255", "terminated") == "indeterminate"


def test_confirmed_terminated_or_shutting_down_is_lost():
    assert call("classify_instance_state", "0", "terminated") == "lost"
    assert call("classify_instance_state", "0", "shutting-down") == "lost"


def test_aged_out_record_is_lost_only_on_api_success():
    # aws CLI prints the literal string "None" for a missing record.
    assert call("classify_instance_state", "0", "None") == "lost"


def test_empty_state_on_api_success_is_indeterminate():
    # Empty output with exit 0 is not a documented AWS response; the
    # 2026-08-04 incident came from trusting empty strings. Never lost.
    assert call("classify_instance_state", "0", "") == "indeterminate"


def test_running_and_pending_are_alive():
    assert call("classify_instance_state", "0", "running") == "alive"
    assert call("classify_instance_state", "0", "pending") == "alive"


def test_stopped_and_stopping_are_indeterminate():
    # A stopped instance still has its EBS volume: relaunching would risk
    # duplicates; leave the decision to a human.
    assert call("classify_instance_state", "0", "stopped") == "indeterminate"
    assert call("classify_instance_state", "0", "stopping") == "indeterminate"


# --- should_relaunch(classification, has_final, duplicate_count) --------

def test_relaunch_only_when_lost_no_final_no_duplicates():
    assert call("should_relaunch", "lost", "no", "0") == "yes"


def test_no_relaunch_when_indeterminate_or_alive():
    assert call("should_relaunch", "indeterminate", "no", "0") == "no"
    assert call("should_relaunch", "alive", "no", "0") == "no"


def test_no_relaunch_when_job_already_finished():
    assert call("should_relaunch", "lost", "yes", "0") == "no"


def test_no_relaunch_when_duplicate_instances_exist():
    assert call("should_relaunch", "lost", "no", "1") == "no"
    assert call("should_relaunch", "lost", "no", "999") == "no"


# --- relaunch_preseed_verdict(ckpt_exists, presigned_url) ---------------
# Guards the account-agnostic S3 bucket path: if a local checkpoint
# exists but the upload/presign chain failed (e.g. the account-id lookup
# returned empty and the bucket name collapsed to a garbage prefix), a
# relaunch would silently start FROM SCRATCH. That must abort instead.

def test_preseed_fresh_start_allowed_only_without_checkpoint():
    assert call("relaunch_preseed_verdict", "no", "") == "proceed-fresh"


def test_preseed_proceeds_with_checkpoint_and_url():
    assert call("relaunch_preseed_verdict", "yes", "https://bucket/key?sig=x") == "proceed"


def test_preseed_aborts_when_checkpoint_exists_but_url_missing():
    assert call("relaunch_preseed_verdict", "yes", "") == "abort"


# --- valid_account_id(account_id) ---------------------------------------

def test_valid_account_id_accepts_12_digits_only():
    assert call("valid_account_id", "099623380651") == "yes"
    assert call("valid_account_id", "") == "no"
    assert call("valid_account_id", "None") == "no"
    assert call("valid_account_id", "12345") == "no"
    assert call("valid_account_id", "abc123456789") == "no"
