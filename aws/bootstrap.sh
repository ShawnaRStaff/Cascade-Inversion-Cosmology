#!/bin/bash
# bootstrap.sh - runs on a freshly-launched EC2 instance.
#
# Brought up via the `--user-data` field of `aws ec2 run-instances`,
# this script:
#   1. Installs Python + git + rsync
#   2. Clones the project repo
#   3. Sets up a venv and installs requirements
#   4. Schedules a safety shutdown so the instance can't run forever
#      if our sweep code hangs (caps cost)
#   5. Launches the FSS sweep in the background, logging to disk
#   6. Writes a status flag (~ubuntu/SWEEP_STATUS) that the rsync
#      step on the local side reads to know when to pull results
#
# User-data scripts run as root with `/` as cwd. We do most work as
# the `ubuntu` user.

set -euxo pipefail

# Bootstrap parameters (sed'd in by launch.sh before upload)
REPO_URL="@REPO_URL@"
BRANCH="@BRANCH@"
SWEEP_FLAGS="@SWEEP_FLAGS@"
SAFETY_HOURS="@SAFETY_HOURS@"

# Update + install minimal deps
export DEBIAN_FRONTEND=noninteractive
apt-get update
# Ubuntu 24.04 dropped awscli from main repos; we don't need it on
# the instance anyway (rsync handles result sync from the local side).
apt-get install -y python3-venv python3-pip git rsync

# Safety: shut down the instance after SAFETY_HOURS, no matter what.
# This is the hard cost cap. If the sweep hangs, the bill stops here.
shutdown -P +$(( SAFETY_HOURS * 60 )) &

# Clone, set up, run.
# Note: heredoc uses ${VAR} for outer-shell substitution (BRANCH,
# REPO_URL, SWEEP_FLAGS were sed'd in above and are normal shell
# variables). Variables we want the INNER shell (sudo -u ubuntu bash)
# to evaluate at runtime must use \${VAR} to escape outer expansion.
sudo -u ubuntu bash <<UBUNTU_EOF
set -euxo pipefail
cd /home/ubuntu
git clone --depth 1 --branch "${BRANCH}" "${REPO_URL}" repo
cd repo
python3 -m venv .venv
.venv/bin/pip install --upgrade pip
.venv/bin/pip install -r requirements.txt
mkdir -p /home/ubuntu/logs

echo "starting" > /home/ubuntu/SWEEP_STATUS

# Write a sequential runner script: run sweep, then mark done, then
# schedule self-termination. Putting these in one chain (vs. trying
# to 'wait' on a sibling process from a backgrounded subshell, which
# silently fails because wait only works for direct children) is the
# only reliable way to do completion-detection here.
cat > /home/ubuntu/sweep_runner.sh <<RUNNER_EOF
#!/bin/bash
cd /home/ubuntu/repo
.venv/bin/python -u scripts/run_milestone6_fss_sweep.py ${SWEEP_FLAGS} \
    >/home/ubuntu/logs/sweep_stdout.log \
    2>/home/ubuntu/logs/sweep_stderr.log
echo "done" > /home/ubuntu/SWEEP_STATUS
# Give 5 min for any final sync, then power off. The instance has
# --instance-initiated-shutdown-behavior=terminate so this terminates.
sudo shutdown -P +5
RUNNER_EOF
chmod +x /home/ubuntu/sweep_runner.sh

nohup /home/ubuntu/sweep_runner.sh > /home/ubuntu/logs/runner.log 2>&1 &
echo \$! > /home/ubuntu/SWEEP_PID
UBUNTU_EOF

echo "bootstrap.sh finished"
