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
apt-get install -y python3-venv python3-pip git rsync awscli

# Safety: shut down the instance after SAFETY_HOURS, no matter what.
# This is the hard cost cap. If the sweep hangs, the bill stops here.
shutdown -P +$(( SAFETY_HOURS * 60 )) &

# Clone, set up, run.
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

# Run the sweep in the background. Output goes to /home/ubuntu/logs/.
nohup .venv/bin/python -u scripts/run_milestone6_fss_sweep.py \
    ${SWEEP_FLAGS} \
    >/home/ubuntu/logs/sweep_stdout.log \
    2>/home/ubuntu/logs/sweep_stderr.log &
SWEEP_PID=\$!
echo "${SWEEP_PID}" > /home/ubuntu/SWEEP_PID

# Wait for sweep, then mark done.
( wait \${SWEEP_PID}; echo "done" > /home/ubuntu/SWEEP_STATUS ) &
UBUNTU_EOF

echo "bootstrap.sh finished"
