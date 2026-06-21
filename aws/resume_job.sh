#!/bin/bash
# resume_job.sh - manually resume a specific interrupted (L, seed) job.
#
# Uploads the local checkpoint to S3, generates a 24-hour presigned URL,
# and launches an on-demand m7i-flex.large instance that will download the
# checkpoint before starting the sweep driver. The new instance is appended
# to .aws_fleet_state so monitor_and_recover.sh continues tracking it.
#
# Usage:
#   ./aws/resume_job.sh <L> <SEED> <SWEEP_SUBDIR>
#
# Example (resume the interrupted s22801 job):
#   ./aws/resume_job.sh 128 22801 fss_sweep_20260521_031056
#
# Note: the monitor uses .aws_fleet_sweep_dir (single value) to scope its
# auto-recovery checks. resume_job.sh writes that file to point at the
# given SWEEP_SUBDIR. If you are already tracking a different batch with
# the monitor, run this BEFORE launching a new batch — or accept that the
# resumed job runs on-demand without auto-recovery (safe, no spot risk).

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

L="${1:?Usage: resume_job.sh <L> <SEED> <SWEEP_SUBDIR>}"
SEED="${2:?Usage: resume_job.sh <L> <SEED> <SWEEP_SUBDIR>}"
SWEEP_SUBDIR="${3:?Usage: resume_job.sh <L> <SEED> <SWEEP_SUBDIR>}"

AWS_REGION="${AWS_REGION:-us-west-2}"
REPO_URL="https://github.com/ShawnaRStaff/Cascade-Inversion-Cosmology.git"
BRANCH="main"
SAFETY_HOURS=120
S3_BUCKET="cascade-cosmo-ckpts-099623380651"

CKPT_LOCAL="data/outputs/${SWEEP_SUBDIR}/L${L}_s${SEED}_ckpt.npz"
if [ ! -f "${CKPT_LOCAL}" ]; then
  echo "ERROR: no checkpoint at ${CKPT_LOCAL}"
  exit 1
fi

echo "Uploading checkpoint to S3..."
S3_KEY="checkpoints/${SWEEP_SUBDIR}/L${L}_s${SEED}_ckpt.npz"
aws s3 mb "s3://${S3_BUCKET}" --region "${AWS_REGION}" 2>/dev/null || true
aws s3 cp "${CKPT_LOCAL}" "s3://${S3_BUCKET}/${S3_KEY}" --region "${AWS_REGION}"

echo "Generating 24-hour presigned URL..."
PRESIGNED=$(aws s3 presign "s3://${S3_BUCKET}/${S3_KEY}" \
  --region "${AWS_REGION}" --expires-in 86400)
CHECKPOINT_URL_B64=$(echo "${PRESIGNED}" | base64 -w 0)

echo "Resolving latest Ubuntu 24.04 AMI..."
AMI_ID=$(aws ec2 describe-images \
  --region "${AWS_REGION}" \
  --owners 099720109477 \
  --filters \
    "Name=name,Values=ubuntu/images/hvm-ssd-gp3/ubuntu-noble-24.04-amd64-server-*" \
    "Name=state,Values=available" \
  --query 'Images | sort_by(@, &CreationDate)[-1].ImageId' \
  --output text)
echo "AMI: ${AMI_ID}"

JOB_FLAGS="--single-job ${L}:${SEED} --sweep-dir data/outputs/${SWEEP_SUBDIR}"
USER_DATA=$(mktemp)
sed \
  -e "s|@REPO_URL@|${REPO_URL}|g" \
  -e "s|@BRANCH@|${BRANCH}|g" \
  -e "s|@SWEEP_FLAGS@|${JOB_FLAGS}|g" \
  -e "s|@SAFETY_HOURS@|${SAFETY_HOURS}|g" \
  -e "s|@CHECKPOINT_URL_B64@|${CHECKPOINT_URL_B64}|g" \
  -e "s|@CHECKPOINT_DEST@|${CKPT_LOCAL}|g" \
  aws/bootstrap.sh > "${USER_DATA}"

echo "Launching on-demand m7i-flex.large for L=${L} seed=${SEED}..."
INSTANCE_ID=$(aws ec2 run-instances \
  --region "${AWS_REGION}" \
  --image-id "${AMI_ID}" \
  --instance-type m7i-flex.large \
  --key-name cosmology-sim-key \
  --security-groups cosmology-sim-sg \
  --user-data "file://${USER_DATA}" \
  --instance-initiated-shutdown-behavior terminate \
  --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=cosmology-sim-resume-L${L}-s${SEED}},{Key=Project,Value=cascade-inversion-cosmology},{Key=L,Value=${L}},{Key=Seed,Value=${SEED}},{Key=Pricing,Value=on-demand-manual-resume}]" \
  --query 'Instances[0].InstanceId' \
  --output text)
rm -f "${USER_DATA}"

echo "Waiting for ${INSTANCE_ID} to reach running state..."
aws ec2 wait instance-running --region "${AWS_REGION}" --instance-ids "${INSTANCE_ID}"
DNS=$(aws ec2 describe-instances \
  --region "${AWS_REGION}" \
  --instance-ids "${INSTANCE_ID}" \
  --query 'Reservations[0].Instances[0].PublicDnsName' \
  --output text)

echo "${INSTANCE_ID}|${DNS}|${L}|${SEED}|$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> .aws_fleet_state
echo "${SWEEP_SUBDIR}" > .aws_fleet_sweep_dir

echo
echo "Resumed."
echo "  Instance:  ${INSTANCE_ID} (${DNS})"
echo "  Job:       L=${L} seed=${SEED}"
echo "  Sweep dir: data/outputs/${SWEEP_SUBDIR}"
echo "  Checkpoint: ${CKPT_LOCAL}"
echo
echo "Monitor: ./aws/monitor_and_recover.sh"
echo "Sync:    ./aws/sync_fleet.sh"
