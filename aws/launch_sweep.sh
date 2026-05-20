#!/bin/bash
# launch_sweep.sh - launch a spot EC2 instance and start the FSS sweep.
#
# This is the local-side script. It:
#   1. Renders the bootstrap.sh template with the chosen parameters
#   2. Spins up a spot EC2 instance with that user-data
#   3. Waits for the instance to come up, prints its public DNS
#   4. Saves the instance ID + DNS to .aws_state for sync/teardown
#
# Prereqs:
#   - aws cli installed and `aws configure` done with IAM credentials
#   - SSH key pair `cosmology-sim-key` exists in your AWS region
#   - A security group exists with port 22 inbound from your IP
#
# Usage:
#   ./aws/launch_sweep.sh [--shakedown]
#
# --shakedown launches a tiny c6i.large instance and runs only the
# L=32 sweep (~3 hours, ~$0.50). Use this once before paying for the
# full sweep, to confirm the pipeline works end-to-end on AWS.

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

# ---- Configuration (edit per session) ----

AWS_REGION="${AWS_REGION:-us-west-2}"
KEY_NAME="${KEY_NAME:-cosmology-sim-key}"
SECURITY_GROUP="${SECURITY_GROUP:-cosmology-sim-sg}"
REPO_URL="${REPO_URL:-https://github.com/ShawnaRStaff/Cascade-Inversion-Cosmology.git}"
BRANCH="${BRANCH:-main}"

# Mode selection
MODE="${1:-full}"
case "${MODE}" in
  --shakedown|shakedown)
    INSTANCE_TYPE="${INSTANCE_TYPE:-c6i.large}"
    SWEEP_FLAGS="--l-list 32"
    SAFETY_HOURS=4
    ;;
  --full|full|"")
    INSTANCE_TYPE="${INSTANCE_TYPE:-c6i.16xlarge}"
    SWEEP_FLAGS=""
    SAFETY_HOURS=144   # 6 days hard cap
    ;;
  *)
    echo "Usage: $0 [--shakedown|--full]"
    exit 1
    ;;
esac

echo "Mode:          ${MODE}"
echo "Region:        ${AWS_REGION}"
echo "Instance type: ${INSTANCE_TYPE}"
echo "Safety hours:  ${SAFETY_HOURS}"
echo "Sweep flags:   '${SWEEP_FLAGS}'"
echo

# ---- Latest Ubuntu 24.04 AMI for the chosen region ----
AMI_ID=$(aws ec2 describe-images \
  --region "${AWS_REGION}" \
  --owners 099720109477 \
  --filters \
    "Name=name,Values=ubuntu/images/hvm-ssd-gp3/ubuntu-noble-24.04-amd64-server-*" \
    "Name=state,Values=available" \
  --query 'Images | sort_by(@, &CreationDate)[-1].ImageId' \
  --output text)

if [ -z "${AMI_ID}" ] || [ "${AMI_ID}" = "None" ]; then
  echo "ERROR: could not resolve Ubuntu 24.04 AMI in ${AWS_REGION}"
  exit 1
fi
echo "AMI: ${AMI_ID}"

# ---- Render bootstrap.sh template ----
USER_DATA=$(mktemp)
trap 'rm -f "${USER_DATA}"' EXIT
sed \
  -e "s|@REPO_URL@|${REPO_URL}|g" \
  -e "s|@BRANCH@|${BRANCH}|g" \
  -e "s|@SWEEP_FLAGS@|${SWEEP_FLAGS}|g" \
  -e "s|@SAFETY_HOURS@|${SAFETY_HOURS}|g" \
  aws/bootstrap.sh > "${USER_DATA}"
echo "rendered user-data to ${USER_DATA}"

# ---- Spin up the spot instance ----
echo "requesting spot ${INSTANCE_TYPE}..."
RUN_OUTPUT=$(aws ec2 run-instances \
  --region "${AWS_REGION}" \
  --image-id "${AMI_ID}" \
  --instance-type "${INSTANCE_TYPE}" \
  --key-name "${KEY_NAME}" \
  --security-groups "${SECURITY_GROUP}" \
  --user-data "file://${USER_DATA}" \
  --instance-initiated-shutdown-behavior terminate \
  --instance-market-options 'MarketType=spot,SpotOptions={SpotInstanceType=one-time,InstanceInterruptionBehavior=terminate}' \
  --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=cosmology-sim-fss-sweep},{Key=Project,Value=cascade-inversion-cosmology}]' \
  --query 'Instances[0].InstanceId' \
  --output text)

INSTANCE_ID="${RUN_OUTPUT}"
echo "instance id: ${INSTANCE_ID}"

echo "waiting for instance to be running..."
aws ec2 wait instance-running --region "${AWS_REGION}" --instance-ids "${INSTANCE_ID}"

PUBLIC_DNS=$(aws ec2 describe-instances \
  --region "${AWS_REGION}" \
  --instance-ids "${INSTANCE_ID}" \
  --query 'Reservations[0].Instances[0].PublicDnsName' \
  --output text)

# Persist state so sync/teardown scripts can find this instance later.
cat > .aws_state <<EOF
INSTANCE_ID=${INSTANCE_ID}
PUBLIC_DNS=${PUBLIC_DNS}
AWS_REGION=${AWS_REGION}
MODE=${MODE}
LAUNCHED_AT=$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF

echo
echo "===================================="
echo "INSTANCE LAUNCHED"
echo "  ID:         ${INSTANCE_ID}"
echo "  Public DNS: ${PUBLIC_DNS}"
echo "  Region:     ${AWS_REGION}"
echo "  Safety cap: ${SAFETY_HOURS} hours"
echo
echo "State saved to .aws_state"
echo
echo "Next steps:"
echo "  - Wait for bootstrap to finish (~3-5 min)"
echo "  - Check progress: ssh -i cosmology-sim-key.pem ubuntu@${PUBLIC_DNS} 'tail -f logs/sweep_stdout.log'"
echo "  - Pull results:   ./aws/sync_results.sh"
echo "  - Tear down:      ./aws/teardown.sh"
echo "===================================="
