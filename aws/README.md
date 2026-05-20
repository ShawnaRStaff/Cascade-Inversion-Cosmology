# AWS automation for the FSS sweep

Scripts to launch, monitor, sync, and tear down a spot EC2 instance
that runs the finite-size scaling sweep (`run_milestone6_fss_sweep.py`).

## Prerequisites (one-time setup)

1. **AWS CLI installed and configured:**
   ```bash
   aws configure
   ```
   Use the IAM access keys downloaded into `cosmology-sim_accessKeys.csv`.
   Set default region to match the one you picked in the console.

2. **SSH key pair in your AWS region:**
   The key `cosmology-sim-key` should already exist in EC2 → Key Pairs
   (downloaded during initial setup as `cosmology-sim-key.pem`). Make
   sure permissions are tight:
   ```bash
   chmod 600 cosmology-sim-key.pem
   ```

3. **Security group permitting SSH from your IP:**
   ```bash
   aws ec2 create-security-group \
     --group-name cosmology-sim-sg \
     --description "SSH access for cosmology sim spot instances"

   MY_IP=$(curl -s https://api.ipify.org)
   aws ec2 authorize-security-group-ingress \
     --group-name cosmology-sim-sg \
     --protocol tcp \
     --port 22 \
     --cidr ${MY_IP}/32
   ```

## Workflow

### Shakedown run (~$0.50, ~3 hours)

Validates the whole pipeline end-to-end with a tiny L=32 sweep on
a cheap `c6i.large` instance. Do this *once* before paying for the
full sweep.

```bash
./aws/launch_sweep.sh --shakedown
# wait ~3 hours, or check status:
ssh -i cosmology-sim-key.pem ubuntu@<PUBLIC_DNS> 'tail -f logs/sweep_stdout.log'

./aws/sync_results.sh
./aws/teardown.sh
```

### Full sweep (~$95, ~5 days wall)

Launches `c6i.16xlarge` (64 vCPU) and runs the full FSS sweep
across L ∈ {32, 48, 64, 96, 128}.

```bash
./aws/launch_sweep.sh --full
# instance has a 6-day safety shutdown built in
# pull progress whenever you want (idempotent):
./aws/sync_results.sh
# when sweep status reads "done":
./aws/teardown.sh
```

### Emergency teardown

If anything seems wrong:

```bash
./aws/teardown.sh --force
```

If you've lost `.aws_state` somehow, find orphan instances:

```bash
aws ec2 describe-instances \
  --filters Name=tag:Project,Values=cascade-inversion-cosmology \
            Name=instance-state-name,Values=running
```

## Cost safety nets

Three layers of cost protection are built in:

1. **Spot pricing** — 60-70% cheaper than on-demand.
2. **Instance-initiated shutdown** — the instance terminates itself
   when its hard safety timer expires (6 hours for shakedown,
   144 hours / 6 days for full). Hard cap.
3. **`InstanceInterruptionBehavior=terminate`** — if AWS reclaims the
   spot, the instance terminates instead of hibernating (which we
   can't easily resume from). Combined with our checkpoint-resume
   code, this means losing at most ~50k drops of progress per L per
   interruption.

## What does the instance do?

When launched:
1. `bootstrap.sh` runs as root via user-data.
2. Installs Python, git, rsync, awscli.
3. Schedules `shutdown -P +N_minutes` (safety cap).
4. Clones the repo, sets up venv, installs requirements.
5. Launches the FSS sweep in the background.
6. Writes `/home/ubuntu/SWEEP_STATUS` (`starting` → `done`) so
   the local sync script knows when to stop polling.

Sweep outputs land at `/home/ubuntu/repo/data/outputs/fss_sweep_*/`.
Logs at `/home/ubuntu/logs/`.

## Files

| file | purpose |
|---|---|
| `bootstrap.sh` | Runs on the EC2 instance via user-data |
| `launch_sweep.sh` | Local: spins up a spot instance |
| `sync_results.sh` | Local: rsync results back |
| `teardown.sh` | Local: terminate the instance |

## After the sweep

Result `.npz` files land in `data/outputs/fss_sweep_{timestamp}/`.
Each one named `L{L}_s{seed}_final.npz`. Plus a `summary.json` with
per-job peak event size, final p, and wall time.

Next step is analysis: fit peak-event ratio vs L scaling law, write
up M6 results. (To be built; see ROADMAP / task list.)
