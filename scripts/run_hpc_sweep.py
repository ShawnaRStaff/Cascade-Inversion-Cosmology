"""hpc parameter sweep: what does the stress->heat fraction actually control?

hpc (heat-per-crack) is the free knob that converts breakdown stress to local heat.
The tip is said to be "robust" across it -- but does it tip at all values?
Is there a minimum hpc below which the cold regime locks in forever?
And does hpc change HOW violent the ignition is (peak temperature, speed)?

Grid: hpc in [0.01, 0.05, 0.1, 0.2, 0.5, 1.0] x drive_amount in [0.5, 1.0, 2.0].
For each: run 3000 steps on L=80, report tip/no-tip, tip_step, max_temp.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from void_cascade.cascade_breakdown import BreakdownParams, run_buildup_tip  # noqa: E402

L     = 60    # moderate size: fast enough for a grid search
STEPS = 3000
SEED  = 0

HPC_VALUES    = [0.01, 0.05, 0.10, 0.20, 0.50, 1.00]
DRIVE_VALUES  = [0.5,  1.0,  2.0]

print(f"hpc sweep  L={L}  steps={STEPS}")
print(f"{'hpc':>6}  {'drive':>6}  {'tip':>8}  {'maxT':>7}  {'resid':>10}")
print("-" * 50)

results = {}
for drive in DRIVE_VALUES:
    for hpc in HPC_VALUES:
        p = BreakdownParams(hpc=hpc, drive_amount=drive, drive_sites=4)
        r = run_buildup_tip(L=L, steps=STEPS, params=p, seed=SEED)
        tip  = r["tip_step"]
        maxT = r["max_temp"]
        res  = r["energy_residual"]
        label = f"tip={tip:5d}" if tip is not None else "no-tip"
        print(f"{hpc:6.2f}  {drive:6.2f}  {label:>8}  {maxT:7.2f}  {res:10.2e}")
        results[(hpc, drive)] = {"tip_step": tip, "max_temp": maxT, "residual": res}

# --- phase diagram plot ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

hpc_arr   = np.array(HPC_VALUES)
drive_arr = np.array(DRIVE_VALUES)

for ax_idx, drive in enumerate(DRIVE_VALUES):
    tips   = [results[(h, drive)]["tip_step"] for h in HPC_VALUES]
    maxTs  = [results[(h, drive)]["max_temp"]  for h in HPC_VALUES]
    tipped = [t is not None for t in tips]
    tip_steps = [t if t is not None else STEPS for t in tips]

    axes[0].plot(hpc_arr, tip_steps, "o-", label=f"drive={drive}")
    axes[1].plot(hpc_arr, maxTs,     "o-", label=f"drive={drive}")

axes[0].set_xlabel("hpc (stress -> heat fraction)")
axes[0].set_ylabel(f"tip step (→{STEPS} = no tip)")
axes[0].set_title("When does ignition happen?")
axes[0].legend(); axes[0].set_xscale("log")
axes[0].axhline(STEPS, color="k", lw=0.8, ls="--", alpha=0.4)

axes[1].set_xlabel("hpc (stress -> heat fraction)")
axes[1].set_ylabel("peak temperature")
axes[1].set_title("How hot does it get?")
axes[1].axhline(2.5, color="k", lw=0.8, ls="--", alpha=0.4, label="e_ign")
axes[1].legend(); axes[1].set_xscale("log")

fig.suptitle(f"hpc sweep  L={L}  steps={STEPS}  drive_sites=4")
fig.tight_layout()

out = REPO / "data" / "outputs" / "hpc_sweep"
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "hpc_phase_diagram.png", dpi=150)
print(f"\nPlot saved: {out}/hpc_phase_diagram.png")
