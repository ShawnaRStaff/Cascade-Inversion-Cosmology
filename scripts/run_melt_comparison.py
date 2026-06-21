"""Melt-gated rigidity: does fall-in intensification accelerate or dampen the cascade?

Sweeps melt_frac in [0, 0.25, 0.5, 1.0]. For each value:
  - Cold substrate, same drive, same seed.
  - Records when the tip happens and how fast n_hot grows afterward.
  - Plots the n_hot curves aligned to their own tip_step so shapes are comparable.

We don't predict which is faster. The model decides.
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

from void_cascade.cascade_breakdown import BreakdownParams, run_melt_arc  # noqa: E402

MELT_FRACS = [0.0, 0.25, 0.5, 1.0]
COLORS     = ["steelblue", "seagreen", "darkorange", "firebrick"]
L     = 60
STEPS = 3200
SEED  = 0
PARAMS = BreakdownParams(hpc=0.1, drive_amount=1.0, drive_sites=4)

print(f"=== Melt-gated intensification sweep ===")
print(f"L={L}  steps={STEPS}  seed={SEED}")
print(f"{'melt_frac':>10}  {'tip_step':>10}  {'max_temp':>9}  {'resid':>10}")

results = {}
for mf in MELT_FRACS:
    r = run_melt_arc(L=L, steps=STEPS, params=PARAMS, seed=SEED, melt_frac=mf,
                     sample_every=20)
    results[mf] = r
    print(f"{mf:>10.2f}  {str(r['tip_step']):>10}  {r['max_temp']:>9.2f}  "
          f"{r['energy_residual']:>10.2e}")

# --- plot 1: absolute n_hot over time ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
for mf, col in zip(MELT_FRACS, COLORS):
    r = results[mf]
    t = np.array(r["t_axis"])
    n = np.array(r["n_hot_trace"])
    ax.plot(t, n, color=col, lw=1.5, label=f"melt={mf}")
    if r["tip_step"] is not None:
        ax.axvline(r["tip_step"], color=col, lw=0.8, ls=":", alpha=0.5)
ax.set_xlabel("step"); ax.set_ylabel("n_hot (cells above ignition temp)")
ax.set_title("n_hot over time (dotted = tip)")
ax.legend(fontsize=9)

# --- plot 2: n_hot relative to each run's own tip_step ---
ax = axes[1]
WINDOW = 600   # steps after tip to show
for mf, col in zip(MELT_FRACS, COLORS):
    r = results[mf]
    if r["tip_step"] is None:
        continue
    t  = np.array(r["t_axis"])
    n  = np.array(r["n_hot_trace"])
    dt = t - r["tip_step"]
    mask = (dt >= 0) & (dt <= WINDOW)
    if mask.sum() < 2:
        continue
    ax.plot(dt[mask], n[mask], color=col, lw=1.5, label=f"melt={mf}")
ax.set_xlabel("steps after tip"); ax.set_ylabel("n_hot")
ax.set_title(f"Cascade shape after tip  (first {WINDOW} steps)")
ax.legend(fontsize=9)

fig.suptitle("Melt-gated rigidity: cold cells absorb hot neighbours' KE as load")
fig.tight_layout()

out = REPO / "data" / "outputs" / "melt_comparison"
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "melt_comparison.png", dpi=150)
print(f"\nPlot saved: {out}/melt_comparison.png")

# --- plain-language summary ---
print("\n=== SUMMARY ===")
tip_steps = {mf: results[mf]["tip_step"] for mf in MELT_FRACS}
tips_valid = {mf: t for mf, t in tip_steps.items() if t is not None}
if len(tips_valid) > 1:
    earliest = min(tips_valid, key=tips_valid.get)
    latest   = max(tips_valid, key=tips_valid.get)
    print(f"  Earliest tip: melt_frac={earliest} at step {tips_valid[earliest]}")
    print(f"  Latest tip:   melt_frac={latest}  at step {tips_valid[latest]}")
    delta = tips_valid[latest] - tips_valid[earliest]
    print(f"  Spread: {delta} steps")
    if delta < 50:
        print("  -> Melt fraction has little effect on TIP TIMING")
    else:
        print("  -> Melt fraction shifts tip timing measurably")
else:
    print("  Not all runs tipped.")
