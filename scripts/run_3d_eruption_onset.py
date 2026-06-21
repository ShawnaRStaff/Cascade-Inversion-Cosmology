"""3D: two questions in one run.

Q1 (eruption front): At default params 3D tips marginally (max_temp=2.47).
   Does a stronger drive or higher hpc produce a SUSTAINED eruption?
   Sweeps hpc=[0.1,0.2,0.3] x drive_amount=[1.0,2.0] at L=20, same per-cell
   drive rate. Reports tip_step and max_temp for each combination.

Q2 (onset ordering): Does the 2D sequence (corr_length -> load_std -> ignition)
   hold in 3D? Is the silent gap longer or shorter with 6 neighbours?
   Runs onset measurement at the strongest-drive combo that tips cleanly.

We don't predict outcomes. Both questions are answered by what the model shows.
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

from void_cascade.cascade_breakdown_3d import (   # noqa: E402
    BreakdownParams3d,
    run_buildup_tip_3d,
    run_onset_measurement_3d,
)

L      = 20
STEPS  = 3200
SEED   = 0
# Per-cell drive rate matched to 2D baseline (4/6400 = 0.000625)
# L=20 cube = 8000 cells -> drive_sites = 5
BASE_DRIVE_SITES = 5

HPC_VALS    = [0.1, 0.2, 0.3]
DRIVE_AMTS  = [1.0, 2.0]

print("=== Q1: 3D eruption front sweep ===")
print(f"L={L}  steps={STEPS}  seed={SEED}")
print(f"{'hpc':>6}  {'drive':>6}  {'tip':>8}  {'maxT':>7}  {'resid':>10}")

results_q1 = {}
for hpc in HPC_VALS:
    for da in DRIVE_AMTS:
        p = BreakdownParams3d(hpc=hpc, drive_amount=da, drive_sites=BASE_DRIVE_SITES)
        r = run_buildup_tip_3d(L=L, steps=STEPS, params=p, seed=SEED)
        results_q1[(hpc, da)] = r
        tip = str(r["tip_step"]) if r["tip_step"] is not None else "none"
        print(f"{hpc:>6.1f}  {da:>6.1f}  {tip:>8}  {r['max_temp']:>7.2f}  "
              f"{r['energy_residual']:>10.2e}")

# Find combos that produced a clear tip (max_temp > e_ign + 0.5 = 3.0)
e_ign = 2.5
clear_tips = {k: v for k, v in results_q1.items()
              if v["tip_step"] is not None and v["max_temp"] > e_ign + 0.5}

print(f"\n  {len(clear_tips)}/{len(results_q1)} combos reached sustained ignition "
      f"(max_temp > {e_ign + 0.5})")

# --- Q2: onset ordering ---
# Use strongest drive that gave a clear tip; fall back to hpc=0.3, drive=2.0
if clear_tips:
    best_k = max(clear_tips, key=lambda k: clear_tips[k]["max_temp"])
    best_hpc, best_da = best_k
else:
    best_hpc, best_da = 0.3, 2.0
    print("  No clear tip found; running onset at hpc=0.3, drive=2.0 anyway.")

print(f"\n=== Q2: 3D onset ordering  (hpc={best_hpc}, drive={best_da}) ===")
p_onset = BreakdownParams3d(hpc=best_hpc, drive_amount=best_da,
                             drive_sites=BASE_DRIVE_SITES)
r_onset = run_onset_measurement_3d(L=L, steps=STEPS, params=p_onset, seed=SEED,
                                    sample_every=20)

t  = np.array(r_onset["t_axis"])
ls = np.array(r_onset["load_std"])
cl = np.array(r_onset["corr_lengths"])
nh = np.array(r_onset["n_hot"])

print(f"  tip_step: {r_onset['tip_step']}")
print(f"  max load_std:  {ls.max():.4f}  at t={t[ls.argmax()]}")
print(f"  max corr_len:  {cl.max():.2f}   at t={t[cl.argmax()]}")
print(f"  n_hot final:   {nh[-1]}")

thresh_10 = {}
for name, series in [("load_std", ls), ("corr_length", cl), ("n_hot", nh)]:
    mx = series.max()
    if mx > 0:
        cross = np.where(series >= 0.1 * mx)[0]
        thresh_10[name] = int(t[cross[0]]) if len(cross) else None
    else:
        thresh_10[name] = None

print("\n  First 10%-of-max crossings:")
for name, step in thresh_10.items():
    print(f"    {name:30s}: step {step}")

ordering = sorted([k for k, v in thresh_10.items() if v is not None],
                  key=lambda k: thresh_10[k])
print(f"\n  3D ordering: {' -> '.join(ordering)}")

# Compare with 2D baseline
print("\n  2D baseline (L=80, hpc=0.1, drive=1.0, seed=0):")
print("    corr_length(0) -> load_std(20) -> ignition(2038)")
print("    silent gap: ~2000 steps between first structure and tip")

if r_onset["tip_step"] is not None and thresh_10.get("corr_length") is not None:
    gap_3d = r_onset["tip_step"] - thresh_10["corr_length"]
    print(f"  3D silent gap: {gap_3d} steps  "
          f"({'longer' if gap_3d > 2000 else 'shorter'} than 2D ~2000)")

# --- plots ---
fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# Q1 heatmaps: tip_step and max_temp
tip_grid  = np.full((len(HPC_VALS), len(DRIVE_AMTS)), np.nan)
temp_grid = np.full((len(HPC_VALS), len(DRIVE_AMTS)), np.nan)
for i, hpc in enumerate(HPC_VALS):
    for j, da in enumerate(DRIVE_AMTS):
        r = results_q1[(hpc, da)]
        if r["tip_step"] is not None:
            tip_grid[i, j]  = r["tip_step"]
        temp_grid[i, j] = r["max_temp"]

ax = axes[0, 0]
im = ax.imshow(tip_grid, aspect="auto", cmap="viridis_r", origin="lower")
ax.set_xticks(range(len(DRIVE_AMTS))); ax.set_xticklabels(DRIVE_AMTS)
ax.set_yticks(range(len(HPC_VALS)));   ax.set_yticklabels(HPC_VALS)
ax.set_xlabel("drive_amount"); ax.set_ylabel("hpc")
ax.set_title("tip_step  (darker = earlier)")
plt.colorbar(im, ax=ax)

ax = axes[0, 1]
im = ax.imshow(temp_grid, aspect="auto", cmap="hot", origin="lower")
ax.set_xticks(range(len(DRIVE_AMTS))); ax.set_xticklabels(DRIVE_AMTS)
ax.set_yticks(range(len(HPC_VALS)));   ax.set_yticklabels(HPC_VALS)
ax.set_xlabel("drive_amount"); ax.set_ylabel("hpc")
ax.set_title("max_temp  (brighter = hotter)")
plt.colorbar(im, ax=ax)

# Q2 onset traces
ax = axes[1, 0]
ax2 = ax.twinx()
ax.plot(t, ls / (ls.max() or 1), color="steelblue",  lw=1.5, label="load_std (norm)")
ax.plot(t, cl / (cl.max() or 1), color="seagreen",   lw=1.5, label="corr_len (norm)")
ax2.plot(t, nh, color="firebrick", lw=1.2, ls="--", label="n_hot")
if r_onset["tip_step"] is not None:
    ax.axvline(r_onset["tip_step"], color="k", lw=1, ls=":", alpha=0.5)
ax.set_xlabel("step"); ax.set_ylabel("normalised value")
ax2.set_ylabel("n_hot", color="firebrick")
ax.set_title(f"3D onset  hpc={best_hpc} drive={best_da}")
ax.legend(loc="upper left", fontsize=8)

# Q2: max_temp traces for all combos
ax = axes[1, 1]
colors = plt.cm.tab10(np.linspace(0, 1, len(results_q1)))
for (hpc, da), col in zip(sorted(results_q1), colors):
    r = results_q1[(hpc, da)]
    t_r = np.array(r["t_axis"])
    T_r = np.array(r["temp_trace"])
    ax.plot(t_r, T_r, color=col, lw=1.2,
            label=f"hpc={hpc} da={da}")
    if r["tip_step"] is not None:
        ax.axvline(r["tip_step"], color=col, lw=0.6, ls=":", alpha=0.4)
ax.axhline(e_ign, color="k", lw=0.8, ls="--", alpha=0.3)
ax.set_xlabel("step"); ax.set_ylabel("max temp")
ax.set_title("3D max_temp traces (dotted = tip)")
ax.legend(fontsize=7, ncol=2)

fig.suptitle(f"3D eruption front + onset  L={L}")
fig.tight_layout()

out = REPO / "data" / "outputs" / "3d_eruption_onset"
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "3d_eruption_onset.png", dpi=150)
print(f"\nPlot saved: {out}/3d_eruption_onset.png")
