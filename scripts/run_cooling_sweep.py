"""Cooling sweep: does a stable cold regime exist when heat bleeds away?

At each step, cells above the ambient baseline lose `cooling` fraction of
their excess heat to the void. Sweeps cooling=[0, 0.01, 0.05, 0.1, 0.3, 0.5, 1.0].

Questions:
  1. Is there a critical cooling rate above which the substrate NEVER tips?
  2. Below that rate, how does cooling change TIP TIMING?
  3. Full energy accounting: E + load + total_cooled == E0 + driven (residual~0).

We don't predict the critical rate. The model decides.
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

from void_cascade.cascade_breakdown import BreakdownParams, run_cooling_arc  # noqa: E402

COOLING_VALS = [0.0, 0.01, 0.05, 0.1, 0.3, 0.5, 1.0]
L     = 60
STEPS = 3500
SEED  = 0
PARAMS = BreakdownParams(hpc=0.1, drive_amount=1.0, drive_sites=4)

print(f"=== Cooling sweep ===")
print(f"L={L}  steps={STEPS}  seed={SEED}")
print(f"{'cooling':>9}  {'tip':>8}  {'maxT':>7}  {'cooled':>12}  {'resid':>10}")

results = {}
for c in COOLING_VALS:
    r = run_cooling_arc(L=L, steps=STEPS, params=PARAMS, seed=SEED, cooling=c)
    results[c] = r
    tip = str(r["tip_step"]) if r["tip_step"] is not None else "none"
    print(f"{c:>9.2f}  {tip:>8}  {r['max_temp']:>7.2f}  "
          f"{r['total_cooled']:>12.2f}  {r['energy_residual']:>10.2e}")

# Find critical cooling rate
tipped    = [c for c in COOLING_VALS if results[c]["tip_step"] is not None]
not_tipped = [c for c in COOLING_VALS if results[c]["tip_step"] is None]

print("\n=== SUMMARY ===")
if tipped and not_tipped:
    crit_lo = max(tipped)
    crit_hi = min(not_tipped)
    print(f"  Critical cooling rate: between {crit_lo} and {crit_hi}")
    print(f"  Below {crit_hi}: substrate tips (heat accumulates faster than it bleeds)")
    print(f"  At/above {crit_hi}: substrate stays cold (cooling wins)")
elif not not_tipped:
    print(f"  All cooling rates tipped -- stable cold regime not found in this range")
elif not tipped:
    print(f"  No cooling rate tipped -- even cooling=0 suppressed the tip "
          f"(unexpected; check params)")

# --- plot ---
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(COOLING_VALS)))

ax = axes[0]
for c, col in zip(COOLING_VALS, colors):
    r = results[c]
    t = np.array(r["t_axis"])
    T = np.array(r["temp_trace"])
    ax.plot(t, T, color=col, lw=1.3, label=f"cool={c}")
    if r["tip_step"] is not None:
        ax.axvline(r["tip_step"], color=col, lw=0.7, ls=":", alpha=0.5)
ax.axhline(2.5, color="k", lw=0.7, ls="--", alpha=0.3)
ax.set_xlabel("step"); ax.set_ylabel("max temp")
ax.set_title("Max temperature over time (dotted = tip)")
ax.legend(fontsize=8, ncol=2)

ax = axes[1]
tip_steps = [results[c]["tip_step"] if results[c]["tip_step"] is not None else STEPS
             for c in COOLING_VALS]
bar_colors = ["firebrick" if results[c]["tip_step"] is not None else "steelblue"
              for c in COOLING_VALS]
bars = ax.bar(range(len(COOLING_VALS)), tip_steps, color=bar_colors, edgecolor="k", lw=0.5)
ax.set_xticks(range(len(COOLING_VALS)))
ax.set_xticklabels([str(c) for c in COOLING_VALS], fontsize=9)
ax.set_xlabel("cooling rate"); ax.set_ylabel("tip step (bar=steps if no tip)")
ax.set_title("Tip step by cooling rate\n(red=tipped, blue=no tip)")
ax.axhline(STEPS, color="steelblue", lw=0.8, ls="--", alpha=0.4)

fig.suptitle(f"Cooling sweep  L={L}  hpc=0.1  drive=1.0")
fig.tight_layout()

out = REPO / "data" / "outputs" / "cooling_sweep"
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "cooling_sweep.png", dpi=150)
print(f"\nPlot saved: {out}/cooling_sweep.png")
