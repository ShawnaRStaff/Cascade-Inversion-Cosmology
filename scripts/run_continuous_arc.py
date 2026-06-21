"""Continuous arc: cold substrate -> buildup -> tip -> eruption front.

Three comparisons in one run:
  (a) No rigidity (rigidity=0)      -- baseline; cold cells are normal gas
  (b) Half rigidity (rigidity=0.5)  -- cold cells lend half their load as pressure
  (c) Full rigidity (rigidity=1.0)  -- cold cells lend ALL load as pressure (stiffest)

For each: does the tip happen, when, how fast does the front grow, and does
rigidity make the interface dynamics more violent (higher peak speed)?

The model leads: we report what the numbers say, not what we hoped.
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

from void_cascade.cascade_breakdown import BreakdownParams, run_full_arc  # noqa: E402

L = 80
BUILDUP = 2500   # steps to run before measuring eruption phase
AFTER   = 700    # steps to continue after buildup (captures the post-tip arc)
SEED    = 0
SAMPLE  = 30     # snapshot every N steps

params = BreakdownParams()  # default: drive_sites=4, drive_amount=1.0, hpc=0.1, thr=2.0

CONFIGS = [
    ("no rigidity",   0.0),
    ("rigid=0.5",     0.5),
    ("rigid=1.0 max", 1.0),
]


def run_config(label, rigidity):
    print(f"\n--- {label} ---")
    r = run_full_arc(
        L=L, steps_buildup=BUILDUP, steps_after=AFTER, params=params,
        seed=SEED, sample_every=SAMPLE, rigidity=rigidity,
    )
    tip = r["tip_step"]
    resid = r["energy_residual"]
    snaps = r["snapshots"]
    post = [s for s in snaps if tip is not None and s["t"] >= tip]
    first_n = post[0]["n_hot"]  if post else 0
    last_n  = post[-1]["n_hot"] if post else 0
    peak_speed = max((s["max_speed"] for s in post), default=0.0)
    print(f"  tip step       : {tip}")
    print(f"  energy residual: {resid:.2e}")
    print(f"  n_hot at tip   : {first_n}  (out of {L*L} cells = {first_n/L**2*100:.1f}%)")
    print(f"  n_hot at end   : {last_n}  ({last_n/L**2*100:.1f}%)")
    print(f"  peak speed (post-tip): {peak_speed:.3f}")
    print(f"\n  Time trace (t, n_hot, max_temp, max_speed):")
    for s in snaps:
        marker = " <-- TIP" if s["phase"] == "tip" else ""
        print(f"    t={s['t']:5d}  n_hot={s['n_hot']:5d}  maxT={s['max_temp']:.2f}  "
              f"spd={s['max_speed']:.3f}  [{s['phase']}]{marker}")
    return r


results = {}
for label, rig in CONFIGS:
    results[label] = run_config(label, rig)

# --- summary table ---
print("\n\n=== SUMMARY ===")
print(f"{'config':<20}  {'tip':>6}  {'n_hot_end':>9}  {'peak_spd':>9}  {'resid':>10}")
for label, rig in CONFIGS:
    r  = results[label]
    tip = r["tip_step"] or -1
    snaps = r["snapshots"]
    post = [s for s in snaps if r["tip_step"] is not None and s["t"] >= r["tip_step"]]
    last_n = post[-1]["n_hot"] if post else 0
    peak_speed = max((s["max_speed"] for s in post), default=0.0)
    print(f"{label:<20}  {tip:>6}  {last_n:>9}  {peak_speed:>9.3f}  {r['energy_residual']:>10.2e}")

# --- plot ---
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
colors = ["steelblue", "darkorange", "firebrick"]

for ax_idx, (label, rig) in enumerate(CONFIGS):
    r = results[label]
    snaps = r["snapshots"]
    t_arr   = [s["t"]         for s in snaps]
    nhot    = [s["n_hot"]     for s in snaps]
    maxtemp = [s["max_temp"]  for s in snaps]
    speed   = [s["max_speed"] for s in snaps]
    c = colors[ax_idx]
    tip = r["tip_step"]

    ax = axes[ax_idx]
    ax2 = ax.twinx()
    ax.plot(t_arr, nhot,    color=c, lw=1.8, label="n_hot cells")
    ax2.plot(t_arr, speed, color=c, lw=1.2, ls="--", alpha=0.7, label="max speed")
    if tip is not None:
        ax.axvline(tip, color="k", lw=1, ls=":", alpha=0.6)
        ax.text(tip + 20, max(nhot) * 0.6, f"tip\n{tip}", fontsize=7, color="k")
    ax.set_title(label, fontsize=10)
    ax.set_xlabel("step")
    ax.set_ylabel("hot cells (n_hot)", color=c)
    ax2.set_ylabel("max speed", color="gray")
    ax2.tick_params(axis="y", colors="gray")

fig.suptitle(
    f"Continuous arc  L={L}  buildup={BUILDUP}  after={AFTER}\n"
    f"Left axis: hot cells (solid).  Right axis: peak speed (dashed)."
)
fig.tight_layout()

out = REPO / "data" / "outputs" / "continuous_arc"
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "arc_comparison.png", dpi=150)
print(f"\nPlot saved: {out}/arc_comparison.png")
