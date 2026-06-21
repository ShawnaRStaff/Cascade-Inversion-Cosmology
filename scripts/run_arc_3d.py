"""3D buildup -> tip: does the arc survive three dimensions?

Runs the 3D engine at L=20 and compares directly with the 2D engine at L=80
(same per-cell drive rate: drive_sites / L^d = constant). Reports:

  - Does the 3D substrate tip? When?
  - Is energy conserved?
  - Does the tip step shift in 3D (different neighbourhood = different criticality)?

Honest scope: L=20 in 3D (8000 cells, comparable to L=90 in 2D). The 3D Manna
critical exponent is τ≈1.47 vs 2D τ=1.27, so larger cascades are expected;
we don't know if the tip step scales the same way.
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

from void_cascade.cascade_breakdown import BreakdownParams, run_buildup_tip          # noqa: E402
from void_cascade.cascade_breakdown_3d import BreakdownParams3d, run_buildup_tip_3d  # noqa: E402

# Per-cell drive rate: 4 sites on L=80 2D = 4/6400 = 0.000625 per cell per step.
# For L=20 3D: 8000 cells. Same per-cell rate => drive_sites = 0.000625 * 8000 = 5.
L_3D  = 20
L_2D  = 80
STEPS = 3200
SEED  = 0

rate_2d = 4 / (L_2D ** 2)              # 0.000625
drive_sites_3d = round(rate_2d * (L_3D ** 3))   # 5 sites on L=20 cube

params_2d = BreakdownParams(drive_sites=4,           drive_amount=1.0, hpc=0.1)
params_3d = BreakdownParams3d(drive_sites=drive_sites_3d, drive_amount=1.0, hpc=0.1)

print(f"=== 3D arc comparison ===")
print(f"2D: L={L_2D}  drive_sites=4  per-cell={4/L_2D**2:.6f}")
print(f"3D: L={L_3D}  drive_sites={drive_sites_3d}  per-cell={drive_sites_3d/L_3D**3:.6f}")

print(f"\n[2D L={L_2D}]...")
r2 = run_buildup_tip(L=L_2D, steps=STEPS, params=params_2d, seed=SEED)
print(f"  tip_step={r2['tip_step']}  max_temp={r2['max_temp']:.2f}  "
      f"resid={r2['energy_residual']:.2e}")

print(f"\n[3D L={L_3D}]...")
r3 = run_buildup_tip_3d(L=L_3D, steps=STEPS, params=params_3d, seed=SEED)
print(f"  tip_step={r3['tip_step']}  max_temp={r3['max_temp']:.2f}  "
      f"resid={r3['energy_residual']:.2e}")

print("\n=== SUMMARY ===")
print(f"{'dim':>5}  {'L':>4}  {'cells':>8}  {'tip':>8}  {'maxT':>7}  {'resid':>10}")
print(f"{'2D':>5}  {L_2D:>4}  {L_2D**2:>8}  "
      f"{str(r2['tip_step']):>8}  {r2['max_temp']:>7.2f}  {r2['energy_residual']:>10.2e}")
print(f"{'3D':>5}  {L_3D:>4}  {L_3D**3:>8}  "
      f"{str(r3['tip_step']):>8}  {r3['max_temp']:>7.2f}  {r3['energy_residual']:>10.2e}")

# --- plot ---
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
for ax, r, label, color in [(axes[0], r2, f"2D L={L_2D}", "steelblue"),
                              (axes[1], r3, f"3D L={L_3D}", "firebrick")]:
    t  = np.array(r["t_axis"])
    T  = np.array(r["temp_trace"])
    ld = np.array(r["load_trace"])
    ax.plot(t, T,  color=color,   lw=1.5, label="max temp")
    ax.plot(t, ld, color="gray",  lw=1.0, ls="--", label="mean load")
    if r["tip_step"] is not None:
        ax.axvline(r["tip_step"], color="k", lw=1, ls=":", alpha=0.6)
        ax.text(r["tip_step"] + 30, 1.6, f"tip\n{r['tip_step']}", fontsize=8)
    ax.axhline(2.5, color="k", lw=0.6, ls="--", alpha=0.3)
    ax.set_title(label); ax.set_xlabel("step"); ax.set_ylabel("temperature / mean load")
    ax.legend(fontsize=8)

fig.suptitle("3D vs 2D buildup -> tip  (same per-cell drive rate)")
fig.tight_layout()

out = REPO / "data" / "outputs" / "arc_3d"
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "arc_3d_vs_2d.png", dpi=150)
print(f"\nPlot saved: {out}/arc_3d_vs_2d.png")
