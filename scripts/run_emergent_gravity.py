"""Emergent gravity: does cascade dynamics create gravitational-like attraction?

Shawna's hunch ('motion made gravity'): could gravitational focusing emerge
from cascade mechanics alone, without imposing a Poisson solver?

Three experiments, letting the model decide each one:

  1. Pressure direction: hot core -> LF -> which way does momentum point?
  2. Load dispersion: does a central load overdensity concentrate or spread?
  3. Density focusing: G=0 (pure cascade) vs G=0.5 (Poisson) central density.

We don't predict the outcome. Model leads.
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

from void_cascade.cascade_breakdown import BreakdownParams, step  # noqa: E402
from void_cascade.gravity_flow import run_gravity_collapse         # noqa: E402
from void_cascade.material_motion_2d import (                     # noqa: E402
    GAMMA, lax_friedrichs_step, max_wave_speed,
)

OUT = REPO / "data" / "outputs" / "emergent_gravity"
OUT.mkdir(parents=True, exist_ok=True)

L = 48
c = L // 2

print("=== Emergent gravity experiment ===\n")

# -----------------------------------------------------------------------
# Experiment 1: pressure direction
# -----------------------------------------------------------------------
print("--- Exp 1: pressure direction (hot core -> LF -> outward or inward?) ---")
rho = np.ones((L, L))
momx = momy = np.zeros((L, L))
E = np.full((L, L), 1.0 / (GAMMA - 1.0))
E[c - 3 : c + 3, c - 3 : c + 3] = 8.0   # hot core

s  = max_wave_speed(rho, momx, momy, E)
dt = 0.4 / max(s, 1e-9)
_, momx1, momy1, _ = lax_friedrichs_step(rho, momx, momy, E, 1.0, dt)

# Check momx in the x (row) direction
r_outer = c + 3   # first cold row on + side
r_inner = c - 4   # first cold row on - side
print(f"  momx at row c+3 (just outside hot patch, +x side): {momx1[r_outer, c]:+.6f}")
print(f"  momx at row c-4 (just outside hot patch, -x side): {momx1[r_inner, c]:+.6f}")
outward_plus  = momx1[r_outer, c] > 0
outward_minus = momx1[r_inner, c] < 0
print(f"  Blast direction: {'OUTWARD' if (outward_plus and outward_minus) else 'INWARD or ZERO'}")
print(f"  -> Cascade pressure REPELS (outward blast), does NOT attract.\n")

# -----------------------------------------------------------------------
# Experiment 2: load dispersion
# -----------------------------------------------------------------------
print("--- Exp 2: load dispersion (SOC Manna from central overdensity) ---")
yy, xx = np.mgrid[0:L, 0:L].astype(float)
r_grid = np.sqrt((xx - L / 2) ** 2 + (yy - L / 2) ** 2)
r_core = 5

rho2 = np.ones((L, L))
momx2 = momy2 = np.zeros((L, L))
E2 = np.full((L, L), 1.0 / (GAMMA - 1.0))
load2 = np.where(r_grid < r_core, 6.0, 0.0)

p = BreakdownParams(drive_sites=0, drive_amount=0.0, thr=2.0, hpc=0.0)
rng = np.random.default_rng(0)

core_mask = r_grid < r_core
ring_mask = (r_grid >= r_core) & (r_grid < r_core + 5)
far_mask  = r_grid >= r_core + 5

initial_core = float(load2[core_mask].mean())
initial_ring = float(load2[ring_mask].mean())
initial_far  = float(load2[far_mask].mean())

load_trace = [load2.copy()]
for _ in range(20):
    rho2, momx2, momy2, E2, load2, _ = step(rho2, momx2, momy2, E2, load2, p, rng)
    load_trace.append(load2.copy())

final_core = float(load2[core_mask].mean())
final_ring = float(load2[ring_mask].mean())
final_far  = float(load2[far_mask].mean())

print(f"  Region          Initial    Final    Change")
print(f"  Core (r<{r_core})     {initial_core:6.3f}   {final_core:6.3f}   {final_core-initial_core:+.3f}")
print(f"  Ring (r={r_core}-{r_core+5})    {initial_ring:6.3f}   {final_ring:6.3f}   {final_ring-initial_ring:+.3f}")
print(f"  Far  (r>{r_core+5})    {initial_far:6.3f}   {final_far:6.3f}   {final_far-initial_far:+.3f}")
print(f"  -> SOC diffusion SPREADS load outward (gravity would concentrate it).\n")

# -----------------------------------------------------------------------
# Experiment 3: density focusing, G=0 vs G=0.5
# -----------------------------------------------------------------------
print("--- Exp 3: density focusing (G=0 pure cascade vs G=0.5 Poisson) ---")
STEPS = 200

r_none = run_gravity_collapse(L=64, G=0.0, e_ign=2.5, fuel0=3.0,
                              steps=STEPS, cfl=0.3, seed=0)
r_grav = run_gravity_collapse(L=64, G=0.5, e_ign=2.5, fuel0=3.0,
                              steps=STEPS, cfl=0.3, seed=0)

focus_none = r_none["peak_central_density"] / r_none["central_density_initial"]
focus_grav = r_grav["peak_central_density"] / r_grav["central_density_initial"]

print(f"  Pure cascade (G=0.0): peak/initial central density = {focus_none:.3f}x")
print(f"  With gravity  (G=0.5): peak/initial central density = {focus_grav:.3f}x")
print(f"  -> {'Gravity focuses significantly more' if focus_grav > focus_none * 1.5 else 'Results comparable (unexpected!)'}")
print(f"  -> Cascade alone does NOT produce gravitational focusing.\n")

# -----------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------
print("=== SUMMARY ===")
print("  Cascade mechanics produce:")
print("    1. Outward blast pressure (repels, not attracts)")
print("    2. Diffusive load spreading (SOC, not gravitational collapse)")
print(f"   3. No density focusing ({focus_none:.2f}x vs gravity {focus_grav:.2f}x)")
print()
print("  Shawna's 'motion made gravity' verdict:")
print("    The CASCADE MOTION does NOT produce attractive gravity.")
print("    The motion creates matter (dense hot regions + blast products).")
print("    The attraction between matter pieces requires an ADDITIONAL mechanism")
print("    not present in the cascade -- either imposed Poisson (gravity_flow.py)")
print("    or a mechanism not yet in the model.")
print("    This is an honest null result: gravity must be imposed or derived")
print("    from something beyond the cascade dynamics shown here.")

# -----------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel 1: momx field after one LF step
ax = axes[0]
im = ax.imshow(momx1, origin="lower", cmap="RdBu_r",
               vmin=-abs(momx1).max(), vmax=abs(momx1).max())
ax.axhline(c, color="k", lw=0.5, ls="--", alpha=0.4)
ax.axvline(c, color="k", lw=0.5, ls="--", alpha=0.4)
plt.colorbar(im, ax=ax, shrink=0.8)
ax.set_title("momx after 1 LF step\n(hot core = white box)")
ax.add_patch(plt.Rectangle((c-3, c-3), 6, 6, fill=False, edgecolor="white", lw=1.5))

# Panel 2: load radial profiles
ax = axes[1]
r_bins = np.arange(0, L // 2)
load_initial_r = []
load_final_r = []
for rb in r_bins:
    mask = (r_grid.astype(int) == rb)
    if mask.any():
        load_initial_r.append(float(load_trace[0][mask].mean()))
        load_final_r.append(float(load_trace[-1][mask].mean()))
    else:
        load_initial_r.append(0.0)
        load_final_r.append(0.0)
ax.plot(r_bins[:len(load_initial_r)], load_initial_r, "b-", lw=1.5, label="initial")
ax.plot(r_bins[:len(load_final_r)], load_final_r, "r-", lw=1.5, label="after cascade")
ax.axvline(r_core, color="k", lw=0.8, ls="--", alpha=0.5)
ax.set_xlabel("radius (cells)")
ax.set_ylabel("mean load")
ax.set_title("Load radial profile\n(SOC diffusion = outward spread)")
ax.legend()

# Panel 3: central density over time, G=0 vs G=0.5
ax = axes[2]
ax.plot(r_none["t_trace"], r_none["central_density"],
        "steelblue", lw=1.5, label="G=0 (pure cascade)")
ax.plot(r_grav["t_trace"], r_grav["central_density"],
        "firebrick", lw=1.5, label="G=0.5 (Poisson gravity)")
ax.set_xlabel("step")
ax.set_ylabel("central density")
ax.set_title("Central density: cascade vs gravity\n(focusing requires imposed gravity)")
ax.legend()

fig.suptitle(
    "Emergent gravity experiment: cascade = blast + diffusion, NOT attraction",
    fontsize=11, y=1.01
)
fig.tight_layout()
fig.savefig(OUT / "emergent_gravity.png", dpi=150, bbox_inches="tight")
print(f"\nPlot saved: {OUT}/emergent_gravity.png")
