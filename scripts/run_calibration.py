"""Dimensional calibration v2: two time anchors, measured event durations.

Supersedes the 2026-06-19 illustrative check, which converted per-SWEEP
front speeds with the per-DROP time anchor (2.6 Myr) and concluded
"everything is 5-16x c". With the timescales separated (see
src/void_cascade/calibration.py) the question becomes: does ANY anchor
set satisfy all constraints simultaneously? This script works the
constraint algebra with measured model quantities:

  - fastest material speed: u0 ~ 1.5 cells/sweep (implosion fluid)
  - peak-event durations: T = 1319/2265/4658/8039 sweeps at L=48/64/96/128
    (measured from M6 finals; T ~ L^1.84)
  - drops to the saturated era: ~3e6 (order of magnitude, L=96-128 runs)
  - M4 length anchor: 8.5 h^-1 Mpc/cell from the xi(r) r_0 match
  - energy anchor candidate: z=0.616 grains/cell <-> vacuum energy density
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from void_cascade.calibration import (
    RHO_VACUUM_J_M3,
    consistency_report,
    event_energy_j,
    j_per_grain_from_vacuum_density,
    min_myr_per_sweep,
    myr_per_drop_from_epoch,
)

MPC_PER_CELL_M4 = 8.5        # xi(r) match (a fit, epoch p=0.65)
U0_CELLS_PER_SWEEP = 1.5     # fastest MATERIAL speed (implosion fluid)
T_PEAK = {48: 1319, 64: 2265, 96: 4658, 128: 8039}   # sweeps, measured
TOPPLES_PEAK_L96 = 8.88e5    # measured mean
N_DROPS_TO_SATURATED_ERA = 3e6   # order of magnitude (L=96-128)
AGE_GYR = 13.8


def main() -> None:
    print("=== Dimensional calibration v2 (two time anchors) ===\n")

    # --- 1. The causal floor on sweep time under the M4 length anchor ---
    t_sweep_min = min_myr_per_sweep(MPC_PER_CELL_M4, U0_CELLS_PER_SWEEP)
    print(f"Causal minimum sweep time at {MPC_PER_CELL_M4} Mpc/cell: "
          f"{t_sweep_min:.1f} Myr/sweep")

    # --- 2. Peak events under that floor: do they fit inside cosmic time? ---
    print(f"\nPeak-event durations if material causality binds avalanche fronts:")
    for L, T in T_PEAK.items():
        dur_gyr = T * t_sweep_min / 1e3
        print(f"  L={L:>3}: T={T} sweeps -> {dur_gyr:,.0f} Gyr "
              f"({dur_gyr / AGE_GYR:.0f}x age of universe)")

    # --- 3. SOC separation + epoch match: what sweep time is allowed? ---
    t_drop = myr_per_drop_from_epoch(int(N_DROPS_TO_SATURATED_ERA), AGE_GYR)
    print(f"\nDrive-time anchor from epoch match: {AGE_GYR} Gyr / "
          f"{N_DROPS_TO_SATURATED_ERA:.0e} drops = {t_drop * 1e6:,.0f} yr/drop")
    t_sweep_soc_max = t_drop / T_PEAK[96]
    print(f"SOC separation (drop interval >> avalanche duration) at L=96 demands "
          f"sweep time << {t_sweep_soc_max * 1e6:.2f} yr")
    lmax_pc = (t_sweep_soc_max * 1e6) * 0.3066 / U0_CELLS_PER_SWEEP  # ly->pc via c*yr
    print(f"With material causality that caps the cell at ~{lmax_pc:.2f} pc "
          f"-- {MPC_PER_CELL_M4 * 1e6 / lmax_pc:.0e}x below the M4 length anchor.")

    print("""
CONCLUSION (the honest one): the four requirements
  (a) material causality on avalanche fronts,
  (b) SOC timescale separation,
  (c) saturated era reached within 13.8 Gyr,
  (d) cell scale = 8.5 Mpc from the xi(r) match
are mutually exclusive by ~7 orders of magnitude. Dropping exactly one
of them gives three self-consistent readings:

  DROP (a): avalanche fronts are ACTIVITY/PATTERN fronts (phase-velocity
    -like, no material transport) and may be superluminal; material
    bounds apply only to the fluid-layer experiments. Then (b)+(c)+(d)
    coexist: 8.5 Mpc cells, ~4,600 yr/drop, sweep time free below ~1 yr.
    This makes front superluminality a REQUIRED framework claim.
  DROP (c): the saturated era lies far beyond 13.8 Gyr -- today's
    universe is PRE-saturation (early in the loading curve), and the
    permanent catastrophe regime is the far future, not the present.
    This contradicts the framework's post-bang reading.
  DROP (d): cells are sub-parsec and the lattice models microphysics;
    the xi(r) match at 8.5 Mpc/cell was a coincidence of shape. The
    L<=192 boxes then span <100 pc and say nothing about LSS.
""")

    # --- 4. Energy anchor (independent of the time tangle) ---
    eps = j_per_grain_from_vacuum_density(RHO_VACUUM_J_M3, MPC_PER_CELL_M4, 0.616)
    e_peak = event_energy_j(int(TOPPLES_PEAK_L96), eps)
    print(f"Energy anchor (z=0.616 <-> rho_vac, {MPC_PER_CELL_M4} Mpc cells):")
    print(f"  1 grain = {eps:.2e} J")
    print(f"  L=96 peak event moves {e_peak:.2e} J "
          f"(~{e_peak / 1.35e70 * 100:.2f}% of observable-universe mass-energy)")
    print(f"  L=128 box side = {MPC_PER_CELL_M4 * 128 / 1e3:.2f} Gpc "
          f"(observable universe ~14 Gpc radius needs L~{int(14e3 / MPC_PER_CELL_M4)})")

    # --- 5. Machine-checkable report for the preferred (activity-front) reading ---
    r = consistency_report(
        mpc_per_cell=MPC_PER_CELL_M4,
        myr_per_sweep=1e-7,          # << 1 yr: activity front, c not binding
        myr_per_drop=t_drop,
        j_per_grain=eps,
        material_speed_cells_per_sweep=0.0,  # no material claim at lattice level
        peak_event_sweeps=T_PEAK[96],
        age_gyr=AGE_GYR,
    )
    print(f"\nActivity-front reading violations: {r['violations'] or 'none'}")


if __name__ == "__main__":
    main()
