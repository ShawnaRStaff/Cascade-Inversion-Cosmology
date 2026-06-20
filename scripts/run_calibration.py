"""Dimensional calibration: what do the model's numbers imply in physical units?

The model is dimensionless. To reach physical numbers we need three ANCHORS:
  - length: how many Mpc per lattice cell
  - time:   how many years per step
  - energy: how many Joules (or K) per model energy unit
These are FREE until pinned by real observations. This script takes the
*illustrative* anchors that M4 used (themselves free parameters from matching
galaxy xi(r)) and works out what the model's speeds become -- and whether that
is even self-consistent. No claim of truth; this is the grounding check the
project's docs call the rate-limiter.
"""

from __future__ import annotations

# Physical constants
C_KMS = 299_792.458          # speed of light, km/s
KM_PER_MPC = 3.0857e19
S_PER_MYR = 3.1557e13

# --- Illustrative anchors (FREE PARAMETERS, from M4 matching; not derived) ---
MPC_PER_CELL = 8.5           # from xi(r) r_0 ~ 5 h^-1 Mpc match (M4)
MYR_PER_STEP = 2.6           # from M4 time-mapping (1 drop ~ 2.6 Myr)
# energy anchor: NONE pinned yet.

CELL_STEP_TO_KMS = MPC_PER_CELL * KM_PER_MPC / (MYR_PER_STEP * S_PER_MYR)


def to_kms(cells_per_step: float) -> float:
    return cells_per_step * CELL_STEP_TO_KMS


def main() -> None:
    print("=== Dimensional calibration (ILLUSTRATIVE -- anchors are free) ===\n")
    print(f"Anchors (free, from M4 matching): {MPC_PER_CELL} Mpc/cell, {MYR_PER_STEP} Myr/step")
    print(f"=> 1 cell/step = {CELL_STEP_TO_KMS:.3e} km/s = {CELL_STEP_TO_KMS / C_KMS:.2f} x c\n")

    speeds = {
        "conserved-energy front (~0.54 cells/step)": 0.54,
        "detonation front (~1.0 cells/step)": 1.0,
        "implosion material velocity (u0 ~1.5 cells/step)": 1.5,
    }
    print(f"{'quantity':<48}{'km/s':>14}{'x c':>10}")
    for name, v in speeds.items():
        kms = to_kms(v)
        print(f"{name:<48}{kms:>14.3e}{kms / C_KMS:>10.2f}")

    print("\nFinding: under these (free) M4 anchors, EVERY model speed of order")
    print("1 cell/step is ~5-11x the speed of light. Implications:")
    print(" - An ACTIVITY/pattern front (no material moving) may be superluminal")
    print("   (like a phase velocity) -- consistent with our 'activity wavefront' label.")
    print(" - But MATERIAL motion (the implosion fluid) must be sub-light -- so the")
    print("   M4 anchors are NOT mutually consistent with material velocities: to make")
    print(f"   u0 sub-light you'd need ~{to_kms(1.5)/C_KMS:.0f}x larger time/step or smaller Mpc/cell.")
    print(" - The ENERGY anchor is entirely unpinned (no temperature/energy-density match yet).")
    print("\nSo calibration is UNDERDETERMINED and the existing free anchors are")
    print("over-constrained for material motion. Pinning length+time+energy from real")
    print("data (e.g. a temperature for the plasma spike, an energy density for z=0.616)")
    print("is the rate-limiter before any quantitative cosmological claim.")


if __name__ == "__main__":
    main()
