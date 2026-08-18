"""Dimensional calibration: lattice units -> physical units, with two time anchors.

The model is dimensionless. Physical claims require anchors, and the
2026-06-19 check ("all model speeds are 5-16x c") conflated two distinct
model timescales when it converted per-sweep front speeds with the
per-drop time anchor. SOC dynamics is built on a separation of
timescales:

  - SWEEP time: one parallel relaxation update inside an avalanche.
    Fronts (detonation, implosion material) move in cells/sweep, so
    causality (v <= c) bounds the sweep time from BELOW.
  - DROP time: the interval between successive grain drops (the drive).
    The SOC limit is drop time >> avalanche duration; cosmologically it
    is pinned by matching a model epoch to a real age, not by causality.

Anchors:
  mpc_per_cell   length  (M4's xi(r) match gives ~8.5 h^-1 Mpc, a fit)
  myr_per_sweep  avalanche-internal time (>= causal minimum)
  myr_per_drop   drive time (free until an epoch match pins it)
  j_per_grain    energy (candidate: z=0.616 grains/cell <-> vacuum
                 energy density rho_Lambda)

All functions are pure arithmetic; the caller supplies anchor values and
model measurements. No hidden state, no I/O.
"""

from __future__ import annotations

# Physical constants (SI-derived, Planck 2018 where cosmological)
C_KMS = 299_792.458        # speed of light [km/s]
MPC_IN_KM = 3.0857e19      # [km/Mpc]
MYR_IN_S = 3.1557e13       # [s/Myr]
RHO_VACUUM_J_M3 = 5.3e-10  # vacuum (dark) energy density [J/m^3], Planck 2018


def min_myr_per_sweep(mpc_per_cell: float, max_speed_cells_per_sweep: float) -> float:
    """Causal lower bound on the sweep time.

    The fastest material front in the model moves max_speed cells per
    sweep. Requiring that to be <= c gives the minimum physical duration
    of one relaxation sweep for a given lattice spacing.
    """
    km_per_sweep = max_speed_cells_per_sweep * mpc_per_cell * MPC_IN_KM
    return km_per_sweep / C_KMS / MYR_IN_S


def front_speed_over_c(
    cells_per_sweep: float, mpc_per_cell: float, myr_per_sweep: float
) -> float:
    """A model front speed expressed as a fraction of c under given anchors."""
    kms = cells_per_sweep * mpc_per_cell * MPC_IN_KM / (myr_per_sweep * MYR_IN_S)
    return kms / C_KMS


def j_per_grain_from_vacuum_density(
    rho_vac_j_m3: float, mpc_per_cell: float, z_sat: float
) -> float:
    """Energy anchor: identify the saturated grain density with vacuum energy.

    The saturated substrate holds z_sat grains per cell (z = 0.616,
    universal across L). If that standing grain density IS the vacuum
    energy density, one grain carries rho * V_cell / z_sat.
    """
    m_per_cell = mpc_per_cell * MPC_IN_KM * 1e3
    return rho_vac_j_m3 * m_per_cell**3 / z_sat


def event_energy_j(n_topples: int, j_per_grain: float) -> float:
    """Energy moved by an avalanche: each topple relocates 2 grains."""
    return 2.0 * n_topples * j_per_grain


def event_duration_myr(n_sweeps: int, myr_per_sweep: float) -> float:
    """Physical duration of an avalanche of n_sweeps parallel updates."""
    return n_sweeps * myr_per_sweep


def myr_per_drop_from_epoch(n_drops_to_now: int, age_gyr: float) -> float:
    """Drive-time anchor from matching a model epoch to a physical age."""
    return age_gyr * 1e3 / n_drops_to_now


def consistency_report(
    mpc_per_cell: float,
    myr_per_sweep: float,
    myr_per_drop: float,
    j_per_grain: float,
    material_speed_cells_per_sweep: float,
    peak_event_sweeps: int | None = None,
    age_gyr: float = 13.8,
) -> dict:
    """Evaluate an anchor set: implied physical quantities + violated constraints.

    Violations checked:
      material_superluminal  — fastest material front exceeds c
      event_outlives_universe — a single avalanche's physical duration
        exceeds the age of the universe (only if peak_event_sweeps given)
      timescale_inversion     — drop interval shorter than one sweep
        (breaks the SOC separation of timescales)
    """
    violations = []
    v_over_c = front_speed_over_c(
        material_speed_cells_per_sweep, mpc_per_cell, myr_per_sweep
    )
    if v_over_c > 1.0:
        violations.append("material_superluminal")

    peak_duration_myr = None
    if peak_event_sweeps is not None:
        peak_duration_myr = event_duration_myr(peak_event_sweeps, myr_per_sweep)
        if peak_duration_myr > age_gyr * 1e3:
            violations.append("event_outlives_universe")

    if myr_per_drop < myr_per_sweep:
        violations.append("timescale_inversion")

    return {
        "material_speed_over_c": v_over_c,
        "min_myr_per_sweep": min_myr_per_sweep(
            mpc_per_cell, material_speed_cells_per_sweep
        ),
        "peak_event_duration_myr": peak_duration_myr,
        "j_per_grain": j_per_grain,
        "violations": violations,
    }
