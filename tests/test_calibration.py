"""Tests for dimensional calibration (lattice units -> physical units).

The model is dimensionless; physical claims need anchors. The 2026-06-19
calibration check found "16x c material motion" — but it converted
per-SWEEP front speeds with the per-DROP time anchor. SOC dynamics has
two separated timescales: avalanche-internal time (sweeps) and drive
time (interval between grain drops). The calibration module keeps them
distinct:

  - mpc_per_cell        (length anchor, from structure matching)
  - myr_per_sweep       (avalanche-internal time; bounded below by
                         causality: fronts must be sub-luminal)
  - myr_per_drop        (drive time; free until pinned by an epoch match)
  - j_per_grain         (energy anchor, e.g. from vacuum energy density)

All functions are pure arithmetic: values in, values out.
"""

import pytest

from void_cascade.calibration import (
    MPC_IN_KM,
    MYR_IN_S,
    C_KMS,
    min_myr_per_sweep,
    front_speed_over_c,
    j_per_grain_from_vacuum_density,
    event_energy_j,
    event_duration_myr,
    myr_per_drop_from_epoch,
    consistency_report,
)


def test_min_sweep_time_enforces_causality_exactly():
    # A front at v cells/sweep moves v * mpc_per_cell per sweep. Demanding
    # <= c gives myr_per_sweep >= v * mpc_per_cell * (Mpc in km) / c / (Myr in s).
    got = min_myr_per_sweep(mpc_per_cell=1.0, max_speed_cells_per_sweep=1.0)
    expected = MPC_IN_KM / C_KMS / MYR_IN_S
    assert got == pytest.approx(expected)
    # Scales linearly in both arguments.
    assert min_myr_per_sweep(8.5, 1.5) == pytest.approx(expected * 8.5 * 1.5)


def test_front_speed_at_exactly_c_when_sweep_time_is_minimal():
    t = min_myr_per_sweep(mpc_per_cell=8.5, max_speed_cells_per_sweep=1.5)
    assert front_speed_over_c(1.5, mpc_per_cell=8.5, myr_per_sweep=t) == pytest.approx(1.0)
    # Half the sweep time -> twice the speed.
    assert front_speed_over_c(1.5, 8.5, t / 2) == pytest.approx(2.0)


def test_energy_anchor_from_vacuum_density():
    # z grains/cell spread over a cell volume must reproduce rho_vac:
    # j_per_grain = rho_vac * cell_volume / z.
    # With a 1 m^3 "cell" the arithmetic is transparent: use mpc_per_cell
    # such that cell volume is exactly (1 Mpc)^3.
    rho = 5.3e-10  # J/m^3
    eps = j_per_grain_from_vacuum_density(rho_vac_j_m3=rho, mpc_per_cell=1.0, z_sat=0.616)
    m_per_mpc = MPC_IN_KM * 1e3
    assert eps == pytest.approx(rho * m_per_mpc**3 / 0.616)


def test_event_energy_counts_two_grains_per_topple():
    assert event_energy_j(n_topples=10, j_per_grain=2.0) == pytest.approx(40.0)


def test_event_duration_is_sweeps_times_sweep_time():
    assert event_duration_myr(n_sweeps=100, myr_per_sweep=0.5) == pytest.approx(50.0)


def test_drop_time_from_epoch_match():
    # If the model reaches "today" after N drops and that spans age_gyr:
    assert myr_per_drop_from_epoch(n_drops_to_now=1000, age_gyr=1.0) == pytest.approx(1.0)


def test_consistency_report_flags_superluminal_material():
    r = consistency_report(
        mpc_per_cell=8.5,
        myr_per_sweep=2.6,   # the old conflated anchor: should violate causality
        myr_per_drop=2.6,
        j_per_grain=1e60,
        material_speed_cells_per_sweep=1.5,
    )
    assert r["material_speed_over_c"] > 1.0
    assert "material_superluminal" in r["violations"]


def test_consistency_report_clean_when_causal():
    t = min_myr_per_sweep(8.5, 1.5)
    r = consistency_report(
        mpc_per_cell=8.5,
        myr_per_sweep=t * 1.01,
        myr_per_drop=100.0,
        j_per_grain=1e60,
        material_speed_cells_per_sweep=1.5,
    )
    assert r["violations"] == []
    assert r["material_speed_over_c"] < 1.0
