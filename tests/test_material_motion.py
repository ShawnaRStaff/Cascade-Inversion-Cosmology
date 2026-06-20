"""Tests for the material-motion (compressible-fluid implosion) model.

Cornerstone: mass and energy conserved (Lax-Friedrichs is conservative; periodic
BC). Plus: pressure/density stay positive (no vacuum blow-up), and the implosion
genuinely compresses + heats the centre and then rebounds.
"""

import numpy as np

from void_cascade.material_motion import (
    converging_initial_conditions,
    internal_energy,
    lax_friedrichs_step,
    primitives,
    run_implosion,
)


def test_mass_and_energy_conserved():
    r = run_implosion(N=120, u0=1.0, steps=200)
    assert abs(r["mass_residual"]) < 1e-6
    assert abs(r["energy_residual"]) < 1e-6


def test_momentum_stays_zero_by_symmetry():
    r = run_implosion(N=120, u0=1.0, steps=200)
    assert abs(r["momentum_residual"]) < 1e-6


def test_pressure_and_density_stay_positive():
    r = run_implosion(N=120, u0=1.0, steps=200)
    assert r["min_density"] > 0.0
    assert r["min_pressure"] > 0.0


def test_implosion_compresses_heats_then_rebounds():
    r = run_implosion(N=160, u0=1.2, steps=300)
    # compresses: centre density rises above its start
    assert r["peak_central_density"] > 1.5 * r["central_density_initial"]
    # heats: internal energy at the compression peak exceeds the start
    assert r["central_heat_at_peak"] > r["central_heat_initial"]
    # rebounds: density falls back after the peak
    assert r["rebounded"]


def test_internal_energy_is_total_minus_kinetic():
    rho, mom, E = converging_initial_conditions(20, u0=1.0)
    ie = internal_energy(rho, mom, E)
    u, _ = primitives(rho, mom, E)
    assert np.allclose(ie, E - 0.5 * rho * u * u)
