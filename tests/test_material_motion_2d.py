"""Tests for the 2D radial implosion.

Cornerstone: mass and energy conserved (Lax-Friedrichs, periodic). Plus:
no vacuum blow-up, and the radial implosion compresses + heats the centre to
a peak and then rebounds.
"""

import numpy as np

from void_cascade.material_motion_2d import (
    internal_energy,
    primitives,
    radial_converging_ic,
    run_radial_implosion,
)


def test_mass_and_energy_conserved():
    r = run_radial_implosion(L=80, u0=1.5, R0=25, steps=200)
    assert abs(r["mass_residual"]) < 1e-6
    assert abs(r["energy_residual"]) < 1e-6


def test_momentum_stays_near_zero_by_symmetry():
    r = run_radial_implosion(L=80, u0=1.5, R0=25, steps=200)
    assert r["momentum_residual"] < 1e-6


def test_pressure_and_density_stay_positive():
    r = run_radial_implosion(L=80, u0=1.5, R0=25, steps=200)
    assert r["min_density"] > 0.0
    assert r["min_pressure"] > 0.0


def test_radial_implosion_compresses_heats_then_rebounds():
    # LF is diffusive, so density compression is modest (~1.3-1.6x); the clear
    # signature is the heat (kinetic->internal) spike, which is strong.
    r = run_radial_implosion(L=100, u0=2.5, R0=35, steps=300, band=3)
    assert r["peak_central_density"] > 1.2 * r["central_density_initial"]  # real compression
    assert r["central_heat_at_peak"] > 1.5 * r["central_heat_initial"]     # plasma-like heating
    assert r["rebounded"]


def test_internal_energy_is_total_minus_kinetic():
    rho, mx, my, E = radial_converging_ic(20, u0=1.0, R0=8)
    ie = internal_energy(rho, mx, my, E)
    u, v, _ = primitives(rho, mx, my, E)
    assert np.allclose(ie, E - 0.5 * rho * (u * u + v * v))
