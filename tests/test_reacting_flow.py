"""Tests for the reacting-flow detonation (motion + self-feeding front unified).

Cornerstone: energy conserved (E_fluid + fuel = const, periodic). And the key
honesty check: with a SMALL ignition the front is genuinely FUEL-DRIVEN -- it
fizzles without enough fuel and self-sustains above a critical fuel level. So
it is a real detonation, not the ignition blast diffusing.
"""

from void_cascade.reacting_flow import run_detonation


def test_energy_conserved():
    r = run_detonation(N=300, fuel0=3.0, e_ign=2.5, ignite_E=3.0, steps=300)
    assert abs(r["energy_residual"]) < 1e-6


def test_pressure_stays_positive():
    r = run_detonation(N=300, fuel0=3.0, e_ign=2.5, ignite_E=3.0, steps=300)
    assert r["min_pressure"] > 0.0


def test_detonation_propagates_with_enough_fuel():
    r = run_detonation(N=400, fuel0=3.0, e_ign=2.5, ignite_E=3.0, steps=400)
    assert r["propagated"]
    assert r["final_burned_fraction"] > 0.5


def test_fizzles_without_fuel():
    # small ignition + no fuel -> the blast dies; not a diffusion artifact.
    r = run_detonation(N=400, fuel0=0.0, e_ign=2.5, ignite_E=3.0, steps=400)
    assert not r["propagated"]
    assert r["final_burned_fraction"] < 0.05


def test_fuel_driven_critical_threshold():
    # below critical fuel -> fizzle; above -> self-sustaining detonation.
    low = run_detonation(N=400, fuel0=0.5, e_ign=2.5, ignite_E=3.0, steps=400)
    high = run_detonation(N=400, fuel0=2.0, e_ign=2.5, ignite_E=3.0, steps=400)
    assert low["final_burned_fraction"] < 0.3
    assert high["final_burned_fraction"] > 0.8
