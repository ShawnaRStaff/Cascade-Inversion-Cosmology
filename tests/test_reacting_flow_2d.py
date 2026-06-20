"""Tests for the end-to-end 2D implosion -> detonation (stages 2+3 in one run).

Cornerstone: energy conserved (E_fluid + fuel = const). And the end-to-end
claim: the IMPLOSION causes the ignition -- with a converging flow it lights
the detonation and the front sweeps the fuel; with no implosion (same fuel),
nothing heats it, so it never ignites.
"""

from void_cascade.reacting_flow_2d import run_end_to_end


def test_energy_conserved():
    r = run_end_to_end(L=90, u0=2.5, R0=28, fuel0=3.0, e_ign=2.5, steps=200)
    assert abs(r["energy_residual"]) < 1e-6


def test_implosion_ignites_detonation_end_to_end():
    r = run_end_to_end(L=90, u0=2.5, R0=28, fuel0=3.0, e_ign=2.5, steps=250)
    assert r["ignite_step"] is not None        # the implosion lit it
    assert r["fuel_burned_fraction"] > 0.5     # detonation swept the fuel
    assert r["max_burn_radius"] > 10           # a front actually expanded


def test_no_implosion_no_ignition():
    # Same loaded fuel, but no converging flow -> nothing heats it -> no
    # ignition. Proves the implosion *causes* the expansion.
    r = run_end_to_end(L=90, u0=0.0, R0=28, fuel0=3.0, e_ign=2.5, steps=250)
    assert r["ignite_step"] is None
    assert r["fuel_burned_fraction"] < 0.05
