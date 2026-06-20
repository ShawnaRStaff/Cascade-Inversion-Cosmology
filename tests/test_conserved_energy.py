"""Tests for the conserved-energy front.

Cornerstone: energy is conserved -- energy_in == potential + kinetic +
boundary_lost, always. Plus: combustion converts potential to kinetic (not
invents it), kinetic stays non-negative, and the front can BOTH sustain (when
the substrate is loaded near threshold) and fizzle (when it isn't) -- so it is
not rigged.
"""

import numpy as np

from void_cascade.conserved_energy import (
    EnergyParams,
    EnergyState,
    combust,
    conservation_residual,
    diffuse_kinetic,
    drive,
    ignite,
    initialize_2d,
    step,
)


def p(**ov):
    d = dict(flip_threshold=1.0, release_fraction=1.0, diffuse=0.5,
             drive_amount=1.0, n_drive_sites=3)
    d.update(ov)
    return EnergyParams(**d)


def test_energy_conserved_over_a_run():
    st = initialize_2d(40)
    rng = np.random.default_rng(0)
    for _ in range(500):
        step(st, p(), rng)
    assert abs(conservation_residual(st)) < 1e-6


def test_energy_conserved_with_ignition():
    st = initialize_2d(40)
    rng = np.random.default_rng(0)
    for _ in range(100):
        step(st, p(), rng)
    ignite(st, p())
    for _ in range(200):
        step(st, p(), rng)
    assert abs(conservation_residual(st)) < 1e-6


def test_combust_converts_potential_to_kinetic_not_invents():
    st = initialize_2d(3)
    st.potential[1, 1] = 1.5  # over threshold 1.0
    before = st.potential.sum() + st.kinetic.sum()
    n = combust(st, p(release_fraction=1.0))
    after = st.potential.sum() + st.kinetic.sum()
    assert n == 1
    assert st.flipped[1, 1]
    assert st.kinetic[1, 1] > 0 and st.potential[1, 1] == 0.0
    assert abs(after - before) < 1e-9  # total energy unchanged by combustion


def test_kinetic_never_negative():
    st = initialize_2d(30)
    rng = np.random.default_rng(1)
    for _ in range(300):
        step(st, p(), rng)
    assert st.kinetic.min() >= 0.0


def test_diffusion_conserves_minus_boundary():
    st = initialize_2d(10)
    st.kinetic[5, 5] = 4.0
    before = st.kinetic.sum()
    diffuse_kinetic(st, p(diffuse=0.5))
    after = st.kinetic.sum()
    assert abs(before - (after + st.boundary_lost)) < 1e-9


def test_front_fizzles_when_substrate_not_loaded():
    # Empty substrate (no potential): ignite -> nothing to prime -> no spread.
    st = initialize_2d(40)
    rng = np.random.default_rng(2)
    ignite(st, p())
    for _ in range(150):
        # no drive: nothing reloads; front has no fuel ahead
        combust(st, p()); diffuse_kinetic(st, p())
    assert st.flipped.mean() < 0.05  # stayed local, fizzled


def test_front_sustains_when_substrate_loaded():
    # Pre-load potential just below threshold everywhere -> ignition primes
    # neighbours over the line -> front spreads.
    st = initialize_2d(40)
    st.potential[:] = 0.95  # threshold 1.0; a nudge tips it
    rng = np.random.default_rng(3)
    ignite(st, p())
    for _ in range(120):
        combust(st, p()); diffuse_kinetic(st, p())
    assert st.flipped.mean() > 0.5  # the front swept the loaded substrate
