"""Tests for conserved energy on the edgeless growing substrate.

Cornerstone: TRUE conservation -- residual ~ 0 AND boundary_lost ~ 0 (no
cooling sink, no boundary sink). Plus: padding fresh substrate accounts its
energy, and the front sustains and grows the domain without ever touching an
edge.
"""

import numpy as np

from void_cascade.conserved_energy import EnergyParams, conservation_residual
from void_cascade.conserved_energy_grow import (
    initialize_2d, pad_energy_domain, run_grow_energy,
)


def p(**ov):
    d = dict(flip_threshold=1.0, release_fraction=1.0, diffuse=0.5,
             drive_amount=1.0, n_drive_sites=1)
    d.update(ov)
    return EnergyParams(**d)


def test_pad_accounts_energy_keeps_conservation():
    st = initialize_2d(8)
    st.potential[:] = 0.8
    st.energy_in = float(st.potential.sum())
    assert abs(conservation_residual(st)) < 1e-9
    pad_energy_domain(st, chunk=4, rng=np.random.default_rng(0))
    # after padding fresh loaded substrate, conservation still exact
    assert abs(conservation_residual(st)) < 1e-9
    assert st.potential.shape == (16, 16)


def test_true_conservation_no_sink():
    r = run_grow_energy(L0=40, propagate_steps=300, p=p(), seed=0,
                        margin=5, chunk=16, max_size=200)
    assert abs(r["conservation_residual"]) < 1e-6   # energy conserved
    assert r["boundary_lost"] < 1e-6                # AND nothing leaked: true conservation


def test_front_sustains_and_grows_edge_free():
    r = run_grow_energy(L0=40, propagate_steps=400, p=p(), seed=1,
                        margin=5, chunk=16, max_size=220)
    assert r["grow_events"] >= 1                 # the front advanced past the start box
    assert r["final_flipped_fraction"] > 0.2     # it actually swept substrate
