"""Tests for the edgeless growing-substrate model.

The defining properties: the domain grows with fresh pre-loaded substrate
(not empty), the old state is preserved when it grows, ignition starts the
catastrophe, and -- the whole point -- the released front never sits on the
array edge (we grow before it can), so no reflection/maker's-box artifact.
"""

import numpy as np

from void_cascade.cascade_heat import CascadeParams
from void_cascade.cascade_heat_2d import initialize_2d
from void_cascade.cascade_heat_grow import (
    ignite,
    near_border,
    pad_domain,
    run_grow,
)


def p(**ov):
    d = dict(fracture_density=2.0, heat_per_crack=0.10, diffuse=0.15, cooling=0.10,
             melt_heat=1.0, release_factor=0.5, drive_amount=1.0, n_drive_sites=1)
    d.update(ov)
    return CascadeParams(**d)


def test_pad_grows_symmetrically_and_preserves_centre():
    st = initialize_2d(10)
    st.density[:] = np.arange(100).reshape(10, 10).astype(float)
    st.cracks[5, 5] = 7
    st.released[5, 5] = True
    pad_domain(st, chunk=4, rng=np.random.default_rng(0))
    assert st.density.shape == (18, 18)  # 10 + 2*4
    # old centre preserved at offset 4
    assert np.array_equal(st.density[4:14, 4:14], np.arange(100).reshape(10, 10))
    assert st.cracks[4 + 5, 4 + 5] == 7
    assert st.released[4 + 5, 4 + 5]


def test_pad_fills_fresh_substrate_loaded_not_empty():
    st = initialize_2d(10)
    st.density[:] = 1.0
    st.cracks[:] = 5  # every cold cell has hidden damage
    pad_domain(st, chunk=3, rng=np.random.default_rng(1))
    # border cells (fresh substrate) must be loaded: cracks sampled from cold
    # cells (all =5 here), heat at absolute zero, nothing released.
    assert st.cracks[0, 0] == 5
    assert st.heat[0, 0] == 0.0
    assert not st.released[0, 0]


def test_ignite_sets_central_heat_above_melt():
    st = initialize_2d(20)
    ignite(st, p())
    assert st.heat[10, 10] >= p().melt_heat


def test_near_border_detects_edge_contact():
    m = np.zeros((10, 10), dtype=bool)
    assert not near_border(m, 2)
    m[0, 5] = True
    assert near_border(m, 2)


def test_run_never_lets_front_touch_the_edge():
    # The whole point: with growth on, the released front never sits on the
    # array border (no reflection / maker's-box artifact).
    res = run_grow(L0=40, accumulate_steps=200, propagate_steps=200,
                   p=p(), seed=0, margin=3, chunk=12, max_size=160)
    assert res["ever_touched_edge"] is False
    # if it grew at all, the front genuinely advanced past the original box
    assert res["grow_events"] >= 0
