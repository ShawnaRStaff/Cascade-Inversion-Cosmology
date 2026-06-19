"""Tests for damage geometries used in the substrate-resilience scrutiny.

Two ways to remove a fraction of the substrate:
- random_damage_mask: cells destroyed uniformly at random (what the old
  resilience experiment did).
- connected_damage_mask: a single connected blob of destroyed cells (the
  proxy for a real crack/event, which spreads through neighbors).

We test the math: exact counts, reproducibility, and the defining
property that the connected damage really is one 6-connected piece while
random damage is not.
"""

import numpy as np
from scipy import ndimage

from void_cascade.damage import (
    connected_damage_mask,
    damage_count,
    random_damage_mask,
)

SHAPE = (12, 12, 12)
_STRUCT_6 = ndimage.generate_binary_structure(rank=3, connectivity=1)


def test_damage_count_is_floor_of_fraction_times_volume():
    assert damage_count((10, 10, 10), 0.0) == 0
    assert damage_count((10, 10, 10), 0.5) == 500
    assert damage_count((10, 10, 10), 1.0) == 1000


def test_random_damage_mask_has_exact_count():
    m = random_damage_mask(SHAPE, 0.25, np.random.default_rng(0))
    assert m.shape == SHAPE
    assert m.dtype == bool
    assert int(m.sum()) == damage_count(SHAPE, 0.25)


def test_random_damage_zero_fraction_is_empty():
    m = random_damage_mask(SHAPE, 0.0, np.random.default_rng(1))
    assert not m.any()


def test_random_damage_reproducible_with_seed():
    a = random_damage_mask(SHAPE, 0.3, np.random.default_rng(7))
    b = random_damage_mask(SHAPE, 0.3, np.random.default_rng(7))
    assert np.array_equal(a, b)


def test_connected_damage_mask_has_exact_count():
    m = connected_damage_mask(SHAPE, 0.25, np.random.default_rng(0))
    assert m.shape == SHAPE
    assert m.dtype == bool
    assert int(m.sum()) == damage_count(SHAPE, 0.25)


def test_connected_damage_is_one_6connected_blob():
    m = connected_damage_mask(SHAPE, 0.3, np.random.default_rng(3))
    _labels, n_components = ndimage.label(m, structure=_STRUCT_6)
    assert n_components == 1


def test_connected_damage_reproducible_with_seed():
    a = connected_damage_mask(SHAPE, 0.4, np.random.default_rng(9))
    b = connected_damage_mask(SHAPE, 0.4, np.random.default_rng(9))
    assert np.array_equal(a, b)


def test_random_damage_is_more_fragmented_than_connected():
    # At the same moderate fraction, scattered damage lands in many
    # disconnected pieces; a grown blob stays a single piece. This is the
    # whole point of separating the two geometries.
    r = random_damage_mask(SHAPE, 0.3, np.random.default_rng(11))
    c = connected_damage_mask(SHAPE, 0.3, np.random.default_rng(11))
    _r_labels, n_r = ndimage.label(r, structure=_STRUCT_6)
    _c_labels, n_c = ndimage.label(c, structure=_STRUCT_6)
    assert n_c == 1
    assert n_r > 1
