"""Tests for the 2D damage-coupled heat model.

2D Manna (threshold 2, sheds 2 grains to random neighbours of 4) is the SOC
class this project validated in M2 (tau~1.27) -- so the 2D version both
grounds the foundation and gives trustworthy heat numbers (1D Manna was
degenerate). Same honest guardrails: no heat source -> stays cold; a real
heat source -> runs away.

Sizes tiny (L=16): heat-off sits at criticality (costlier), so few steps.
"""

import numpy as np

from void_cascade.cascade_heat import CascadeParams
from void_cascade.cascade_heat_2d import run_2d


def base(**overrides) -> CascadeParams:
    p = dict(
        fracture_density=2.0,
        heat_per_crack=0.15,
        diffuse=0.15,
        cooling=0.05,
        melt_heat=1.0,
        release_factor=0.5,
        drive_amount=1.0,
        n_drive_sites=1,
    )
    p.update(overrides)
    return CascadeParams(**p)


def test_heat_never_below_absolute_zero():
    r = run_2d(L=16, n_steps=800, p=base(), seed=0)
    assert r["min_heat_overall"] >= 0.0


def test_stays_cold_with_no_heat_source():
    r = run_2d(L=16, n_steps=400, p=base(heat_per_crack=0.0), seed=0)
    assert r["ran_away"] is False
    assert r["peak_heat_overall"] == 0.0


def test_runs_away_with_a_real_heat_source():
    r = run_2d(L=16, n_steps=800, p=base(cooling=0.05), seed=0)
    assert r["ran_away"] is True


def test_cascades_happen():
    r = run_2d(L=16, n_steps=400, p=base(heat_per_crack=0.0), seed=0)
    assert r["overall_max_cascade"] > 5


def test_cracks_accumulate_repeatedly():
    r = run_2d(L=16, n_steps=400, p=base(heat_per_crack=0.0), seed=0)
    assert r["max_cracks_on_a_cell"] > 1
