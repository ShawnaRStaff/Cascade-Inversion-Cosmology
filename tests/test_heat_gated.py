"""Tests for the heat-gated substrate model (cold-frozen -> heat runaway).

The model must support BOTH outcomes honestly:
- with strong cooling it stays cold forever (never runs away) -- the guardrail,
- with no cooling it eventually runs away -- proof the mechanism works.
Plus the defining rules: heat never goes below absolute zero, cracks pile up
while frozen, and a cell only releases its cracks once it's hot enough.
"""

import numpy as np

from void_cascade.heat_gated import (
    HeatParams,
    HeatState,
    initialize,
    release,
    run,
)


def base_params(**overrides) -> HeatParams:
    p = dict(
        fracture_density=1.0,
        fracture_relief=1.0,
        drive_rate=0.02,
        heat_per_crack=0.15,
        diffuse=0.2,
        cooling=0.02,
        melt_heat=1.0,
        release_factor=0.5,
    )
    p.update(overrides)
    return HeatParams(**p)


def test_heat_never_below_absolute_zero():
    res = run(L=200, n_steps=2000, p=base_params(), seed=0)
    assert res["min_heat_overall"] >= 0.0


def test_strong_cooling_stays_cold():
    # The guardrail: enough cooling and it must NEVER run away.
    res = run(L=200, n_steps=5000, p=base_params(cooling=0.5), seed=0)
    assert res["ran_away"] is False
    assert res["peak_heat_overall"] < base_params().melt_heat


def test_no_cooling_runs_away():
    # The mechanism must be ABLE to run away, or it's a dead model.
    res = run(L=200, n_steps=5000, p=base_params(cooling=0.0), seed=0)
    assert res["ran_away"] is True


def test_cracks_accumulate_while_frozen():
    # Strong cooling keeps it frozen: cracks pile up (repeated fracture) but
    # nothing releases -- hidden damage that "looks intact".
    res = run(L=200, n_steps=5000, p=base_params(cooling=0.5), seed=0)
    assert res["max_cracks_on_a_cell"] > 1
    assert res["fraction_released_final"] == 0.0


def test_release_requires_melt_heat():
    L = 5
    state = HeatState(
        density=np.zeros(L),
        cracks=np.array([0, 3, 0, 0, 0], dtype=np.int64),
        heat=np.zeros(L),
        released=np.zeros(L, dtype=bool),
    )
    p = base_params(melt_heat=1.0, release_factor=0.5)

    # Below melt: nothing releases.
    n, freed = release(state, p)
    assert n == 0 and freed == 0
    assert state.cracks[1] == 3

    # Push that cell over the melt point: it releases its cracks once.
    state.heat[1] = 1.5
    n, freed = release(state, p)
    assert n == 1 and freed == 3
    assert state.cracks[1] == 0
    assert state.released[1]
    assert state.heat[1] > 1.5  # got a heat burst
