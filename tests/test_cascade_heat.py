"""Tests for the damage-coupled (scattered -> cascading) heat model.

Fracturing uses the validated Manna rule (threshold 2, sheds 2 grains to
random neighbours), so a single received grain only fractures a cell already
near threshold -> bounded, self-organising avalanches. Early (sparse)
avalanches are tiny and scattered; as the line fills toward critical they
grow and cluster, concentrating the heat.

Key finding this version exposes: with clustering, the dial that controls
tipping is the HEAT-MAKING strength, not cooling -- concentrated avalanche
heat can outrun even heavy cooling. So the honest guardrail here is "no heat
source -> stays cold" (tipping must come from heat, not bookkeeping), and "a
real heat source -> runs away, even under heavy cooling."

Sizes kept small (L=30): heat-on cases tip then go quiet (fast); heat-off
cases sit at criticality (costlier), so they use few steps -- enough to show
the early scattered phase grow into cascades.
"""

import numpy as np

from void_cascade.cascade_heat import CascadeParams, run


def base(**overrides) -> CascadeParams:
    p = dict(
        fracture_density=2.0,
        heat_per_crack=0.15,
        diffuse=0.2,
        cooling=0.05,
        melt_heat=1.0,
        release_factor=0.5,
        drive_amount=1.0,
        n_drive_sites=1,
    )
    p.update(overrides)
    return CascadeParams(**p)


def test_heat_never_below_absolute_zero():
    r = run(L=30, n_steps=1500, p=base(), seed=0)
    assert r["min_heat_overall"] >= 0.0


def test_stays_cold_with_no_heat_source():
    # Guardrail: tipping must come from heat. No heat made -> never tips,
    # even though the avalanches still happen.
    r = run(L=30, n_steps=300, p=base(heat_per_crack=0.0), seed=0)
    assert r["ran_away"] is False
    assert r["peak_heat_overall"] == 0.0


def test_runs_away_with_a_real_heat_source():
    r = run(L=30, n_steps=1500, p=base(cooling=0.05), seed=0)
    assert r["ran_away"] is True


def test_clustered_heat_tips_even_under_heavy_cooling():
    # The headline: avalanche clustering concentrates heat enough that even
    # strong cooling cannot stop it (the scattered model needed zero cooling).
    r = run(L=30, n_steps=1500, p=base(cooling=0.3), seed=0)
    assert r["ran_away"] is True


def test_fracturing_goes_scattered_then_cascading():
    # Heat off so it stays cold -> we see the pure avalanche dynamics grow
    # from the sparse early lattice into clustered cascades.
    r = run(L=30, n_steps=300, p=base(heat_per_crack=0.0), seed=0)
    cs = np.asarray(r["cascade_sizes"])
    early = cs[: len(cs) // 10].mean()   # sparse start
    late = cs[len(cs) // 2:].mean()      # filled / critical
    assert r["overall_max_cascade"] > 5  # real cascades, not all size-1
    assert late > early


def test_cracks_accumulate_repeatedly():
    r = run(L=30, n_steps=300, p=base(heat_per_crack=0.0), seed=0)
    assert r["max_cracks_on_a_cell"] > 1
