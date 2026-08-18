"""Tests for per-event mask geometry (M7 asymptote discriminator).

The physics question: is a single avalanche confined to the carrier
network (cells at z=1 before the event), or does it recruit sink
cells (z=0) mid-event? Confinement supports the 62%-carrier-cap
asymptote; growing sink participation kills it.

Contract: classify_event(carriers_pre, mask) -> flat counts;
run_peak_event_geometry(L, n_drops, seed) -> record-breaking event
series + peak-event geometry, deterministic given seed.
"""

import numpy as np
import pytest

from void_cascade.event_geometry import classify_event, run_peak_event_geometry


def test_classify_event_counts_carriers_and_sinks_exactly():
    # 2x2x2 lattice: mark two cells as carriers, rest sinks.
    carriers_pre = np.zeros((2, 2, 2), dtype=bool)
    carriers_pre[0, 0, 0] = True
    carriers_pre[1, 1, 1] = True
    # Event touched one carrier and one sink.
    mask = np.zeros((2, 2, 2), dtype=bool)
    mask[0, 0, 0] = True   # carrier
    mask[0, 1, 1] = True   # sink
    g = classify_event(carriers_pre, mask)
    assert g["unique"] == 2
    assert g["toppled_carriers"] == 1
    assert g["toppled_sinks"] == 1
    assert g["sink_participation"] == pytest.approx(0.5)
    assert g["carrier_coverage"] == pytest.approx(0.5)  # 1 of 2 carriers


def test_classify_event_empty_mask_gives_zero_ratios():
    carriers_pre = np.ones((2, 2, 2), dtype=bool)
    mask = np.zeros((2, 2, 2), dtype=bool)
    g = classify_event(carriers_pre, mask)
    assert g["unique"] == 0
    assert g["toppled_carriers"] == 0
    assert g["toppled_sinks"] == 0
    assert g["sink_participation"] == 0.0
    assert g["carrier_coverage"] == 0.0


def test_classify_event_partition_is_exact():
    # toppled_carriers + toppled_sinks must equal unique for any inputs.
    rng = np.random.default_rng(7)
    carriers_pre = rng.random((5, 5, 5)) < 0.6
    mask = rng.random((5, 5, 5)) < 0.3
    g = classify_event(carriers_pre, mask)
    assert g["toppled_carriers"] + g["toppled_sinks"] == g["unique"]
    assert g["unique"] == int(mask.sum())


def test_classify_event_shape_mismatch_raises():
    with pytest.raises(ValueError):
        classify_event(np.zeros((2, 2, 2), dtype=bool), np.zeros((3, 3, 3), dtype=bool))


def test_run_records_are_increasing_and_consistent():
    out = run_peak_event_geometry(L=8, n_drops=2000, seed=42)
    rec = out["records"]
    # At least one avalanche must have occurred in 2000 drops on L=8.
    assert len(rec["drop"]) >= 1
    uniq = rec["unique"]
    # Record-breaking series: strictly increasing unique size, increasing drop index.
    assert np.all(np.diff(uniq) > 0)
    assert np.all(np.diff(rec["drop"]) > 0)
    # Partition holds for every record.
    assert np.all(rec["toppled_carriers"] + rec["toppled_sinks"] == uniq)
    # The last record IS the peak of the whole run.
    assert uniq[-1] == out["unique_sizes"].max()
    assert rec["drop"][-1] == int(np.argmax(out["unique_sizes"]))


def test_run_peak_mask_matches_reported_counts():
    out = run_peak_event_geometry(L=8, n_drops=2000, seed=42)
    g = classify_event(out["peak_carriers_pre"], out["peak_mask"])
    rec = out["records"]
    assert g["unique"] == rec["unique"][-1]
    assert g["toppled_carriers"] == rec["toppled_carriers"][-1]
    assert g["toppled_sinks"] == rec["toppled_sinks"][-1]


def test_run_is_deterministic_given_seed():
    a = run_peak_event_geometry(L=8, n_drops=1000, seed=3)
    b = run_peak_event_geometry(L=8, n_drops=1000, seed=3)
    assert np.array_equal(a["unique_sizes"], b["unique_sizes"])
    assert np.array_equal(a["records"]["toppled_sinks"], b["records"]["toppled_sinks"])
    assert np.array_equal(a["peak_mask"], b["peak_mask"])


def test_run_matches_plain_run_trajectory():
    # Instrumentation must not perturb the dynamics: same seed must give
    # the same avalanche sizes as the uninstrumented driver.
    from void_cascade.sandpile_3d import run_with_ever_toppled

    out = run_peak_event_geometry(L=8, n_drops=500, seed=11)
    _, sizes, _, _ = run_with_ever_toppled(L=8, n_drops=500, seed=11)
    assert np.array_equal(out["sizes"], sizes)
