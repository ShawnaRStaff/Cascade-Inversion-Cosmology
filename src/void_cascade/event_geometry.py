"""Per-event mask geometry: is a single avalanche confined to the carrier network?

M6 established that the peak avalanche touches a growing fraction of the
lattice (36% -> 43% unique cells for L=48 -> 128) but that the L -> infinity
asymptote cannot be decided by comparing mean peak sizes: the candidate
models (capped at the 62% carrier fraction, capped near 80-90%, or full
coverage) predict means that differ by less than the seed-to-seed scatter
at any affordable L.

The structural discriminator, instead of a statistical one: classify every
cell the avalanche topples by its state immediately BEFORE the event.
In a stable Manna configuration each cell holds z = 0 ("sink") or z = 1
("carrier"). A carrier topples after receiving one grain; a sink must
receive two grains mid-event before it can fire. If peak events are
confined to the carrier network (sink participation ~ 0, or shrinking
with L), the 62%-carrier-cap asymptote is the physical picture. If sink
participation is substantial and grows with event size, the cap story
dies and the asymptote lies above the carrier fraction.

We record the classification for every record-breaking event (each new
maximum of unique cells touched), which gives the size-progression of
sink participation within a single run at negligible cost.
"""

from __future__ import annotations

import numpy as np

from void_cascade.sandpile_3d import drive, initialize, relax


def classify_event(carriers_pre: np.ndarray, mask: np.ndarray) -> dict:
    """Classify an avalanche's toppled cells by pre-event carrier/sink state.

    Parameters
    ----------
    carriers_pre : (L, L, L) bool
        True where z == 1 immediately before the event (the stable
        configuration before the drive grain landed). Every False cell
        was a sink (z == 0); stable Manna configurations have no other
        values.
    mask : (L, L, L) bool
        True for every cell that toppled at least once during the event.

    Returns
    -------
    dict with integer counts (unique, toppled_carriers, toppled_sinks)
    and ratios (sink_participation = toppled_sinks / unique,
    carrier_coverage = toppled_carriers / total carriers). Ratios are
    0.0 when their denominator is zero.
    """
    if carriers_pre.shape != mask.shape:
        raise ValueError(
            f"shape mismatch: carriers_pre {carriers_pre.shape} vs mask {mask.shape}"
        )
    unique = int(mask.sum())
    toppled_carriers = int((mask & carriers_pre).sum())
    toppled_sinks = unique - toppled_carriers
    n_carriers = int(carriers_pre.sum())
    return {
        "unique": unique,
        "toppled_carriers": toppled_carriers,
        "toppled_sinks": toppled_sinks,
        "sink_participation": toppled_sinks / unique if unique else 0.0,
        "carrier_coverage": toppled_carriers / n_carriers if n_carriers else 0.0,
    }


def run_peak_event_geometry(L: int, n_drops: int, seed: int | None = None) -> dict:
    """Drive an empty lattice, recording geometry of record-breaking events.

    Identical dynamics and RNG stream as run_with_ever_toppled (drive one
    grain, relax with support tracking), so a given seed reproduces the
    exact avalanche trajectory of the uninstrumented M6 runs. The only
    additions are a per-drop boolean snapshot of the carrier field taken
    before the drive grain lands, and a classification whenever an event
    sets a new record for unique cells touched.

    Returns
    -------
    dict with:
      sizes, unique_sizes : (n_drops,) int64 arrays
      records : dict of parallel arrays, one entry per record-breaking
        event — drop, unique, toppled_carriers, toppled_sinks,
        sink_participation, carrier_coverage, n_carriers_pre
      peak_mask, peak_carriers_pre : (L, L, L) bool for the peak event
      final_z : (L, L, L) int64; grains_lost : int
    """
    rng = np.random.default_rng(seed)
    state = initialize(L)
    sizes = np.zeros(n_drops, dtype=np.int64)
    unique_sizes = np.zeros(n_drops, dtype=np.int64)

    rec_drop: list[int] = []
    rec_geom: list[dict] = []
    rec_ncar: list[int] = []
    peak_unique = 0
    peak_mask = np.zeros((L, L, L), dtype=bool)
    peak_carriers_pre = np.zeros((L, L, L), dtype=bool)

    for t in range(n_drops):
        # Stable configurations hold only z in {0, 1}: carriers are z == 1.
        carriers_pre = state.z == 1
        drive(state, rng)
        s, _, mask = relax(state, rng, track_support=True)
        assert mask is not None  # track_support=True always yields a mask
        sizes[t] = s
        u = int(mask.sum())
        unique_sizes[t] = u
        if u > peak_unique:
            peak_unique = u
            peak_mask = mask
            peak_carriers_pre = carriers_pre
            rec_drop.append(t)
            rec_geom.append(classify_event(carriers_pre, mask))
            rec_ncar.append(int(carriers_pre.sum()))

    records = {
        "drop": np.array(rec_drop, dtype=np.int64),
        "unique": np.array([g["unique"] for g in rec_geom], dtype=np.int64),
        "toppled_carriers": np.array(
            [g["toppled_carriers"] for g in rec_geom], dtype=np.int64
        ),
        "toppled_sinks": np.array([g["toppled_sinks"] for g in rec_geom], dtype=np.int64),
        "sink_participation": np.array(
            [g["sink_participation"] for g in rec_geom], dtype=np.float64
        ),
        "carrier_coverage": np.array(
            [g["carrier_coverage"] for g in rec_geom], dtype=np.float64
        ),
        "n_carriers_pre": np.array(rec_ncar, dtype=np.int64),
    }
    return {
        "sizes": sizes,
        "unique_sizes": unique_sizes,
        "records": records,
        "peak_mask": peak_mask,
        "peak_carriers_pre": peak_carriers_pre,
        "final_z": state.z,
        "grains_lost": state.grains_lost,
    }
