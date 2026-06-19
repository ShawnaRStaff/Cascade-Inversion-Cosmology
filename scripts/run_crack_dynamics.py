"""Experiment 2: does a crack heal, or feed itself?

Carve one connected crack-shaped hole into the saturated substrate (set
those cells empty), then let the EXISTING dynamics run with input
continuing. No new physics, no "pull" added. We just watch what the
current rules do to the hole over time:

  - heal  = the crack refills, the biggest empty region shrinks back
            toward the normal background of small voids.
  - feed  = the empty region grows beyond the original crack, or the
            emptiness spreads across the box on its own.

We also do a no-input check first: with input off, relax once and count
spontaneous topples. (Expected zero — emptying cells can never push a
cell over threshold. Recorded so the claim is measured, not assumed.)
"""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.damage import connected_damage_mask  # noqa: E402
from void_cascade.percolation import check_spanning  # noqa: E402
from void_cascade.sandpile_3d import MannaState3D, drive, relax  # noqa: E402

CRACK_FRACTION = 0.30
N_DROPS = 60_000
SAMPLE_EVERY = 500
_STRUCT_6 = ndimage.generate_binary_structure(rank=3, connectivity=1)


def load_saturated_field() -> np.ndarray:
    path = REPO_ROOT / "data" / "outputs" / "fss_sweep_20260521_031056" / "L48_s14800_final.npz"
    return np.asarray(np.load(path, allow_pickle=True)["final_z"]).astype(np.int64)


def largest_empty_region(z: np.ndarray) -> int:
    """Size of the biggest connected (face-neighbour) blob of empty cells."""
    labels, n = ndimage.label(z == 0, structure=_STRUCT_6)
    if n == 0:
        return 0
    return int(np.bincount(labels.ravel())[1:].max())


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"crack_dynamics_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)

    z0 = load_saturated_field()
    shape = z0.shape
    print(f"Loaded L={shape[0]} field, z_avg={z0.mean():.4f}")
    print(f"Background biggest empty region (before any crack): {largest_empty_region(z0)} cells")

    crack = connected_damage_mask(shape, CRACK_FRACTION, np.random.default_rng(0))
    crack_size = int(crack.sum())
    print(f"Carving a connected crack: {crack_size} cells ({CRACK_FRACTION:.0%})")

    # --- No-input check ---
    z = z0.copy()
    z[crack] = 0
    state_off = MannaState3D(z=z.copy(), grains_lost=0)
    s_off, _, _ = relax(state_off, np.random.default_rng(1), track_support=False)
    print(f"[input OFF] spontaneous topples after carving: {s_off}")

    # --- Input ON: watch the hole over time ---
    state = MannaState3D(z=z.copy(), grains_lost=0)
    rng = np.random.default_rng(2)
    crack_flat = crack.ravel()

    drops_axis, crack_refill, biggest_empty, empty_total, empty_spans = [], [], [], [], []

    def sample(t: int) -> None:
        zf = state.z.ravel()
        drops_axis.append(t)
        crack_refill.append(float((zf[crack_flat] >= 1).mean()))
        biggest_empty.append(largest_empty_region(state.z))
        empty_total.append(int((state.z == 0).sum()))
        empty_spans.append(bool(check_spanning(state.z == 0).percolates))

    sample(0)
    for t in range(1, N_DROPS + 1):
        drive(state, rng)
        relax(state, rng, track_support=False)
        if t % SAMPLE_EVERY == 0:
            sample(t)

    summary = {
        "crack_fraction": CRACK_FRACTION,
        "crack_size_cells": crack_size,
        "n_drops": N_DROPS,
        "input_off_spontaneous_topples": int(s_off),
        "background_biggest_empty_before_crack": largest_empty_region(z0),
        "crack_refill_start": crack_refill[0],
        "crack_refill_end": crack_refill[-1],
        "biggest_empty_start": biggest_empty[0],
        "biggest_empty_end": biggest_empty[-1],
        "empty_spans_start": empty_spans[0],
        "empty_spans_end": empty_spans[-1],
    }
    print("\n=== summary ===")
    for k, v in summary.items():
        print(f"  {k}: {v}")

    with open(out_dir / "results.json", "w") as f:
        json.dump({"summary": summary,
                   "traces": {"drops": drops_axis, "crack_refill": crack_refill,
                              "biggest_empty": biggest_empty, "empty_total": empty_total}},
                  f, indent=2)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    axes[0].plot(drops_axis, crack_refill, "b-")
    axes[0].axhline(0.62, color="k", ls=":", alpha=0.6, label="equilibrium ~62%")
    axes[0].set_xlabel("drops after carving the crack")
    axes[0].set_ylabel("fraction of original crack cells refilled (z>=1)")
    axes[0].set_title("Does the crack refill? (healing)")
    axes[0].legend(); axes[0].grid(True, alpha=0.3)

    axes[1].plot(drops_axis, biggest_empty, "r-", label="biggest empty region")
    axes[1].axhline(crack_size, color="orange", ls="--", alpha=0.7, label="original crack size")
    axes[1].axhline(summary["background_biggest_empty_before_crack"], color="green", ls=":", alpha=0.7,
                    label="background (no crack)")
    axes[1].set_xlabel("drops after carving the crack")
    axes[1].set_ylabel("cells in biggest connected empty region")
    axes[1].set_title("Does the hole grow (feed) or shrink (heal)?")
    axes[1].legend(); axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_dir / "crack_dynamics.png", dpi=150)
    print(f"\nResults: {out_dir}/results.json\nPlot: {out_dir}/crack_dynamics.png")


if __name__ == "__main__":
    main()
