"""Experiment 1: the cracked-window test.

Re-asks the resilience question the RIGHT way. The old experiment counted
how many cells still held grains (a density). This one asks whether the
leftover substrate is still one connected piece reaching edge-to-edge, or
whether it has broken into separate islands.

We remove a fraction of the substrate two ways:
  - random  (scattered, as the old experiment did)
  - connected (one grown blob, the proxy for a real spreading crack)
and for each we report, on the LEFTOVER (every cell the event did not
remove):
  - does it still span the box? (joined-up edge-to-edge)
  - how many separate pieces is it in?
  - what fraction of the leftover is in its single biggest piece?

We also print the OLD density number (carrier fraction) next to it, so we
can see exactly where the old conclusion came from and how it differs.

Pure geometry; the grain values do not affect connectivity. Runs in
seconds on this machine. No new physics.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.damage import connected_damage_mask, random_damage_mask  # noqa: E402
from void_cascade.percolation import check_spanning  # noqa: E402

FRACTIONS = [0.30, 0.50, 0.70, 0.91, 0.99]
SEEDS = [0, 1, 2, 3, 4]
GEOMETRIES = {"random": random_damage_mask, "connected": connected_damage_mask}


def load_saturated_field() -> np.ndarray:
    path = REPO_ROOT / "data" / "outputs" / "fss_sweep_20260521_031056" / "L48_s14800_final.npz"
    z = np.asarray(np.load(path, allow_pickle=True)["final_z"]).astype(np.int64)
    print(f"Loaded L={z.shape[0]} saturated field, z_avg={z.mean():.4f}")
    return z


def measure(z: np.ndarray, damage: np.ndarray) -> dict:
    """Connectivity of the leftover + the old density number, after damage."""
    intact = ~damage
    res = check_spanning(intact)
    n_intact = int(intact.sum())
    # Old test's number: set damaged cells to 0, measure carrier fraction.
    z_after = z.copy()
    z_after[damage] = 0
    carrier_fraction = float((z_after >= 1).mean())
    return {
        "leftover_spans": bool(res.percolates),
        "leftover_pieces": int(res.n_clusters),
        "leftover_biggest_piece_frac": (res.largest_cluster_size / n_intact) if n_intact else 0.0,
        "old_carrier_fraction": carrier_fraction,
    }


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"resilience_connectivity_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    z = load_saturated_field()
    shape = z.shape

    print(f"\n=== Cracked-window test ===\nOutput: {out_dir}\n")
    print(f"{'geometry':<11}{'removed':>8}{'spans?':>8}{'pieces':>8}{'biggest%':>10}{'old carrier%':>14}")

    results = []
    for geom_name, geom_fn in GEOMETRIES.items():
        for fraction in FRACTIONS:
            per_seed = []
            for seed in SEEDS:
                damage = geom_fn(shape, fraction, np.random.default_rng(seed))
                per_seed.append(measure(z, damage))
            spans_rate = float(np.mean([r["leftover_spans"] for r in per_seed]))
            pieces = float(np.mean([r["leftover_pieces"] for r in per_seed]))
            biggest = float(np.mean([r["leftover_biggest_piece_frac"] for r in per_seed]))
            carrier = float(np.mean([r["old_carrier_fraction"] for r in per_seed]))
            print(f"{geom_name:<11}{fraction:>7.0%}{spans_rate:>8.0%}{pieces:>8.0f}{biggest:>9.1%}{carrier:>13.1%}")
            results.append({
                "geometry": geom_name,
                "fraction_removed": fraction,
                "n_seeds": len(SEEDS),
                "leftover_spans_rate": spans_rate,
                "leftover_pieces_mean": pieces,
                "leftover_biggest_piece_frac_mean": biggest,
                "old_carrier_fraction_mean": carrier,
            })

    with open(out_dir / "results.json", "w") as f:
        json.dump(results, f, indent=2)

    # Plot: biggest-piece fraction vs removed, both geometries, plus the
    # old density number for comparison.
    fig, ax = plt.subplots(figsize=(9, 6))
    for geom_name in GEOMETRIES:
        rows = [r for r in results if r["geometry"] == geom_name]
        x = [r["fraction_removed"] for r in rows]
        ax.plot(x, [r["leftover_biggest_piece_frac_mean"] for r in rows], "o-",
                label=f"leftover biggest piece ({geom_name})")
    rows = [r for r in results if r["geometry"] == "random"]
    ax.plot([r["fraction_removed"] for r in rows],
            [r["old_carrier_fraction_mean"] for r in rows], "x--", color="gray",
            label="OLD test number (carrier density)")
    ax.set_xlabel("fraction of substrate removed")
    ax.set_ylabel("fraction of leftover in its single biggest piece")
    ax.set_title("Does the leftover hold together? (connected vs random damage)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "cracked_window.png", dpi=150)
    print(f"\nResults: {out_dir}/results.json\nPlot: {out_dir}/cracked_window.png")


if __name__ == "__main__":
    main()
