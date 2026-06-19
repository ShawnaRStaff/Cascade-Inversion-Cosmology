"""Finish out M5: does the permanent catastrophe plateau survive without healing?

Milestone 5's headline: the substrate enters a permanent, never-ending
regime of huge events (sustained 50-70% of L^3 at L=48, z_avg=0.616 held
indefinitely). That permanence is fed by self-healing — cells empty, refill
and re-fire forever.

Here we re-run the same lattice but with IRREVERSIBLE fracture: a cell
cracks at most once. Two forms:
  - sink: fractured cells still soak up grains, never fire again
  - hole: grains routed into fractured cells are lost

We measure, we do not presuppose. M5 (healing) reference at L=48:
  peak single event ~71% of L^3, big events sustained indefinitely.
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

from void_cascade.irreversible import run  # noqa: E402

L = 48
N_DROPS = 400_000
SEED = 4800
M5_HEALING_PEAK_PCT = 71.0  # documented L=48 healing peak, for reference


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"m5_irreversible_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== M5 closure: irreversible fracture at L={L} ===\nOutput: {out_dir}\n")

    summaries = {}
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    for mode in ("sink", "hole"):
        print(f"--- mode = {mode} ---")
        res = run(L=L, n_drops=N_DROPS, mode=mode, seed=SEED)
        sizes = np.asarray(res.pop("sizes"), dtype=np.int64)
        np.savez_compressed(out_dir / f"{mode}_sizes.npz", sizes=sizes)
        summaries[mode] = res
        for key in ("drops_done", "total_topplings", "max_event_pct_of_volume",
                    "final_fractured_fraction", "fully_fractured_at_drop",
                    "events_after_full_fracture_max", "mean_event_last_tenth"):
            print(f"    {key}: {res[key]}")
        # conservation sanity
        print(f"    conservation ok: {res['grains_accounted'] == res['grains_in']}")
        print()

        drops = np.arange(sizes.size)
        axes[0].plot(drops, sizes / res["volume"] * 100.0, lw=0.5, label=f"{mode}")
        axes[1].plot(res["frac_trace_drops"], res["frac_trace_values"], label=f"{mode}")

    axes[0].axhline(M5_HEALING_PEAK_PCT, color="r", ls="--", alpha=0.6,
                    label=f"M5 healing peak ~{M5_HEALING_PEAK_PCT:.0f}%")
    axes[0].set_xlabel("drop")
    axes[0].set_ylabel("event size (% of L^3)")
    axes[0].set_title("Event sizes: does a catastrophe plateau appear without healing?")
    axes[0].legend(); axes[0].grid(True, alpha=0.3)

    axes[1].set_xlabel("drop")
    axes[1].set_ylabel("fraction of cells fractured")
    axes[1].set_title("How much of the substrate has fractured")
    axes[1].legend(); axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_dir / "m5_irreversible.png", dpi=150)

    with open(out_dir / "results.json", "w") as f:
        json.dump(summaries, f, indent=2)
    print(f"Results: {out_dir}/results.json\nPlot: {out_dir}/m5_irreversible.png")


if __name__ == "__main__":
    main()
