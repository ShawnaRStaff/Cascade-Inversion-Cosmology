"""Edgeless growing substrate: does the ignited front sustain, with no edge?

Loads + accumulates hidden damage, ignites one spot, then propagates by
heat-release cascade while the domain GROWS fresh substrate ahead of the
front -- so the answer can't be confused by an edge, a reflection, or a wrap.

Honest labels: distances are local-grid (flat-local) cell counts; we impose
NO edge and NO global geometry, only local connectivity.
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

from void_cascade.cascade_heat import CascadeParams  # noqa: E402
from void_cascade.cascade_heat_grow import run_grow  # noqa: E402

SEEDS = [0, 1, 2]


def params():
    return CascadeParams(fracture_density=2.0, heat_per_crack=0.10, diffuse=0.15,
                         cooling=0.10, melt_heat=1.0, release_factor=0.5,
                         drive_amount=1.0, n_drive_sites=1)


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"cascade_heat_grow_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== Edgeless growing substrate ===\nOutput: {out_dir}\n")

    runs = []
    fig, ax = plt.subplots(figsize=(9, 6))
    for s in SEEDS:
        r = run_grow(L0=60, accumulate_steps=800, propagate_steps=1500,
                     p=params(), seed=s, margin=4, chunk=20, max_size=320)
        runs.append(r)
        print(f"seed {s}: cracks@ignition~{r['mean_cracks_at_ignition']}  "
              f"grow_events={r['grow_events']}  final_size={r['final_size']}  "
              f"max_front_radius={r['max_front_radius']:.1f}  "
              f"sustained_and_grew={r['front_sustained_and_grew']}  "
              f"touched_edge={r['ever_touched_edge']}  "
              f"capped={r['capped_at_max_size']}")
        t = [tt for tt, _, _, _ in r["trace"]]
        fr = [ff for _, _, ff, _ in r["trace"]]
        ax.plot(t, fr, label=f"seed {s} (grew x{r['grow_events']}, final {r['final_size']})")

    ax.set_xlabel("steps after ignition")
    ax.set_ylabel("front radius (cells from centre)")
    ax.set_title("Edgeless front: does it keep advancing into fresh substrate?")
    ax.legend(); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(out_dir / "grow_front.png", dpi=150)

    summary = [{k: r[k] for k in ("mean_cracks_at_ignition", "grow_events",
                                  "final_size", "max_front_radius",
                                  "front_sustained_and_grew", "ever_touched_edge",
                                  "capped_at_max_size", "propagate_steps_run")}
               for r in runs]
    with open(out_dir / "results.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nResults: {out_dir}/results.json\nPlot: {out_dir}/grow_front.png")


if __name__ == "__main__":
    main()
