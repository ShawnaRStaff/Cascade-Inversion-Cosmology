"""The full picture: conserved energy on an edgeless growing substrate.

Reports the cornerstone (true conservation: residual ~0 AND boundary_lost ~0)
and the front behaviour (constant speed? sustains?). Flat-local lattice units;
ignited; activity wavefront (not material motion) -- all flagged.
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

from void_cascade.conserved_energy import EnergyParams  # noqa: E402
from void_cascade.conserved_energy_grow import run_grow_energy  # noqa: E402

SEEDS = [0, 1, 2]


def params():
    return EnergyParams(flip_threshold=1.0, release_fraction=1.0, diffuse=0.5,
                        drive_amount=1.0, n_drive_sites=1)


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"conserved_energy_grow_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== Conserved energy, edgeless growing substrate ===\nOutput: {out_dir}\n")

    runs = []
    fig, ax = plt.subplots(figsize=(9, 6))
    for s in SEEDS:
        r = run_grow_energy(L0=60, propagate_steps=1500, p=params(), seed=s,
                            margin=5, chunk=20, max_size=320)
        runs.append(r)
        # front speed via linear fit of radius vs time (pre-cap)
        t = np.array([tt for tt, _, _, _ in r["trace"]], float)
        fr = np.array([ff for _, _, ff, _ in r["trace"]], float)
        speed = float(np.polyfit(t, fr, 1)[0]) if t.size > 2 else 0.0
        r["front_speed"] = speed
        print(f"seed {s}: grow_events={r['grow_events']} final_size={r['final_size']} "
              f"max_front_radius={r['max_front_radius']:.1f} speed={speed:.3f} cells/step "
              f"sustained={r['front_sustained_and_grew']}  "
              f"CONSERVATION residual={r['conservation_residual']:.2e} "
              f"boundary_lost={r['boundary_lost']:.2e}")
        ax.plot(t, fr, label=f"seed {s} ({speed:.2f} cells/step, grew x{r['grow_events']})")

    ax.set_xlabel("steps after ignition"); ax.set_ylabel("front radius (cells)")
    ax.set_title("Conserved-energy front on edgeless substrate (energy truly conserved)")
    ax.legend(); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(out_dir / "conserved_energy_grow.png", dpi=150)

    summary = [{k: r[k] for k in ("grow_events", "final_size", "max_front_radius",
                                  "front_speed", "front_sustained_and_grew",
                                  "conservation_residual", "boundary_lost", "energy_in",
                                  "final_flipped_fraction")} for r in runs]
    with open(out_dir / "results.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nResults: {out_dir}/results.json\nPlot: {out_dir}/conserved_energy_grow.png")


if __name__ == "__main__":
    main()
