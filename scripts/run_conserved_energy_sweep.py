"""Conserved-energy front: where does it sustain vs fizzle?

The model-led, theory-relevant question: HOW LOADED must the pre-existing
substrate be (relative to threshold) for the ignited catastrophe to
propagate -- and how does that trade off against how much energy each
combustion releases?

Fast fixed box (sustain/fizzle only needs "does the front sweep the box or
stay local", not edgelessness). For each (substrate loading, release
fraction): pre-load potential to that level, ignite the centre, run, and
record the final flipped fraction. sustain ~ swept; fizzle ~ stayed local.

Flat-local lattice units; ignited; activity wavefront. Conservation still
holds (energy only changes form + tracked boundary loss).
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

from void_cascade.conserved_energy import (  # noqa: E402
    EnergyParams, combust, conservation_residual, diffuse_kinetic, ignite, initialize_2d,
)

L = 80
STEPS = 500
THRESH = 1.0
LOADS = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
RELEASE_FRACTIONS = [0.5, 0.75, 1.0]
DIFFUSE = 0.5


def run_one(load, rf, seed):
    p = EnergyParams(flip_threshold=THRESH, release_fraction=rf, diffuse=DIFFUSE,
                     drive_amount=1.0, n_drive_sites=1)
    rng = np.random.default_rng(seed)
    st = initialize_2d(L)
    st.potential[:] = np.clip(load * THRESH + rng.uniform(-0.02, 0.02, (L, L)), 0.0, None)
    st.energy_in = float(st.potential.sum())
    ignite(st, p)
    for _ in range(STEPS):
        combust(st, p)
        diffuse_kinetic(st, p)
    return float(st.flipped.mean()), float(conservation_residual(st) - (-st.boundary_lost))


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"conserved_energy_sweep_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== Conserved-energy sustain/fizzle sweep ===\nOutput: {out_dir}\n")
    print(f"{'load':>6}" + "".join(f"  rf={rf:<5}" for rf in RELEASE_FRACTIONS))

    grid = {}
    for load in LOADS:
        row = []
        for rf in RELEASE_FRACTIONS:
            frac, _resid = run_one(load, rf, seed=0)
            grid[(load, rf)] = frac
            row.append(frac)
        print(f"{load:>6.2f}" + "".join(f"  {f:>6.1%}" for f in row))

    # critical load per release fraction = first load with flipped fraction > 0.5
    crit = {}
    for rf in RELEASE_FRACTIONS:
        c = next((load for load in LOADS if grid[(load, rf)] > 0.5), None)
        crit[rf] = c
    print("\nCritical loading (first load that sweeps the box) per release fraction:")
    for rf in RELEASE_FRACTIONS:
        print(f"  release_fraction={rf}: critical load ~ {crit[rf]}")

    summary = {"L": L, "steps": STEPS, "threshold": THRESH, "diffuse": DIFFUSE,
               "grid": {f"{load}_{rf}": grid[(load, rf)] for load in LOADS for rf in RELEASE_FRACTIONS},
               "critical_load_by_release_fraction": {str(rf): crit[rf] for rf in RELEASE_FRACTIONS}}
    with open(out_dir / "results.json", "w") as f:
        json.dump(summary, f, indent=2)

    fig, ax = plt.subplots(figsize=(9, 6))
    for rf in RELEASE_FRACTIONS:
        ax.plot(LOADS, [grid[(load, rf)] for load in LOADS], "o-",
                label=f"release fraction {rf}")
    ax.axhline(0.5, color="gray", ls=":", alpha=0.6, label="sweep/fizzle line")
    ax.set_xlabel("substrate loading (mean potential / threshold)")
    ax.set_ylabel("final flipped fraction")
    ax.set_title("Conserved-energy front: how loaded must the substrate be to propagate?")
    ax.legend(); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(out_dir / "sweep.png", dpi=150)
    print(f"\nResults: {out_dir}/results.json\nPlot: {out_dir}/sweep.png")


if __name__ == "__main__":
    main()
