"""Reacting-flow detonation: the unification (motion + self-feeding front).

Part A: one detonation through loaded substrate -- front advances, fuel-driven.
Part B: fuel sweep -- find the critical fuel below which it fizzles and above
which the front self-sustains (proof it is fuel-driven, not the ignition blast
diffusing).

Honest scope: 1D reacting Euler (Lax-Friedrichs); abstract units; the exact
front speed is LF-diffusion-limited (the existence of the fuel-driven self-
sustaining front is the result).
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

from void_cascade.reacting_flow import run_detonation  # noqa: E402

FUELS = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"reacting_flow_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== Reacting-flow detonation (motion + self-feeding front) ===\nOutput: {out_dir}\n")

    # Part A: a propagating detonation.
    a = run_detonation(N=500, fuel0=3.0, e_ign=2.5, ignite_E=3.0, steps=500)
    print(f"[A] fuel=3.0: propagated={a['propagated']} burned={a['final_burned_fraction']:.2f} "
          f"front_speed={a['front_speed']:.3f} energy_residual={a['energy_residual']:.2e}")

    # Part B: fuel sweep -> critical fuel.
    print("\n[B] fuel sweep (small ignition, is propagation fuel-driven?):")
    sweep = []
    for f in FUELS:
        r = run_detonation(N=400, fuel0=f, e_ign=2.5, ignite_E=3.0, steps=400)
        sweep.append((f, r["final_burned_fraction"], r["propagated"]))
        print(f"   fuel={f}: burned={r['final_burned_fraction']:.2f} propagated={r['propagated']}")
    crit = next((f for f, b, p in sweep if p), None)
    print(f"\nCritical fuel for self-sustaining detonation ~ {crit}")

    summary = {"part_A": {k: a[k] for k in ("fuel0", "e_ign", "propagated",
                          "final_burned_fraction", "front_speed", "energy_residual")},
               "fuel_sweep": [{"fuel": f, "burned": b, "propagated": p} for f, b, p in sweep],
               "critical_fuel": crit}
    with open(out_dir / "results.json", "w") as f:
        json.dump(summary, f, indent=2)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].plot(a["t_trace"], a["front"]); axes[0].set_xlabel("step")
    axes[0].set_ylabel("detonation front position (cells from centre)")
    axes[0].set_title(f"Detonation propagates (fuel=3.0, speed~{a['front_speed']:.2f})")
    axes[0].grid(True, alpha=0.3)
    axes[1].plot([f for f, _, _ in sweep], [b for _, b, _ in sweep], "o-")
    axes[1].axhline(0.5, color="gray", ls=":", alpha=0.6)
    axes[1].set_xlabel("fuel (substrate loading)"); axes[1].set_ylabel("final burned fraction")
    axes[1].set_title("Fuel-driven: critical fuel for self-sustaining front")
    axes[1].grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(out_dir / "reacting_flow.png", dpi=150)
    print(f"\nResults: {out_dir}/results.json\nPlot: {out_dir}/reacting_flow.png")


if __name__ == "__main__":
    main()
