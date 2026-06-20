"""End-to-end (stages 2+3): 2D implosion -> ignition -> detonation expansion.

One run: a ring implodes, collapses to a plasma spike, IGNITES the loaded
fuel, and a detonation front expands outward. Inversion -> expansion in one
continuous model, the implosion causing the ignition. Energy conserved.

Honest scope: 2D reacting Euler (LF), abstract units; stages 2+3 (buildup
folded in later if warranted).
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

from void_cascade.reacting_flow_2d import run_end_to_end  # noqa: E402

SNAPS = (0, 40, 80, 160, 300, 499)


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"reacting_flow_2d_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== End-to-end: 2D implosion -> ignition -> detonation ===\nOutput: {out_dir}\n")

    r = run_end_to_end(L=180, u0=2.5, R0=55, fuel0=3.0, e_ign=2.5, steps=500, band=4, snap_steps=SNAPS)
    for k in ("energy_residual", "ignite_step", "peak_density_step",
              "peak_central_density", "fuel_burned_fraction", "max_burn_radius"):
        print(f"  {k}: {r[k]}")
    print(f"\n  implosion peaks ~step {r['peak_density_step']}, ignites ~step {r['ignite_step']}, "
          f"detonation burns {r['fuel_burned_fraction']:.0%} of the fuel out to radius {r['max_burn_radius']:.0f}")

    summary = {k: r[k] for k in r if k not in ("t_trace", "central_density", "central_heat",
                                               "burned_fraction", "burn_radius", "snapshots")}
    with open(out_dir / "results.json", "w") as f:
        json.dump(summary, f, indent=2)

    snaps = r["snapshots"]
    ks = sorted(snaps)
    fig, axes = plt.subplots(1, len(ks), figsize=(3.0 * len(ks), 3.3))
    for ax, k in zip(axes, ks):
        ax.imshow(snaps[k], cmap="inferno", interpolation="nearest")
        ax.set_title(f"step {k}"); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Density: implode (ring) -> collapse/ignite -> detonation front out")
    fig.tight_layout(); fig.savefig(out_dir / "end_to_end_snaps.png", dpi=150)

    fig2, ax2 = plt.subplots(1, 2, figsize=(13, 5))
    ax2[0].plot(r["t_trace"], r["central_heat"], color="orange")
    if r["ignite_step"] is not None:
        ax2[0].axvline(r["ignite_step"], color="r", ls=":", label="ignition")
    ax2[0].set_title("central heat (implosion -> plasma -> ignition)"); ax2[0].set_xlabel("step")
    ax2[0].legend(); ax2[0].grid(True, alpha=0.3)
    ax2[1].plot(r["t_trace"], r["burn_radius"], color="green")
    ax2[1].set_title("detonation front radius (the expansion)"); ax2[1].set_xlabel("step")
    ax2[1].grid(True, alpha=0.3)
    fig2.tight_layout(); fig2.savefig(out_dir / "end_to_end_traces.png", dpi=150)

    print(f"\nResults: {out_dir}/results.json")
    print(f"Plots: {out_dir}/end_to_end_snaps.png , {out_dir}/end_to_end_traces.png")


if __name__ == "__main__":
    main()
