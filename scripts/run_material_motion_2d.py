"""The synthesis run: 2D radial implosion -> plasma -> rebound shell.

A ring of substrate rushes inward, collapses toward the centre (compresses +
heats to a plasma-like spike), and rebounds as an outward-expanding shell --
the inversion AND the launch of the expansion, in one run. Real material
motion; energy conserved to machine precision.

Honest scope: generic 2D compressible fluid; flat-local grid; abstract units
(no calibration). Density compression is modest (Lax-Friedrichs is diffusive);
the clear signature is the heat spike and the rebound shell.
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

from void_cascade.material_motion_2d import run_radial_implosion  # noqa: E402

L, STEPS = 180, 600
SNAPS = (0, 120, 240, 360, 480, 599)


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"material_motion_2d_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== 2D radial implosion -> plasma -> rebound shell ===\nOutput: {out_dir}\n")

    r = run_radial_implosion(L=L, u0=2.5, R0=55, steps=STEPS, band=3, snap_steps=SNAPS)
    for k in ("mass_residual", "energy_residual", "momentum_residual", "min_pressure",
              "peak_central_density", "central_density_initial", "peak_density_step",
              "central_heat_initial", "central_heat_at_peak", "rebounded"):
        print(f"  {k}: {r[k]}")
    print(f"\n  compression {r['peak_central_density']/r['central_density_initial']:.2f}x  "
          f"heat spike {r['central_heat_at_peak']/r['central_heat_initial']:.2f}x")

    summary = {k: r[k] for k in r if k not in ("t_trace", "central_density", "central_heat", "snapshots")}
    with open(out_dir / "results.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Density snapshots: converging ring -> collapsed core -> rebound shell.
    snaps = r["snapshots"]
    ks = sorted(snaps)
    fig, axes = plt.subplots(1, len(ks), figsize=(3.1 * len(ks), 3.4))
    for ax, k in zip(axes, ks):
        ax.imshow(snaps[k], cmap="inferno", interpolation="nearest")
        ax.set_title(f"step {k}"); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Density: implosion (ring in) -> collapse -> rebound shell (out)")
    fig.tight_layout(); fig.savefig(out_dir / "implosion_2d_snaps.png", dpi=150)

    # Central density + heat over time.
    fig2, ax2 = plt.subplots(1, 2, figsize=(13, 5))
    ax2[0].plot(r["t_trace"], r["central_density"]); ax2[0].axvline(r["peak_density_step"], color="r", ls=":")
    ax2[0].set_title("central density (collapse -> peak -> rebound)"); ax2[0].set_xlabel("step"); ax2[0].grid(True, alpha=0.3)
    ax2[1].plot(r["t_trace"], r["central_heat"], color="orange"); ax2[1].axvline(r["peak_density_step"], color="r", ls=":")
    ax2[1].set_title("central heat / plasma spike (internal energy)"); ax2[1].set_xlabel("step"); ax2[1].grid(True, alpha=0.3)
    fig2.tight_layout(); fig2.savefig(out_dir / "implosion_2d_central.png", dpi=150)

    print(f"\nResults: {out_dir}/results.json")
    print(f"Plots: {out_dir}/implosion_2d_snaps.png , {out_dir}/implosion_2d_central.png")


if __name__ == "__main__":
    main()
