"""Material motion: watch the substrate implode, heat, and rebound.

Converging material (an implosion) -> compression at the centre -> kinetic
energy converts to internal (heat: plasma-like spike) -> pressure rebounds the
material outward. Real momentum; energy conserved. The inversion mechanism the
theory hinges on, now with moving material.

Honest scope: generic 1D compressible fluid (Euler), connected to the theory
as the substrate falling in on itself and rebounding. 1D; radial/2D next.
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

from void_cascade.material_motion import run_implosion  # noqa: E402


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"material_motion_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== Material motion: implosion -> heat -> rebound ===\nOutput: {out_dir}\n")

    r = run_implosion(N=300, u0=1.5, steps=700, cfl=0.4, band=6)
    for k in ("mass_residual", "energy_residual", "momentum_residual", "min_pressure",
              "peak_central_density", "central_density_initial", "peak_density_step",
              "central_heat_initial", "central_heat_at_peak", "rebounded"):
        print(f"  {k}: {r[k]}")
    print(f"\n  compression factor: {r['peak_central_density']/r['central_density_initial']:.2f}x")
    print(f"  central heat rise:  {r['central_heat_at_peak']/r['central_heat_initial']:.2f}x")

    summary = {k: r[k] for k in ("N", "u0", "steps", "mass_residual", "energy_residual",
                                 "momentum_residual", "min_pressure", "peak_central_density",
                                 "central_density_initial", "peak_density_step",
                                 "central_heat_initial", "central_heat_at_peak", "rebounded")}
    with open(out_dir / "results.json", "w") as f:
        json.dump(summary, f, indent=2)

    t = r["t_trace"]
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    axes[0].plot(t, r["central_density"]); axes[0].set_title("central density (implosion -> peak -> rebound)")
    axes[0].set_xlabel("step"); axes[0].set_ylabel("density"); axes[0].grid(True, alpha=0.3)
    axes[0].axvline(r["peak_density_step"], color="r", ls=":", alpha=0.6)
    axes[1].plot(t, r["central_heat"], color="orange"); axes[1].set_title("central heat (internal energy)")
    axes[1].set_xlabel("step"); axes[1].set_ylabel("internal energy"); axes[1].grid(True, alpha=0.3)
    axes[1].axvline(r["peak_density_step"], color="r", ls=":", alpha=0.6)
    axes[2].plot(t, r["central_velocity"], color="green"); axes[2].axhline(0, color="k", lw=0.6)
    axes[2].set_title("central velocity (in <0 ... rebound >0)")
    axes[2].set_xlabel("step"); axes[2].set_ylabel("velocity"); axes[2].grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(out_dir / "implosion.png", dpi=150)
    print(f"\nResults: {out_dir}/results.json\nPlot: {out_dir}/implosion.png")


if __name__ == "__main__":
    main()
