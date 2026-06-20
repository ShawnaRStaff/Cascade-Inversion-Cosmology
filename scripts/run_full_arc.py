"""Full arc: buildup (emergent loaded substrate) -> implosion -> detonation.

Folds stage 1 (the buildup) into the end-to-end by REMOVING the imposed
uniform fuel: the loaded substrate now EMERGES from the SOC Manna buildup
(slow drive + threshold fracturing on the same grid), giving a spatially
structured fuel field. That emergent fuel feeds the implosion->detonation.

Stages: (1) SOC buildup loads + fractures the substrate -> emergent fuel;
(2) implosion collapses + heats -> ignites the emergent fuel; (3) detonation
front expands. Energy conserved in the fluid stage.

HONEST remaining imposition: the implosion TRIGGER (the converging initial
flow) is still imposed -- making *that* emerge needs a chosen collapse driver
(gravity / rigidity-loss / fracture-voids), a real theoretical decision we
should not impose. This run folds in the buildup (fuel emerges); the trigger
is the next, separate decision.
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
from void_cascade.cascade_heat_2d import avalanche_2d, drive_2d, initialize_2d  # noqa: E402
from void_cascade.reacting_flow_2d import run_end_to_end  # noqa: E402

L = 100
BUILD_STEPS = 2500
FUEL_SCALE = 4.0


def buildup(L, steps, seed):
    """SOC Manna loading (heat off): slow drive + threshold fracturing ->
    a loaded, spatially structured substrate (the emergent fuel)."""
    p = CascadeParams(fracture_density=2.0, heat_per_crack=0.0, diffuse=0.15,
                      cooling=0.05, melt_heat=1.0, release_factor=0.5,
                      drive_amount=1.0, n_drive_sites=5)
    rng = np.random.default_rng(seed)
    st = initialize_2d(L)
    maxsw = 50 * L * L
    for _ in range(steps):
        drive_2d(st, p, rng)
        avalanche_2d(st, p, rng, maxsw)
    return st.density.copy(), st.cracks.copy()


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"full_arc_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== Full arc: buildup -> implosion -> detonation ===\nOutput: {out_dir}\n")

    print("[1] SOC buildup (emergent loaded substrate)...")
    density, cracks = buildup(L, BUILD_STEPS, seed=0)
    fuel_field = density * FUEL_SCALE
    print(f"    emergent fuel: mean={fuel_field.mean():.2f} max={fuel_field.max():.2f} "
          f"loaded-fraction(>1)={float((fuel_field > 1).mean()):.2f} "
          f"mean_cracks={cracks.mean():.1f} (spatially structured, NOT uniform)")

    print("\n[2-3] implosion -> ignite emergent fuel -> detonation...")
    r = run_end_to_end(L=L, u0=2.5, R0=32, e_ign=2.5, steps=400, band=3,
                       fuel_field=fuel_field, snap_steps=(0, 40, 120, 399))
    print(f"    ignite_step={r['ignite_step']} fuel_burned={r['fuel_burned_fraction']:.2f} "
          f"max_burn_radius={r['max_burn_radius']:.0f} energy_residual={r['energy_residual']:.1e}")

    summary = {"buildup": {"L": L, "build_steps": BUILD_STEPS, "fuel_mean": float(fuel_field.mean()),
                           "fuel_max": float(fuel_field.max()),
                           "loaded_fraction": float((fuel_field > 1).mean()),
                           "mean_cracks": float(cracks.mean())},
               "arc": {k: r[k] for k in ("ignite_step", "fuel_burned_fraction",
                                         "max_burn_radius", "energy_residual",
                                         "peak_density_step")}}
    with open(out_dir / "results.json", "w") as f:
        json.dump(summary, f, indent=2)

    fig, axes = plt.subplots(1, 5, figsize=(16, 3.4))
    axes[0].imshow(fuel_field, cmap="viridis"); axes[0].set_title("emergent fuel\n(from SOC buildup)")
    axes[0].set_xticks([]); axes[0].set_yticks([])
    for ax, k in zip(axes[1:], sorted(r["snapshots"])):
        ax.imshow(r["snapshots"][k], cmap="inferno"); ax.set_title(f"density step {k}")
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Full arc: emergent substrate -> implosion -> detonation through structured fuel")
    fig.tight_layout(); fig.savefig(out_dir / "full_arc.png", dpi=150)
    print(f"\nResults: {out_dir}/results.json\nPlot: {out_dir}/full_arc.png")


if __name__ == "__main__":
    main()
