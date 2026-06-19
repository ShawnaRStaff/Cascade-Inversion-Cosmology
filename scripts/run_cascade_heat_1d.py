"""Damage-coupled heat model, 1D: does avalanche clustering let it tip under
realistic cooling -- unlike the scattered model?

The earlier scattered heat_gated model only tipped at essentially zero
cooling (critical cooling ~0.003). Here fracturing is sandpile-coupled, so
breaks come in clustered avalanches that dump heat concentrated in space and
time. We sweep the cooling rate and find how high it can go before the
substrate stops tipping -- the headline comparison.

Also captures one heat-off run to show the fracturing itself going from
scattered (sparse early lattice) to cascading (filled, critical).

1D, local. Tipping runs go quiet fast; only the stay-cold runs sit at
criticality, so sizes are kept modest.
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

from void_cascade.cascade_heat import CascadeParams, run  # noqa: E402

L = 60
N_STEPS = 2500
SEED = 7
HEAT_PER_CRACK = 0.15
COOLINGS = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
SCATTERED_MODEL_CRITICAL_COOLING = 0.003  # from heat_gated_1d run, for reference


def params(cooling: float, heat_per_crack: float = HEAT_PER_CRACK) -> CascadeParams:
    return CascadeParams(
        fracture_density=2.0,
        heat_per_crack=heat_per_crack,
        diffuse=0.2,
        cooling=cooling,
        melt_heat=1.0,
        release_factor=0.5,
        drive_amount=1.0,
        n_drive_sites=1,
    )


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"cascade_heat_1d_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== Damage-coupled heat, 1D: tipping vs cooling ===\nOutput: {out_dir}\n")

    # 1) Cooling sweep (heat on).
    print(f"{'cooling':>9}{'ran away?':>11}{'tip step':>10}{'peak heat':>11}")
    rows = []
    heat_traces = {}
    for c in COOLINGS:
        r = run(L=L, n_steps=N_STEPS, p=params(c), seed=SEED)
        rows.append(r)
        heat_traces[c] = (r["steps_axis"], r["peak_heat_trace"])
        print(f"{c:>9.3f}{str(r['ran_away']):>11}{str(r['first_release_step']):>10}{r['peak_heat_overall']:>11.3f}")

    ran = [r["cooling"] for r in rows if r["ran_away"]]
    cold = [r["cooling"] for r in rows if not r["ran_away"]]
    critical = (max(ran), min(cold)) if ran and cold else None
    print()
    if ran:
        print(f"Clustered model tips up to cooling = {max(ran)} "
              f"(scattered model only tipped below ~{SCATTERED_MODEL_CRITICAL_COOLING}).")
    if critical:
        print(f"Critical cooling between {critical[0]} (tips) and {critical[1]} (stays cold).")

    # 2) Heat-off run to show scattered -> cascading.
    rg = run(L=50, n_steps=400, p=params(0.05, heat_per_crack=0.0), seed=1)
    cs = np.asarray(rg["cascade_sizes"])

    summary = {
        "L": L, "n_steps": N_STEPS, "heat_per_crack": HEAT_PER_CRACK, "seed": SEED,
        "scattered_model_critical_cooling": SCATTERED_MODEL_CRITICAL_COOLING,
        "clustered_tips_up_to_cooling": max(ran) if ran else None,
        "critical_between": critical,
        "sweep": [{k: r[k] for k in ("cooling", "ran_away", "first_release_step",
                                     "peak_heat_overall")} for r in rows],
        "growth_demo": {"early_10pct_mean_avalanche": float(cs[: len(cs)//10].mean()),
                        "late_50pct_mean_avalanche": float(cs[len(cs)//2:].mean()),
                        "max_avalanche": int(cs.max())},
    }
    with open(out_dir / "results.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Plots.
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    axes[0].plot([r["cooling"] for r in rows], [1 if r["ran_away"] else 0 for r in rows], "o-")
    axes[0].axvline(SCATTERED_MODEL_CRITICAL_COOLING, color="gray", ls=":",
                    label=f"scattered model died at ~{SCATTERED_MODEL_CRITICAL_COOLING}")
    axes[0].set_xlabel("cooling rate"); axes[0].set_ylabel("ran away? (1=yes)")
    axes[0].set_title("Clustering lets it tip at FAR higher cooling")
    axes[0].legend(); axes[0].grid(True, alpha=0.3)

    axes[1].plot(np.arange(cs.size), cs, lw=0.6)
    axes[1].set_xlabel("step"); axes[1].set_ylabel("avalanche size")
    axes[1].set_title("Fracturing: scattered (early) -> cascading (filled)")
    axes[1].grid(True, alpha=0.3)

    tip_c = next((r["cooling"] for r in rows if r["ran_away"]), None)
    if tip_c is not None:
        s, h = heat_traces[tip_c]
        axes[2].plot(s, h, label=f"cooling={tip_c}")
    axes[2].axhline(1.0, color="r", ls="--", alpha=0.6, label="melt point")
    axes[2].set_xlabel("step"); axes[2].set_ylabel("hottest cell's heat")
    axes[2].set_title("Heat over time for a tipping case")
    axes[2].legend(); axes[2].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_dir / "cascade_heat_1d.png", dpi=150)
    print(f"\nResults: {out_dir}/results.json\nPlot: {out_dir}/cascade_heat_1d.png")


if __name__ == "__main__":
    main()
