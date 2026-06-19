"""Heat-gated substrate, 1D: does the cold buildup tip into a runaway?

Sweeps the cooling rate (how fast heat escapes to the surrounding cold) and
asks, for each value: does the substrate stay frozen forever, or does
breaking outpace cooling and run away? Finds the critical line between the
two, and shows the heat-over-time shape (long quiet, then sudden tip?).

Everything else is held fixed and principled; cooling is the guardrail dial.
1D, fast, local. No tuning to force a runaway -- strong cooling stays cold.
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

from void_cascade.heat_gated import HeatParams, run  # noqa: E402

L = 400
N_STEPS = 10_000
SEED = 7
COOLINGS = [0.0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5]


def params(cooling: float) -> HeatParams:
    return HeatParams(
        fracture_density=1.0,
        fracture_relief=1.0,
        drive_rate=0.02,
        heat_per_crack=0.15,
        diffuse=0.2,
        cooling=cooling,
        melt_heat=1.0,
        release_factor=0.5,
    )


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"heat_gated_1d_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== Heat-gated 1D: stay cold or run away? ===\nOutput: {out_dir}\n")
    print(f"{'cooling':>9}{'ran away?':>11}{'tip step':>10}{'peak heat':>11}{'% released':>12}{'hidden cracks@tip':>19}")

    rows = []
    traces = {}
    for c in COOLINGS:
        res = run(L=L, n_steps=N_STEPS, p=params(c), seed=SEED)
        rows.append(res)
        traces[c] = {"steps": res["steps_axis"], "peak_heat": res["peak_heat_trace"]}
        tip = res["first_release_step"]
        print(f"{c:>9.3f}{str(res['ran_away']):>11}{str(tip):>10}"
              f"{res['peak_heat_overall']:>11.3f}{res['fraction_released_final']*100:>11.1f}%"
              f"{res['hidden_cracks_at_first_release']:>19}")

    # Critical cooling: lowest cooling at which it stays cold.
    cold = [r["cooling"] for r in rows if not r["ran_away"]]
    ranaway = [r["cooling"] for r in rows if r["ran_away"]]
    critical = (min(cold), max(ranaway)) if cold and ranaway else None
    print(f"\nRan away for cooling <= {max(ranaway) if ranaway else 'none'}; "
          f"stayed cold for cooling >= {min(cold) if cold else 'none'}.")
    if critical:
        print(f"Critical line sits between cooling = {critical[1]} (runs away) and {critical[0]} (stays cold).")

    summary = {
        "params_fixed": {k: getattr(params(0.0), k) for k in
                         ("fracture_density", "fracture_relief", "drive_rate",
                          "heat_per_crack", "diffuse", "melt_heat", "release_factor")},
        "L": L, "n_steps": N_STEPS, "seed": SEED,
        "critical_between": critical,
        "rows": [{k: r[k] for k in ("cooling", "ran_away", "first_release_step",
                                    "peak_heat_overall", "fraction_released_final",
                                    "hidden_cracks_at_first_release")} for r in rows],
    }
    with open(out_dir / "results.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Plot. Left: did it run away vs cooling. Right: heat-over-time for a
    # runaway case and a stayed-cold case (the tip shape).
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    axes[0].plot([r["cooling"] for r in rows], [1 if r["ran_away"] else 0 for r in rows], "o-")
    axes[0].set_xlabel("cooling rate (heat lost per step)")
    axes[0].set_ylabel("ran away?  (1 = yes, 0 = stayed cold)")
    axes[0].set_title("The critical line: where the cold stops tipping")
    axes[0].grid(True, alpha=0.3)

    runaway_case = next((r["cooling"] for r in rows if r["ran_away"]), None)
    cold_case = next((r["cooling"] for r in rows if not r["ran_away"]), None)
    for c, label in [(runaway_case, "runs away"), (cold_case, "stays cold")]:
        if c is not None:
            axes[1].plot(traces[c]["steps"], traces[c]["peak_heat"], label=f"cooling={c} ({label})")
    axes[1].axhline(1.0, color="r", ls="--", alpha=0.6, label="melt point")
    axes[1].set_xlabel("step (time)")
    axes[1].set_ylabel("hottest cell's heat")
    axes[1].set_title("Long quiet, then a sudden tip?")
    axes[1].legend(); axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_dir / "heat_gated_1d.png", dpi=150)
    print(f"\nResults: {out_dir}/results.json\nPlot: {out_dir}/heat_gated_1d.png")


if __name__ == "__main__":
    main()
