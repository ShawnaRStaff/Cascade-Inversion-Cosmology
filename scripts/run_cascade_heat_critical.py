"""Damage-coupled heat, 1D: find the critical heat-making, and the
long-quiet-then-sudden tip.

The cooling sweep showed clustering lets it tip across all realistic cooling
when heat-making is strong (0.15) -- but then it tips instantly, with no
quiet buildup. The real dial is heat-making strength. Here we FIX cooling at
a realistic value and sweep the heat made per fracture, to find:
  - the critical heat-making (boundary: stays cold vs runs away), and
  - near that edge, a long quiet phase that then SUDDENLY tips -- the
    narrative shape (eons frozen, then catastrophe).

Built so weak heat-making stays cold (guardrail). 1D, local. Stay-cold runs
sit at criticality (costlier), so L and steps are modest.
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

L = 50
N_STEPS = 4000
SEED = 7
COOLING = 0.1  # fixed, realistic
HEATS = [0.0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.15]


def params(heat_per_crack: float) -> CascadeParams:
    return CascadeParams(
        fracture_density=2.0,
        heat_per_crack=heat_per_crack,
        diffuse=0.2,
        cooling=COOLING,
        melt_heat=1.0,
        release_factor=0.5,
        drive_amount=1.0,
        n_drive_sites=1,
    )


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"cascade_heat_critical_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== Critical heat-making at cooling={COOLING} ===\nOutput: {out_dir}\n")
    print(f"{'heat/crack':>11}{'ran away?':>11}{'tip step':>10}{'peak heat':>11}")

    rows = []
    traces = {}
    for h in HEATS:
        r = run(L=L, n_steps=N_STEPS, p=params(h), seed=SEED)
        r["heat_per_crack"] = h  # run() doesn't echo it back; attach for reporting
        rows.append(r)
        traces[h] = (r["steps_axis"], r["peak_heat_trace"])
        print(f"{h:>11.3f}{str(r['ran_away']):>11}{str(r['first_release_step']):>10}{r['peak_heat_overall']:>11.3f}")

    ran = [r for r in rows if r["ran_away"]]
    cold = [r for r in rows if not r["ran_away"]]
    crit_lo = max((r["heat_per_crack"] for r in cold), default=None)
    crit_hi = min((r["heat_per_crack"] for r in ran), default=None)
    print()
    if crit_hi is not None:
        print(f"Critical heat-making between {crit_lo} (stays cold) and {crit_hi} (runs away).")
    # The longest-quiet tipping case = the smallest heat that still tips.
    near_crit = min((r for r in ran), key=lambda r: r["heat_per_crack"], default=None)
    if near_crit is not None:
        print(f"Nearest-critical tip: heat={near_crit['heat_per_crack']} tips at step "
              f"{near_crit['first_release_step']} of {N_STEPS} "
              f"({near_crit['first_release_step']/N_STEPS*100:.0f}% in -- long quiet then sudden).")

    summary = {
        "L": L, "n_steps": N_STEPS, "cooling": COOLING, "seed": SEED,
        "critical_between": [crit_lo, crit_hi],
        "sweep": [{"heat_per_crack": r["heat_per_crack"] if "heat_per_crack" in r else None,
                   "ran_away": r["ran_away"], "first_release_step": r["first_release_step"],
                   "peak_heat_overall": r["peak_heat_overall"]} for r in rows],
    }
    with open(out_dir / "results.json", "w") as f:
        json.dump(summary, f, indent=2)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    hs = [r["heat_per_crack"] if "heat_per_crack" in r else None for r in rows]
    axes[0].plot(HEATS, [1 if r["ran_away"] else 0 for r in rows], "o-", label="ran away?")
    tip_steps = [(r["first_release_step"] if r["first_release_step"] is not None else N_STEPS) for r in rows]
    ax0b = axes[0].twinx()
    ax0b.plot(HEATS, tip_steps, "s--", color="orange", alpha=0.7, label="tip step")
    ax0b.set_ylabel("tip step (later = longer quiet)")
    axes[0].set_xlabel("heat made per fracture"); axes[0].set_ylabel("ran away? (1=yes)")
    axes[0].set_title(f"Critical heat-making at cooling={COOLING}")
    axes[0].grid(True, alpha=0.3)

    if near_crit is not None:
        s, h = traces[near_crit["heat_per_crack"]]
        axes[1].plot(s, h, label=f"near-critical heat={near_crit['heat_per_crack']} (tips)")
    if cold:
        cc = max(cold, key=lambda r: r["heat_per_crack"])
        s, h = traces[cc["heat_per_crack"]]
        axes[1].plot(s, h, label=f"just-below heat={cc['heat_per_crack']} (stays cold)")
    axes[1].axhline(1.0, color="r", ls="--", alpha=0.6, label="melt point")
    axes[1].set_xlabel("step (time)"); axes[1].set_ylabel("hottest cell's heat")
    axes[1].set_title("Long quiet, then a sudden tip")
    axes[1].legend(); axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_dir / "cascade_heat_critical.png", dpi=150)
    print(f"\nResults: {out_dir}/results.json\nPlot: {out_dir}/cascade_heat_critical.png")


if __name__ == "__main__":
    main()
