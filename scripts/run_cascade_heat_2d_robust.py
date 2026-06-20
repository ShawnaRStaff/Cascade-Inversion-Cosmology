"""2D heat-cascade: (1) robustness across seeds, (2) front speed.

(1) Does it reliably tip and spread as a front, not just on one lucky seed?
    - A fast batch (strong heat: tips quickly) for many seeds.
    - A slow batch (weak heat) to check the LONG-quiet-then-sudden shape is
      robust, not a fluke.
(2) Does the released region grow at a roughly CONSTANT rate? A constant
    front speed (linear size-vs-time) is the prerequisite for a Hubble-like
    distance-time relation (calibration to real H0 is deferred).

We track, each step after the tip, the released count and the released
region's radius (radius of gyration), and break shortly after so we don't
pay for the long quiet tail more than once.
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
from void_cascade.cascade_heat_2d import initialize_2d, step_2d  # noqa: E402

TRACE_STEPS = 50


def params(cooling, heat):
    return CascadeParams(fracture_density=2.0, heat_per_crack=heat, diffuse=0.15,
                         cooling=cooling, melt_heat=1.0, release_factor=0.5,
                         drive_amount=1.0, n_drive_sites=1)


def released_stats(mask: np.ndarray) -> tuple[int, float]:
    ys, xs = np.nonzero(mask)
    n = xs.size
    if n == 0:
        return 0, 0.0
    rg = float(np.sqrt(((xs - xs.mean()) ** 2 + (ys - ys.mean()) ** 2).mean()))
    return n, rg


def run_with_front_trace(L, n_steps, p, seed):
    rng = np.random.default_rng(seed)
    st = initialize_2d(L)
    maxsw = 50 * L * L
    tip = None
    trace = []
    for t in range(n_steps):
        _, n_rel = step_2d(st, p, rng, maxsw)
        if tip is None and n_rel > 0:
            tip = t
        if tip is not None:
            dt = t - tip
            if dt <= TRACE_STEPS:
                n, rg = released_stats(st.released)
                trace.append((dt, n, rg))
            else:
                break
    return {"seed": seed, "tip_step": tip, "ran_away": bool(st.released.any()),
            "released_after_trace": float(st.released.mean()), "trace": trace}


def front_grows(trace) -> bool:
    """A front = released region starts small and grows over several steps."""
    if len(trace) < 4:
        return False
    n0 = trace[0][1]
    nlate = trace[min(len(trace) - 1, 8)][1]
    return n0 < 40 and nlate > 4 * max(n0, 1)


def fit_speed(trace):
    """Fit radius vs time-after-tip; slope ~ front speed, plus R^2."""
    t = np.array([d for d, _, _ in trace], float)
    rg = np.array([r for _, _, r in trace], float)
    sel = (t >= 1)
    if sel.sum() < 3:
        return None, None
    slope, intercept = np.polyfit(t[sel], rg[sel], 1)
    pred = slope * t[sel] + intercept
    ss_res = ((rg[sel] - pred) ** 2).sum()
    ss_tot = ((rg[sel] - rg[sel].mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return float(slope), float(r2)


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"cascade_heat_2d_robust_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== 2D robustness + front speed ===\nOutput: {out_dir}\n")

    L = 40
    fast = []
    print("[1a] fast batch (heat=0.15), 8 seeds:")
    for s in range(8):
        r = run_with_front_trace(L, 1500, params(0.1, 0.15), seed=100 + s)
        spd, r2 = fit_speed(r["trace"])
        r["front"] = front_grows(r["trace"]); r["speed"] = spd; r["speed_r2"] = r2
        fast.append(r)
        print(f"   seed {100+s}: tip={r['tip_step']} front={r['front']} "
              f"speed={spd:.3f}/step R2={r2:.3f}" if spd else
              f"   seed {100+s}: tip={r['tip_step']} front={r['front']} (no speed fit)")

    slow = []
    print("\n[1b] slow batch (heat=0.10) — long-quiet robustness, 3 seeds:")
    for s in range(3):
        r = run_with_front_trace(L, 5000, params(0.1, 0.10), seed=200 + s)
        r["front"] = front_grows(r["trace"])
        slow.append(r)
        print(f"   seed {200+s}: tip={r['tip_step']} of 5000 front={r['front']}")

    tips_fast = [r["tip_step"] for r in fast if r["tip_step"] is not None]
    speeds = [r["speed"] for r in fast if r["speed"] is not None]
    summary = {
        "fast_batch": {"n": len(fast), "all_ran_away": all(r["ran_away"] for r in fast),
                       "all_front": all(r["front"] for r in fast),
                       "tip_step_mean": float(np.mean(tips_fast)) if tips_fast else None,
                       "tip_step_std": float(np.std(tips_fast)) if tips_fast else None,
                       "speed_mean": float(np.mean(speeds)) if speeds else None,
                       "speed_std": float(np.std(speeds)) if speeds else None,
                       "speed_r2_mean": float(np.mean([r["speed_r2"] for r in fast if r["speed_r2"] is not None]))},
        "slow_batch": {"n": len(slow), "all_ran_away": all(r["ran_away"] for r in slow),
                       "all_front": all(r["front"] for r in slow),
                       "tip_steps": [r["tip_step"] for r in slow]},
    }
    print("\n=== summary ===")
    print(json.dumps(summary, indent=2))
    with open(out_dir / "results.json", "w") as f:
        json.dump(summary, f, indent=2)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    axes[0].bar(range(len(tips_fast)), tips_fast, alpha=0.7, label="fast (heat=0.15)")
    axes[0].bar(range(len(tips_fast), len(tips_fast) + len(slow)),
                [r["tip_step"] for r in slow], alpha=0.7, color="orange", label="slow (heat=0.10)")
    axes[0].set_xlabel("seed"); axes[0].set_ylabel("tip step")
    axes[0].set_title("Robustness: every seed tips (long quiet for weak heat)")
    axes[0].legend(); axes[0].grid(True, alpha=0.3)

    for r in fast:
        if r["trace"]:
            t = [d for d, _, _ in r["trace"]]; rg = [g for _, _, g in r["trace"]]
            axes[1].plot(t, rg, alpha=0.6)
    axes[1].set_xlabel("steps after tip"); axes[1].set_ylabel("released region radius")
    axes[1].set_title(f"Front speed: radius vs time (mean {summary['fast_batch']['speed_mean']:.2f}/step)"
                      if speeds else "Front growth")
    axes[1].grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(out_dir / "robust_front.png", dpi=150)
    print(f"\nResults: {out_dir}/results.json\nPlot: {out_dir}/robust_front.png")


if __name__ == "__main__":
    main()
