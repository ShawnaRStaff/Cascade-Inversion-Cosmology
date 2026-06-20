"""2D damage-coupled heat: (d) ground it as proper SOC, and (a) does the
catastrophe spread as a front?

Part A (grounding): heat off, measure the avalanche-size exponent tau. 2D
Manna's validated value is ~1.27 (M2). If we recover that, the foundation is
genuinely grounded (unlike 1D's degenerate ~0.78).

Part B (front): heat on at realistic cooling. When it tips, snapshot the
released region over the next steps -- does it grow outward from a seed (a
front) or appear everywhere at once?
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
from void_cascade.cascade_heat_2d import run_2d  # noqa: E402


def params(cooling, heat_per_crack):
    return CascadeParams(fracture_density=2.0, heat_per_crack=heat_per_crack,
                         diffuse=0.15, cooling=cooling, melt_heat=1.0,
                         release_factor=0.5, drive_amount=1.0, n_drive_sites=1)


def fit_tau(sizes: np.ndarray) -> tuple[float, np.ndarray, np.ndarray, tuple]:
    smax = sizes.max()
    bins = np.unique(np.round(np.logspace(0, np.log10(smax), 22)).astype(int))
    counts, edges = np.histogram(sizes, bins=bins)
    centers = np.sqrt(edges[:-1] * edges[1:])
    pdf = counts / np.diff(edges) / sizes.size
    ok = counts > 0
    centers, pdf = centers[ok], pdf[ok]
    lo, hi = 2.0, smax / 5.0
    sel = (centers >= lo) & (centers <= hi)
    slope, intercept = np.polyfit(np.log10(centers[sel]), np.log10(pdf[sel]), 1)
    return -slope, centers, pdf, (lo, hi, intercept, slope)


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"cascade_heat_2d_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== 2D cascade heat: grounding + front ===\nOutput: {out_dir}\n")

    # Part A: SOC grounding in 2D.
    print("[A] SOC grounding (heat off, L=20)...")
    ra = run_2d(L=20, n_steps=5000, p=params(0.05, 0.0), seed=0)
    sizes = np.asarray(ra["cascade_sizes"])[1000:]
    sizes = sizes[sizes > 0]
    tau, centers, pdf, (lo, hi, intercept, slope) = fit_tau(sizes)
    print(f"    avalanches={sizes.size} max={sizes.max()} mean={sizes.mean():.1f}")
    print(f"    tau = {tau:.3f}  (2D Manna reference ~1.27; 1D was degenerate ~0.78)")

    # Part B: front.
    print("\n[B] Front (heat on, L=48, cooling=0.1)...")
    rb = run_2d(L=48, n_steps=3000, p=params(0.10, 0.10), seed=1)
    print(f"    ran_away={rb['ran_away']} tip_step={rb['first_release_step']} "
          f"snapshots={len(rb['snapshots'])} final_released={rb['fraction_released_final']:.2%}")

    summary = {
        "soc_grounding": {"L": 20, "tau": float(tau), "reference_2D_Manna": 1.27,
                          "n_avalanches": int(sizes.size), "max": int(sizes.max())},
        "front": {"L": 48, "cooling": 0.10, "heat_per_crack": 0.10,
                  "ran_away": rb["ran_away"], "tip_step": rb["first_release_step"],
                  "final_released_fraction": rb["fraction_released_final"],
                  "released_at_snapshots": [{"after_tip": s["after_tip"],
                                             "n_released": int(s["released"].sum())}
                                            for s in rb["snapshots"]]},
    }
    with open(out_dir / "results.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Plot A: tau power law.
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.loglog(centers, pdf, "o", label="2D avalanche PDF")
    xs = np.array([lo, hi])
    ax.loglog(xs, 10 ** intercept * xs ** slope, "r-", label=f"fit tau={tau:.2f} (ref 1.27)")
    ax.set_xlabel("avalanche size"); ax.set_ylabel("P(s)")
    ax.set_title("2D SOC grounding: avalanche power law")
    ax.legend(); ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout(); fig.savefig(out_dir / "soc_2d.png", dpi=150)

    # Plot B: released-region snapshots (the front).
    snaps = rb["snapshots"]
    if snaps:
        ncol = len(snaps)
        fig2, axes = plt.subplots(1, ncol, figsize=(3 * ncol, 3.4))
        if ncol == 1:
            axes = [axes]
        for ax, s in zip(axes, snaps):
            ax.imshow(s["released"], cmap="hot", interpolation="nearest")
            ax.set_title(f"+{s['after_tip']} steps\n{int(s['released'].sum())} released")
            ax.set_xticks([]); ax.set_yticks([])
        fig2.suptitle("Does the catastrophe spread as a front? (released cells over time)")
        fig2.tight_layout(); fig2.savefig(out_dir / "front_2d.png", dpi=150)

    print(f"\nResults: {out_dir}/results.json")
    print(f"Plots: {out_dir}/soc_2d.png , {out_dir}/front_2d.png")


if __name__ == "__main__":
    main()
