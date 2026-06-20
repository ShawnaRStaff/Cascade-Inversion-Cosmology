"""Grounding check: is cascade_heat's fracturing genuine SOC?

Self-organized criticality has one unmistakable signature: a power-law
distribution of avalanche sizes over several decades. The cosmological-SOC
literature (Bak-Tang-Wiesenfeld; Manna; Moffat's SOC universe;
Carfora-Marzuoli) rests on exactly this. If our damage-coupled model is to
claim the SOC foundation, its avalanches must show that power law.

We run it heat-off (so it sits at criticality, no releases), collect
avalanche sizes, and fit the slope tau in log-log. Reference points from
this project's own validated sandpiles: 1D Oslo tau~1.55, 2D Manna ~1.27,
3D Manna ~1.35. Our model is a 1D continuous-density Manna variant, so we
report tau and compare honestly -- the headline is "clean power law or not".
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
N_STEPS = 12_000
TRANSIENT = 2_000  # drop while it fills to criticality
SEED = 0


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"cascade_soc_check_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    p = CascadeParams(fracture_density=2.0, heat_per_crack=0.0, diffuse=0.2,
                      cooling=0.05, melt_heat=1.0, release_factor=0.5,
                      drive_amount=1.0, n_drive_sites=1)
    print(f"=== SOC grounding check (cascade_heat, heat off) ===\nOutput: {out_dir}")
    r = run(L=L, n_steps=N_STEPS, p=p, seed=SEED)
    sizes = np.asarray(r["cascade_sizes"])[TRANSIENT:]
    sizes = sizes[sizes > 0]
    print(f"avalanches: {sizes.size}, max size {sizes.max()}, mean {sizes.mean():.1f}")

    # Log-binned PDF.
    smax = sizes.max()
    bins = np.unique(np.round(np.logspace(0, np.log10(smax), 25)).astype(int))
    counts, edges = np.histogram(sizes, bins=bins)
    widths = np.diff(edges)
    centers = np.sqrt(edges[:-1] * edges[1:])
    pdf = counts / widths / sizes.size
    ok = (counts > 0) & (pdf > 0)
    centers, pdf = centers[ok], pdf[ok]

    # Fit tau over the scaling region (drop the smallest bin and the last
    # couple, which feel the lower cutoff and the finite-size bump).
    fit_lo, fit_hi = 2.0, smax / 5.0
    sel = (centers >= fit_lo) & (centers <= fit_hi)
    slope, intercept = np.polyfit(np.log10(centers[sel]), np.log10(pdf[sel]), 1)
    tau = -slope
    print(f"\nFitted avalanche-size exponent tau = {tau:.3f} "
          f"(fit range size {fit_lo:.0f}-{fit_hi:.0f}, {sel.sum()} bins)")
    print("Reference (this project's validated sandpiles): 1D Oslo ~1.55, "
          "2D Manna ~1.27, 3D Manna ~1.35.")

    summary = {"L": L, "n_steps": N_STEPS, "n_avalanches": int(sizes.size),
               "max_size": int(smax), "tau": float(tau),
               "fit_range": [fit_lo, float(fit_hi)], "n_fit_bins": int(sel.sum())}
    with open(out_dir / "results.json", "w") as f:
        json.dump(summary, f, indent=2)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.loglog(centers, pdf, "o", label="avalanche-size PDF")
    xs = np.array([fit_lo, fit_hi])
    ax.loglog(xs, 10 ** intercept * xs ** slope, "r-", label=f"power-law fit, tau={tau:.2f}")
    ax.set_xlabel("avalanche size s"); ax.set_ylabel("P(s)")
    ax.set_title("SOC fingerprint: is the avalanche distribution a power law?")
    ax.legend(); ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "soc_check.png", dpi=150)
    print(f"\nResults: {out_dir}/results.json\nPlot: {out_dir}/soc_check.png")


if __name__ == "__main__":
    main()
