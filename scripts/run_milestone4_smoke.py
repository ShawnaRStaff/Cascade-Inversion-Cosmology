"""Milestone 4 smoke: xi(r) on the saved 3D Manna ever-toppled mask.

Loads the L=24 mask written by run_milestone3_smoke.py, computes the
two-point correlation function, and overlays a published galaxy ξ(r)
reference form so we can see, at a glance, whether the shape of the
Manna cluster correlation resembles galaxy clustering at all.

This is a *smoke* test for the M4 pipeline, not a real cosmological
comparison. L=24 is too small to extract a meaningful r range, and the
boundary effects from open BC vs the FFT estimator's periodic
assumption are larger here than at L=128. The point is to confirm the
xi(r) machinery runs end-to-end and produces a curve with sensible
shape before we commit to a real run.

Run:
    .venv/bin/python scripts/run_milestone4_smoke.py [path/to/smoke.npz]
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.correlation import (  # noqa: E402
    power_law_galaxy_xi,
    two_point_correlation_3d,
)
from void_cascade.percolation import check_spanning  # noqa: E402


def latest_smoke_npz() -> Path:
    outdir = REPO_ROOT / "data" / "outputs"
    candidates = sorted(outdir.glob("manna_3d_smoke_*.npz"))
    if not candidates:
        raise FileNotFoundError("No manna_3d_smoke_*.npz in data/outputs/")
    return candidates[-1]


def main() -> None:
    if len(sys.argv) > 1:
        path = Path(sys.argv[1])
    else:
        path = latest_smoke_npz()
    print(f"Loading {path}")
    data = np.load(path)
    ever = data["ever_toppled"]
    L = int(data["L"])
    mean_occ = float(ever.mean())
    print(f"L={L}, ever_toppled fraction = {mean_occ:.4f}, "
          f"n_toppled = {int(ever.sum())}")

    span = check_spanning(ever)
    print(f"Spanning: x={span.spans_x} y={span.spans_y} z={span.spans_z}; "
          f"largest cluster {span.largest_cluster_size}")

    print()
    print("Computing two-point correlation xi(r)...")
    r, xi, n_pairs = two_point_correlation_3d(
        ever, r_max=L / 2.0, n_bins=15, log_bins=True
    )
    print(f"{'r (lattice)':>12}  {'xi(r)':>10}  {'n_pairs':>10}")
    for ri, xii, ni in zip(r, xi, n_pairs):
        if np.isnan(xii):
            print(f"{ri:>12.3f}  {'n/a':>10}  {ni:>10}")
        else:
            print(f"{ri:>12.3f}  {xii:>10.4f}  {ni:>10}")

    # Filter for plotting: bins with enough pair samples
    valid = ~np.isnan(xi) & (n_pairs > 50) & (xi > 0)
    r_plot = r[valid]
    xi_plot = xi[valid]

    # Fit a power law xi(r) = A * r^(-gamma) over the small-r range
    # where xi is positive. Restrict to r <= L/4 to avoid the
    # boundary-artifact region.
    fit_sel = r_plot <= L / 4.0
    fit_result = None
    if fit_sel.sum() >= 3:
        log_r = np.log10(r_plot[fit_sel])
        log_xi = np.log10(xi_plot[fit_sel])
        coeffs = np.polyfit(log_r, log_xi, 1)
        gamma_fit = -coeffs[0]
        log_A = coeffs[1]
        # Solve for the implied r0: xi(r0) = 1 means A * r0^-gamma = 1
        # => r0 = A^(1/gamma)
        r0_fit = 10.0 ** (log_A / gamma_fit)
        fit_result = (gamma_fit, r0_fit, log_A)
        print()
        print(f"Power-law fit on r in [{r_plot[fit_sel].min():.2f}, "
              f"{r_plot[fit_sel].max():.2f}] (r <= L/4):")
        print(f"  xi(r) = (r / r0)^(-gamma) with")
        print(f"    gamma = {gamma_fit:.3f}  (galaxy ref ~ 1.8)")
        print(f"    r0    = {r0_fit:.3f} lattice units")
        print(f"  This r0 sets the physical scale: at the galaxy-ξ(r) fit")
        print(f"  r_0 ~ 5 h^-1 Mpc, our lattice spacing maps to")
        print(f"  ~ {5.0 / r0_fit:.2f} h^-1 Mpc / lattice cell.")

    # Plot
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(r_plot, xi_plot, "o", color="C0", label="Manna ever-toppled $\\xi(r)$")

    # Galaxy reference, scaled so its r0 matches our fitted r0 (if any),
    # otherwise placed at r0 = L/8 just for visual comparison.
    r_ref = np.logspace(np.log10(r_plot.min()), np.log10(r_plot.max()), 50)
    if fit_result is not None:
        ax.loglog(r_ref, power_law_galaxy_xi(r_ref, r0=fit_result[1], gamma=1.8),
                  "--", color="C1", alpha=0.7,
                  label=r"galaxy ref $(r/r_0)^{-1.8}$, $r_0=$ fitted")
        ax.loglog(r_ref, fit_result[2] is not None
                  and 10 ** (fit_result[2] - fit_result[0] * np.log10(r_ref)),
                  "r:", alpha=0.6,
                  label=rf"power-law fit $\gamma={fit_result[0]:.2f}$")
    else:
        ax.loglog(r_ref, power_law_galaxy_xi(r_ref, r0=L / 8.0, gamma=1.8),
                  "--", color="C1", alpha=0.7,
                  label=r"galaxy ref $(r/r_0)^{-1.8}$, $r_0 = L/8$")
    ax.axhline(0, color="0.7", linestyle=":", alpha=0.5)
    ax.set_xlabel(r"separation $r$ (lattice units)")
    ax.set_ylabel(r"$\xi(r)$")
    ax.set_title(f"3D Manna two-point correlation (smoke, L={L})")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)

    outdir = REPO_ROOT / "data" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    fig_path = outdir / f"manna_3d_xi_smoke_{stamp}.png"
    fig.tight_layout()
    fig.savefig(fig_path, dpi=150)
    print()
    print(f"Wrote {fig_path}")


if __name__ == "__main__":
    main()
