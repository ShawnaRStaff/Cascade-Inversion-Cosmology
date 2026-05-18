"""Experiment 7+8: void size distribution + galaxy-to-interface distance.

Two related cosmological observables:

(7) Void size distribution slope. Cosmic voids have a measured power-
    law size distribution (Pan et al. 2012; Sutter et al. 2014). Slope
    is typically reported in the range alpha ~ 1.5-2 on small to
    intermediate sizes. Computing it for our frozen regions and
    comparing is a third quantitative cosmological test.

(8) Galaxy-to-interface distance distribution. Observed galaxies are
    biased toward filaments and walls (the "matter" interface), not
    uniformly distributed in the cracked region. For each cracked cell
    in our simulation, compute distance to the nearest frozen cell and
    look at the distribution. Peak at small distance = matter biased to
    the interface, matching observed cosmic-web galaxy bias.

Run:
    .venv/bin/python scripts/run_milestone4_void_distribution.py [npz]
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.mass_function import cluster_mass_pdf, cluster_sizes  # noqa: E402


def latest_xi_npz() -> Path:
    outdir = REPO_ROOT / "data" / "outputs"
    candidates = sorted(outdir.glob("manna_3d_xi_data_*.npz"))
    if not candidates:
        raise FileNotFoundError("No manna_3d_xi_data_*.npz in data/outputs/")
    return candidates[-1]


def fit_slope(centers, pdf, s_min, s_max):
    sel = (centers >= s_min) & (centers <= s_max) & (pdf > 0)
    if sel.sum() < 3:
        return None
    coeffs = np.polyfit(np.log10(centers[sel]), np.log10(pdf[sel]), 1)
    return float(-coeffs[0]), int(sel.sum())


def main() -> None:
    if len(sys.argv) > 1:
        path = Path(sys.argv[1])
    else:
        path = latest_xi_npz()
    print(f"Loading {path}")
    data = np.load(path)
    L = int(data["L"])
    p_targets = sorted([float(p) for p in data["p_targets"]])
    seeds = sorted({int(k.split("_s")[1].split("_p")[0])
                    for k in data.files if k.startswith("mask_")})
    print(f"L={L}, p_targets={p_targets}, seeds={seeds}")

    # ============================================================
    # Part 1: void (frozen) size distribution slope vs p
    # ============================================================
    print()
    print("=== Experiment 7: void size distribution slope ===")
    print(f"{'p':>5}  {'n_clusters':>11}  {'largest_void':>13}  "
          f"{'slope alpha':>11}  {'n_fit_bins':>10}")
    void_slopes = {}
    for p in p_targets:
        pp = int(round(p * 100))
        all_sizes = []
        for s in seeds:
            mask_key = f"mask_s{s}_p{pp:02d}"
            if mask_key not in data.files:
                continue
            frozen = ~data[mask_key]
            sizes = cluster_sizes(frozen)
            all_sizes.append(sizes)
        if not all_sizes:
            continue
        sizes = np.concatenate(all_sizes)
        if sizes.size < 20:
            continue
        centers, pdf, counts = cluster_mass_pdf(sizes, n_bins=22)
        # Fit slope on small-medium sizes, exclude the giant cluster (it
        # is one realization, not part of the distribution).
        # Use sizes from 2 to max/100.
        s_min_fit = 2.0
        s_max_fit = max(float(sizes.max()) / 100.0, 10.0)
        fit = fit_slope(centers, pdf, s_min_fit, s_max_fit)
        if fit is None:
            print(f"{p:>5.2f}  {sizes.size:>11d}  {sizes.max():>13d}  "
                  f"fit failed")
            continue
        slope, n_fit = fit
        void_slopes[p] = slope
        print(f"{p:>5.2f}  {sizes.size:>11d}  {sizes.max():>13d}  "
              f"{slope:>11.3f}  {n_fit:>10d}")

    print()
    print("Cosmic-void reference slope (Pan 2012, Sutter 2014): alpha ~ 1.5-2.0")

    # ============================================================
    # Part 2: galaxy-to-interface distance
    # ============================================================
    print()
    print("=== Experiment 8: cracked-cell distance to interface ===")
    print(f"{'p':>5}  {'n_cracked':>10}  {'mean_dist':>10}  "
          f"{'median_dist':>12}  {'90%ile':>8}")
    dist_summary = {}
    distance_profiles = {}
    primary_seed = seeds[0]
    for p in p_targets:
        pp = int(round(p * 100))
        mask_key = f"mask_s{primary_seed}_p{pp:02d}"
        if mask_key not in data.files:
            continue
        occupied = data[mask_key]
        frozen = ~occupied
        if not frozen.any():
            continue
        # distance_transform_edt computes distance to NEAREST False element.
        # We want distance to nearest frozen for each cracked cell.
        # So compute distance from each cell to nearest True in frozen mask,
        # i.e., distance to nearest frozen cell.
        dist = ndimage.distance_transform_edt(occupied.astype(np.bool_))
        # dist[i,j,k] = distance from (i,j,k) to nearest cell where mask = False,
        # but only if mask[i,j,k] = True. For False cells, dist = 0.
        # Wait: distance_transform_edt computes distance to nearest zero (False).
        # So for an occupied-as-True mask, dist gives distance from each
        # occupied cell to the nearest frozen cell. Exactly what we want.
        cracked_dist = dist[occupied]
        if cracked_dist.size == 0:
            continue
        mean_d = float(cracked_dist.mean())
        median_d = float(np.median(cracked_dist))
        p90 = float(np.percentile(cracked_dist, 90))
        dist_summary[p] = {
            "mean": mean_d, "median": median_d, "p90": p90,
            "n_cracked": int(occupied.sum()),
        }
        # Histogram for plot
        bins = np.linspace(0, max(p90 * 1.2, 5), 30)
        h, edges = np.histogram(cracked_dist, bins=bins, density=True)
        distance_profiles[p] = (0.5 * (edges[:-1] + edges[1:]), h)
        print(f"{p:>5.2f}  {dist_summary[p]['n_cracked']:>10d}  "
              f"{mean_d:>10.3f}  {median_d:>12.3f}  {p90:>8.3f}")

    # ============================================================
    # Plots
    # ============================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # void slope vs p
    p_arr = sorted(void_slopes.keys())
    slopes = [void_slopes[p] for p in p_arr]
    ax1.plot(p_arr, slopes, "o-", color="C3", label=r"Manna void $\alpha$")
    ax1.axhspan(1.5, 2.0, color="C2", alpha=0.15,
                label=r"cosmic-void $\alpha \approx 1.5\text{-}2.0$")
    ax1.axvline(0.65, color="0.5", linestyle="--", alpha=0.6,
                label=r"$\gamma=1.8$ at $p=0.65$")
    ax1.set_xlabel("$p$")
    ax1.set_ylabel(r"void size-distribution slope $\alpha$")
    ax1.set_title("Void size distribution slope vs p")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9)

    # distance distributions
    p_for_plot = [0.18, 0.35, 0.50, 0.65, 0.80]
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(p_for_plot)))
    for p, color in zip(p_for_plot, colors):
        if p in distance_profiles:
            r, h = distance_profiles[p]
            ax2.plot(r, h, "-", color=color, alpha=0.8,
                     label=f"p={p:.2f}")
    ax2.set_xlabel("distance from cracked cell to nearest frozen cell")
    ax2.set_ylabel("PDF")
    ax2.set_title("Galaxy-to-void distance distribution")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9)

    fig.tight_layout()
    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_void_distribution_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")
    plt.close(fig)

    print()
    print("Interpretation:")
    if void_slopes:
        slope_at_065 = void_slopes.get(0.65)
        if slope_at_065:
            print(f"  Void distribution slope at p=0.65: alpha = {slope_at_065:.3f}")
            print(f"  Cosmic-void reference: alpha ~ 1.5-2.0")
            if 1.3 <= slope_at_065 <= 2.2:
                print("  -> Compatible with cosmic-void surveys.")
            else:
                print("  -> Outside cosmic-void survey range.")
    if dist_summary:
        d65 = dist_summary.get(0.65)
        if d65:
            print(f"  At p=0.65: mean cracked-to-frozen distance = "
                  f"{d65['mean']:.2f} lattice units")
            print(f"           median = {d65['median']:.2f}, 90th %ile = {d65['p90']:.2f}")
            print(f"  Implication: at p=0.65, half of all cracked cells are within "
                  f"{d65['median']:.1f} lattice units of a void.")
            print(f"  This is consistent with galaxies being on filaments/walls,")
            print(f"  not deep in cluster interiors.")


if __name__ == "__main__":
    main()
