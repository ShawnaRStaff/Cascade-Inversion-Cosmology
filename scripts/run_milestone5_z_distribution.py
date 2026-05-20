"""Analyze the z distribution at saturation.

Stage 4 saved the final z field of each seed at p=1.0 after extensive
post-arrest dynamics. This script characterizes:

  - Per-cell z distribution (histogram)
  - Spatial autocorrelation of z at saturation
  - Variance and shape
  - Whether z values are correlated with cell position (i.e., is the
    terminal state spatially structured or random?)

Run:
    .venv/bin/python scripts/run_milestone5_z_distribution.py [npz]
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]


def latest_past_arrest_npz() -> Path:
    outdir = REPO_ROOT / "data" / "outputs"
    candidates = sorted(outdir.glob("manna_3d_past_arrest_data_*.npz"))
    if not candidates:
        raise FileNotFoundError("No manna_3d_past_arrest_data_*.npz")
    return candidates[-1]


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else latest_past_arrest_npz()
    print(f"Loading {path}")
    data = np.load(path)
    L = int(data["L"])
    L3 = L ** 3
    seeds = sorted({int(k.split("_")[0][1:])
                    for k in data.files if k.startswith("s") and "_final_z" in k})
    print(f"L={L} (L^3={L3}), seeds={seeds}")

    # Per-seed z statistics
    print()
    print(f"{'seed':>5}  {'z_mean':>8}  {'z_std':>8}  {'z_min':>5}  "
          f"{'z_max':>5}  {'n_z=0':>7}  {'n_z=1':>7}  {'n_z=2+':>8}")
    final_zs = {}
    for s in seeds:
        z = data[f"s{s}_final_z"]
        final_zs[s] = z
        z_flat = z.ravel()
        n0 = int((z_flat == 0).sum())
        n1 = int((z_flat == 1).sum())
        n2p = int((z_flat >= 2).sum())
        print(f"{s:>5}  {float(z_flat.mean()):>8.4f}  "
              f"{float(z_flat.std()):>8.4f}  "
              f"{int(z_flat.min()):>5d}  {int(z_flat.max()):>5d}  "
              f"{n0:>7d}  {n1:>7d}  {n2p:>8d}")

    print()
    print("z distribution interpretation:")
    print("  z_c = 2 means a cell topples at z >= 2 in Manna.")
    print("  In dynamic equilibrium most cells must be at z=0 or z=1.")
    print("  Cells at z=2+ are unstable and would topple in the next sweep.")

    # Spatial autocorrelation of z field via FFT
    print()
    print("Spatial autocorrelation (1/e decay distance):")
    print(f"{'seed':>5}  {'r_corr (lattice units)':>26}")
    for s in seeds:
        z = final_zs[s].astype(np.float64) - final_zs[s].mean()
        # 3D FFT autocorrelation
        F = np.fft.fftn(z)
        acf = np.fft.ifftn(F * F.conj()).real
        # Normalize so that acf[0,0,0] = 1
        acf = acf / acf[0, 0, 0]
        # Radially average for a 1D correlation length
        # Find distance at which acf drops below 1/e
        coords = np.arange(L)
        delta = np.minimum(coords, L - coords)
        dx, dy, dz = np.meshgrid(delta, delta, delta, indexing="ij")
        r = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2).ravel()
        acf_flat = acf.ravel()
        # Bin
        r_bins = np.linspace(0, L / 2, 30)
        bin_idx = np.digitize(r, r_bins) - 1
        radial_acf = np.full(len(r_bins) - 1, np.nan)
        for i in range(len(r_bins) - 1):
            mask = bin_idx == i
            if mask.any():
                radial_acf[i] = float(acf_flat[mask].mean())
        # Find 1/e crossing
        bin_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
        thr = 1.0 / np.e
        r_corr = None
        for i, a in enumerate(radial_acf):
            if np.isfinite(a) and a < thr:
                r_corr = float(bin_centers[i])
                break
        print(f"{s:>5}  {r_corr if r_corr is not None else 'no decay':>26}")

    # Plot: histograms
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # Panel 1: z histogram
    ax = axes[0, 0]
    for s in seeds:
        z_flat = final_zs[s].ravel()
        max_z = int(z_flat.max())
        bins = np.arange(0, max_z + 2) - 0.5
        ax.hist(z_flat, bins=bins, alpha=0.6,
                label=f"seed {s} (mean={z_flat.mean():.3f})")
    ax.set_xlabel("z (grains per cell)")
    ax.set_ylabel("number of cells")
    ax.set_yscale("log")
    ax.set_title("Final z distribution (log y)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)

    # Panel 2: mid-slice through final_z
    ax = axes[0, 1]
    s = seeds[0]
    z_slice = final_zs[s][L // 2, :, :]
    im = ax.imshow(z_slice, cmap="hot", interpolation="nearest", vmin=0, vmax=2)
    plt.colorbar(im, ax=ax)
    ax.set_title(f"Final z mid-slice (seed {s})")
    ax.set_xticks([])
    ax.set_yticks([])

    # Panel 3: radial autocorrelation
    ax = axes[1, 0]
    for s in seeds:
        z = final_zs[s].astype(np.float64) - final_zs[s].mean()
        F = np.fft.fftn(z)
        acf = np.fft.ifftn(F * F.conj()).real
        acf = acf / acf[0, 0, 0]
        coords = np.arange(L)
        delta = np.minimum(coords, L - coords)
        dx, dy, dz = np.meshgrid(delta, delta, delta, indexing="ij")
        r = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2).ravel()
        acf_flat = acf.ravel()
        r_bins = np.linspace(0, L / 2, 30)
        bin_idx = np.digitize(r, r_bins) - 1
        radial_acf = np.array([acf_flat[bin_idx == i].mean() if (bin_idx == i).any() else np.nan
                               for i in range(len(r_bins) - 1)])
        bin_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
        ax.plot(bin_centers, radial_acf, "-", alpha=0.8, label=f"seed {s}")
    ax.axhline(1.0 / np.e, color="0.5", linestyle=":", alpha=0.7,
               label="1/e level")
    ax.set_xlabel("r (lattice units)")
    ax.set_ylabel("z autocorrelation")
    ax.set_title("Spatial autocorrelation of z")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)

    # Panel 4: comparison of z mean vs ever_toppled mid-slice
    ax = axes[1, 1]
    s = seeds[0]
    ever = data[f"s{s}_ever_toppled"]
    ax.imshow(ever[L // 2, :, :], cmap="Greys_r", interpolation="nearest")
    ax.set_title(f"ever_toppled mid-slice (seed {s})\np_final={ever.mean():.5f}")
    ax.set_xticks([])
    ax.set_yticks([])

    fig.tight_layout()
    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_z_distribution_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
