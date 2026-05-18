"""Refit the void size distribution slope with various minimum-size cuts.

The morning analysis fit the void slope over s in [2, max/100], giving
alpha=3.0 at p=0.65 vs the cosmic-void reference 1.5-2.0. The likely
reason: pixel-scale "voids" (singletons, pairs) dominate the fit, while
real cosmic-void surveys exclude voids smaller than ~5 Mpc.

This script refits with a progressively larger minimum cluster size and
sees whether the cosmic-void range is reached. If it is, the original
disagreement was a counting artifact. If it is not, there is a real
shape disagreement and we need bigger statistics.

Cut levels tested: 2, 5, 10, 20, 50, 100 cells. With our scale-mapping
of ~8.5 Mpc per cell, the cuts correspond approximately to:
  2 cells  ~ 17 Mpc (below survey cuts)
  5 cells  ~ 43 Mpc (around survey lower bound)
  10 cells ~ 85 Mpc (typical "interesting void" size)
  20 cells ~ 170 Mpc (large voids only)
  50 cells ~ 425 Mpc (very large, sparse)
  100 cells ~ 850 Mpc (only the spanning void at low p)

Run:
    .venv/bin/python scripts/run_milestone4_void_slope_recut.py [npz]
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.mass_function import cluster_mass_pdf, cluster_sizes  # noqa: E402


CUTS = [2, 5, 10, 20, 50, 100]


def latest_xi_npz() -> Path:
    outdir = REPO_ROOT / "data" / "outputs"
    candidates = sorted(outdir.glob("manna_3d_xi_data_*.npz"))
    if not candidates:
        raise FileNotFoundError("No manna_3d_xi_data_*.npz in data/outputs/")
    return candidates[-1]


def fit_slope(sizes: np.ndarray, s_min: float, s_max: float,
              n_bins: int = 20) -> tuple[float, int] | None:
    """Log-log slope fit on cluster sizes >= s_min, <= s_max."""
    if sizes.size == 0:
        return None
    sizes = sizes[(sizes >= s_min) & (sizes <= s_max)]
    if sizes.size < 30:
        return None
    centers, pdf, counts = cluster_mass_pdf(sizes, n_bins=n_bins, min_size=int(s_min))
    sel = (centers > 0) & (pdf > 0)
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

    # For each p, pool void cluster sizes across seeds; then fit slope at
    # each min-size cut. Cap upper end at max/10 so the spanning giant
    # void doesn't drive the fit.
    print()
    header = "p   " + "  ".join(f"cut>={c:>4}" for c in CUTS)
    print(header)
    table: dict[float, dict[int, float | None]] = {}
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
        table[p] = {}
        row = f"{p:>4.2f}"
        for cut in CUTS:
            # cap at max/10 to leave room for power-law decay
            s_max = max(float(sizes.max()) / 10.0, 5.0 * cut)
            fit = fit_slope(sizes, float(cut), s_max)
            if fit is None:
                table[p][cut] = None
                row += "        --"
            else:
                slope, nbins = fit
                table[p][cut] = slope
                row += f"     {slope:5.2f}"
        print(row)

    print()
    print("Cosmic-void reference (Pan 2012, Sutter 2014): alpha ~ 1.5-2.0")

    # --- find the cut at which p=0.65 enters cosmic-void range ---
    print()
    if 0.65 in table:
        print(f"At p=0.65:")
        for cut in CUTS:
            v = table[0.65][cut]
            if v is None:
                print(f"  cut>={cut}: insufficient data")
            else:
                in_range = "  <-- in cosmic-void range" if 1.3 <= v <= 2.2 else ""
                print(f"  cut>={cut} cells (~{cut * 8.5:.0f} Mpc): alpha = {v:.3f}{in_range}")

    # --- plot ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # alpha vs p, one curve per cut
    cuts_sorted = sorted(CUTS)
    colors = plt.cm.plasma(np.linspace(0.15, 0.85, len(cuts_sorted)))
    for cut, color in zip(cuts_sorted, colors):
        p_arr = []
        slopes = []
        for p in p_targets:
            v = table.get(p, {}).get(cut)
            if v is not None:
                p_arr.append(p)
                slopes.append(v)
        if p_arr:
            ax1.plot(p_arr, slopes, "o-", color=color, alpha=0.8,
                     label=fr"min size $\geq$ {cut} cells")
    ax1.axhspan(1.5, 2.0, color="C2", alpha=0.15,
                label=r"cosmic-void $\alpha \approx 1.5\text{-}2.0$")
    ax1.axvline(0.65, color="0.5", linestyle="--", alpha=0.6,
                label=r"$\gamma=1.8$ crossing at $p=0.65$")
    ax1.set_xlabel("$p$")
    ax1.set_ylabel(r"void size-distribution slope $\alpha$")
    ax1.set_title("Void slope vs p, by minimum-size cut")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8)

    # alpha vs cut, one curve per p value
    plot_p = [0.50, 0.55, 0.60, 0.65, 0.70, 0.80]
    colors2 = plt.cm.viridis(np.linspace(0.15, 0.85, len(plot_p)))
    for p, color in zip(plot_p, colors2):
        cut_arr = []
        slopes = []
        for cut in cuts_sorted:
            v = table.get(p, {}).get(cut)
            if v is not None:
                cut_arr.append(cut)
                slopes.append(v)
        if cut_arr:
            ax2.semilogx(cut_arr, slopes, "o-", color=color, alpha=0.8,
                         label=f"p={p:.2f}")
    ax2.axhspan(1.5, 2.0, color="C2", alpha=0.15,
                label=r"cosmic-void $\alpha$")
    ax2.set_xlabel("minimum cluster size (cells)")
    ax2.set_ylabel(r"slope $\alpha$")
    ax2.set_title("Slope vs minimum-size cut, at fixed p")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_void_slope_recut_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
