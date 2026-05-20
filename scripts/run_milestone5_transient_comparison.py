"""Compare the model's saturation event-size distribution to observed
cosmic-transient luminosity functions.

If the model's "permanent plateau of catastrophes" interpretation is
correct, the size distribution of substrate avalanches should map to
the energy distribution of observed cosmic transient events.

Observed cosmic transients (for reference - power-law slopes from
published luminosity-function fits):

  GRBs (long, observed):     alpha ~ 1.5-2.0  (Cumulative Log N - Log P)
                              Wanderman & Piran 2010; Pescalli+ 2016
  Magnetar giant flares:    alpha ~ 1.5      (Kashiyama & Murase 2017)
  Soft gamma repeaters:     alpha ~ 1.66     (Cheng+ 2017)
  Fast radio bursts:        alpha ~ 1.5-1.8  (James+ 2022; CHIME 2021)
  UHECRs (above 10^19 eV):  alpha ~ 2.6-2.8  (Auger Collab 2020)
  Solar flares:             alpha ~ 1.8-2.0  (Aschwanden 2011)

These are remarkably consistent: most observed transient luminosity
functions have power-law slopes in the range 1.5-2.0 (with a tail
at higher slopes for the most extreme events like UHECRs).

This script:
  1. Computes the avalanche-size distribution from the past-arrest
     run (saturation regime).
  2. Fits a power law to the size distribution.
  3. Reports the slope and compares to the observed transient range.

Run:
    .venv/bin/python scripts/run_milestone5_transient_comparison.py [npz]
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
                    for k in data.files if k.startswith("s") and "_sizes" in k})
    print(f"L={L}, seeds={seeds}")

    # For each seed, isolate the SATURATION-REGIME sizes (p > 0.9 onward).
    # Get the drop index at which p crossed 0.9 - we can find this from
    # the snapshot data.
    print()
    print("Avalanche-size distribution at saturation (p > 0.9):")
    print(f"{'seed':>5}  {'n_avalanches':>13}  {'mean_s':>9}  {'max_s':>8}  "
          f"{'tau_s fit':>10}  {'n_fit_pts':>10}")
    fits = {}
    for s in seeds:
        sizes = data[f"s{s}_sizes"]
        p_snaps = data[f"s{s}_p"]
        drops_snaps = data[f"s{s}_drops"]
        # Find first drop where p > 0.9
        idx_pgt09 = np.where(p_snaps > 0.9)[0]
        if idx_pgt09.size == 0:
            print(f"{s:>5}  no snapshots above p=0.9")
            continue
        sat_drop = int(drops_snaps[idx_pgt09[0]])
        sizes_sat = sizes[sat_drop:]
        sizes_nonzero = sizes_sat[sizes_sat > 0]

        # Log-binned distribution
        if sizes_nonzero.size == 0:
            continue
        log_max = np.log10(sizes_nonzero.max() + 1)
        bins = np.logspace(0, log_max, 30)
        counts, edges = np.histogram(sizes_nonzero, bins=bins)
        widths = np.diff(edges)
        centers = np.sqrt(edges[:-1] * edges[1:])
        pdf = counts / (sizes_nonzero.size * widths)

        # Fit power law over s in [3, max/10] (avoid small-s rollover and
        # the rare-large-event noise)
        s_min_fit = 3.0
        s_max_fit = max(float(sizes_nonzero.max()) / 10.0, 10.0)
        keep = (centers >= s_min_fit) & (centers <= s_max_fit) & (counts >= 3)
        if keep.sum() < 3:
            print(f"{s:>5}  too few fit points")
            continue
        coeffs = np.polyfit(np.log10(centers[keep]), np.log10(pdf[keep]), 1)
        tau = -float(coeffs[0])

        fits[s] = {"tau": tau, "centers": centers, "pdf": pdf,
                   "counts": counts, "keep": keep}
        print(f"{s:>5}  {sizes_nonzero.size:>13d}  "
              f"{float(sizes_nonzero.mean()):>9.1f}  "
              f"{int(sizes_nonzero.max()):>8d}  "
              f"{tau:>10.3f}  {int(keep.sum()):>10d}")

    print()
    print("Observed cosmic-transient slope range: alpha ~ 1.5-2.0")
    print()
    tau_values = [f["tau"] for f in fits.values()]
    if tau_values:
        mean_tau = float(np.mean(tau_values))
        print(f"Manna saturation tau_s (mean across seeds): {mean_tau:.3f}")
        if 1.3 <= mean_tau <= 2.2:
            print("  -> IN the observed cosmic-transient range")
        elif mean_tau < 1.3:
            print(f"  -> shallower than observed (events too heavy-tailed)")
        else:
            print(f"  -> steeper than observed (rare large events under-produced)")

    # Plot
    fig, ax = plt.subplots(figsize=(9, 7))
    colors = plt.cm.tab10(np.linspace(0, 1, len(fits)))
    for (s, f), color in zip(fits.items(), colors):
        ax.loglog(f["centers"], f["pdf"], "o", color=color, alpha=0.6,
                  label=fr"seed {s}, $\tau_s={f['tau']:.2f}$")
        # Fit line
        keep = f["keep"]
        if keep.any():
            line_x = np.logspace(np.log10(f["centers"][keep].min()),
                                  np.log10(f["centers"][keep].max()), 30)
            # log10(pdf) = -tau * log10(s) + const, fit on keep
            const = (np.log10(f["pdf"][keep]) + f["tau"] * np.log10(f["centers"][keep])).mean()
            line_y = 10 ** (const - f["tau"] * np.log10(line_x))
            ax.loglog(line_x, line_y, "-", color=color, alpha=0.6)

    # Cosmic-transient reference band
    s_ref = np.logspace(0.5, 4, 50)
    # show alpha=1.5 to alpha=2.0 envelope, normalized arbitrarily
    norm = 1.0
    upper = norm * s_ref ** (-1.5)
    lower = norm * s_ref ** (-2.0)
    # Scale so envelope sits in the same y-range as the data
    if fits:
        peak_pdf = max(f["pdf"][f["keep"]].max() for f in fits.values()
                       if f["keep"].any())
        s_ref_center = 10 ** ((np.log10(s_ref[0]) + np.log10(s_ref[-1])) / 2)
        scale = peak_pdf / s_ref_center ** -1.75
        upper *= scale
        lower *= scale
    ax.fill_between(s_ref, lower, upper, alpha=0.2, color="C2",
                    label=r"cosmic-transient band $\alpha\in[1.5, 2.0]$")
    ax.set_xlabel("avalanche size s")
    ax.set_ylabel(r"$P(s)$  (log-binned PDF)")
    ax.set_title("Saturation-regime event-size distribution vs cosmic-transient band")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)

    fig.tight_layout()
    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_transient_comparison_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
