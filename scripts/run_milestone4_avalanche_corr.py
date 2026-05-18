"""Avalanche temporal correlations.

Question: in each regime (pre-percolation, eye, storm-building), are
the avalanche sizes a sequence of independent draws, or are large
avalanches clustered in time?

If avalanches are independent, the autocorrelation of the size
sequence drops to zero after lag 1. If they cluster, the
autocorrelation persists across many lags (memory in the dynamics).

Three regimes:
  pre-percolation:  drops [0, drop@p_c]            ; p ~ 0    - 0.18
  eye:              drops [drop@p_c, drop@p=0.80]  ; p ~ 0.18 - 0.80
  storm-building:   drops [drop@p=0.80, end]       ; p ~ 0.80 - 0.90

For each, compute:
  - autocorrelation of avalanche size at lags 1, 2, ..., 100
  - autocorrelation of the binary indicator "size > L = 96" (rare-event
    clustering)
  - characteristic decorrelation lag

Cosmological tie-in: if observable cosmic structure-formation events
(galaxy mergers, supernovae, cluster mergers) cluster in time, that
matches an interpretation where we live in a non-trivially-correlated
regime of the substrate dynamics.

Run:
    .venv/bin/python scripts/run_milestone4_avalanche_corr.py [npz]
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))


def latest_xi_npz() -> Path:
    outdir = REPO_ROOT / "data" / "outputs"
    candidates = sorted(outdir.glob("manna_3d_xi_data_*.npz"))
    if not candidates:
        raise FileNotFoundError("No manna_3d_xi_data_*.npz in data/outputs/")
    return candidates[-1]


def autocorr(x: np.ndarray, max_lag: int) -> np.ndarray:
    """Compute autocorrelation up to max_lag using FFT."""
    x = x - x.mean()
    if x.std() == 0:
        return np.zeros(max_lag + 1)
    n = x.size
    # Pad to next power of two for FFT efficiency
    pad = 1
    while pad < 2 * n:
        pad *= 2
    F = np.fft.fft(x, n=pad)
    acf_full = np.fft.ifft(F * F.conj()).real
    acf = acf_full[: max_lag + 1] / (acf_full[0])
    return acf


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else latest_xi_npz()
    print(f"Loading {path}")
    data = np.load(path)
    L = int(data["L"])
    seeds = sorted({int(k.split("_s")[1])
                    for k in data.files if k.startswith("sizes_s")})
    primary_seed = seeds[0]
    sizes = data[f"sizes_s{primary_seed}"]
    print(f"L={L}, primary seed={primary_seed}, n_drops={sizes.size}")

    # Get the drop indices that bound each regime
    drop_at_p_c = int(data[f"drop_s{primary_seed}_p18"])
    drop_at_p80 = int(data[f"drop_s{primary_seed}_p80"])
    drop_end = sizes.size

    regimes = {
        "pre-percolation (p<0.18)": (0, drop_at_p_c),
        "eye (p in [0.18, 0.80])": (drop_at_p_c, drop_at_p80),
        "storm building (p in [0.80, 0.90])": (drop_at_p80, drop_end),
    }

    print()
    print("Regime sizes:")
    for name, (lo, hi) in regimes.items():
        print(f"  {name}: {hi - lo} drops")
    print()

    max_lag = 200
    print(f"Computing autocorrelation up to lag {max_lag}...")

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))

    # Row 1: autocorrelation of size series
    # Row 2: autocorrelation of "large avalanche" indicator (size > L)
    for col, (name, (lo, hi)) in enumerate(regimes.items()):
        win_sizes = sizes[lo:hi].astype(np.float64)
        if win_sizes.size < max_lag + 100:
            continue

        # ACF of raw sizes
        acf = autocorr(win_sizes, max_lag)
        # ACF of "size > L" indicator
        large = (win_sizes > L).astype(np.float64)
        if large.std() > 0:
            acf_large = autocorr(large, max_lag)
        else:
            acf_large = np.zeros(max_lag + 1)

        # Decorrelation lag: first lag where ACF drops below 0.1
        decor_lag_size = None
        for lag in range(1, max_lag + 1):
            if abs(acf[lag]) < 0.1:
                decor_lag_size = lag
                break

        decor_lag_large = None
        for lag in range(1, max_lag + 1):
            if abs(acf_large[lag]) < 0.1:
                decor_lag_large = lag
                break

        # Stats for printing
        large_frac = float(large.mean())
        print(f"\n{name}")
        print(f"  drops: {win_sizes.size}")
        print(f"  fraction with size > L={L}: {large_frac:.4f}")
        print(f"  size-ACF decorrelation lag (where |ACF|<0.1): "
              f"{decor_lag_size if decor_lag_size else f'>{max_lag}'}")
        print(f"  large-event-ACF decorrelation lag: "
              f"{decor_lag_large if decor_lag_large else f'>{max_lag}'}")
        print(f"  size-ACF at lag 1: {acf[1]:.3f}")
        print(f"  size-ACF at lag 10: {acf[10]:.3f}")
        print(f"  size-ACF at lag 100: {acf[100]:.3f}")

        # Plot top row: size ACF
        ax = axes[0, col]
        ax.plot(np.arange(max_lag + 1), acf, "-", color="C0")
        ax.axhline(0, color="0.6", linestyle="--", linewidth=0.5)
        ax.axhline(0.1, color="0.7", linestyle=":", linewidth=0.5)
        ax.axhline(-0.1, color="0.7", linestyle=":", linewidth=0.5)
        ax.set_xlabel("lag")
        ax.set_ylabel("ACF of size")
        ax.set_title(name, fontsize=10)
        ax.set_xscale("symlog", linthresh=1)
        ax.grid(True, alpha=0.3)

        # Plot bottom row: large-event ACF
        ax = axes[1, col]
        ax.plot(np.arange(max_lag + 1), acf_large, "-", color="C3")
        ax.axhline(0, color="0.6", linestyle="--", linewidth=0.5)
        ax.axhline(0.1, color="0.7", linestyle=":", linewidth=0.5)
        ax.axhline(-0.1, color="0.7", linestyle=":", linewidth=0.5)
        ax.set_xlabel("lag")
        ax.set_ylabel(f"ACF of (size > {L}) indicator")
        ax.set_xscale("symlog", linthresh=1)
        ax.grid(True, alpha=0.3)

    fig.suptitle("Avalanche temporal correlations across regimes")
    fig.tight_layout()
    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_avalanche_corr_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
