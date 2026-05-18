"""Cluster mass function n(s) for the 3D Manna ever-toppled snapshots.

Loads the saved npz of multi-snapshot ever-toppled masks, identifies
connected components in each, and reports the cluster size distribution
at each p value. Compares to the standard percolation-theory reference
(tau_s ~ 2.189 in 3D at p_c for random percolation) and the
cosmological Press-Schechter halo mass function form
(n(M) ~ M^-2 exp(-(M/M_star)^beta)).

Run:
    .venv/bin/python scripts/run_milestone4_mass_function.py [path/to/xi_data.npz]
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.mass_function import (  # noqa: E402
    cluster_mass_pdf,
    cluster_sizes,
)


def latest_xi_npz() -> Path:
    outdir = REPO_ROOT / "data" / "outputs"
    candidates = sorted(outdir.glob("manna_3d_xi_data_*.npz"))
    if not candidates:
        raise FileNotFoundError("No manna_3d_xi_data_*.npz in data/outputs/")
    return candidates[-1]


def fit_power_law(centers: np.ndarray, pdf: np.ndarray,
                  s_min: float, s_max: float
                  ) -> tuple[float, float, int] | None:
    sel = (centers >= s_min) & (centers <= s_max) & (pdf > 0)
    if sel.sum() < 3:
        return None
    log_s = np.log10(centers[sel])
    log_p = np.log10(pdf[sel])
    coeffs = np.polyfit(log_s, log_p, 1)
    tau = -float(coeffs[0])
    log_A = float(coeffs[1])
    return tau, log_A, int(sel.sum())


def main() -> None:
    if len(sys.argv) > 1:
        path = Path(sys.argv[1])
    else:
        path = latest_xi_npz()
    print(f"Loading {path}")
    data = np.load(path)
    L = int(data["L"])
    p_targets = list(data["p_targets"])
    p_keys = [int(round(p * 100)) for p in p_targets]
    print(f"L = {L}, p_targets = {p_targets}")

    # Pool cluster sizes across seeds for each p
    per_p_sizes: dict[float, np.ndarray] = {}
    for p, pp in zip(p_targets, p_keys):
        chunks = []
        for k in data.files:
            if k.startswith(f"mask_") and k.endswith(f"_p{pp:02d}"):
                mask = data[k]
                sizes = cluster_sizes(mask)
                chunks.append(sizes)
        if chunks:
            per_p_sizes[p] = np.concatenate(chunks)

    if not per_p_sizes:
        print("No mask data found in npz; cannot compute mass function.")
        return

    print()
    print(f"{'p':>6}  {'n_seeds':>7}  {'n_clusters':>10}  {'largest':>9}  "
          f"{'mean':>8}  {'median':>8}")
    p_values = sorted(per_p_sizes.keys())
    summary = {}
    for p in p_values:
        sizes = per_p_sizes[p]
        n_seeds = 0
        for k in data.files:
            if k.startswith("mask_") and k.endswith(f"_p{int(round(p * 100)):02d}"):
                n_seeds += 1
        summary[p] = {
            "sizes": sizes,
            "n_seeds": n_seeds,
            "n_clusters": int(sizes.size),
            "largest": int(sizes.max()),
            "mean": float(sizes.mean()),
            "median": float(np.median(sizes)),
        }
        print(f"{p:>6.2f}  {n_seeds:>7d}  {summary[p]['n_clusters']:>10d}  "
              f"{summary[p]['largest']:>9d}  {summary[p]['mean']:>8.1f}  "
              f"{summary[p]['median']:>8.1f}")

    # Power-law fit on each p's mass function over a scaling range
    print()
    print("Power-law fits n(s) = A * s^-tau_s over s in [3, largest/3]:")
    print(f"{'p':>6}  {'tau_s':>8}  {'logA':>8}  {'n_fit':>6}")
    fits = {}
    for p in p_values:
        sizes = summary[p]["sizes"]
        if sizes.size < 30:
            print(f"{p:>6.2f}  too few clusters for fit")
            continue
        centers, pdf, _ = cluster_mass_pdf(sizes, n_bins=18)
        s_min_fit = 3.0
        s_max_fit = max(summary[p]["largest"] / 3.0, 10.0)
        fit = fit_power_law(centers, pdf, s_min_fit, s_max_fit)
        if fit is None:
            print(f"{p:>6.2f}  fit failed")
            continue
        tau, log_A, n_fit = fit
        fits[p] = (tau, log_A, centers, pdf)
        print(f"{p:>6.2f}  {tau:>8.3f}  {log_A:>8.3f}  {n_fit:>6d}")

    # Plot
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(p_values)))
    fig, ax = plt.subplots(figsize=(8, 6))
    for p, color in zip(p_values, colors):
        sizes = summary[p]["sizes"]
        if sizes.size < 10:
            continue
        centers, pdf, _ = cluster_mass_pdf(sizes, n_bins=18)
        ax.loglog(centers, pdf, "o-", color=color, alpha=0.7,
                  label=f"p={p:.2f}")
        if p in fits:
            tau, log_A, _, _ = fits[p]
            r_line = np.linspace(centers.min(), centers.max(), 30)
            ax.loglog(r_line, 10 ** (log_A) * r_line ** (-tau),
                      "--", color=color, alpha=0.4)

    # Press-Schechter reference (alpha=2 small-mass)
    s_ref = np.logspace(0.5, 3.5, 50)
    ps_ref = s_ref ** (-2.0) * np.exp(-(s_ref / 300.0))
    ps_ref *= 0.01 / ps_ref.max()
    ax.loglog(s_ref, ps_ref, ":", color="0.4", alpha=0.7,
              label=r"Press-Schechter $\alpha=2$ ref shape")

    ax.set_xlabel("cluster size $s$  (number of toppled sites)")
    ax.set_ylabel(r"$n(s)$  (cluster-size PDF, log-binned)")
    ax.set_title(f"3D Manna cluster mass function across snapshots (L={L})")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_mass_function_{stamp}.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
