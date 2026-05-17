"""Milestone 2 demo: 2D Manna sandpile with FSS and cluster geometry.

Runs the Manna model for several L, extracts (tau, D) and (alpha, z) via
the moment method, and fits the fractal dimension D_f of avalanche
support clusters from the (area, R_g) population at the largest L.

Reference values for 2D Manna (Lubeck 2000; Chessa, Marinari, Vespignani,
Zapperi 1999):
  tau    ~ 1.275
  D      ~ 2.76
  alpha  ~ 1.51
  z      ~ 1.55
  D_f    ~ 2.0 (compact support; cluster fills its bounding box)

The conservation relation D*(2-tau) = 1 holds in Manna because the bulk
toppling rule conserves grains exactly (only the boundary loses them).

Run with the repo venv:
    .venv/bin/python scripts/run_milestone2.py

Outputs:
- oslo_2d_*.npz: raw steady-state samples + cluster geometry data
- oslo_2d_size_fss_*.png: P(s) overlay + data collapse
- oslo_2d_size_moments_*.png: log<s^k> vs log L
- oslo_2d_duration_fss_*.png, oslo_2d_duration_moments_*.png: same for T
- oslo_2d_cluster_fractal_*.png: area vs R_g log-log with D_f fit
"""

from __future__ import annotations

import sys
import time
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.analysis import log_binned_pdf  # noqa: E402
from void_cascade.cluster_geometry import fit_fractal_dimension  # noqa: E402
from void_cascade.sandpile_2d import run, run_with_clusters  # noqa: E402
from void_cascade.scaling import (  # noqa: E402
    collapsed_distribution,
    moment_scaling,
)


# Manna 2D fills to mean density ~ 0.72, so the initial fill needs
# roughly 0.72 * L^2 drives; we double it for the transient. n_drops
# leaves at least ~30k steady-state samples per L.
RUN_PLAN = [
    # (L, n_drops, seed, track_clusters_in_steady_state)
    # n_drops = 2*L^2 transient + ~30k-60k steady-state samples. The
    # bottleneck is L=128 (~20 min); the other sizes finish quickly.
    (16, 25_000, 200, False),
    (32, 40_000, 201, False),
    (64, 50_000, 202, False),
    (128, 90_000, 203, True),
]


def simulate_all() -> dict:
    results = {}
    for L, n_drops, seed, track in RUN_PLAN:
        transient = 2 * L * L
        if n_drops <= transient + 5_000:
            raise ValueError(
                f"L={L}: n_drops={n_drops} leaves too few steady-state samples "
                f"after transient {transient}."
            )
        t0 = time.time()
        if track:
            state, sizes, durations, areas, rgyr = run_with_clusters(
                L=L, n_drops=n_drops, seed=seed, burn_in=transient
            )
        else:
            state, sizes, durations = run(L=L, n_drops=n_drops, seed=seed)
            areas = np.zeros_like(sizes)
            rgyr = np.full(n_drops, np.nan)
        elapsed = time.time() - t0
        results[L] = {
            "sizes_ss": sizes[transient:],
            "durations_ss": durations[transient:],
            "areas_ss": areas[transient:],
            "rgyr_ss": rgyr[transient:],
            "mean_density": float(state.z.mean()),
            "grains_lost": int(state.grains_lost),
            "transient": transient,
            "n_drops": n_drops,
            "elapsed_seconds": elapsed,
            "tracked": track,
        }
        s_ss = results[L]["sizes_ss"]
        print(
            f"  L={L:3d}  n={n_drops:6d}  transient={transient:6d}  "
            f"steady={s_ss.size:6d}  max_s={int(s_ss.max()):7d}  "
            f"<rho>={state.z.mean():.3f}  t={elapsed:.1f}s"
        )
    return results


def fit_and_summarize(results: dict) -> tuple:
    sizes_by_L = {L: r["sizes_ss"] for L, r in results.items()}
    durations_by_L = {L: r["durations_ss"] for L, r in results.items()}

    size_fit = moment_scaling(sizes_by_L, k_values=[1, 2, 3, 4])
    duration_fit = moment_scaling(durations_by_L, k_values=[1, 2, 3, 4])

    tau = size_fit.tail_exponent
    D = size_fit.cutoff_exponent
    print()
    print("Avalanche size scaling (2D Manna reference: tau~1.275, D~2.76)")
    print(f"  tau = {tau:.3f} +/- {size_fit.tail_exponent_err:.3f}")
    print(f"  D   = {D:.3f} +/- {size_fit.cutoff_exponent_err:.3f}")
    print("  sigma_k slopes (predicted D*(1+k-tau)):")
    for k, s, se in zip(size_fit.k_values, size_fit.sigma_k, size_fit.sigma_k_err):
        pred = D * (1.0 + k - tau)
        print(f"    k={int(k)}: sigma={s:.3f} +/- {se:.3f}  (predicted {pred:.3f})")
    print(f"  conservation check D*(2-tau) ~ 1: {D * (2 - tau):.3f}")

    alpha = duration_fit.tail_exponent
    z = duration_fit.cutoff_exponent
    print()
    print("Avalanche duration scaling (2D Manna reference: alpha~1.51, z~1.55)")
    print(f"  alpha = {alpha:.3f} +/- {duration_fit.tail_exponent_err:.3f}")
    print(f"  z     = {z:.3f} +/- {duration_fit.cutoff_exponent_err:.3f}")

    return size_fit, duration_fit


def make_dist_plots(results, fit, observable: str, outdir: Path, stamp: str):
    """Generic raw-distribution + data-collapse plot."""
    Ls = sorted(results.keys())
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(Ls)))
    key = "sizes_ss" if observable == "size" else "durations_ss"
    label_x = "$s$" if observable == "size" else "$T$"
    a_sym = r"\tau" if observable == "size" else r"\alpha"
    B_sym = "D" if observable == "size" else "z"

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    n_bins = 40 if observable == "size" else 30
    for L, color in zip(Ls, colors):
        x = results[L][key]
        centers, pdf, _ = log_binned_pdf(x, n_bins=n_bins)
        ax.loglog(centers, pdf, "o", color=color, alpha=0.7, label=f"L={L}")
    ax.set_xlabel(f"avalanche {observable} {label_x}")
    ax.set_ylabel(f"$P({label_x[1:-1]})$")
    ax.set_title(f"Raw {observable} distributions across L")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    ax = axes[1]
    a = fit.tail_exponent
    B = fit.cutoff_exponent
    for L, color in zip(Ls, colors):
        x_s, y_s = collapsed_distribution(
            results[L][key], float(L), a, B, n_bins=n_bins
        )
        ax.loglog(x_s, y_s, "o-", color=color, alpha=0.7, label=f"L={L}")
    ax.set_xlabel(f"${label_x[1:-1]} / L^{{{B_sym}}}$")
    ax.set_ylabel(f"${label_x[1:-1]}^{{{a_sym}}} P({label_x[1:-1]})$")
    ax.set_title(rf"Data collapse: ${a_sym}={a:.3f}$, ${B_sym}={B:.3f}$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    p = outdir / f"manna_2d_{observable}_fss_{stamp}.png"
    fig.savefig(p, dpi=150)
    print(f"  wrote {p}")
    plt.close(fig)


def make_moment_plot(fit, observable: str, outdir: Path, stamp: str):
    fig, ax = plt.subplots(figsize=(6, 5))
    L_arr = fit.L_values
    for i, k in enumerate(fit.k_values):
        ax.loglog(L_arr, fit.moments[i], "o-", label=rf"$k={int(k)}$")
    ax.set_xlabel("$L$")
    sym = "s" if observable == "size" else "T"
    ax.set_ylabel(rf"$\langle {sym}^k \rangle$")
    ax.set_title(f"Moment scaling for avalanche {observable}s")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    p = outdir / f"manna_2d_{observable}_moments_{stamp}.png"
    fig.savefig(p, dpi=150)
    print(f"  wrote {p}")
    plt.close(fig)


def make_fractal_plot(results, outdir: Path, stamp: str) -> tuple[float, float] | None:
    """Fit D_f from area vs R_g at the largest tracked L."""
    tracked_Ls = [L for L, r in results.items() if r["tracked"]]
    if not tracked_Ls:
        print("  no L tracked clusters, skipping fractal-dimension plot")
        return None
    L = max(tracked_Ls)
    areas = results[L]["areas_ss"]
    rgyr = results[L]["rgyr_ss"]

    # Restrict the fit to R_g well inside the lattice; the upper knee is
    # where R_g ~ L/4 due to the finite-size cutoff of the avalanche extent.
    r_max = L / 4.0
    D_f, D_f_err, centers, mean_area = fit_fractal_dimension(
        areas, rgyr, r_min=2.0, r_max=r_max, n_bins=18
    )
    print()
    print(f"Cluster fractal dimension (from L={L} cluster geometry, "
          f"reference D_f ~ 2.0 for compact 2D Manna)")
    print(f"  D_f = {D_f:.3f} +/- {D_f_err:.3f}  (fit R_g in [2, L/4]={r_max:.1f})")

    fig, ax = plt.subplots(figsize=(7, 5))
    valid = np.isfinite(rgyr) & (areas > 0)
    ax.loglog(rgyr[valid], areas[valid], ".", color="0.6", alpha=0.25,
              markersize=2, label="raw avalanches")
    ax.loglog(centers, mean_area, "o", color="C1", label="binned mean")
    rline = np.array([centers.min(), centers.max()])
    a_pred = mean_area[0] * (rline / centers[0]) ** D_f
    ax.loglog(rline, a_pred, "r-",
              label=rf"fit $a\propto R_g^{{{D_f:.3f}}}$")
    ax.set_xlabel(r"radius of gyration $R_g$")
    ax.set_ylabel(r"avalanche area $a$")
    ax.set_title(f"2D Manna cluster fractal dimension, L={L}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    p = outdir / f"manna_2d_cluster_fractal_{stamp}.png"
    fig.savefig(p, dpi=150)
    print(f"  wrote {p}")
    plt.close(fig)
    return D_f, D_f_err


def save_raw(results, outdir: Path, stamp: str) -> Path:
    save_kwargs: dict = {}
    for L, r in results.items():
        save_kwargs[f"sizes_L{L}"] = r["sizes_ss"]
        save_kwargs[f"durations_L{L}"] = r["durations_ss"]
        if r["tracked"]:
            save_kwargs[f"areas_L{L}"] = r["areas_ss"]
            save_kwargs[f"rgyr_L{L}"] = r["rgyr_ss"]
    path = outdir / f"manna_2d_fss_data_{stamp}.npz"
    np.savez_compressed(path, **save_kwargs)
    return path


def main() -> None:
    print("Simulating 2D Manna across multiple L:")
    t_start = time.time()
    results = simulate_all()
    sim_elapsed = time.time() - t_start
    print(f"Total simulation time: {sim_elapsed:.1f}s")

    size_fit, duration_fit = fit_and_summarize(results)

    outdir = REPO_ROOT / "data" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    npz_path = save_raw(results, outdir, stamp)
    print(f"\nSaved raw steady-state samples: {npz_path}")

    print("\nWriting plots:")
    make_dist_plots(results, size_fit, "size", outdir, stamp)
    make_moment_plot(size_fit, "size", outdir, stamp)
    make_dist_plots(results, duration_fit, "duration", outdir, stamp)
    make_moment_plot(duration_fit, "duration", outdir, stamp)
    make_fractal_plot(results, outdir, stamp)


if __name__ == "__main__":
    main()
