"""Milestone 1 finite-size scaling: 1D Oslo across several L.

Runs the Oslo model for L in {32, 64, 128, 256}, extracts the avalanche size
and duration distributions in steady state, and fits (tau, D) and (alpha, z)
by the moment method. Produces:

  * Overlay of P(s) curves for each L
  * Data collapse plot s^tau P(s) vs s/L^D
  * Moment scaling plot log<s^k> vs log L for k in [1, 2, 3, 4]
  * Same three plots for durations T
  * A summary printout with (tau, D, alpha, z) and the consistency check
    D*(2-tau) ~ 1 (the conservation relation for boundary-driven SOC).

Reference values for 1D Oslo (Pruessner & Jensen 2003):
  tau = 1.557(7), D = 2.245(10), alpha ~ 1.76, z = 1.422(7).

Raw simulation data is saved as .npz so the analysis can be rerun without
resimulating.

Run with the repo venv:
    .venv/bin/python scripts/run_milestone1_fss.py
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
from void_cascade.sandpile_1d import run  # noqa: E402
from void_cascade.scaling import (  # noqa: E402
    collapsed_distribution,
    moment_scaling,
)


# Per-L number of drops. Transient is set inside the run (2*L^2) and the
# n_drops choice has to leave at least ~3e4 steady-state samples.
RUN_PLAN = [
    # (L, n_drops, seed)
    (32, 60_000, 101),
    (64, 80_000, 102),
    (128, 120_000, 103),
    (256, 250_000, 104),
]


def simulate_all() -> dict:
    results = {}
    for L, n_drops, seed in RUN_PLAN:
        transient = 2 * L * L
        if n_drops <= transient + 1000:
            raise ValueError(
                f"L={L}: n_drops={n_drops} leaves no steady-state samples "
                f"after transient {transient}."
            )
        t0 = time.time()
        state, sizes, durations = run(L=L, n_drops=n_drops, seed=seed)
        elapsed = time.time() - t0
        s_ss = sizes[transient:]
        T_ss = durations[transient:]
        results[L] = {
            "sizes_full": sizes,
            "durations_full": durations,
            "sizes_ss": s_ss,
            "durations_ss": T_ss,
            "mean_slope": float(state.z.mean()),
            "grains_lost": int(state.grains_lost),
            "transient": transient,
            "n_drops": n_drops,
            "elapsed_seconds": elapsed,
        }
        print(
            f"  L={L:3d}  n={n_drops:6d}  transient={transient:6d}  "
            f"steady={s_ss.size:6d}  max_s={int(s_ss.max()):7d}  "
            f"<slope>={state.z.mean():.3f}  t={elapsed:.1f}s"
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
    print("Avalanche size scaling (1D Oslo reference: tau=1.557(7), D=2.245(10))")
    print(f"  tau = {tau:.3f} +/- {size_fit.tail_exponent_err:.3f}")
    print(f"  D   = {D:.3f} +/- {size_fit.cutoff_exponent_err:.3f}")
    print("  sigma_k slopes (D*(1+k-tau)):")
    for k, s, se in zip(
        size_fit.k_values, size_fit.sigma_k, size_fit.sigma_k_err
    ):
        predicted = D * (1.0 + k - tau)
        print(f"    k={int(k)}: sigma={s:.3f} +/- {se:.3f}  (predicted {predicted:.3f})")
    print(f"  conservation check D*(2-tau) ~ 1: {D * (2 - tau):.3f}")

    alpha = duration_fit.tail_exponent
    z = duration_fit.cutoff_exponent
    print()
    print("Avalanche duration scaling (1D Oslo reference: alpha~1.76, z=1.422(7))")
    print(f"  alpha = {alpha:.3f} +/- {duration_fit.tail_exponent_err:.3f}")
    print(f"  z     = {z:.3f} +/- {duration_fit.cutoff_exponent_err:.3f}")
    print("  sigma_k slopes (z*(1+k-alpha)):")
    for k, s, se in zip(
        duration_fit.k_values, duration_fit.sigma_k, duration_fit.sigma_k_err
    ):
        predicted = z * (1.0 + k - alpha)
        print(f"    k={int(k)}: sigma={s:.3f} +/- {se:.3f}  (predicted {predicted:.3f})")

    return size_fit, duration_fit


def make_plots(results: dict, size_fit, duration_fit, outdir: Path, stamp: str) -> None:
    Ls = sorted(results.keys())
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(Ls)))

    # --- Sizes: raw distributions ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    for L, color in zip(Ls, colors):
        s_ss = results[L]["sizes_ss"]
        centers, pdf, _ = log_binned_pdf(s_ss, n_bins=40)
        ax.loglog(centers, pdf, "o", color=color, alpha=0.7, label=f"L={L}")
    ax.set_xlabel("avalanche size $s$")
    ax.set_ylabel("$P(s)$")
    ax.set_title("Raw distributions across L")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    # --- Sizes: collapse ---
    ax = axes[1]
    tau = size_fit.tail_exponent
    D = size_fit.cutoff_exponent
    for L, color in zip(Ls, colors):
        x_s, y_s = collapsed_distribution(
            results[L]["sizes_ss"], float(L), tau, D, n_bins=40
        )
        ax.loglog(x_s, y_s, "o-", color=color, alpha=0.7, label=f"L={L}")
    ax.set_xlabel(r"$s / L^D$")
    ax.set_ylabel(r"$s^{\tau} P(s)$")
    ax.set_title(rf"Data collapse: $\tau={tau:.3f}$, $D={D:.3f}$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    p = outdir / f"oslo_1d_size_fss_{stamp}.png"
    fig.savefig(p, dpi=150)
    print(f"  wrote {p}")
    plt.close(fig)

    # --- Sizes: moment scaling ---
    fig, ax = plt.subplots(figsize=(6, 5))
    L_arr = size_fit.L_values
    for i, k in enumerate(size_fit.k_values):
        ax.loglog(L_arr, size_fit.moments[i], "o-", label=rf"$k={int(k)}$")
        # plot fitted slope
        ax.loglog(
            L_arr,
            10 ** (
                size_fit.sigma_k[i] * np.log10(L_arr)
                + np.log10(size_fit.moments[i, 0])
                - size_fit.sigma_k[i] * np.log10(L_arr[0])
            ),
            "--",
            alpha=0.4,
        )
    ax.set_xlabel("$L$")
    ax.set_ylabel(r"$\langle s^k \rangle$")
    ax.set_title("Moment scaling for avalanche sizes")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    p = outdir / f"oslo_1d_size_moments_{stamp}.png"
    fig.savefig(p, dpi=150)
    print(f"  wrote {p}")
    plt.close(fig)

    # --- Durations: collapse + raw ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    for L, color in zip(Ls, colors):
        T_ss = results[L]["durations_ss"]
        centers, pdf, _ = log_binned_pdf(T_ss, n_bins=30)
        ax.loglog(centers, pdf, "o", color=color, alpha=0.7, label=f"L={L}")
    ax.set_xlabel("avalanche duration $T$")
    ax.set_ylabel("$P(T)$")
    ax.set_title("Raw distributions across L")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    ax = axes[1]
    alpha = duration_fit.tail_exponent
    z = duration_fit.cutoff_exponent
    for L, color in zip(Ls, colors):
        x_s, y_s = collapsed_distribution(
            results[L]["durations_ss"], float(L), alpha, z, n_bins=30
        )
        ax.loglog(x_s, y_s, "o-", color=color, alpha=0.7, label=f"L={L}")
    ax.set_xlabel(r"$T / L^z$")
    ax.set_ylabel(r"$T^{\alpha} P(T)$")
    ax.set_title(rf"Data collapse: $\alpha={alpha:.3f}$, $z={z:.3f}$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    p = outdir / f"oslo_1d_duration_fss_{stamp}.png"
    fig.savefig(p, dpi=150)
    print(f"  wrote {p}")
    plt.close(fig)

    # --- Durations: moment scaling ---
    fig, ax = plt.subplots(figsize=(6, 5))
    L_arr = duration_fit.L_values
    for i, k in enumerate(duration_fit.k_values):
        ax.loglog(L_arr, duration_fit.moments[i], "o-", label=rf"$k={int(k)}$")
    ax.set_xlabel("$L$")
    ax.set_ylabel(r"$\langle T^k \rangle$")
    ax.set_title("Moment scaling for avalanche durations")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    p = outdir / f"oslo_1d_duration_moments_{stamp}.png"
    fig.savefig(p, dpi=150)
    print(f"  wrote {p}")
    plt.close(fig)


def main() -> None:
    print("Simulating 1D Oslo across multiple L:")
    t_start = time.time()
    results = simulate_all()
    sim_elapsed = time.time() - t_start
    print(f"Total simulation time: {sim_elapsed:.1f}s")

    size_fit, duration_fit = fit_and_summarize(results)

    outdir = REPO_ROOT / "data" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Save raw steady-state data so the analysis can be rerun later.
    save_kwargs: dict = {}
    for L, r in results.items():
        save_kwargs[f"sizes_L{L}"] = r["sizes_ss"]
        save_kwargs[f"durations_L{L}"] = r["durations_ss"]
    npz_path = outdir / f"oslo_1d_fss_data_{stamp}.npz"
    np.savez_compressed(npz_path, **save_kwargs)
    print(f"\nSaved raw steady-state samples: {npz_path}")

    print("\nWriting plots:")
    make_plots(results, size_fit, duration_fit, outdir, stamp)


if __name__ == "__main__":
    main()
