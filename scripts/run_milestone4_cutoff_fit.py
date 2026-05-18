"""M4 Speculation #1: fit power-law-with-cutoff to xi(r) at every p.

Tests the hypothesis that the apparent gamma we measured is a fitting
artifact, and the underlying correlation function has a fixed
universal alpha plus a p-dependent xi_corr.

Loads the most recent xi(r) data, refits each p with the form
    xi(r) = A * r^(-alpha) * exp(-r / xi_corr)
and reports whether alpha is constant across p.

Decisive question for the design plan:
- If alpha is constant (within fit error) across all p:
    -> Speculation #1 is supported. Skip #2. Do #3.
- If alpha varies systematically with p:
    -> Speculation #1 fails. Need crossover physics (#2) first.

Run:
    .venv/bin/python scripts/run_milestone4_cutoff_fit.py [path/to/xi_data.npz]
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
    fit_power_law_with_cutoff,
    power_law_with_cutoff,
)


def latest_xi_npz() -> Path:
    outdir = REPO_ROOT / "data" / "outputs"
    candidates = sorted(outdir.glob("manna_3d_xi_data_*.npz"))
    if not candidates:
        raise FileNotFoundError("No manna_3d_xi_data_*.npz in data/outputs/")
    return candidates[-1]


def main() -> None:
    if len(sys.argv) > 1:
        path = Path(sys.argv[1])
    else:
        path = latest_xi_npz()
    print(f"Loading {path}")
    data = np.load(path)
    L = int(data["L"])
    p_targets = sorted([float(p) for p in data["p_targets"]])
    print(f"L = {L}, p_targets = {p_targets}")

    print()
    print(f"{'p':>6}  {'alpha':>10}  {'alpha_err':>10}  {'xi_corr':>10}  "
          f"{'xi_corr_err':>12}  {'A':>10}  {'n_fit':>5}")
    fits: dict[float, tuple[float, float, float, np.ndarray]] = {}

    # Fit range: same as the pure-power-law fits used before, so the
    # two analyses are directly comparable.
    r_min_fit = 1.0
    r_max_fit = float(L) / 2.0

    for p in p_targets:
        key_r = f"xi_p{int(round(p * 100)):02d}_r"
        key_mean = f"xi_p{int(round(p * 100)):02d}_mean"
        key_std = f"xi_p{int(round(p * 100)):02d}_std"
        key_npairs = f"xi_p{int(round(p * 100)):02d}_npairs"
        if key_r not in data.files:
            continue
        r = data[key_r]
        xi = data[key_mean]
        sigma = data[key_std]
        n_pairs = data[key_npairs]

        # Use the same selection as before (positive xi, enough pairs)
        valid = np.isfinite(xi) & (xi > 0) & (n_pairs > 50)
        if valid.sum() < 4:
            print(f"{p:>6.2f}   too few valid points")
            continue

        # For curve_fit, weight by stderr where it's finite and > 0
        sig_for_fit = np.where(sigma > 0, sigma, np.nan)

        result = fit_power_law_with_cutoff(
            r[valid], xi[valid],
            sigma=sig_for_fit[valid] if np.isfinite(sig_for_fit[valid]).all() else None,
            r_min=r_min_fit, r_max=r_max_fit,
        )
        if result is None:
            print(f"{p:>6.2f}  fit failed")
            continue
        A, alpha, xi_corr, perr = result
        fits[p] = (A, alpha, xi_corr, perr)
        print(f"{p:>6.2f}  {alpha:>10.3f}  {perr[1]:>10.3f}  "
              f"{xi_corr:>10.3f}  {perr[2]:>12.3f}  {A:>10.3f}  {int(valid.sum()):>5}")

    if not fits:
        print("No fits succeeded.")
        return

    # --- decision: is alpha constant? ---
    print()
    alpha_values = np.array([fits[p][1] for p in fits])
    alpha_errs = np.array([fits[p][3][1] for p in fits])
    weights = 1.0 / np.maximum(alpha_errs, 1e-6) ** 2
    alpha_weighted_mean = float((alpha_values * weights).sum() / weights.sum())
    alpha_weighted_se = float(1.0 / np.sqrt(weights.sum()))
    chi2 = float(((alpha_values - alpha_weighted_mean) ** 2 * weights).sum())
    dof = max(len(alpha_values) - 1, 1)
    print(f"alpha (weighted mean over all p) = {alpha_weighted_mean:.3f} +/- {alpha_weighted_se:.3f}")
    print(f"chi^2 / dof for 'alpha constant' = {chi2 / dof:.2f}  (dof = {dof})")
    if chi2 / dof < 3.0:
        print("  -> alpha is consistent with being CONSTANT across p.")
        print("     Speculation #1 SUPPORTED. Drop #2. Proceed to #3.")
    else:
        print("  -> alpha varies systematically with p.")
        print("     Speculation #1 NOT supported. Run #2 (crossover regimes).")

    # --- xi_corr(p) trajectory ---
    print()
    print("xi_corr(p) trajectory:")
    for p in sorted(fits):
        _A, _alpha, xi_corr, perr = fits[p]
        print(f"  p={p:.2f}  xi_corr = {xi_corr:7.3f} +/- {perr[2]:.3f} lattice units")

    # --- plot ---
    p_values = sorted(fits.keys())
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(p_values)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Panel 1: alpha vs p with weighted-mean band
    alphas = [fits[p][1] for p in p_values]
    alpha_errs = [fits[p][3][1] for p in p_values]
    ax1.errorbar(p_values, alphas, yerr=alpha_errs, fmt="o-", color="C0",
                 capsize=4, label="fitted alpha per p")
    ax1.axhline(alpha_weighted_mean, color="C1", linestyle="-",
                label=fr"weighted mean $\alpha = {alpha_weighted_mean:.3f}$"
                      fr" ($\chi^2/\nu = {chi2/dof:.2f}$)")
    ax1.fill_between([min(p_values), max(p_values)],
                     alpha_weighted_mean - alpha_weighted_se,
                     alpha_weighted_mean + alpha_weighted_se,
                     color="C1", alpha=0.2)
    ax1.set_xlabel("$p$  (fractured fraction)")
    ax1.set_ylabel(r"underlying $\alpha$ from $\xi(r) = A r^{-\alpha} e^{-r/\xi_{\rm corr}}$")
    ax1.set_title("Speculation #1 test: is alpha constant?")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9)

    # Panel 2: xi_corr vs p (log y)
    xi_corrs = [fits[p][2] for p in p_values]
    xi_corr_errs = [fits[p][3][2] for p in p_values]
    ax2.errorbar(p_values, xi_corrs, yerr=xi_corr_errs, fmt="o-", color="C2",
                 capsize=4)
    ax2.set_xlabel("$p$")
    ax2.set_ylabel(r"$\xi_{\rm corr}$ (lattice units)")
    ax2.set_title(r"Correlation length $\xi_{\rm corr}(p)$")
    ax2.set_yscale("log")
    ax2.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_xi_cutoff_fit_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")

    # Panel 3 (separate fig): xi(r) curves with fitted form overlaid
    fig2, ax = plt.subplots(figsize=(8, 6))
    for p, color in zip(p_values, colors):
        key_r = f"xi_p{int(round(p * 100)):02d}_r"
        key_mean = f"xi_p{int(round(p * 100)):02d}_mean"
        key_npairs = f"xi_p{int(round(p * 100)):02d}_npairs"
        r = data[key_r]
        xi = data[key_mean]
        n_pairs = data[key_npairs]
        valid = np.isfinite(xi) & (xi > 0) & (n_pairs > 50)
        ax.loglog(r[valid], xi[valid], "o", color=color, alpha=0.6, markersize=4)
        if p in fits:
            A, alpha, xi_corr, _ = fits[p]
            r_line = np.logspace(np.log10(r[valid].min()),
                                 np.log10(r[valid].max()), 80)
            ax.loglog(r_line, power_law_with_cutoff(r_line, A, alpha, xi_corr),
                      "-", color=color, alpha=0.8,
                      label=fr"p={p:.2f}: $\alpha={alpha:.2f}$, $\xi_c={xi_corr:.2f}$")
    ax.set_xlabel("$r$ (lattice units)")
    ax.set_ylabel(r"$\xi(r)$")
    ax.set_title("Manna ξ(r) with power-law-with-cutoff fits per p")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    fig2.tight_layout()
    out_path2 = outdir / f"manna_3d_xi_cutoff_curves_{stamp}.png"
    fig2.savefig(out_path2, dpi=150)
    print(f"Wrote {out_path2}")


if __name__ == "__main__":
    main()
