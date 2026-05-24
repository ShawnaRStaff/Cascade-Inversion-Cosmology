"""M6 sub-track: refine the asymptote extrapolation with multiple functional forms.

The linear extrapolation A(L) = A_inf + c/L gave A_inf = 0.50 with large
finite-size coefficient. With only 4 L values and accelerating growth,
the true asymptote is uncertain.

This script tries multiple functional forms:
  (a) Linear in 1/L:    A_inf + c/L                  (default earlier)
  (b) Power-law in 1/L: A_inf + c/L^p                (fit p)
  (c) Logarithmic:      A_inf - c*ln(L_ref/L) / something
  (d) Exponential:      A_inf * (1 - exp(-L/L0))     (saturating)
  (e) Stretched exp:    A_inf * (1 - exp(-(L/L0)^q))

Reports the range of plausible asymptote values across all fits to
honestly characterize the uncertainty.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

REPO_ROOT = Path(__file__).resolve().parents[1]
SWEEP_DIR = REPO_ROOT / "data" / "outputs" / "fss_sweep_20260521_031056"


def load_data():
    """Get (L, max_area_pct) pairs across all runs, averaged per L."""
    by_L: dict[int, list[float]] = {}
    for f in sorted(SWEEP_DIR.glob("L*_*_final.npz")):
        d = np.load(f, allow_pickle=True)
        L = int(d["L"])
        if "unique_sizes" not in d.files:
            continue
        unique = np.asarray(d["unique_sizes"])
        max_area = int(unique.max())
        by_L.setdefault(L, []).append(max_area / (L ** 3))
    Ls = sorted(by_L.keys())
    means = [float(np.mean(by_L[L])) for L in Ls]
    stds = [float(np.std(by_L[L])) for L in Ls]
    return np.array(Ls), np.array(means), np.array(stds)


def fit_linear_1_over_L(Ls, ys):
    """A(L) = A_inf + c/L"""
    def f(L, A_inf, c):
        return A_inf + c / L
    popt, _ = curve_fit(f, Ls, ys)
    return popt[0], {"c": popt[1]}


def fit_power_law_1_over_L(Ls, ys):
    """A(L) = A_inf - c/L^p"""
    def f(L, A_inf, c, p):
        return A_inf - c / L ** p
    try:
        popt, _ = curve_fit(f, Ls, ys, p0=[0.6, 5.0, 1.0], maxfev=5000)
        return popt[0], {"c": popt[1], "p": popt[2]}
    except Exception as e:
        return None, {"error": str(e)}


def fit_exponential(Ls, ys):
    """A(L) = A_inf * (1 - exp(-L/L0))"""
    def f(L, A_inf, L0):
        return A_inf * (1 - np.exp(-L / L0))
    try:
        popt, _ = curve_fit(f, Ls, ys, p0=[0.6, 100.0], maxfev=5000)
        return popt[0], {"L0": popt[1]}
    except Exception as e:
        return None, {"error": str(e)}


def fit_stretched_exp(Ls, ys):
    """A(L) = A_inf * (1 - exp(-(L/L0)^q))"""
    def f(L, A_inf, L0, q):
        return A_inf * (1 - np.exp(-(L / L0) ** q))
    try:
        popt, _ = curve_fit(f, Ls, ys, p0=[0.7, 80.0, 0.7], maxfev=5000)
        return popt[0], {"L0": popt[1], "q": popt[2]}
    except Exception as e:
        return None, {"error": str(e)}


def fit_logarithmic(Ls, ys):
    """A(L) = A_inf - c/ln(L)"""
    def f(L, A_inf, c):
        return A_inf - c / np.log(L)
    try:
        popt, _ = curve_fit(f, Ls, ys, p0=[0.7, 1.0], maxfev=5000)
        return popt[0], {"c": popt[1]}
    except Exception as e:
        return None, {"error": str(e)}


def main():
    Ls, means, stds = load_data()
    out_dir = REPO_ROOT / "data" / "outputs" / f"asymptote_refinement_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== Data points ===")
    for L, m, s in zip(Ls, means, stds):
        print(f"  L={L:>3}: max_area/L^3 = {m:.4f} +/- {s:.4f}")

    fits = {}
    print(f"\n=== Asymptote fits (different functional forms) ===")

    for name, fitter in [
        ("linear_1/L", fit_linear_1_over_L),
        ("power_law_1/L^p", fit_power_law_1_over_L),
        ("exponential_saturation", fit_exponential),
        ("stretched_exp", fit_stretched_exp),
        ("logarithmic", fit_logarithmic),
    ]:
        result = fitter(Ls, means)
        if isinstance(result, tuple) and result[0] is not None:
            A_inf, params = result
            fits[name] = {"A_inf": float(A_inf), "params": {k: float(v) for k, v in params.items()}}
            print(f"  {name:>22}: A_inf = {A_inf:.4f}  params = {fits[name]['params']}")
        else:
            fits[name] = {"error": "fit failed"}
            print(f"  {name:>22}: FIT FAILED")

    # Range across fits
    valid_A_infs = [f["A_inf"] for f in fits.values() if "A_inf" in f]
    if valid_A_infs:
        print(f"\n=== Asymptote range across fits ===")
        print(f"  min A_inf = {min(valid_A_infs):.3f}")
        print(f"  max A_inf = {max(valid_A_infs):.3f}")
        print(f"  median    = {float(np.median(valid_A_infs)):.3f}")
        print(f"  mean      = {float(np.mean(valid_A_infs)):.3f}")
        print()
        print(f"  Carrier fraction reference: 0.620")
        print(f"  Full lattice reference:    1.000")

    # Plot all fits
    fig, ax = plt.subplots(figsize=(10, 7))
    L_plot = np.linspace(Ls.min(), max(Ls.max() * 4, 1024), 200)
    ax.errorbar(Ls, means, yerr=stds, fmt="o", markersize=10, capsize=5,
                color="black", label="data")
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
    fit_funcs = {
        "linear_1/L": lambda L, p: p["A_inf"] + p["params"]["c"] / L,
        "power_law_1/L^p": lambda L, p: p["A_inf"] - p["params"]["c"] / L ** p["params"]["p"],
        "exponential_saturation": lambda L, p: p["A_inf"] * (1 - np.exp(-L / p["params"]["L0"])),
        "stretched_exp": lambda L, p: p["A_inf"] * (1 - np.exp(-(L / p["params"]["L0"]) ** p["params"]["q"])),
        "logarithmic": lambda L, p: p["A_inf"] - p["params"]["c"] / np.log(L),
    }
    for (name, fit_info), color in zip(fits.items(), colors):
        if "A_inf" not in fit_info:
            continue
        y_plot = fit_funcs[name](L_plot, fit_info)
        ax.plot(L_plot, y_plot, "-", color=color, alpha=0.7,
                label=f"{name}: A_inf={fit_info['A_inf']:.3f}")
        # Mark asymptote
        ax.axhline(fit_info["A_inf"], color=color, linestyle=":", alpha=0.4)

    ax.axhline(0.62, color="g", linestyle="--", alpha=0.5, label="carrier fraction (0.62)")
    ax.axhline(1.0, color="k", linestyle="--", alpha=0.3, label="full lattice (1.0)")
    ax.set_xlabel("L")
    ax.set_ylabel(r"max area / $L^3$")
    ax.set_title("Asymptote extrapolation: multiple functional forms")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(True, alpha=0.3)
    ax.set_xscale("log")

    fig.tight_layout()
    fig.savefig(out_dir / "asymptote_fits.png", dpi=150)
    print(f"\nPlot: {out_dir / 'asymptote_fits.png'}")

    with open(out_dir / "results.json", "w") as f:
        json.dump({
            "data": {"Ls": Ls.tolist(), "means": means.tolist(), "stds": stds.tolist()},
            "fits": fits,
            "A_inf_range": {"min": min(valid_A_infs), "max": max(valid_A_infs),
                            "median": float(np.median(valid_A_infs))} if valid_A_infs else None,
        }, f, indent=2)
    print(f"JSON: {out_dir / 'results.json'}")


if __name__ == "__main__":
    main()
