"""M6 sub-track: moment-method analysis of avalanche area distribution.

The H&P methodology for extracting tau_a and D_a from finite-size data:

  <a^n> ~ L^(mu_n)   with   mu_n = D_a * (1 + n - tau_a)

For each n, log<a^n> vs log L is linear with slope mu_n. Then mu_n itself
is linear in n:
  mu_n = D_a + D_a * n - D_a * tau_a
       = (1 - tau_a) * D_a + n * D_a

So plotting mu_n vs n gives D_a as slope and (1 - tau_a) * D_a as
intercept. From those: D_a = slope, tau_a = 1 - intercept / D_a.

This is a STRONGER analysis than our earlier log-binned histogram fit
because it uses the full bulk of the distribution, not just the linear
scaling regime. The literature canonical values for 3D Abelian Manna
on a cubic lattice (Huynh & Pruessner 2012):
  D_a = 3.003 +/- 0.014
  tau_a = 1.442 +/- 0.012

If our parallel-update implementation gives consistent values upon
finite-L extrapolation, we are in the canonical universality class.

Inputs: existing FSS sweep data (no new compute)
Outputs: data/outputs/moment_analysis_<stamp>/
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
SWEEP_DIR = REPO_ROOT / "data" / "outputs" / "fss_sweep_20260521_031056"

D_A_LITERATURE = 3.003
D_A_LIT_ERR = 0.014
TAU_A_LITERATURE = 1.442
TAU_A_LIT_ERR = 0.012

MOMENTS = [1, 2, 3, 4, 5]


def find_saturation_start(snapshots: list[dict]) -> int:
    for s in snapshots:
        if s.get("p", 0) >= 0.95:
            return int(s["drop"])
    return -1


def compute_moments(unique_sizes: np.ndarray, n_list: list[int]) -> dict[int, float]:
    """Compute <a^n> for n in n_list, using nonzero events only."""
    nonzero = unique_sizes[unique_sizes > 0].astype(np.float64)
    out = {}
    for n in n_list:
        if len(nonzero) == 0:
            out[n] = float("nan")
        else:
            out[n] = float(np.mean(nonzero ** n))
    return out


def main():
    out_dir = REPO_ROOT / "data" / "outputs" / f"moment_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load and compute moments per run
    by_L: dict[int, list[dict]] = {}
    for f in sorted(SWEEP_DIR.glob("L*_*_final.npz")):
        d = np.load(f, allow_pickle=True)
        if "unique_sizes" not in d.files:
            continue
        L = int(d["L"])
        seed = int(d["seed"])
        unique = np.asarray(d["unique_sizes"])
        snapshots = list(d["snapshots"]) if "snapshots" in d.files else []
        sat_start = find_saturation_start(snapshots)
        if sat_start < 0:
            sat_start = len(unique) // 2
        sat = unique[sat_start:]
        moments = compute_moments(sat, MOMENTS)
        by_L.setdefault(L, []).append({"seed": seed, "moments": moments,
                                       "n_events": int((sat > 0).sum())})

    # Aggregate per L
    Ls = sorted(by_L.keys())
    print(f"=== Moments <a^n> per L (averaged over seeds) ===\n")
    print(f"  {'L':>4}  {'n_seeds':>7}  ", end="")
    for n in MOMENTS:
        print(f"{'<a^' + str(n) + '>':>15}", end="  ")
    print()

    L_moments: dict[int, dict[int, float]] = {}
    for L in Ls:
        runs = by_L[L]
        L_moments[L] = {}
        print(f"  {L:>4}  {len(runs):>7}  ", end="")
        for n in MOMENTS:
            vals = [r["moments"][n] for r in runs]
            mean = float(np.mean(vals))
            L_moments[L][n] = mean
            print(f"{mean:>15.4g}", end="  ")
        print()

    # For each n, fit <a^n> ~ L^mu_n
    print(f"\n=== Per-moment scaling fits ===")
    print(f"  log<a^n> = mu_n * log(L) + const")
    print(f"  {'n':>3}  {'mu_n':>8}  {'sigma_mu':>8}")
    mu_n_values: dict[int, float] = {}
    mu_n_errors: dict[int, float] = {}
    for n in MOMENTS:
        ys = np.array([L_moments[L][n] for L in Ls])
        xs = np.array(Ls, dtype=np.float64)
        log_y = np.log10(ys)
        log_x = np.log10(xs)
        slope, intercept = np.polyfit(log_x, log_y, 1)
        # Residuals -> err
        residuals = log_y - (slope * log_x + intercept)
        sigma = float(np.std(residuals)) if len(residuals) > 2 else 0.0
        # Standard error of slope
        sx2 = float(np.sum((log_x - log_x.mean()) ** 2))
        slope_err = sigma / np.sqrt(sx2) if sx2 > 0 else 0.0
        mu_n_values[n] = float(slope)
        mu_n_errors[n] = slope_err
        print(f"  {n:>3}  {slope:>8.4f}  {slope_err:>8.4f}")

    # Fit mu_n vs n: mu_n = (1 - tau_a) * D_a + n * D_a
    # slope = D_a, intercept = (1 - tau_a) * D_a, so tau_a = 1 - intercept / D_a
    ns_arr = np.array(MOMENTS, dtype=np.float64)
    mu_arr = np.array([mu_n_values[n] for n in MOMENTS])
    mu_err = np.array([mu_n_errors[n] for n in MOMENTS])

    # Weighted fit
    weights = 1.0 / np.maximum(mu_err, 1e-6) ** 2
    p, V = np.polyfit(ns_arr, mu_arr, 1, w=weights, cov=True)
    D_a_fit = float(p[0])
    intercept_fit = float(p[1])
    tau_a_fit = 1.0 - intercept_fit / D_a_fit
    D_a_err = float(np.sqrt(V[0, 0]))
    intercept_err = float(np.sqrt(V[1, 1]))
    # Error propagation for tau_a
    tau_a_err = np.sqrt(
        (intercept_err / D_a_fit) ** 2 +
        (intercept_fit * D_a_err / D_a_fit ** 2) ** 2
    )

    print(f"\n=== Extracted exponents (from mu_n vs n fit) ===")
    print(f"  D_a  = {D_a_fit:.4f} +/- {D_a_err:.4f}")
    print(f"  tau_a = {tau_a_fit:.4f} +/- {tau_a_err:.4f}")
    print(f"\n  Literature (Huynh & Pruessner 2012):")
    print(f"  D_a  = {D_A_LITERATURE} +/- {D_A_LIT_ERR}")
    print(f"  tau_a = {TAU_A_LITERATURE} +/- {TAU_A_LIT_ERR}")
    print()
    da_gap = abs(D_a_fit - D_A_LITERATURE) / max(D_a_err, 1e-6)
    ta_gap = abs(tau_a_fit - TAU_A_LITERATURE) / max(tau_a_err, 1e-6)
    print(f"  D_a vs literature:   {da_gap:.1f} sigma")
    print(f"  tau_a vs literature: {ta_gap:.1f} sigma")

    if da_gap < 2 and ta_gap < 2:
        verdict = "BOTH consistent with literature within 2 sigma. Our implementation reproduces canonical 3D Manna exponents."
    elif da_gap < 2:
        verdict = "D_a consistent with literature; tau_a deviates. Could be finite-L effect."
    elif ta_gap < 2:
        verdict = "tau_a consistent with literature; D_a deviates."
    else:
        verdict = "Both exponents deviate from literature. Investigate."
    print(f"\n  Verdict: {verdict}")

    # Plots
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    for n, color in zip(MOMENTS, ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]):
        ys = np.array([L_moments[L][n] for L in Ls])
        ax.loglog(Ls, ys, "o-", color=color, label=f"n={n}, μ={mu_n_values[n]:.3f}")
    ax.set_xlabel("L")
    ax.set_ylabel(r"$\langle a^n \rangle$")
    ax.set_title("Avalanche area moments vs L")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)

    ax = axes[1]
    ax.errorbar(ns_arr, mu_arr, yerr=mu_err, fmt="o", markersize=10, capsize=5,
                color="black", label="data")
    n_plot = np.linspace(0, max(MOMENTS) + 1, 50)
    ax.plot(n_plot, D_a_fit * n_plot + intercept_fit, "-", color="tab:blue",
            label=f"fit: D_a={D_a_fit:.3f}, τ_a={tau_a_fit:.3f}")
    # Literature line
    intercept_lit = (1 - TAU_A_LITERATURE) * D_A_LITERATURE
    ax.plot(n_plot, D_A_LITERATURE * n_plot + intercept_lit, "--", color="tab:red", alpha=0.6,
            label=f"literature: D_a={D_A_LITERATURE}, τ_a={TAU_A_LITERATURE}")
    ax.set_xlabel("n")
    ax.set_ylabel(r"$\mu_n$  (slope of $\log\langle a^n\rangle$ vs $\log L$)")
    ax.set_title(r"Moment exponents $\mu_n = D_a \cdot (1 + n - \tau_a)$")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_dir / "moment_analysis.png", dpi=150)
    print(f"\nPlot: {out_dir / 'moment_analysis.png'}")

    # JSON
    out = {
        "literature": {"D_a": D_A_LITERATURE, "D_a_err": D_A_LIT_ERR,
                       "tau_a": TAU_A_LITERATURE, "tau_a_err": TAU_A_LIT_ERR},
        "our_fit": {"D_a": D_a_fit, "D_a_err": D_a_err,
                    "tau_a": tau_a_fit, "tau_a_err": tau_a_err},
        "gap_sigma": {"D_a": float(da_gap), "tau_a": float(ta_gap)},
        "verdict": verdict,
        "mu_n": {str(n): {"value": mu_n_values[n], "err": mu_n_errors[n]} for n in MOMENTS},
        "moments_per_L": {str(L): {str(n): L_moments[L][n] for n in MOMENTS} for L in Ls},
    }
    with open(out_dir / "results.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"JSON: {out_dir / 'results.json'}")


if __name__ == "__main__":
    main()
