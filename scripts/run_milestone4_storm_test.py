"""Eye-of-the-storm test: avalanche distribution shape and divergence.

Tests three predictions of the eye-of-the-storm reading:
  (1) Avalanche distribution shape (tau_s, s_max from a power-law-with-
      cutoff fit) is approximately stationary in p ~ 0.18-0.80 - the
      'eye' where the universe lives quietly between Big Bang
      (percolation) and inversion event.
  (2) The shape departs sharply from stationary at some p > 0.80 -
      'leaving the eye'.
  (3) The divergence of s_max (or other moments) toward some critical
      value p* < 1.0 - the model's predicted inversion point.

All three measurements use the per-drop avalanche size data already
saved in the most recent xi npz file.

Run:
    .venv/bin/python scripts/run_milestone4_storm_test.py [npz]
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))


def latest_xi_npz() -> Path:
    outdir = REPO_ROOT / "data" / "outputs"
    candidates = sorted(outdir.glob("manna_3d_xi_data_*.npz"))
    if not candidates:
        raise FileNotFoundError("No manna_3d_xi_data_*.npz in data/outputs/")
    return candidates[-1]


def log_binned_pdf(sizes: np.ndarray, n_bins: int = 30
                   ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Log-binned PDF of cluster sizes, dropping size-0."""
    sizes = sizes[sizes > 0]
    if sizes.size == 0:
        return np.array([]), np.array([]), np.array([])
    s_max = float(sizes.max())
    edges = np.logspace(0, np.log10(s_max) + 1e-9, n_bins + 1)
    counts, _ = np.histogram(sizes, bins=edges)
    widths = np.diff(edges)
    centers = np.sqrt(edges[:-1] * edges[1:])
    pdf = counts / (sizes.size * widths)
    keep = counts > 0
    return centers[keep], pdf[keep], counts[keep]


def fit_powerlaw_cutoff(centers, pdf, counts, min_s=2.0
                        ) -> tuple[float, float, float] | None:
    """Fit P(s) = A s^(-tau) exp(-s/s_c)."""
    sel = (centers >= min_s) & (pdf > 0) & (counts >= 3)
    if sel.sum() < 4:
        return None

    def model(s, A, tau, s_c):
        return A * s ** (-tau) * np.exp(-s / s_c)

    # Initial guesses
    A0 = pdf[sel][0] * centers[sel][0] ** 1.5
    p0 = [A0, 1.5, centers[sel].max()]
    try:
        popt, _ = curve_fit(model, centers[sel], pdf[sel], p0=p0,
                            sigma=np.maximum(pdf[sel], 1e-12) / np.sqrt(counts[sel]),
                            maxfev=10_000,
                            bounds=([1e-30, 0.0, 1.0], [np.inf, 5.0, 1e10]))
    except Exception:
        return None
    A, tau, s_c = float(popt[0]), float(popt[1]), float(popt[2])
    if s_c < 1 or tau < 0 or tau > 5:
        return None
    return tau, s_c, A


def main() -> None:
    if len(sys.argv) > 1:
        path = Path(sys.argv[1])
    else:
        path = latest_xi_npz()
    print(f"Loading {path}")
    data = np.load(path)
    L = int(data["L"])
    p_targets = sorted([float(p) for p in data["p_targets"]])
    seeds = sorted({int(k.split("_s")[1])
                    for k in data.files if k.startswith("sizes_s")})
    print(f"L={L}, seeds={seeds}, p_targets={p_targets}")
    print()

    # Pool avalanche sizes across seeds within each p-window
    windows: list[dict] = []
    for i, p_hi in enumerate(p_targets):
        p_lo = p_targets[i - 1] if i > 0 else 0.0
        all_sizes = []
        for s in seeds:
            sizes = data[f"sizes_s{s}"]
            drop_hi = int(data[f"drop_s{s}_p{int(round(p_hi * 100)):02d}"])
            drop_lo = (int(data[f"drop_s{s}_p{int(round(p_lo * 100)):02d}"])
                       if p_lo > 0 else 0)
            all_sizes.append(sizes[drop_lo:drop_hi])
        pooled = np.concatenate(all_sizes)
        nonzero = pooled[pooled > 0]
        windows.append({
            "p_lo": p_lo, "p_hi": p_hi,
            "n_drops": int(pooled.size),
            "n_avalanches": int(nonzero.size),
            "mean_size": float(nonzero.mean()) if nonzero.size else 0.0,
            "max_size": int(nonzero.max()) if nonzero.size else 0,
            "sizes": nonzero,
        })

    # Fit power-law-with-cutoff to each window
    print(f"{'window':>14}  {'n_aval':>10}  {'mean_s':>9}  {'max_s':>9}  "
          f"{'tau_s':>7}  {'s_cutoff':>10}")
    fits = []
    for w in windows:
        if w["sizes"].size < 200:
            print(f"{w['p_lo']:.2f}->{w['p_hi']:.2f}    "
                  f"{w['n_avalanches']:>10d}  {w['mean_size']:>9.2f}  "
                  f"{w['max_size']:>9d}     too few avalanches")
            continue
        centers, pdf, counts = log_binned_pdf(w["sizes"], n_bins=22)
        fit = fit_powerlaw_cutoff(centers, pdf, counts)
        if fit is None:
            tau = s_c = float("nan")
        else:
            tau, s_c, _A = fit
        fits.append({**w, "tau": tau, "s_c": s_c})
        print(f"{w['p_lo']:.2f}->{w['p_hi']:.2f}    "
              f"{w['n_avalanches']:>10d}  {w['mean_size']:>9.2f}  "
              f"{w['max_size']:>9d}  {tau:>7.3f}  {s_c:>10.1f}")

    # Question (1): is tau_s stationary in the eye (p ~ 0.18-0.80)?
    print()
    eye_fits = [f for f in fits
                if not np.isnan(f["tau"]) and 0.18 <= f["p_hi"] <= 0.80]
    if len(eye_fits) >= 3:
        taus_eye = np.array([f["tau"] for f in eye_fits])
        print(f"Eye-region (p in [0.18, 0.80]) tau_s values: "
              f"{[f'{t:.3f}' for t in taus_eye]}")
        print(f"  mean = {taus_eye.mean():.3f}, std = {taus_eye.std():.3f}, "
              f"range = [{taus_eye.min():.3f}, {taus_eye.max():.3f}]")
        if taus_eye.std() < 0.15 and (taus_eye.max() - taus_eye.min()) < 0.4:
            print("  -> tau_s is approximately stationary in the eye")
        else:
            print("  -> tau_s varies through the eye - shape NOT stationary")

    # Question (2): where does the shape change?
    print()
    if len(fits) >= 3:
        print("Shape change indicator (s_cutoff growth rate):")
        for f in fits:
            if np.isnan(f["s_c"]):
                continue
            print(f"  p<={f['p_hi']:.2f}  s_c = {f['s_c']:>9.1f}")

    # Question (3): extrapolate s_max divergence to find inversion point p*
    print()
    print("Extrapolating max(s) divergence ~ (p* - p)^(-alpha):")
    p_arr = np.array([w["p_hi"] for w in windows if w["max_size"] > 0])
    s_max_arr = np.array([w["max_size"] for w in windows if w["max_size"] > 0])
    # Fit log(max_s) = -alpha * log(p* - p) + C, looking for p* in (max(p), 1.5]
    def divergence_model(p, p_star, alpha, C):
        return C - alpha * np.log10(np.maximum(p_star - p, 1e-9))
    try:
        popt, pcov = curve_fit(
            divergence_model, p_arr, np.log10(s_max_arr),
            p0=[1.0, 2.0, 0.0],
            bounds=([p_arr.max() + 0.001, 0.1, -10.0],
                    [2.0, 10.0, 20.0]),
            maxfev=20_000,
        )
        p_star, alpha, C = popt
        perr = np.sqrt(np.diag(pcov))
        print(f"  p* (inversion point) = {p_star:.3f} +/- {perr[0]:.3f}")
        print(f"  alpha (divergence exponent) = {alpha:.3f} +/- {perr[1]:.3f}")
        # Residuals
        pred = divergence_model(p_arr, *popt)
        residuals = np.log10(s_max_arr) - pred
        rms = float(np.sqrt(np.mean(residuals ** 2)))
        print(f"  RMS residual in log10(max_s) = {rms:.3f}")
        if p_star < 1.0:
            print(f"  -> implies inversion event at p={p_star:.3f}")
        else:
            print(f"  -> p* > 1; saturation rather than finite-time divergence")
    except Exception as e:
        print(f"  divergence fit failed: {e}")
        p_star, alpha, C = None, None, None

    # --- plot ---
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))

    # Panel 1: tau_s vs p, with eye region shaded
    tau_p = []
    tau_v = []
    sc_p = []
    sc_v = []
    for f in fits:
        if not np.isnan(f["tau"]):
            tau_p.append(0.5 * (f["p_lo"] + f["p_hi"]))
            tau_v.append(f["tau"])
        if not np.isnan(f["s_c"]):
            sc_p.append(0.5 * (f["p_lo"] + f["p_hi"]))
            sc_v.append(f["s_c"])

    ax = axes[0]
    ax.plot(tau_p, tau_v, "o-", color="C0", label=r"$\tau_s$")
    ax.axvspan(0.18, 0.80, color="C2", alpha=0.15, label="proposed eye")
    ax.axhline(2.189, color="0.5", linestyle=":", label="3D percolation 2.189")
    ax.axhline(2.0, color="0.6", linestyle=":", label="Press-Schechter 2.0")
    ax.set_xlabel("$p$ (window midpoint)")
    ax.set_ylabel(r"$\tau_s$ (avalanche-size slope)")
    ax.set_title("Avalanche size distribution slope vs p")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    # Panel 2: s_cutoff vs p, log-scale
    ax = axes[1]
    ax.semilogy(sc_p, sc_v, "o-", color="C3")
    ax.axvspan(0.18, 0.80, color="C2", alpha=0.15)
    ax.set_xlabel("$p$ (window midpoint)")
    ax.set_ylabel(r"$s_c$ (cutoff scale)")
    ax.set_title("Cutoff scale (log) vs p — divergence?")
    ax.grid(True, which="both", alpha=0.3)

    # Panel 3: max_s vs p with divergence fit
    ax = axes[2]
    ax.semilogy(p_arr, s_max_arr, "o-", color="C4", label="data max(s)")
    if p_star is not None:
        p_plot = np.linspace(0.05, min(p_star - 0.01, 0.99), 200)
        s_pred = 10 ** divergence_model(p_plot, p_star, alpha, C)
        ax.semilogy(p_plot, s_pred, "r--",
                    label=fr"fit: $p^* = {p_star:.3f}, \alpha = {alpha:.2f}$")
        ax.axvline(p_star, color="r", linestyle=":", alpha=0.5,
                   label=fr"$p^* = {p_star:.3f}$")
    ax.axvspan(0.18, 0.80, color="C2", alpha=0.15, label="proposed eye")
    ax.set_xlabel("$p$")
    ax.set_ylabel("max avalanche size in window")
    ax.set_title("max(s) divergence")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    fig.tight_layout()
    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_storm_test_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
