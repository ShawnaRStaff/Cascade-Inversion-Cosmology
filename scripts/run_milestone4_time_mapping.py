"""Time-mapping experiment.

The 'linear time' assumption (1 sim drop = 1 cosmic time unit) makes
gamma(t) vary sharply because gamma(p) varies sharply. The
eye-of-the-storm reading suggests cosmic time should be stretched in
the eye (where dp/drops is small) and compressed near the inversion
(where dp/drops is large).

This script tests several candidate time mappings against the eye
hypothesis. For each mapping, we compute gamma vs cosmic-time and ask:
is gamma approximately flat through the eye?

Time mappings tested:
  (A) Linear in p (the morning's 'free time-mapping parameter' fit)
  (B) Linear in drops (one drop = one cosmic-time unit)
  (C) Activity-rate: dt_cosmic ~ 1 / <s> (small avalanches = large dt)
  (D) Inverse cracking rate: dt_cosmic ~ 1 / (dp/d_drops)
  (E) Cracking rate: dt_cosmic ~ dp/d_drops (the OPPOSITE - cosmic time
      proportional to substrate activity)
  (F) Avalanche-volume-weighted: dt_cosmic ~ s (cosmic time proportional
      to total fracture activity)

For each mapping, normalize cosmic time to [0, 1] and replot gamma(t).
A mapping that flattens gamma in the eye supports the eye hypothesis;
a mapping where gamma still varies sharply doesn't.

Run:
    .venv/bin/python scripts/run_milestone4_time_mapping.py [npz]
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


def fit_gamma(r, xi, npairs, r_min=2.0, r_max=24.0):
    sel = (r >= r_min) & (r <= r_max) & np.isfinite(xi) & (xi > 0) & (npairs > 50)
    if sel.sum() < 3:
        return None
    coeffs = np.polyfit(np.log10(r[sel]), np.log10(xi[sel]), 1)
    return -float(coeffs[0])


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else latest_xi_npz()
    print(f"Loading {path}")
    data = np.load(path)
    L = int(data["L"])
    p_targets = sorted([float(p) for p in data["p_targets"]])
    seeds = sorted({int(k.split("_s")[1])
                    for k in data.files if k.startswith("sizes_s")})

    # Build per-window data: drops_at_p, gamma_at_p, mean_avalanche_size_in_window
    p_arr = np.array(p_targets)
    primary_seed = seeds[0]
    drops_at_p = np.array(
        [int(data[f"drop_s{primary_seed}_p{int(round(p * 100)):02d}"])
         for p in p_targets]
    )

    sizes = data[f"sizes_s{primary_seed}"]

    # Avalanche stats per window (between consecutive snapshots)
    mean_s_in_window = []
    sum_s_in_window = []
    for i in range(len(p_targets)):
        lo = drops_at_p[i - 1] if i > 0 else 0
        hi = drops_at_p[i]
        win = sizes[lo:hi]
        win_nz = win[win > 0]
        mean_s_in_window.append(float(win_nz.mean()) if win_nz.size else 0.0)
        sum_s_in_window.append(float(win.sum()))
    mean_s_in_window = np.array(mean_s_in_window)
    sum_s_in_window = np.array(sum_s_in_window)

    # Gamma at each snapshot
    gammas = []
    for p in p_targets:
        pp = int(round(p * 100))
        r = data[f"xi_p{pp:02d}_r"]
        xi = data[f"xi_p{pp:02d}_mean"]
        npairs = data[f"xi_p{pp:02d}_npairs"]
        g = fit_gamma(r, xi, npairs)
        gammas.append(g if g is not None else np.nan)
    gammas = np.array(gammas)

    # dp/d_drops per window
    drops_diffs = np.diff(drops_at_p, prepend=0)
    p_diffs = np.diff(p_arr, prepend=0)
    dp_per_drop_window = p_diffs / np.maximum(drops_diffs, 1)
    drops_per_dp_window = drops_diffs / np.maximum(p_diffs, 1e-9)

    print()
    print(f"{'p':>5}  {'drops':>9}  {'gamma':>6}  {'dp/drop':>10}  {'drops/dp':>10}  {'<s>_window':>11}")
    for i in range(len(p_targets)):
        print(f"{p_arr[i]:>5.2f}  {drops_at_p[i]:>9d}  {gammas[i]:>6.2f}  "
              f"{dp_per_drop_window[i]:>10.2e}  {drops_per_dp_window[i]:>10.0f}  "
              f"{mean_s_in_window[i]:>11.2f}")

    # Now construct cumulative cosmic time under each mapping.
    # For mapping with dt_cosmic ~ f(p), cumulative t = integral of f.
    # We have discrete windows; integrate accordingly.

    def cumulative_time(rate_per_window):
        """Given dt per window, return cumulative time at each snapshot p."""
        increments = np.array(rate_per_window) * drops_diffs
        return np.cumsum(increments)

    # Compute six mappings
    # (A) Linear in p
    t_A = p_arr.copy()
    # (B) Linear in drops
    t_B = drops_at_p.astype(float).copy()
    # (C) Activity rate dt ~ 1/<s>  (per drop). Per-window rate = 1/<s>_window.
    t_C = cumulative_time(1.0 / np.maximum(mean_s_in_window, 1e-9))
    # (D) Inverse cracking rate dt ~ 1/(dp/drops)
    t_D = cumulative_time(drops_per_dp_window)
    # (E) Cracking rate dt ~ dp/drops
    t_E = cumulative_time(dp_per_drop_window)
    # (F) Avalanche volume: dt ~ sum(s) in window per drop = <s>
    t_F = cumulative_time(mean_s_in_window)

    mappings = {
        "A: linear in p": t_A,
        "B: linear in drops": t_B,
        "C: 1/<s> (slow when active)": t_C,
        "D: 1/(dp/drops) (slow when changing)": t_D,
        "E: dp/drops (fast when changing)": t_E,
        "F: <s> (fast when active)": t_F,
    }

    # For each mapping, compute the fraction of cosmic time spent in the eye
    # (defined as p in [0.18, 0.80]) and the gamma range traversed in that
    # fraction.
    print()
    print("How well each mapping concentrates cosmic time in the eye:")
    print(f"{'mapping':>40}  {'eye time frac':>14}  {'gamma range in eye':>20}")
    eye_mask = (p_arr >= 0.18) & (p_arr <= 0.80)

    for name, t in mappings.items():
        t_norm = (t - t.min()) / (t.max() - t.min()) if t.max() > t.min() else t
        if eye_mask.sum() < 2:
            continue
        t_eye = t_norm[eye_mask]
        eye_frac = t_eye[-1] - t_eye[0]
        # Gamma range in eye
        g_eye = gammas[eye_mask]
        g_range = float(np.nanmax(g_eye) - np.nanmin(g_eye))
        print(f"{name:>40}  {eye_frac:>14.3f}  {g_range:>20.3f}")

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    flat = axes.flatten()
    for i, (name, t) in enumerate(mappings.items()):
        ax = flat[i]
        t_norm = (t - t.min()) / max(t.max() - t.min(), 1e-9)
        ax.plot(t_norm, gammas, "o-", color="C0")
        # Eye shading
        eye_t = t_norm[eye_mask]
        if eye_t.size:
            ax.axvspan(eye_t[0], eye_t[-1], color="C2", alpha=0.15,
                       label="eye p∈[0.18, 0.80]")
        # Galaxy reference band
        ax.axhspan(1.66, 1.91, color="C1", alpha=0.2,
                   label=r"observed $\gamma$ range")
        ax.axhline(1.82, color="C1", linestyle=":", alpha=0.7)
        # Annotate p values
        for ti, pi, gi in zip(t_norm, p_arr, gammas):
            if np.isfinite(gi):
                ax.annotate(f"{pi:.2f}", (ti, gi), fontsize=7, alpha=0.5,
                            xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel(f"normalized cosmic time ({name.split(':')[0]})")
        ax.set_ylabel(r"$\gamma$")
        ax.set_title(name, fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="best")

    fig.suptitle(r"$\gamma$ vs cosmic time under six candidate time mappings",
                 fontsize=12)
    fig.tight_layout()
    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_time_mapping_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")

    # Quantitative test of "eye-flatness": std of gamma values whose
    # normalized cosmic time falls in the middle 80% of the eye, for each
    # mapping. If gamma is flat in the eye, the std should be small.
    print()
    print("Eye-flatness test: std of gamma at snapshots inside the eye")
    print(f"{'mapping':>40}  {'std(gamma) in eye':>20}")
    for name, t in mappings.items():
        if eye_mask.sum() < 3:
            continue
        g_eye = gammas[eye_mask]
        g_eye = g_eye[np.isfinite(g_eye)]
        std = float(g_eye.std())
        print(f"{name:>40}  {std:>20.3f}")
    print()
    print("Note: the eye's gamma range (3.62 down to 1.12 in p=0.18->0.80)")
    print("is intrinsic to gamma(p). NO time mapping can reduce this range.")
    print("What a 'good' mapping does is stretch the eye, putting most")
    print("cosmic time at one end of the gamma range so the OBSERVED")
    print("gamma at most cosmic times is near a single value.")


if __name__ == "__main__":
    main()
