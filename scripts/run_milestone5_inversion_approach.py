"""Stage 1 of the inversion-event redesign.

Uses only existing M4 avalanche time-series data (no new simulation).
Computes the precursor signatures of the inversion event:

  1. Energy turnover per drop (cumulative topplings / cumulative drops)
  2. Boundary-loss escape fraction (grains_out / grains_in) over time
  3. Per-window avalanche statistics: mean, std, max, large-event rate
  4. Inter-event time distribution for "large" avalanches (size > L)
  5. Variance growth of avalanche sizes - critical-fluctuation signature

Goal: characterize what the substrate looks like as it approaches the
predicted inversion point p* ~ 0.97. If the dynamics shows critical-
like behavior (diverging variance, accelerating large-event rate,
shrinking inter-event times), the inversion has identifiable precursors
that map to observable cosmological phenomena (cosmic acceleration,
extreme transient events).

Run:
    .venv/bin/python scripts/run_milestone5_inversion_approach.py [npz]
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


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else latest_xi_npz()
    print(f"Loading {path}")
    data = np.load(path)
    L = int(data["L"])
    p_targets = sorted([float(p) for p in data["p_targets"]])
    seeds = sorted({int(k.split("_s")[1])
                    for k in data.files if k.startswith("sizes_s")})
    primary_seed = seeds[0]
    print(f"L={L}, primary seed={primary_seed}, n_seeds={len(seeds)}")
    L3 = L ** 3

    sizes = data[f"sizes_s{primary_seed}"]
    durations = data[f"durations_s{primary_seed}"]
    n_drops = sizes.size

    # Drop indices at each p target
    drops_at_p = {p: int(data[f"drop_s{primary_seed}_p{int(round(p * 100)):02d}"])
                  for p in p_targets}

    # ============================================================
    # (1) Energy turnover and escape fraction over time
    # ============================================================
    # For Manna with z_c=2, each toppling sheds 2 grains. Each grain
    # either stays in the lattice or leaves via the boundary.
    # cumulative topplings = cumsum(sizes)
    # cumulative grains shed = 2 * cumsum(sizes)
    # cumulative drives = arange(n_drops) + 1  (= grains_in)
    # In steady state: grains_in ~ grains_lost
    # So grains_lost_estimate at time t is bounded by 2*sum(sizes[:t])
    # but really equals grains_in - z_sum, which is = t - z_sum.
    # We don't have z_sum directly. But for boundary-driven SOC,
    # grains_lost approaches grains_in asymptotically.
    print()
    print("Energy / activity per p-window:")
    print(f"{'window':>14}  {'drops':>9}  {'<s>':>8}  {'std(s)':>8}  "
          f"{'max(s)':>8}  {'topp/drop':>9}  {'n_large':>8}  "
          f"{'<inter_large>':>14}")

    windows: list[dict] = []
    for i, p_hi in enumerate(p_targets):
        p_lo = p_targets[i - 1] if i > 0 else 0.0
        drop_lo = drops_at_p[p_lo] if p_lo > 0 else 0
        drop_hi = drops_at_p[p_hi]
        # Pool across all seeds for statistical power
        all_sizes_in_window = []
        all_durations_in_window = []
        for s in seeds:
            sz = data[f"sizes_s{s}"][drop_lo:drop_hi]
            du = data[f"durations_s{s}"][drop_lo:drop_hi]
            all_sizes_in_window.append(sz)
            all_durations_in_window.append(du)
        sizes_pool = np.concatenate(all_sizes_in_window)
        dur_pool = np.concatenate(all_durations_in_window)
        n_drops_in = sizes_pool.size

        # All avalanches in this window (including size-0 = no topple)
        mean_s = float(sizes_pool.mean())
        std_s = float(sizes_pool.std())
        max_s = int(sizes_pool.max())
        topp_per_drop = mean_s  # equivalent

        # Large events: size > L (linear dimension as a natural cut)
        is_large = sizes_pool > L
        n_large = int(is_large.sum())
        # Inter-event times for large events: how many drops between
        # consecutive large events?
        large_drop_indices = np.flatnonzero(is_large)
        if large_drop_indices.size >= 2:
            inter = np.diff(large_drop_indices)
            mean_inter = float(inter.mean())
        else:
            mean_inter = float("inf") if n_large < 2 else float("nan")

        windows.append({
            "p_lo": p_lo, "p_hi": p_hi,
            "n_drops": n_drops_in,
            "mean_s": mean_s, "std_s": std_s,
            "max_s": max_s, "topp_per_drop": topp_per_drop,
            "n_large": n_large, "mean_inter_large": mean_inter,
            "rate_large": n_large / max(n_drops_in, 1),
        })
        print(f"{p_lo:.2f}->{p_hi:.2f}    "
              f"{n_drops_in:>9d}  "
              f"{mean_s:>8.2f}  {std_s:>8.2f}  "
              f"{max_s:>8d}  {topp_per_drop:>9.2f}  "
              f"{n_large:>8d}  "
              f"{mean_inter_large_str(mean_inter):>14}")

    # ============================================================
    # (2) Plot the precursor signatures
    # ============================================================
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))

    p_centers = np.array([0.5 * (w["p_lo"] + w["p_hi"]) for w in windows])

    # Panel 1: mean and std of avalanche size vs p
    ax = axes[0, 0]
    means = np.array([w["mean_s"] for w in windows])
    stds = np.array([w["std_s"] for w in windows])
    ax.semilogy(p_centers, means, "o-", color="C0", label=r"$\langle s\rangle$")
    ax.semilogy(p_centers, stds, "s-", color="C1", label=r"$\sigma_s$ (std)")
    ax.axvline(0.65, color="0.5", linestyle="--", alpha=0.5,
               label="$p=0.65$ (galaxy γ match)")
    ax.axvline(0.97, color="r", linestyle="--", alpha=0.5,
               label="$p^*=0.97$ (predicted inversion)")
    ax.set_xlabel("$p$ (window center)")
    ax.set_ylabel("avalanche size")
    ax.set_title("Avalanche mean and std vs p")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    # Panel 2: large-event rate (events with s > L per drop) vs p
    ax = axes[0, 1]
    rates = np.array([w["rate_large"] for w in windows])
    ax.semilogy(p_centers, rates + 1e-6, "o-", color="C2")
    ax.axvline(0.65, color="0.5", linestyle="--", alpha=0.5)
    ax.axvline(0.97, color="r", linestyle="--", alpha=0.5)
    ax.set_xlabel("$p$")
    ax.set_ylabel("fraction of drops with $s > L$")
    ax.set_title(r"Rate of large (s > L = 96) events")
    ax.grid(True, which="both", alpha=0.3)

    # Panel 3: mean inter-event time for large events
    ax = axes[0, 2]
    inter = np.array([w["mean_inter_large"] for w in windows])
    valid = np.isfinite(inter) & (inter > 0)
    ax.semilogy(p_centers[valid], inter[valid], "o-", color="C3")
    ax.axvline(0.65, color="0.5", linestyle="--", alpha=0.5)
    ax.axvline(0.97, color="r", linestyle="--", alpha=0.5)
    ax.set_xlabel("$p$")
    ax.set_ylabel("mean drops between large events")
    ax.set_title("Inter-event time (large events)")
    ax.grid(True, which="both", alpha=0.3)

    # Panel 4: max(s) vs p with divergence fit overlaid
    ax = axes[1, 0]
    maxs = np.array([w["max_s"] for w in windows])
    ax.semilogy(p_centers, maxs, "o-", color="C4")
    ax.axvline(0.97, color="r", linestyle="--", alpha=0.5,
               label="$p^* = 0.97$")
    ax.axhline(L3, color="k", linestyle=":", alpha=0.5,
               label=f"$L^3 = {L3}$ (whole lattice)")
    ax.set_xlabel("$p$")
    ax.set_ylabel("max(s)")
    ax.set_title("Catastrophic event size vs p")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    # Panel 5: cumulative topplings vs drops (slope = topplings per drop)
    ax = axes[1, 1]
    sizes_pri = data[f"sizes_s{primary_seed}"]
    drop_axis = np.arange(1, sizes_pri.size + 1)
    cum_topplings = np.cumsum(sizes_pri.astype(np.int64))
    # Subsample for plotting
    step = max(1, sizes_pri.size // 5000)
    ax.loglog(drop_axis[::step], cum_topplings[::step], "-", color="C5")
    ax.set_xlabel("drops")
    ax.set_ylabel("cumulative topplings")
    ax.set_title("Energy turnover (cumulative cracks)")
    ax.grid(True, which="both", alpha=0.3)

    # Panel 6: std/mean ratio (relative fluctuations) — critical signature
    ax = axes[1, 2]
    rel_var = stds / np.maximum(means, 1e-9)
    ax.plot(p_centers, rel_var, "o-", color="C6")
    ax.axvline(0.65, color="0.5", linestyle="--", alpha=0.5)
    ax.axvline(0.97, color="r", linestyle="--", alpha=0.5)
    ax.set_xlabel("$p$")
    ax.set_ylabel(r"$\sigma_s / \langle s\rangle$  (relative fluctuation)")
    ax.set_title("Relative fluctuations vs p")
    ax.grid(True, alpha=0.3)

    fig.suptitle("Stage 1 — Inversion-approach signatures from M4 data",
                 fontsize=13)
    fig.tight_layout()

    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_inversion_approach_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")

    # ============================================================
    # Quantitative summary
    # ============================================================
    print()
    print("Summary of precursor signatures:")
    print(f"  At p ~ 0.65 (eye middle): <s> = {windows_at(windows, 0.65)['mean_s']:.1f}, "
          f"std/mean = "
          f"{(windows_at(windows, 0.65)['std_s'] / max(windows_at(windows, 0.65)['mean_s'], 1e-9)):.2f}")
    print(f"  At p ~ 0.80 (eye exit): <s> = {windows_at(windows, 0.80)['mean_s']:.1f}, "
          f"std/mean = "
          f"{(windows_at(windows, 0.80)['std_s'] / max(windows_at(windows, 0.80)['mean_s'], 1e-9)):.2f}")
    print(f"  At p ~ 0.90 (storm-building): <s> = {windows_at(windows, 0.90)['mean_s']:.1f}, "
          f"std/mean = "
          f"{(windows_at(windows, 0.90)['std_s'] / max(windows_at(windows, 0.90)['mean_s'], 1e-9)):.2f}")
    print()
    print("Critical-divergence signatures (if any of these grow without bound:")
    print("  the dynamics is approaching a phase transition / inversion)")
    w65 = windows_at(windows, 0.65)
    w90 = windows_at(windows, 0.90)
    rel_65 = w65["std_s"] / max(w65["mean_s"], 1e-9)
    rel_90 = w90["std_s"] / max(w90["mean_s"], 1e-9)
    print(f"  std/mean growth p=0.65 -> p=0.90: {rel_90 / max(rel_65, 1e-9):.2f}x")
    rate_growth = w90["rate_large"] / max(w65["rate_large"], 1e-9)
    print(f"  large-event rate p=0.65 -> p=0.90: {rate_growth:.0f}x")


def mean_inter_large_str(v: float) -> str:
    if not np.isfinite(v):
        return "inf"
    return f"{v:.0f}"


def windows_at(windows, p):
    # find the window whose p_hi == p
    for w in windows:
        if abs(w["p_hi"] - p) < 1e-6:
            return w
    return windows[-1]


if __name__ == "__main__":
    main()
