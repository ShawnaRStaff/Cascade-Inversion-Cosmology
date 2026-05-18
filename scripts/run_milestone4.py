"""Milestone 4: two-point correlation xi(r) of the ever-toppled set
across multiple snapshots in time, at L=128, multi-seed.

For each seed, drive the 3D Manna lattice from empty and save the
ever_toppled mask at each of several target fractured fractions p:

    p = 0.10  (pre-spanning; sparse fractured set)
    p = 0.18  (around p_c; first global connectivity)
    p = 0.25  (post-spanning; cluster densifying)
    p = 0.35  (late evolution; substantial occupation)

For each saved mask, compute xi(r) via the FFT estimator. Average xi(r)
across seeds at each p; report mean and seed-to-seed std.

Compare to the canonical galaxy two-point correlation reference,
xi(r) ~ (r/r_0)^-gamma with r_0 ~ 5 h^-1 Mpc and gamma ~ 1.8 on small
scales (Peebles 1980; Davis & Peebles 1983; Zehavi et al. 2011).
Fit a power law to the Manna xi(r) over the cosmologically meaningful
range (r small enough to avoid the boundary artifact at r ~ L/2 from
the periodic FFT on an open-boundary simulation; r large enough to
escape lattice-spacing artifacts).

Outputs:
- manna_3d_xi_data_<stamp>.npz : per-(p, seed) xi(r) curves and masks
- manna_3d_xi_curves_<stamp>.png : xi(r) overlay across p targets
- manna_3d_xi_vs_galaxy_<stamp>.png : Manna xi(r) vs galaxy power law

Run:
    .venv/bin/python -u scripts/run_milestone4.py
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

from void_cascade.correlation import (  # noqa: E402
    power_law_galaxy_xi,
    two_point_correlation_3d,
)
from void_cascade.sandpile_3d import run_with_ever_toppled  # noqa: E402


# Run plan. L=96 is set by M3's measured drops-to-span scaling: at L=96
# we know p=0.18 takes ~355k drops, so p=0.35 takes roughly ~760k. At
# L=128 the same targets would take ~22h of compute; the trade-off was
# reviewed and L=96 gives the same ξ(r) shape with manageable wall time.
# n_drops_max is set well above the empirical p=0.35 drops estimate to
# give margin without wasting compute when targets are hit early.
L = 96
SEEDS = [1000, 1001, 1002]
P_TARGETS = [0.10, 0.18, 0.25, 0.35, 0.50, 0.55, 0.60, 0.65, 0.70, 0.80, 0.85, 0.90]
N_DROPS_MAX = 2_000_000
CHECK_EVERY = 100   # check fractured-fraction every 100 drops


def run_one_seed(seed: int) -> dict:
    """Drive a single seed until all target p values are reached.

    Returns dict with: seed, snapshots (dict {p: (mask, drop, p_actual)}),
    wall_seconds, drops_executed.
    """
    targets_sorted = sorted(P_TARGETS)
    next_idx = [0]   # mutable so the closure can advance
    snapshots: dict[float, tuple[np.ndarray, int, float]] = {}

    def cb(t, ever, sizes, durations):
        p_current = float(ever.mean())
        # Drain as many targets as the current p satisfies
        while next_idx[0] < len(targets_sorted) and p_current >= targets_sorted[next_idx[0]]:
            target = targets_sorted[next_idx[0]]
            snapshots[target] = (ever.copy(), t + 1, p_current)
            print(f"    [seed={seed}] hit p={target:.3f} at drop {t+1:6d} "
                  f"(actual p={p_current:.4f})")
            next_idx[0] += 1
        if next_idx[0] >= len(targets_sorted):
            return True
        return False

    t0 = time.time()
    state, sizes, durations, _ever = run_with_ever_toppled(
        L=L, n_drops=N_DROPS_MAX, seed=seed,
        check_every=CHECK_EVERY, percolation_callback=cb,
    )
    wall = time.time() - t0
    print(f"  [seed={seed}] DONE in {wall:.0f}s, drops executed = {sizes.size}, "
          f"snapshots = {len(snapshots)} / {len(P_TARGETS)}")
    return {
        "seed": seed,
        "snapshots": snapshots,
        "wall_seconds": wall,
        "drops_executed": int(sizes.size),
        "sizes": sizes,
        "durations": durations,
    }


def compute_xi_per_snapshot(seed_results: list[dict], n_bins: int) -> dict:
    """Compute xi(r) for every (seed, p) snapshot.

    Returns dict {p: {seed: (r, xi, n_pairs)}}.
    """
    by_p: dict[float, dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]]] = {}
    for sr in seed_results:
        seed = sr["seed"]
        for p, (mask, _drop, _p_actual) in sr["snapshots"].items():
            r, xi, np_arr = two_point_correlation_3d(
                mask, r_max=L / 2.0, n_bins=n_bins, log_bins=True, min_r=1.0
            )
            by_p.setdefault(p, {})[seed] = (r, xi, np_arr)
    return by_p


def aggregate_xi(by_p: dict) -> dict:
    """Average xi(r) across seeds at each p.

    Returns dict {p: {r_centers, xi_mean, xi_std, n_pairs_avg, n_seeds}}.
    """
    agg = {}
    for p, per_seed in by_p.items():
        seeds = sorted(per_seed.keys())
        # All seeds use the same bins (same L, same params), so r_centers
        # should match; take the first one's r.
        r0, xi0, np0 = per_seed[seeds[0]]
        xi_matrix = np.stack([per_seed[s][1] for s in seeds], axis=0)
        np_matrix = np.stack([per_seed[s][2] for s in seeds], axis=0)
        valid = ~np.isnan(xi_matrix)
        # Mask out per-(seed, bin) NaNs; average over seeds that contributed
        with np.errstate(invalid="ignore"):
            xi_mean = np.where(valid.any(axis=0),
                               np.nanmean(xi_matrix, axis=0),
                               np.nan)
            xi_std = np.where(valid.sum(axis=0) > 1,
                              np.nanstd(xi_matrix, axis=0, ddof=1),
                              0.0)
        n_pairs_avg = np_matrix.mean(axis=0)
        agg[p] = {
            "r": r0,
            "xi_mean": xi_mean,
            "xi_std": xi_std,
            "n_pairs_avg": n_pairs_avg,
            "n_seeds": len(seeds),
        }
    return agg


def fit_power_law(r: np.ndarray, xi: np.ndarray,
                  r_min: float, r_max: float
                  ) -> tuple[float, float] | None:
    """Fit xi(r) = (r/r0)^(-gamma) on [r_min, r_max]. Returns (gamma, r0)
    or None if there are too few valid points or xi is non-positive."""
    sel = (r >= r_min) & (r <= r_max) & np.isfinite(xi) & (xi > 0)
    if sel.sum() < 3:
        return None
    log_r = np.log10(r[sel])
    log_xi = np.log10(xi[sel])
    coeffs = np.polyfit(log_r, log_xi, 1)
    gamma = -float(coeffs[0])
    log_A = float(coeffs[1])
    if gamma <= 0:
        return None
    r0 = 10.0 ** (log_A / gamma)
    return gamma, r0


def make_plots(agg: dict, outdir: Path, stamp: str) -> None:
    p_values = sorted(agg.keys())
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(p_values)))

    # --- xi(r) overlay ---
    fig, ax = plt.subplots(figsize=(8, 6))
    for p, color in zip(p_values, colors):
        a = agg[p]
        r = a["r"]
        xi = a["xi_mean"]
        std = a["xi_std"]
        valid = np.isfinite(xi) & (xi > 0) & (a["n_pairs_avg"] > 50)
        ax.loglog(r[valid], xi[valid], "o-", color=color, alpha=0.8,
                  label=f"p = {p:.2f}  (n_seeds={a['n_seeds']})")
        upper = xi[valid] + std[valid]
        lower = xi[valid] - std[valid]
        lower = np.maximum(lower, 1e-6)
        ax.fill_between(r[valid], lower, upper, color=color, alpha=0.15)
    # Galaxy power-law reference: gamma=1.8, r0 free for shape comparison.
    r_ref = np.logspace(0.0, np.log10(L / 2.0), 80)
    for r0_ref, ls in [(5.0, ":"), (10.0, "--"), (20.0, "-.")]:
        ax.loglog(r_ref, power_law_galaxy_xi(r_ref, r0=r0_ref, gamma=1.8),
                  ls, color="0.4", alpha=0.6,
                  label=fr"galaxy ref $\gamma=1.8$, $r_0={r0_ref:.0f}$")
    ax.set_xlabel("separation $r$  (lattice units)")
    ax.set_ylabel(r"$\xi(r)$")
    ax.set_title(f"3D Manna two-point correlation across snapshots (L={L})")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    path = outdir / f"manna_3d_xi_curves_{stamp}.png"
    fig.savefig(path, dpi=150)
    print(f"  wrote {path}")
    plt.close(fig)

    # --- power-law fits per snapshot ---
    fig, ax = plt.subplots(figsize=(8, 6))
    print()
    print("Power-law fits xi(r) = (r/r0)^-gamma on r in [2, L/4]:")
    print(f"{'p':>6}  {'gamma':>8}  {'r0_lat':>8}  {'r_min':>6}  {'r_max':>6}  {'n_fit':>5}")
    fit_results = {}
    for p, color in zip(p_values, colors):
        a = agg[p]
        r = a["r"]
        xi = a["xi_mean"]
        # Fit range: skip the very smallest r (shot noise) and avoid the
        # boundary-artifact region above ~L/4.
        r_min_fit = 2.0
        r_max_fit = float(L) / 4.0
        fit = fit_power_law(r, xi, r_min_fit, r_max_fit)
        sel_for_count = (r >= r_min_fit) & (r <= r_max_fit) & np.isfinite(xi) & (xi > 0)
        n_fit = int(sel_for_count.sum())
        if fit is None:
            print(f"{p:>6.2f}  {'--':>8}  {'--':>8}  {r_min_fit:>6.1f}  "
                  f"{r_max_fit:>6.1f}  {n_fit:>5}  (fit failed)")
            continue
        gamma, r0 = fit
        fit_results[p] = (gamma, r0)
        print(f"{p:>6.2f}  {gamma:>8.3f}  {r0:>8.3f}  {r_min_fit:>6.1f}  "
              f"{r_max_fit:>6.1f}  {n_fit:>5}")

        valid = sel_for_count
        ax.loglog(r[valid], xi[valid], "o", color=color,
                  label=fr"p={p:.2f}  $\gamma={gamma:.2f}, r_0={r0:.2f}$")
        r_line = np.linspace(r_min_fit, r_max_fit, 50)
        ax.loglog(r_line, (r_line / r0) ** (-gamma), "--", color=color, alpha=0.6)

    ax.loglog(r_ref, power_law_galaxy_xi(r_ref, r0=5.0, gamma=1.8),
              ":", color="0.4", alpha=0.7,
              label=r"galaxy ref $\gamma=1.8, r_0=5$")
    ax.set_xlabel("separation $r$  (lattice units)")
    ax.set_ylabel(r"$\xi(r)$")
    ax.set_title(f"Power-law fits in cosmologically-meaningful range (L={L})")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = outdir / f"manna_3d_xi_vs_galaxy_{stamp}.png"
    fig.savefig(path, dpi=150)
    print(f"  wrote {path}")
    plt.close(fig)


def save_raw(seed_results: list[dict], agg: dict, outdir: Path, stamp: str) -> Path:
    path = outdir / f"manna_3d_xi_data_{stamp}.npz"
    kw: dict = {"L": np.int64(L), "p_targets": np.asarray(P_TARGETS)}
    for sr in seed_results:
        for p, (mask, drop, p_actual) in sr["snapshots"].items():
            tag = f"mask_s{sr['seed']}_p{int(round(p * 100)):02d}"
            kw[tag] = mask
            kw[f"drop_s{sr['seed']}_p{int(round(p * 100)):02d}"] = np.int64(drop)
        # Full per-drop avalanche size/duration arrays per seed, for the
        # avalanche-distribution-by-p analysis.
        kw[f"sizes_s{sr['seed']}"] = sr["sizes"].astype(np.int32, copy=False)
        kw[f"durations_s{sr['seed']}"] = sr["durations"].astype(np.int32, copy=False)
    for p, a in agg.items():
        tag = f"xi_p{int(round(p * 100)):02d}"
        kw[f"{tag}_r"] = a["r"]
        kw[f"{tag}_mean"] = a["xi_mean"]
        kw[f"{tag}_std"] = a["xi_std"]
        kw[f"{tag}_npairs"] = a["n_pairs_avg"]
    np.savez_compressed(path, **kw)
    return path


def main() -> None:
    print(f"M4 multi-snapshot xi(r): L={L}, seeds={SEEDS}, p_targets={P_TARGETS}")
    print("=" * 72)
    t_start = time.time()
    seed_results: list[dict] = []
    for seed in SEEDS:
        print(f"  seed={seed}")
        seed_results.append(run_one_seed(seed))
    total = time.time() - t_start
    print(f"All seeds done. Total wall time: {total:.0f}s")

    print()
    print("Computing xi(r) for each snapshot...")
    by_p = compute_xi_per_snapshot(seed_results, n_bins=22)
    agg = aggregate_xi(by_p)

    outdir = REPO_ROOT / "data" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    npz = save_raw(seed_results, agg, outdir, stamp)
    print(f"Saved raw data to {npz}")

    print()
    print("Writing plots:")
    make_plots(agg, outdir, stamp)
    print()
    print("M4 done.")


if __name__ == "__main__":
    main()
