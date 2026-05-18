"""Local-clock observer experiment.

Hypothesis: time = count of local fracture events, not uniform drop
counter. Observers (galaxies, us) sitting in regions of low local event
rate have slow clocks. The flat observed gamma(z) across cosmic history
should fall out naturally if we use observer-local-time instead of
global p as the time axis.

This script tests the hypothesis directly. For several candidate
observer regions:
  - central fixed cube of various sizes
  - "deep interior" cells (cells whose neighborhood was fully cracked
    early)
compute the cumulative local event count as global p evolves, then plot
the model's gamma(p) curve against this observer-local-time axis. If
the gamma curve flattens out across a long stretch of observer-time
(many observer-local-events with only small change in gamma), the
non-linear-time hypothesis is supported. If gamma still varies sharply
in observer-time, the hypothesis is falsified.

Run:
    .venv/bin/python scripts/run_milestone4_local_clock.py [npz]
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


def fit_gamma(r: np.ndarray, xi: np.ndarray, npairs: np.ndarray,
              r_min: float, r_max: float) -> float | None:
    sel = (r >= r_min) & (r <= r_max) & np.isfinite(xi) & (xi > 0) & (npairs > 50)
    if sel.sum() < 3:
        return None
    coeffs = np.polyfit(np.log10(r[sel]), np.log10(xi[sel]), 1)
    return -float(coeffs[0])


def observer_event_history(masks: dict[float, np.ndarray],
                           observer_mask: np.ndarray
                           ) -> dict[float, int]:
    """For each global p, count cells in observer_mask that have toppled."""
    result = {}
    for p, mask_at_p in masks.items():
        # cells in observer region that are cracked at this p
        result[p] = int((mask_at_p & observer_mask).sum())
    return result


def main() -> None:
    if len(sys.argv) > 1:
        path = Path(sys.argv[1])
    else:
        path = latest_xi_npz()
    print(f"Loading {path}")
    data = np.load(path)
    L = int(data["L"])
    p_targets = sorted([float(p) for p in data["p_targets"]])
    seeds = sorted({int(k.split("_s")[1].split("_p")[0])
                    for k in data.files if k.startswith("mask_")})
    print(f"L={L}, p_targets={p_targets}, seeds={seeds}")

    # Compute global gamma at each p (averaged across seeds via the saved
    # mean xi(r) curves).
    print()
    print("Step 1: gamma(p) (already in the npz aggregation)")
    gamma_p = {}
    for p in p_targets:
        pp = int(round(p * 100))
        r = data[f"xi_p{pp:02d}_r"]
        xi = data[f"xi_p{pp:02d}_mean"]
        npairs = data[f"xi_p{pp:02d}_npairs"]
        g = fit_gamma(r, xi, npairs, 2.0, L / 4.0)
        gamma_p[p] = g
        print(f"  p={p:.2f}  gamma={g:.3f}")

    # Load primary-seed masks
    primary_seed = seeds[0]
    masks = {}
    for p in p_targets:
        pp = int(round(p * 100))
        masks[p] = data[f"mask_s{primary_seed}_p{pp:02d}"]
    print()
    print(f"Loaded {len(masks)} masks for seed {primary_seed}")

    # Define several observer regions:
    # (a) Central 32^3 cube
    # (b) Central 16^3 cube
    # (c) Central 8^3 cube
    # (d) Single central cell + 1-cell-radius neighborhood
    # (e) "Deep early-cracked region" — cells that toppled before p=0.35
    observers = {}
    L_arr = np.arange(L)
    ii, jj, kk = np.meshgrid(L_arr, L_arr, L_arr, indexing="ij")
    cx = cy = cz = L / 2.0

    # Central cubes
    for cube_side in (8, 16, 32, 48):
        in_cube = (
            (np.abs(ii - cx) <= cube_side / 2)
            & (np.abs(jj - cy) <= cube_side / 2)
            & (np.abs(kk - cz) <= cube_side / 2)
        )
        observers[f"central_{cube_side}^3"] = in_cube

    # "Deep cracked early" — cells in the cracked region at p=0.35
    if 0.35 in masks:
        observers["cracked_by_0.35"] = masks[0.35].copy()
    # "Cracked by 0.18" (right after percolation)
    if 0.18 in masks:
        observers["cracked_by_0.18"] = masks[0.18].copy()

    print()
    print("Observer regions:")
    for name, mask in observers.items():
        print(f"  {name}: {int(mask.sum())} cells "
              f"({mask.mean():.4f} of lattice)")

    # For each observer, compute local event count vs p
    print()
    print("Step 2: local event count at each p, per observer region")
    obs_events = {}
    for name, obs_mask in observers.items():
        events = observer_event_history(masks, obs_mask)
        obs_events[name] = events
        print(f"\n  Observer: {name}")
        print(f"  {'p':>5}  {'local_events':>13}  {'fraction_of_obs':>16}  {'gamma':>8}")
        for p in p_targets:
            ne = events[p]
            frac = ne / max(obs_mask.sum(), 1)
            g = gamma_p[p] if gamma_p[p] is not None else float("nan")
            print(f"  {p:>5.2f}  {ne:>13d}  {frac:>16.4f}  {g:>8.3f}")

    # Plot: gamma vs observer-local-event-count for each observer
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    flat = axes.flatten()
    cosmic_obs_range = (1.66, 1.91)
    cosmic_obs_mean = 1.82

    for idx, (name, events) in enumerate(obs_events.items()):
        if idx >= len(flat):
            break
        ax = flat[idx]
        ne_arr = np.array([events[p] for p in p_targets])
        g_arr = np.array([gamma_p[p] if gamma_p[p] is not None else np.nan
                          for p in p_targets])
        # Sort by event count
        order = np.argsort(ne_arr)
        ne_arr = ne_arr[order]
        g_arr = g_arr[order]

        ax.plot(ne_arr, g_arr, "o-", color="C0")
        # mark the cosmic observed range
        ax.axhspan(cosmic_obs_range[0], cosmic_obs_range[1], color="C1", alpha=0.2,
                   label=fr"observed $\gamma \in [{cosmic_obs_range[0]}, {cosmic_obs_range[1]}]$")
        ax.axhline(cosmic_obs_mean, color="C1", linestyle=":", alpha=0.7)
        # Annotate which p values these are
        for p in p_targets:
            ax.annotate(f"{p:.2f}", (events[p], gamma_p[p] or 0),
                        fontsize=7, alpha=0.6,
                        xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel("observer-local event count (cumulative)")
        ax.set_ylabel(r"$\gamma$")
        ax.set_title(f"Observer = {name}\n({int(observers[name].sum())} cells)",
                     fontsize=10)
        ax.set_xscale("symlog", linthresh=10)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

    for idx in range(len(obs_events), len(flat)):
        flat[idx].axis("off")

    fig.suptitle(r"$\gamma$ vs observer-local-event-count (non-linear time test)",
                 fontsize=12)
    fig.tight_layout()
    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_local_clock_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")
    plt.close(fig)

    # For each observer, compute how compressed the gamma=1.82 plateau
    # ends up being. Look at the range of observer-local-events for which
    # gamma is in [1.5, 2.0] (the observed cosmic range).
    print()
    print("Step 3: how much observer-local-time does the universe spend at "
          "observed gamma values?")
    print()
    print(f"{'observer':>25}  {'events at γ=2.0':>17}  {'events at γ=1.5':>17}  "
          f"{'ratio':>8}  {'flat span?':>12}")
    for name, events in obs_events.items():
        # Get gamma values sorted by p
        p_sorted = sorted(p_targets)
        g_arr = np.array([gamma_p[p] for p in p_sorted])
        ne_arr = np.array([events[p] for p in p_sorted])

        # Find p where gamma crosses 2.0 and 1.5
        ne_g20 = np.interp(2.0, g_arr[::-1], ne_arr[::-1]) if g_arr.min() < 2.0 else None
        ne_g15 = np.interp(1.5, g_arr[::-1], ne_arr[::-1]) if g_arr.min() < 1.5 else None

        if ne_g20 is not None and ne_g15 is not None and ne_g20 > 0:
            ratio = ne_g15 / ne_g20
            in_cosmic_range = (
                f"{ne_g20:.0f} -> {ne_g15:.0f}"
            )
        else:
            ratio = float("nan")
            in_cosmic_range = "n/a"
        print(f"{name:>25}  {ne_g20 or 'n/a':>17}  {ne_g15 or 'n/a':>17}  "
              f"{ratio:>8.2f}  {in_cosmic_range:>12}")

    print()
    print("Interpretation guide:")
    print("  If 'ratio' (events at γ=1.5 / events at γ=2.0) is LARGE for")
    print("  an observer, that observer's local clock RAPIDLY accumulates")
    print("  time while γ slowly drifts from 2.0 to 1.5 - matching the")
    print("  observed near-flat γ(z) across cosmic history. SMALL ratio")
    print("  means γ moves through 2.0->1.5 quickly in observer time -")
    print("  inconsistent with observation.")


if __name__ == "__main__":
    main()
