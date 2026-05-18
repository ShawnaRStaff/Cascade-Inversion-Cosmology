"""Experiment 1: local event rate map. Test of time-as-events.

Hypothesis: if "time" is the count of fracture events in the substrate
rather than a uniform external parameter, then time should flow at very
different rates in different regions of the substrate. Cells deep in
the already-cracked region see no new events (their local clocks have
stopped). Cells at the cluster boundary see frequent events (their
clocks tick fast). Cells in voids see nothing (their clocks have not
started). The implications:

  - The uniform "one drop = one time unit" mapping we have been using
    is wrong.
  - The cosmological gamma(z) match at p=0.65 may not require a free
    time-mapping parameter; it may be predicted from the dynamics if
    cosmic time is built from local event counts rather than global
    drop counts.

This script tests the hypothesis. For each consecutive pair of saved
masks at (p_i, p_{i+1}), we compute the set of cells that newly toppled
in that window via XOR. We then ask:

  1. Are the new events spatially concentrated (clustered) or uniform?
  2. Does the spatial distribution of new events change qualitatively
     across the trajectory (different morphology at different p)?
  3. What is the variance in local event rate across cells, normalized
     by mean?

Outputs:
  - Per-window event-density slices through the lattice
  - Spatial variance of event rate as a function of p
  - Comparison to "uniform random fill" expectation

Run:
    .venv/bin/python scripts/run_milestone4_event_rate.py [npz]
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
    print()

    # We'll work with seed 1000 (first seed) for spatial diagnostics,
    # and average across seeds for the variance metric.
    primary_seed = seeds[0]

    # Compute per-window new-event masks (cells that toppled in this window)
    # using XOR of consecutive cumulative masks. Each "window" is the
    # interval [p_i, p_{i+1}].
    print(f"{'window':>14}  {'p_lo':>5}  {'p_hi':>5}  {'new_cells':>10}  "
          f"{'frac_of_L^3':>11}  {'spatial_CV':>10}  {'drops_in_win':>12}")
    windows: list[dict] = []
    for i in range(len(p_targets) - 1):
        p_lo = p_targets[i]
        p_hi = p_targets[i + 1]
        pp_lo = int(round(p_lo * 100))
        pp_hi = int(round(p_hi * 100))

        new_cells_per_seed = []
        drops_in_win_per_seed = []
        for s in seeds:
            mask_lo_key = f"mask_s{s}_p{pp_lo:02d}"
            mask_hi_key = f"mask_s{s}_p{pp_hi:02d}"
            drop_lo_key = f"drop_s{s}_p{pp_lo:02d}"
            drop_hi_key = f"drop_s{s}_p{pp_hi:02d}"
            if mask_lo_key not in data.files or mask_hi_key not in data.files:
                continue
            new_mask = data[mask_hi_key] & ~data[mask_lo_key]
            new_cells_per_seed.append(int(new_mask.sum()))
            drops_in_win_per_seed.append(int(data[drop_hi_key] - data[drop_lo_key]))

        if not new_cells_per_seed:
            continue
        new_cells_avg = float(np.mean(new_cells_per_seed))
        drops_in_win = float(np.mean(drops_in_win_per_seed))

        # Spatial coefficient of variation: how non-uniform is the
        # event-rate map? Compute by coarse-graining the new_mask of the
        # primary seed into 8x8x8 blocks and looking at counts per block.
        mask_lo_key = f"mask_s{primary_seed}_p{pp_lo:02d}"
        mask_hi_key = f"mask_s{primary_seed}_p{pp_hi:02d}"
        if mask_lo_key in data.files and mask_hi_key in data.files:
            new_mask = data[mask_hi_key] & ~data[mask_lo_key]
            b = 8
            if L % b == 0:
                tiled = new_mask.reshape(L // b, b, L // b, b, L // b, b)
                per_block = tiled.sum(axis=(1, 3, 5)).ravel()
            else:
                per_block = np.array([new_mask.sum()])
            mean_per_block = float(per_block.mean())
            std_per_block = float(per_block.std())
            cv = std_per_block / mean_per_block if mean_per_block > 0 else 0.0
        else:
            cv = np.nan

        windows.append({
            "p_lo": p_lo, "p_hi": p_hi,
            "new_cells_avg": new_cells_avg,
            "drops_in_win": drops_in_win,
            "events_per_drop": new_cells_avg / max(drops_in_win, 1),
            "spatial_cv": cv,
        })
        print(f"{p_lo:.2f}->{p_hi:.2f}    {p_lo:5.2f}  {p_hi:5.2f}  "
              f"{new_cells_avg:>10.0f}  {new_cells_avg / L**3:>11.4f}  "
              f"{cv:>10.3f}  {drops_in_win:>12.0f}")

    # --- Plot per-window event density slices (primary seed) for visual ---
    fig, axes = plt.subplots(3, 4, figsize=(14, 11))
    flat_axes = axes.flatten()
    for idx, w in enumerate(windows[:12]):
        ax = flat_axes[idx]
        pp_lo = int(round(w["p_lo"] * 100))
        pp_hi = int(round(w["p_hi"] * 100))
        mask_lo_key = f"mask_s{primary_seed}_p{pp_lo:02d}"
        mask_hi_key = f"mask_s{primary_seed}_p{pp_hi:02d}"
        if mask_lo_key in data.files and mask_hi_key in data.files:
            new_mask = data[mask_hi_key] & ~data[mask_lo_key]
            mid = L // 2
            slc = new_mask[mid, :, :]
            ax.imshow(slc, cmap="hot", interpolation="nearest")
            ax.set_title(f"window {w['p_lo']:.2f}->{w['p_hi']:.2f}\n"
                         f"new cells = {int(w['new_cells_avg']):,}",
                         fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
    for idx in range(len(windows), len(flat_axes)):
        flat_axes[idx].axis("off")
    fig.suptitle(f"Local event-density (mid-plane slice through L={L} lattice)\n"
                 f"bright = cells that fractured in this p-window")
    fig.tight_layout()

    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_event_rate_slices_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")
    plt.close(fig)

    # --- Plot spatial variance + events-per-drop vs p ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    p_centers = [0.5 * (w["p_lo"] + w["p_hi"]) for w in windows]
    cvs = [w["spatial_cv"] for w in windows]
    epd = [w["events_per_drop"] for w in windows]
    ax1.plot(p_centers, cvs, "o-", color="C2",
             label="spatial CV of event density")
    # For a uniform random fill into a box, CV across 8^3=512 blocks at
    # rate ~n_new/512 per block: sigma = sqrt(rate); CV = 1/sqrt(rate).
    # For typical n_new~10000 in our windows, rate~20, CV_uniform~0.22.
    ax1.axhline(0.22, color="gray", linestyle=":",
                label="uniform-random-fill CV reference (~0.22)")
    ax1.set_xlabel("$p$  (window midpoint)")
    ax1.set_ylabel(r"spatial coefficient of variation of event density")
    ax1.set_title("Event-density non-uniformity vs p")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9)

    ax2.plot(p_centers, epd, "o-", color="C4")
    ax2.set_xlabel("$p$  (window midpoint)")
    ax2.set_ylabel("new cells fractured per drop")
    ax2.set_title("Local fracture efficiency vs p")
    ax2.grid(True, alpha=0.3)
    ax2.set_yscale("log")

    fig.tight_layout()
    out_path2 = outdir / f"manna_3d_event_rate_metrics_{stamp}.png"
    fig.savefig(out_path2, dpi=150)
    print(f"Wrote {out_path2}")
    plt.close(fig)

    # --- Interpretation summary ---
    print()
    print("Interpretation:")
    if cvs:
        max_cv_idx = int(np.argmax(cvs))
        print(f"  Most spatially clustered events: window "
              f"p={windows[max_cv_idx]['p_lo']:.2f}->{windows[max_cv_idx]['p_hi']:.2f}"
              f", CV={cvs[max_cv_idx]:.3f}")
        min_cv_idx = int(np.argmin(cvs))
        print(f"  Most uniform events: window "
              f"p={windows[min_cv_idx]['p_lo']:.2f}->{windows[min_cv_idx]['p_hi']:.2f}"
              f", CV={cvs[min_cv_idx]:.3f}")
    if epd:
        max_epd_idx = int(np.argmax(epd))
        print(f"  Most efficient fracturing (highest events/drop): window "
              f"p={windows[max_epd_idx]['p_lo']:.2f}->{windows[max_epd_idx]['p_hi']:.2f}"
              f", {epd[max_epd_idx]:.2f} events/drop")


if __name__ == "__main__":
    main()
