"""Experiment 5+6: active interface fractal dimension across p.

The active interface = cracked cells adjacent to frozen cells. This is
where new fracture events happen and where "time" is concentrated in
the event-time framework. The geometry of this interface is directly
comparable to the cosmic-web walls and filaments observed in galaxy
surveys, which have fractal dimension D~2.0-2.5.

For each p:
  1. Identify the interface: cracked cell with at least one frozen
     6-nearest-neighbor.
  2. Box-count the interface fractal dimension.
  3. Also compute surface-to-volume ratio of the cracked cluster
     (interface / cluster size) — peaks where the cluster is most
     porous.

Run:
    .venv/bin/python scripts/run_milestone4_interface.py [npz]
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.cluster_geometry import box_count_dimension  # noqa: E402


_STRUCT_6 = ndimage.generate_binary_structure(rank=3, connectivity=1)


def interface_mask(occupied: np.ndarray) -> np.ndarray:
    """Cracked cells adjacent to at least one frozen cell."""
    frozen = ~occupied
    # Dilate the frozen region by one cell with 6-conn structure, then
    # intersect with the occupied set: gives occupied cells that touch
    # frozen cells.
    dilated_frozen = ndimage.binary_dilation(frozen, structure=_STRUCT_6)
    return occupied & dilated_frozen


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

    primary_seed = seeds[0]

    print(f"{'p':>5}  {'interface':>10}  {'cluster':>10}  {'S/V':>8}  "
          f"{'D_box':>8}  {'D_err':>8}")

    summary = {}
    for p in p_targets:
        pp = int(round(p * 100))
        mask_key = f"mask_s{primary_seed}_p{pp:02d}"
        if mask_key not in data.files:
            continue
        occupied = data[mask_key]
        intf = interface_mask(occupied)
        n_intf = int(intf.sum())
        n_cluster = int(occupied.sum())
        sv = n_intf / max(n_cluster, 1)

        # Box-count the interface (the cosmologically interesting object)
        try:
            D, Derr, _, _ = box_count_dimension(intf)
        except Exception as e:
            D, Derr = np.nan, np.nan

        summary[p] = {
            "n_intf": n_intf,
            "n_cluster": n_cluster,
            "sv_ratio": sv,
            "D_box": D,
            "D_err": Derr,
        }
        print(f"{p:>5.2f}  {n_intf:>10d}  {n_cluster:>10d}  {sv:>8.4f}  "
              f"{D:>8.3f}  {Derr:>8.3f}")

    # --- Plot interface fractal dimension and S/V vs p ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    p_arr = sorted(summary.keys())
    D_arr = [summary[p]["D_box"] for p in p_arr]
    Derr_arr = [summary[p]["D_err"] for p in p_arr]
    sv_arr = [summary[p]["sv_ratio"] for p in p_arr]

    ax1.errorbar(p_arr, D_arr, yerr=Derr_arr, fmt="o-", color="C0", capsize=4,
                 label=r"interface $D_{\rm box}$")
    ax1.axhspan(2.0, 2.5, color="C2", alpha=0.15,
                label=r"cosmic-web filament/wall $D\approx 2.0\text{-}2.5$")
    ax1.axvline(0.65, color="0.5", linestyle="--", alpha=0.6,
                label=r"$\gamma=1.8$ crossing at $p=0.65$")
    ax1.set_xlabel("$p$")
    ax1.set_ylabel(r"$D_{\rm box}$ of active interface")
    ax1.set_title("Interface fractal dimension vs p")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9, loc="best")

    ax2.plot(p_arr, sv_arr, "o-", color="C3")
    ax2.axvline(0.65, color="0.5", linestyle="--", alpha=0.6)
    ax2.set_xlabel("$p$")
    ax2.set_ylabel("interface cells / cracked cells (S/V)")
    ax2.set_title("Cluster porosity (surface-to-volume) vs p")
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_interface_geometry_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")
    plt.close(fig)

    # --- Interface mid-slice visualization ---
    fig, axes = plt.subplots(2, 5, figsize=(15, 6))
    flat_axes = axes.flatten()
    for idx, p in enumerate(p_arr[:10]):
        if idx >= len(flat_axes):
            break
        pp = int(round(p * 100))
        mask_key = f"mask_s{primary_seed}_p{pp:02d}"
        occupied = data[mask_key]
        intf = interface_mask(occupied)
        mid = L // 2
        slc = intf[mid, :, :]
        flat_axes[idx].imshow(slc, cmap="hot", interpolation="nearest")
        flat_axes[idx].set_title(f"p={p:.2f}\nintf={int(intf.sum()):,}",
                                 fontsize=9)
        flat_axes[idx].set_xticks([])
        flat_axes[idx].set_yticks([])
    fig.suptitle(f"Active interface (mid-slice through L={L})\n"
                 "bright = cracked cells touching frozen cells",
                 fontsize=10)
    fig.tight_layout()
    out_path2 = outdir / f"manna_3d_interface_slices_{stamp}.png"
    fig.savefig(out_path2, dpi=150)
    print(f"Wrote {out_path2}")
    plt.close(fig)

    # --- interpretation ---
    print()
    if D_arr:
        D_at_065 = summary.get(0.65, {}).get("D_box", np.nan)
        if not np.isnan(D_at_065):
            print(f"Interface fractal dimension at p=0.65: D = {D_at_065:.3f}")
            print(f"Cosmic-web filament/wall reference: D ~ 2.0-2.5")
            if 1.9 <= D_at_065 <= 2.6:
                print("  -> Compatible with cosmic web filament/wall D.")
            else:
                print("  -> Outside cosmic web reference range.")
        # Peak S/V
        sv_peak_idx = int(np.argmax(sv_arr))
        print(f"Surface-to-volume peaks at p={p_arr[sv_peak_idx]:.2f}, "
              f"S/V = {sv_arr[sv_peak_idx]:.4f}")


if __name__ == "__main__":
    main()
