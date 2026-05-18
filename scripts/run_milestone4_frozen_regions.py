"""Experiment 3: frozen-region analysis.

If time = events, then regions where no events have ever happened are
'time-empty.' The substrate is not uniformly clocked; there are
patches where time has never started. Characterizing these gives a
direct picture of the dynamical heterogeneity of the substrate.

For each p, identify:
  - the always-frozen set (cells that have never toppled up to this p)
  - its connected-component structure
  - its size distribution
  - how it shrinks/fragments as p grows

Compare the spatial structure of frozen voids to observed cosmic voids
(another testable observable for the model: void size distribution).

Run:
    .venv/bin/python scripts/run_milestone4_frozen_regions.py [npz]
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.mass_function import cluster_sizes  # noqa: E402
from void_cascade.percolation import check_spanning  # noqa: E402


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

    # For each (p, seed), compute the FROZEN set (the inverse of
    # ever_toppled), its connected components, and various metrics.
    print(f"{'p':>5}  {'frozen_frac':>11}  {'n_frozen_cl':>12}  {'largest_frozen':>14}  "
          f"{'largest/frozen':>14}  {'frozen_spans':>12}")
    summary = {}
    primary_seed = seeds[0]
    for p in p_targets:
        pp = int(round(p * 100))
        frozen_fracs = []
        largest_fracs = []
        spans = []
        n_clusters_list = []
        for s in seeds:
            mask_key = f"mask_s{s}_p{pp:02d}"
            if mask_key not in data.files:
                continue
            frozen = ~data[mask_key]
            frozen_total = int(frozen.sum())
            if frozen_total == 0:
                continue
            sizes = cluster_sizes(frozen)
            n_clusters = sizes.size
            largest = int(sizes.max()) if sizes.size else 0
            r = check_spanning(frozen)
            frozen_fracs.append(frozen_total / L ** 3)
            largest_fracs.append(largest / frozen_total)
            spans.append(int(r.percolates))
            n_clusters_list.append(n_clusters)
        if not frozen_fracs:
            continue
        summary[p] = {
            "frozen_frac": float(np.mean(frozen_fracs)),
            "n_clusters": float(np.mean(n_clusters_list)),
            "largest_in_frozen_frac": float(np.mean(largest_fracs)),
            "spans_frac": float(np.mean(spans)),
        }
        s = summary[p]
        print(f"{p:>5.2f}  {s['frozen_frac']:>11.4f}  "
              f"{s['n_clusters']:>12.0f}  "
              f"{int(s['largest_in_frozen_frac'] * s['frozen_frac'] * L ** 3):>14d}  "
              f"{s['largest_in_frozen_frac']:>14.4f}  "
              f"{s['spans_frac']:>12.2f}")

    # --- frozen-region size distribution at a few p values ---
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_p_values = [0.10, 0.35, 0.50, 0.65, 0.80]
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(plot_p_values)))
    for p, color in zip(plot_p_values, colors):
        pp = int(round(p * 100))
        all_sizes = []
        for s in seeds:
            mask_key = f"mask_s{s}_p{pp:02d}"
            if mask_key not in data.files:
                continue
            frozen = ~data[mask_key]
            sizes = cluster_sizes(frozen)
            all_sizes.append(sizes)
        if not all_sizes:
            continue
        sizes = np.concatenate(all_sizes)
        if sizes.size == 0:
            continue
        # log-binned histogram
        bins = np.logspace(0, np.log10(sizes.max() + 1), 25)
        counts, edges = np.histogram(sizes, bins=bins)
        centers = np.sqrt(edges[:-1] * edges[1:])
        widths = np.diff(edges)
        keep = counts > 0
        pdf_full = counts / (sizes.size * widths)
        ax.loglog(centers[keep], pdf_full[keep], "o-", color=color, alpha=0.7,
                  label=f"p={p:.2f}, $\\rho_{{frozen}}={1-p:.2f}$")
    ax.set_xlabel("size of frozen cluster (cells)")
    ax.set_ylabel("PDF $n(s)$")
    ax.set_title("Frozen-region (time-empty) size distribution across p")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_frozen_distribution_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")
    plt.close(fig)

    # --- Frozen-region mid-slice at several p values ---
    fig, axes = plt.subplots(2, 5, figsize=(15, 6))
    flat_axes = axes.flatten()
    show_p = p_targets[:10]
    for idx, p in enumerate(show_p):
        if idx >= len(flat_axes):
            break
        pp = int(round(p * 100))
        mask_key = f"mask_s{primary_seed}_p{pp:02d}"
        if mask_key not in data.files:
            flat_axes[idx].axis("off")
            continue
        frozen = ~data[mask_key]
        mid = L // 2
        slc = frozen[mid, :, :]
        flat_axes[idx].imshow(slc, cmap="Greys", interpolation="nearest")
        flat_axes[idx].set_title(f"p={p:.2f}\nfrozen={float(slc.mean()):.3f}",
                                 fontsize=9)
        flat_axes[idx].set_xticks([])
        flat_axes[idx].set_yticks([])
    for idx in range(len(show_p), len(flat_axes)):
        flat_axes[idx].axis("off")
    fig.suptitle(f"Frozen (time-empty) regions — mid-slice through L={L}\n"
                 f"dark = frozen, white = cracked", fontsize=10)
    fig.tight_layout()
    out_path2 = outdir / f"manna_3d_frozen_slices_{stamp}.png"
    fig.savefig(out_path2, dpi=150)
    print(f"Wrote {out_path2}")
    plt.close(fig)

    # --- summary metrics plot ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    p_values = sorted(summary.keys())
    frozen_fracs = [summary[p]["frozen_frac"] for p in p_values]
    largest_fracs = [summary[p]["largest_in_frozen_frac"] for p in p_values]
    spans_fracs = [summary[p]["spans_frac"] for p in p_values]
    n_clusters = [summary[p]["n_clusters"] for p in p_values]

    ax1.plot(p_values, frozen_fracs, "o-", color="C0", label="frozen fraction of lattice")
    ax1.plot(p_values, largest_fracs, "s-", color="C3",
             label="largest-frozen / total-frozen")
    ax1.plot(p_values, spans_fracs, "^-", color="C2",
             label="fraction of seeds where frozen set spans")
    ax1.axvline(0.65, color="0.5", linestyle="--", alpha=0.6,
                label="γ=1.8 crossing at p=0.65")
    ax1.set_xlabel("$p$")
    ax1.set_ylabel("fraction")
    ax1.set_title("Frozen-region structure vs p")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8)

    ax2.plot(p_values, n_clusters, "o-", color="C1")
    ax2.axvline(0.65, color="0.5", linestyle="--", alpha=0.6)
    ax2.set_xlabel("$p$")
    ax2.set_ylabel("number of frozen clusters")
    ax2.set_title("Number of disconnected frozen regions vs p")
    ax2.set_yscale("log")
    ax2.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    out_path3 = outdir / f"manna_3d_frozen_metrics_{stamp}.png"
    fig.savefig(out_path3, dpi=150)
    print(f"Wrote {out_path3}")
    plt.close(fig)

    # --- interpretation ---
    print()
    print("Interpretation:")
    # Find where frozen set starts to fragment fast (n_clusters takes off)
    nc = np.array(n_clusters)
    p_arr = np.array(p_values)
    # The point where n_clusters first exceeds e.g. 100
    fragment_p = None
    for p, n in zip(p_values, n_clusters):
        if n > 100:
            fragment_p = p
            break
    if fragment_p is not None:
        print(f"  Frozen set fragments past p={fragment_p:.2f} "
              f"(n_clusters > 100). Before this, frozen region is essentially "
              f"one big connected void.")

    # When does the largest-frozen / total-frozen ratio drop below 0.99?
    for p, lf in zip(p_values, largest_fracs):
        if lf < 0.99:
            print(f"  Largest frozen cluster drops below 99% of frozen mass at "
                  f"p={p:.2f}. Before this, virtually all 'time-empty' is one "
                  f"connected region.")
            break


if __name__ == "__main__":
    main()
