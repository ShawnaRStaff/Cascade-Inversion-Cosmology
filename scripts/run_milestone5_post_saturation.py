"""Stage 3 of the inversion-event redesign: post-saturation state.

After Stage 2 showed that the dynamics enters a 'plateau of catastrophes'
between p~0.95 and arrest at p~0.9998, the natural question is: what
does the substrate look like at arrest? Is it:
  - thermalized (uniform z, no structure)
  - patchy (residual geometric memory of the cascade history)
  - ready to restart (low-z everywhere, fresh substrate)

These are the cyclic-cosmology-relevant questions. If the post-arrest
state looks like a fresh pristine substrate, the model is compatible
with cyclic cosmology. If it's permanently locked in a wreck, it isn't.

Uses the Stage 2 saved data (no new simulation):
  - ever_toppled mask at arrest
  - grains_lost over time
  - sizes and durations time series

Computes:
  (3a) Energy balance: total grains_in vs grains_lost vs grains_stored
  (3b) Spatial structure of ever_toppled at arrest (fully saturated?)
  (3c) Time history of energy stored in lattice
  (3d) Plateau-of-catastrophes characterization

Run:
    .venv/bin/python scripts/run_milestone5_post_saturation.py [npz]
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


def latest_saturation_npz() -> Path:
    outdir = REPO_ROOT / "data" / "outputs"
    candidates = sorted(outdir.glob("manna_3d_saturation_data_*.npz"))
    if not candidates:
        raise FileNotFoundError("No manna_3d_saturation_data_*.npz")
    return candidates[-1]


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else latest_saturation_npz()
    print(f"Loading {path}")
    data = np.load(path)
    L = int(data["L"])
    L3 = L ** 3
    seeds = sorted({int(k.split("_")[0][1:])
                    for k in data.files if k.startswith("s") and "_" in k})
    print(f"L={L} (L^3={L3}), seeds={seeds}")
    print()

    # ============================================================
    # (3a) Energy balance at arrest
    # ============================================================
    print("(3a) Energy balance at arrest:")
    print(f"{'seed':>5}  {'grains_in':>11}  {'grains_lost':>12}  "
          f"{'grains_stored':>14}  {'fraction_in_lattice':>21}")
    energy_per_seed = {}
    for s in seeds:
        sizes_full = data[f"s{s}_sizes"]
        grains_lost_history = data[f"s{s}_grains_lost"]
        grains_in_total = int(sizes_full.size)
        grains_lost_final = int(grains_lost_history[-1])
        grains_stored = grains_in_total - grains_lost_final
        frac_in = grains_stored / max(grains_in_total, 1)
        energy_per_seed[s] = {
            "in": grains_in_total,
            "lost": grains_lost_final,
            "stored": grains_stored,
            "frac_in_lattice": frac_in,
        }
        print(f"{s:>5}  {grains_in_total:>11d}  {grains_lost_final:>12d}  "
              f"{grains_stored:>14d}  {frac_in:>21.4f}")

    # ============================================================
    # (3b) Spatial structure of ever_toppled at arrest
    # ============================================================
    print()
    print("(3b) Spatial structure of ever_toppled at arrest:")
    print(f"{'seed':>5}  {'p_final':>9}  {'n_frozen':>10}  "
          f"{'n_clusters':>11}  {'largest_void':>13}")
    for s in seeds:
        ever = data[f"s{s}_ever_toppled"]
        frozen = ~ever
        n_frozen = int(frozen.sum())
        if n_frozen > 0:
            struct = ndimage.generate_binary_structure(3, 1)
            labels, n_cl = ndimage.label(frozen, structure=struct)
            cluster_sizes = np.bincount(labels.ravel())[1:]
            largest = int(cluster_sizes.max()) if cluster_sizes.size else 0
        else:
            n_cl, largest = 0, 0
        p = float(ever.mean())
        print(f"{s:>5}  {p:>9.4f}  {n_frozen:>10d}  "
              f"{n_cl:>11d}  {largest:>13d}")

    # ============================================================
    # (3c) Time history of energy stored
    # ============================================================
    # grains stored in lattice at snapshot t = drops_at_t - grains_lost_at_t
    print()
    print("(3c) Energy time-history per seed (final 6 snapshots):")
    for s in seeds:
        drops = data[f"s{s}_drops"]
        gl = data[f"s{s}_grains_lost"]
        p_arr = data[f"s{s}_p"]
        stored = drops - gl
        print(f"  Seed {s}:")
        for i in range(max(0, len(drops) - 6), len(drops)):
            print(f"    drop={drops[i]:>7d}  p={p_arr[i]:.4f}  "
                  f"stored={stored[i]:>7d}  lost={int(gl[i]):>7d}  "
                  f"stored/L^3 = {stored[i]/L3:.4f}")

    # ============================================================
    # (3d) Plateau-of-catastrophes: characterize the late regime
    # ============================================================
    print()
    print("(3d) Plateau of catastrophes (drops with very large avalanches):")
    print(f"{'seed':>5}  {'drops_>L^3/10':>14}  {'drops_>L^3/4':>14}  "
          f"{'drops_>L^3/2':>14}  {'biggest':>9}")
    for s in seeds:
        sizes_full = data[f"s{s}_sizes"]
        n_tenth = int((sizes_full > L3 / 10).sum())
        n_quarter = int((sizes_full > L3 / 4).sum())
        n_half = int((sizes_full > L3 / 2).sum())
        biggest = int(sizes_full.max())
        print(f"{s:>5}  {n_tenth:>14d}  {n_quarter:>14d}  "
              f"{n_half:>14d}  {biggest:>9d}")

    # --- Plots ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))

    # Panel 1: stored vs lost over time
    ax = axes[0, 0]
    for s in seeds:
        drops = data[f"s{s}_drops"]
        stored = drops - data[f"s{s}_grains_lost"]
        ax.plot(drops, stored, "-", alpha=0.7, label=f"seed {s} stored")
    ax.set_xlabel("drops (grains_in)")
    ax.set_ylabel("grains stored in lattice")
    ax.set_title("Energy held by the lattice over time")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    # Panel 2: escape fraction
    ax = axes[0, 1]
    for s in seeds:
        drops = data[f"s{s}_drops"]
        escape_rate = data[f"s{s}_grains_lost"] / np.maximum(drops, 1)
        ax.plot(drops, escape_rate, "-", alpha=0.7, label=f"seed {s}")
    ax.set_xlabel("drops")
    ax.set_ylabel("fraction lost / fraction in")
    ax.set_title("Escape fraction (1 = fully thermalized)")
    ax.axhline(1.0, color="k", linestyle="--", alpha=0.4,
               label="full thermal balance")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    # Panel 3: catastrophic event distribution (cumulative count > threshold)
    ax = axes[1, 0]
    for s in seeds:
        sizes_full = data[f"s{s}_sizes"].astype(np.int64)
        # cumulative count of events larger than s (descending)
        # use logspaced thresholds
        thresholds = np.logspace(1, np.log10(sizes_full.max() + 1), 50)
        counts = np.array([(sizes_full >= t).sum() for t in thresholds])
        ax.loglog(thresholds, counts + 1e-9, "-", alpha=0.7,
                  label=f"seed {s}")
    ax.axvline(L3 / 2, color="r", linestyle=":", alpha=0.5,
               label="L^3/2 (half lattice)")
    ax.axvline(L3, color="k", linestyle=":", alpha=0.5,
               label="L^3 (whole lattice)")
    ax.set_xlabel("event size s")
    ax.set_ylabel("number of events with size >= s")
    ax.set_title("Cumulative event-size distribution")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    # Panel 4: spatial slice of ever_toppled at arrest, primary seed
    ax = axes[1, 1]
    s = seeds[0]
    ever = data[f"s{s}_ever_toppled"]
    mid = L // 2
    ax.imshow(ever[mid, :, :], cmap="Greys_r", interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    p = float(ever.mean())
    ax.set_title(f"ever_toppled mid-slice (seed {s}, p={p:.4f})")

    fig.suptitle("Stage 3 — Post-saturation state characterization",
                 fontsize=13)
    fig.tight_layout()

    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_post_saturation_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")

    # ============================================================
    # Interpretation summary
    # ============================================================
    print()
    print("=" * 50)
    print("Interpretation:")
    print("=" * 50)
    avg_frac_in = float(np.mean([energy_per_seed[s]["frac_in_lattice"]
                                  for s in seeds]))
    print(f"  At saturation, ~{avg_frac_in*100:.1f}% of input grains")
    print(f"  remain in the lattice; the rest escaped through the")
    print(f"  boundary. The substrate is a 'partial heat sink'.")
    print()
    # Are voids near gone?
    n_frozen_avg = float(np.mean([
        int((~data[f"s{s}_ever_toppled"]).sum()) for s in seeds
    ]))
    print(f"  Average frozen cells remaining: {n_frozen_avg:.0f} ({n_frozen_avg/L3*100:.2f}% of L^3)")
    if n_frozen_avg / L3 < 0.001:
        print(f"  -> substrate is essentially fully cracked. Almost no")
        print(f"     pristine 'fresh substrate' remains.")
    else:
        print(f"  -> {n_frozen_avg/L3*100:.2f}% of substrate remains pristine.")


if __name__ == "__main__":
    main()
