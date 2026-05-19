"""Stage 4 of the inversion-event redesign: drive past arrest.

Stage 3 showed that at p ~ 0.9998 the dynamics auto-arrests with:
  - ~53% of input grains still trapped in the lattice
  - Only ~0.02% of cells pristine
  - Stable terminal z density of ~0.61 grains/cell

But auto-arrest just means 'p hasn't changed in 20k drops.' That
doesn't mean the dynamics is *truly* terminal. The lattice may be in a
metastable state where new grains keep being added but mostly fall off
the boundary, until eventually a new catastrophe nucleates.

This script disables auto-arrest and drives for 250k drops total
(roughly twice as long as the run that auto-arrested). We watch what
happens after the apparent terminal state:
  - Does the lattice keep absorbing grains, or do they all escape?
  - Do new large avalanches still fire?
  - Does p ever change again?
  - Is the substrate in a stable steady-state or a long-lived
    metastable state with eventual transitions?

Also saves the final z field this time (Stage 3 couldn't because
Stage 2 didn't save it). This will enable post-arrest z-distribution
analysis.

Run:
    .venv/bin/python -u scripts/run_milestone5_past_arrest.py
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

from void_cascade.sandpile_3d import drive, initialize, relax  # noqa: E402


L = 48
SEEDS = [2000, 2001]  # 2 seeds is enough for this follow-up
N_DROPS_MAX = 250_000  # ~2x what auto-arrest took before
SNAPSHOT_EVERY = 2_000


def run_one_seed(seed: int) -> dict:
    print(f"  [seed={seed}] starting...")
    t0 = time.time()
    rng = np.random.default_rng(seed)
    state = initialize(L)
    ever_toppled = np.zeros((L, L, L), dtype=bool)
    sizes = np.zeros(N_DROPS_MAX, dtype=np.int32)
    durations = np.zeros(N_DROPS_MAX, dtype=np.int32)

    snapshots: list[dict] = []
    t_drop = 0

    while t_drop < N_DROPS_MAX:
        drive(state, rng)
        s, T, mask = relax(state, rng, track_support=True)
        sizes[t_drop] = s
        durations[t_drop] = T
        if mask is not None:
            ever_toppled |= mask
        t_drop += 1

        if t_drop % SNAPSHOT_EVERY == 0:
            p = float(ever_toppled.mean())
            window_sizes = sizes[max(0, t_drop - SNAPSHOT_EVERY):t_drop]
            mean_size_w = float(window_sizes.mean())
            max_size_w = int(window_sizes.max())
            z_sum = int(state.z.sum())
            z_mean = float(state.z.mean())
            z_std = float(state.z.std())
            snapshots.append({
                "drop": t_drop,
                "p": p,
                "mean_size_window": mean_size_w,
                "max_size_window": max_size_w,
                "grains_lost": int(state.grains_lost),
                "z_sum": z_sum,
                "z_mean": z_mean,
                "z_std": z_std,
                "stored_per_L3": z_sum / (L ** 3),
            })

            # NO AUTO-ARREST. We want to see what happens past arrest.
            # Print every 10 snapshots, or every snapshot once p > 0.99
            if t_drop % (SNAPSHOT_EVERY * 10) == 0 or p > 0.99:
                elapsed = time.time() - t0
                rate = t_drop / max(elapsed, 1e-9)
                print(f"    [seed={seed}] drop={t_drop:>7d} p={p:.5f} "
                      f"<s>_win={mean_size_w:>7.1f} max_s={max_size_w:>6d} "
                      f"z_avg={z_mean:.3f} z_std={z_std:.3f} "
                      f"stored/L3={z_sum/(L**3):.4f} "
                      f"rate={rate:.0f}/s", flush=True)

    wall = time.time() - t0
    print(f"  [seed={seed}] DONE in {wall:.0f}s, drops={t_drop}, "
          f"final p={float(ever_toppled.mean()):.5f}, "
          f"max(s) reached = {int(sizes[:t_drop].max())}")
    return {
        "seed": seed,
        "wall_seconds": wall,
        "drops_executed": t_drop,
        "snapshots": snapshots,
        "sizes": sizes[:t_drop],
        "durations": durations[:t_drop],
        "ever_toppled": ever_toppled,
        "final_z": state.z.copy(),
        "final_p": float(ever_toppled.mean()),
        "final_grains_lost": int(state.grains_lost),
    }


def main() -> None:
    print(f"M5 Stage 4: drive past arrest, L={L}")
    print(f"  seeds={SEEDS}, n_drops_max={N_DROPS_MAX:,}")
    print(f"  snapshot every {SNAPSHOT_EVERY} drops")
    print(f"  AUTO-ARREST DISABLED")
    print("=" * 72)

    t_start = time.time()
    results = []
    for seed in SEEDS:
        results.append(run_one_seed(seed))
    total = time.time() - t_start
    print(f"\nAll seeds done. Total wall: {total:.0f}s")

    # Summary
    L3 = L ** 3
    print()
    print(f"{'seed':>5}  {'drops':>9}  {'final p':>10}  {'final z_avg':>13}  "
          f"{'final lost':>11}  {'biggest':>9}")
    for r in results:
        z_avg = float(r["final_z"].mean())
        print(f"{r['seed']:>5}  {r['drops_executed']:>9d}  "
              f"{r['final_p']:>10.5f}  {z_avg:>13.4f}  "
              f"{r['final_grains_lost']:>11d}  {int(r['sizes'].max()):>9d}")

    # Plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))

    # Panel 1: z stored over time
    ax = axes[0, 0]
    for r in results:
        snaps = r["snapshots"]
        d = [s["drop"] for s in snaps]
        z = [s["z_sum"] for s in snaps]
        ax.plot(d, z, "-", alpha=0.7, label=f"seed {r['seed']}")
    ax.set_xlabel("drops")
    ax.set_ylabel("total z in lattice")
    ax.set_title("Energy stored in lattice (past arrest)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    # Panel 2: max(s) over time, log scale
    ax = axes[0, 1]
    for r in results:
        snaps = r["snapshots"]
        d = [s["drop"] for s in snaps]
        m = [s["max_size_window"] for s in snaps]
        ax.semilogy(d, m, "-", alpha=0.7, label=f"seed {r['seed']}")
    ax.axhline(L3, color="k", linestyle=":", alpha=0.5,
               label=f"L^3 = {L3}")
    ax.set_xlabel("drops")
    ax.set_ylabel("max(s) per window")
    ax.set_title("Catastrophic events past arrest")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    # Panel 3: z distribution at final state
    ax = axes[1, 0]
    for r in results:
        z_final = r["final_z"].ravel()
        # Histogram of z values
        max_z = max(int(z_final.max()), 2)
        bins = np.arange(0, max_z + 2) - 0.5
        ax.hist(z_final, bins=bins, alpha=0.5,
                label=f"seed {r['seed']} (mean={z_final.mean():.3f})")
    ax.set_xlabel("z (grains per cell)")
    ax.set_ylabel("number of cells")
    ax.set_title("Final z distribution per cell")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    # Panel 4: stored vs total grains_in to see if balance shifts
    ax = axes[1, 1]
    for r in results:
        snaps = r["snapshots"]
        d = np.array([s["drop"] for s in snaps])
        gl = np.array([s["grains_lost"] for s in snaps])
        stored = d - gl
        # plot stored as fraction of grains_in
        ax.plot(d, stored / np.maximum(d, 1), "-", alpha=0.7,
                label=f"seed {r['seed']}")
    ax.set_xlabel("drops")
    ax.set_ylabel("fraction of grains_in still in lattice")
    ax.set_title("Lattice retention fraction (does it asymptote?)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    fig.tight_layout()
    outdir = REPO_ROOT / "data" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_past_arrest_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")

    # Save raw data
    npz_path = outdir / f"manna_3d_past_arrest_data_{stamp}.npz"
    kw: dict = {"L": np.int64(L)}
    for r in results:
        s = r["seed"]
        snaps = r["snapshots"]
        kw[f"s{s}_drops"] = np.array([sn["drop"] for sn in snaps])
        kw[f"s{s}_p"] = np.array([sn["p"] for sn in snaps])
        kw[f"s{s}_mean_size"] = np.array([sn["mean_size_window"] for sn in snaps])
        kw[f"s{s}_max_size"] = np.array([sn["max_size_window"] for sn in snaps])
        kw[f"s{s}_grains_lost"] = np.array([sn["grains_lost"] for sn in snaps])
        kw[f"s{s}_z_sum"] = np.array([sn["z_sum"] for sn in snaps])
        kw[f"s{s}_z_mean"] = np.array([sn["z_mean"] for sn in snaps])
        kw[f"s{s}_z_std"] = np.array([sn["z_std"] for sn in snaps])
        kw[f"s{s}_sizes"] = r["sizes"]
        kw[f"s{s}_durations"] = r["durations"]
        kw[f"s{s}_ever_toppled"] = r["ever_toppled"]
        kw[f"s{s}_final_z"] = r["final_z"]
        kw[f"s{s}_final_p"] = np.float64(r["final_p"])
    np.savez_compressed(npz_path, **kw)
    print(f"Saved raw data to {npz_path}")


if __name__ == "__main__":
    main()
