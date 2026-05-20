"""M5 follow-up: larger-L saturation test.

L=48 in Stage 4 reached p=1.0 with peak event sizes of 71% of lattice
but never produced a 100% (full-lattice-spanning) event in a single
drop. Question: is the 71% ceiling a finite-size effect that lifts at
larger L, or a genuine asymptotic property?

Larger L = more cells in possible cascade chain. If max(s)/L^3 ratio
stays the same (~0.71), it's intrinsic. If max(s)/L^3 grows with L,
the model approaches 100% events at large L.

L=64 with no auto-arrest. 2 seeds. Drive to natural p=1.0 plus ~25%
margin.

Run:
    .venv/bin/python -u scripts/run_milestone5_larger_L.py
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


L = 64
SEEDS = [3000, 3001]
N_DROPS_MAX = 600_000  # L^3 scale: ~3x what L=48 needed
SNAPSHOT_EVERY = 5_000


def run_one_seed(seed: int) -> dict:
    print(f"  [seed={seed}] L={L} starting...")
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
            z_mean = float(state.z.mean())
            snapshots.append({
                "drop": t_drop, "p": p,
                "mean_size_window": mean_size_w, "max_size_window": max_size_w,
                "z_mean": z_mean, "grains_lost": int(state.grains_lost),
            })
            elapsed = time.time() - t0
            rate = t_drop / max(elapsed, 1e-9)
            # Print every 5 snapshots before p=0.9, every snapshot after
            if t_drop % (SNAPSHOT_EVERY * 5) == 0 or p > 0.9:
                L3 = L ** 3
                pct = max_size_w / L3 * 100
                print(f"    [seed={seed}] drop={t_drop:>7d} p={p:.5f} "
                      f"<s>_win={mean_size_w:>6.1f} max_s={max_size_w:>6d} "
                      f"(max%={pct:>4.1f}) z={z_mean:.3f} rate={rate:.0f}/s "
                      f"elapsed={elapsed:.0f}s", flush=True)

    wall = time.time() - t0
    print(f"  [seed={seed}] DONE in {wall:.0f}s drops={t_drop} "
          f"final p={float(ever_toppled.mean()):.5f} "
          f"max(s) reached = {int(sizes[:t_drop].max())} "
          f"({int(sizes[:t_drop].max()) / L**3 * 100:.1f}% of L^3)")
    return {
        "seed": seed, "wall_seconds": wall, "drops_executed": t_drop,
        "snapshots": snapshots, "sizes": sizes[:t_drop],
        "durations": durations[:t_drop], "ever_toppled": ever_toppled,
        "final_z": state.z.copy(), "final_p": float(ever_toppled.mean()),
    }


def main() -> None:
    print(f"M5 larger-L: L={L} (L^3={L**3}), seeds={SEEDS}")
    print(f"  N_DROPS_MAX = {N_DROPS_MAX:,}")
    print("=" * 72)

    t_start = time.time()
    results = []
    for seed in SEEDS:
        results.append(run_one_seed(seed))
    total = time.time() - t_start
    print(f"\nAll seeds done. Total wall: {total:.0f}s")

    L3 = L ** 3
    print()
    print(f"{'seed':>5}  {'drops':>9}  {'final p':>10}  {'max(s)':>9}  "
          f"{'max%':>7}  {'wall':>7}")
    for r in results:
        max_s = int(r["sizes"].max())
        print(f"{r['seed']:>5}  {r['drops_executed']:>9d}  "
              f"{r['final_p']:>10.5f}  {max_s:>9d}  "
              f"{max_s/L3*100:>7.2f}%  {r['wall_seconds']:>7.0f}s")

    # Plots
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    ax = axes[0, 0]
    for r in results:
        snaps = r["snapshots"]
        d = [s["drop"] for s in snaps]
        p = [s["p"] for s in snaps]
        ax.plot(d, p, "-", alpha=0.7, label=f"seed {r['seed']}")
    ax.set_xlabel("drops")
    ax.set_ylabel("p")
    ax.set_title(f"p(drop), L={L}")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    for r in results:
        snaps = r["snapshots"]
        d = [s["drop"] for s in snaps]
        m = [s["max_size_window"] for s in snaps]
        ax.semilogy(d, m, "-", alpha=0.7, label=f"seed {r['seed']}")
    ax.axhline(L3, color="k", linestyle=":", alpha=0.5, label=f"L^3={L3}")
    ax.axhline(L3 / 2, color="0.5", linestyle=":", alpha=0.5)
    ax.set_xlabel("drops")
    ax.set_ylabel("max(s) per window")
    ax.set_title(f"Catastrophe size, L={L}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    # peak max(s) per p bin to show ceiling
    ax = axes[1, 0]
    for r in results:
        snaps = r["snapshots"]
        p = np.array([s["p"] for s in snaps])
        m = np.array([s["max_size_window"] for s in snaps])
        ax.plot(p, m / L3, "-", alpha=0.7, label=f"seed {r['seed']}")
    ax.axhline(1.0, color="r", linestyle="--", alpha=0.5,
               label="full lattice (100%)")
    ax.set_xlabel("p")
    ax.set_ylabel("max(s) / L^3")
    ax.set_title(f"Catastrophe ceiling, L={L}")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    # z mean over time
    ax = axes[1, 1]
    for r in results:
        snaps = r["snapshots"]
        d = [s["drop"] for s in snaps]
        z = [s["z_mean"] for s in snaps]
        ax.plot(d, z, "-", alpha=0.7, label=f"seed {r['seed']}")
    ax.set_xlabel("drops")
    ax.set_ylabel("z_avg")
    ax.set_title(f"Energy density, L={L}")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    fig.tight_layout()
    outdir = REPO_ROOT / "data" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_larger_L_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print(f"\nWrote {out_path}")

    npz_path = outdir / f"manna_3d_larger_L_data_{stamp}.npz"
    kw: dict = {"L": np.int64(L)}
    for r in results:
        s = r["seed"]
        snaps = r["snapshots"]
        kw[f"s{s}_drops"] = np.array([sn["drop"] for sn in snaps])
        kw[f"s{s}_p"] = np.array([sn["p"] for sn in snaps])
        kw[f"s{s}_max_size"] = np.array([sn["max_size_window"] for sn in snaps])
        kw[f"s{s}_mean_size"] = np.array([sn["mean_size_window"] for sn in snaps])
        kw[f"s{s}_z_mean"] = np.array([sn["z_mean"] for sn in snaps])
        kw[f"s{s}_sizes"] = r["sizes"]
        kw[f"s{s}_ever_toppled"] = r["ever_toppled"]
        kw[f"s{s}_final_z"] = r["final_z"]
    np.savez_compressed(npz_path, **kw)
    print(f"Saved raw data to {npz_path}")


if __name__ == "__main__":
    main()
