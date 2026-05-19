"""Stage 2 of the inversion-event redesign.

Drives simple 3D Manna past p=0.90 toward saturation to observe what
the dynamics naturally produces at the predicted inversion point
p* ~ 0.97.

We use the simple (single-state) Manna from M3 - faster and cleaner
than the multi-state variant. The question is whether catastrophic
spanning avalanches emerge naturally from the dynamics as the
substrate approaches full saturation.

Configuration:
  - L = 48 (small enough for fast iteration; large enough for spatial
            structure)
  - 3 seeds
  - n_drops_max = 500_000 (well past expected saturation)
  - snapshot every 2000 drops
  - track per-drop sizes, durations, and ever_toppled

Goal: capture the moment(s) when single events span huge fractions of
the lattice. Map the trajectory through the predicted inversion point.

Run:
    .venv/bin/python -u scripts/run_milestone5_saturation.py
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
SEEDS = [2000, 2001, 2002]
N_DROPS_MAX = 500_000
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
            snapshots.append({
                "drop": t_drop,
                "p": p,
                "mean_size_window": mean_size_w,
                "max_size_window": max_size_w,
                "grains_lost": int(state.grains_lost),
            })

            elapsed = time.time() - t0
            rate = t_drop / max(elapsed, 1e-9)
            # Print every 10 snapshots only, to avoid flood
            if t_drop % (SNAPSHOT_EVERY * 10) == 0 or p > 0.90:
                print(f"    [seed={seed}] drop={t_drop:>7d} "
                      f"p={p:.4f} <s>_win={mean_size_w:>7.1f} "
                      f"max_s={max_size_w:>6d} rate={rate:.0f}/s "
                      f"elapsed={elapsed:.0f}s", flush=True)

            # Auto-arrest if p has saturated (less than 0.0005 change over 20k drops)
            if len(snapshots) >= 10:
                recent_ps = [s["p"] for s in snapshots[-10:]]
                if max(recent_ps) - min(recent_ps) < 0.0005:
                    print(f"  [seed={seed}] saturated at p={p:.4f}, "
                          f"halting after {t_drop} drops")
                    break

    wall = time.time() - t0
    print(f"  [seed={seed}] DONE in {wall:.0f}s, drops={t_drop}, "
          f"final p={float(ever_toppled.mean()):.4f}, "
          f"max(s) reached = {int(sizes[:t_drop].max())}")
    return {
        "seed": seed,
        "wall_seconds": wall,
        "drops_executed": t_drop,
        "snapshots": snapshots,
        "sizes": sizes[:t_drop],
        "durations": durations[:t_drop],
        "ever_toppled": ever_toppled,
        "final_p": float(ever_toppled.mean()),
        "final_grains_lost": int(state.grains_lost),
    }


def main() -> None:
    print(f"M5 Stage 2: simple Manna driven to saturation, L={L}")
    print(f"  seeds={SEEDS}, n_drops_max={N_DROPS_MAX:,}")
    print(f"  snapshot every {SNAPSHOT_EVERY} drops")
    print("=" * 72)

    t_start = time.time()
    results = []
    for seed in SEEDS:
        results.append(run_one_seed(seed))
    total = time.time() - t_start
    print(f"\nAll seeds done. Total wall: {total:.0f}s")

    # --- Summary ---
    print()
    print(f"{'seed':>5}  {'drops':>9}  {'final p':>9}  {'max(s)':>9}  "
          f"{'wall':>7}")
    L3 = L ** 3
    for r in results:
        max_s = int(r["sizes"].max())
        print(f"{r['seed']:>5}  {r['drops_executed']:>9d}  "
              f"{r['final_p']:>9.4f}  {max_s:>9d}  ({max_s/L3*100:>4.1f}% of L^3)  "
              f"{r['wall_seconds']:>7.0f}s")

    # --- Plots ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))

    # p(drop)
    ax = axes[0, 0]
    for r in results:
        snaps = r["snapshots"]
        d = [s["drop"] for s in snaps]
        p = [s["p"] for s in snaps]
        ax.plot(d, p, "-", alpha=0.7, label=f"seed {r['seed']}")
    ax.axhline(0.97, color="r", linestyle="--", alpha=0.5,
               label=r"predicted $p^* = 0.97$")
    ax.set_xlabel("drops")
    ax.set_ylabel("p (fraction ever-toppled)")
    ax.set_title("Approach to saturation")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    # max(s) per window (the catastrophic-event tracker)
    ax = axes[0, 1]
    for r in results:
        snaps = r["snapshots"]
        d = [s["drop"] for s in snaps]
        ms = [s["max_size_window"] for s in snaps]
        ax.semilogy(d, ms, "-", alpha=0.7, label=f"seed {r['seed']}")
    ax.axhline(L3, color="k", linestyle=":", alpha=0.5,
               label=f"$L^3 = {L3}$ (full lattice)")
    ax.axhline(L3 / 2, color="0.5", linestyle=":", alpha=0.5,
               label="half lattice")
    ax.set_xlabel("drops")
    ax.set_ylabel("max(s) in window")
    ax.set_title("Catastrophic event size vs time")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    # max(s) vs p (where in the trajectory are catastrophic events?)
    ax = axes[1, 0]
    for r in results:
        snaps = r["snapshots"]
        p = np.array([s["p"] for s in snaps])
        ms = np.array([s["max_size_window"] for s in snaps])
        ax.semilogy(p, ms, "-", alpha=0.7, label=f"seed {r['seed']}")
    ax.axhline(L3, color="k", linestyle=":", alpha=0.5,
               label=f"$L^3 = {L3}$")
    ax.axvline(0.97, color="r", linestyle="--", alpha=0.5)
    ax.set_xlabel("p")
    ax.set_ylabel("max(s) in window")
    ax.set_title("Catastrophic event size vs p")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    # cumulative grains lost (boundary dissipation rate)
    ax = axes[1, 1]
    for r in results:
        snaps = r["snapshots"]
        d = [s["drop"] for s in snaps]
        gl = [s["grains_lost"] for s in snaps]
        ax.plot(d, gl, "-", alpha=0.7, label=f"seed {r['seed']}")
    # Reference: grains_in = drops; if all grains escape, this = drops
    drop_line = np.linspace(0, max((s["drop"] for r in results for s in r["snapshots"])), 100)
    ax.plot(drop_line, drop_line, "k:", alpha=0.5,
            label="grains_lost = grains_in (full saturation)")
    ax.set_xlabel("drops")
    ax.set_ylabel("cumulative grains lost")
    ax.set_title("Energy dissipation through boundary")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    fig.tight_layout()
    outdir = REPO_ROOT / "data" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_saturation_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")

    # --- Save raw data ---
    npz_path = outdir / f"manna_3d_saturation_data_{stamp}.npz"
    kw: dict = {"L": np.int64(L)}
    for r in results:
        s = r["seed"]
        snaps = r["snapshots"]
        kw[f"s{s}_drops"] = np.array([sn["drop"] for sn in snaps])
        kw[f"s{s}_p"] = np.array([sn["p"] for sn in snaps])
        kw[f"s{s}_mean_size"] = np.array([sn["mean_size_window"] for sn in snaps])
        kw[f"s{s}_max_size"] = np.array([sn["max_size_window"] for sn in snaps])
        kw[f"s{s}_grains_lost"] = np.array([sn["grains_lost"] for sn in snaps])
        kw[f"s{s}_sizes"] = r["sizes"]
        kw[f"s{s}_durations"] = r["durations"]
        kw[f"s{s}_ever_toppled"] = r["ever_toppled"]
        kw[f"s{s}_final_p"] = np.float64(r["final_p"])
    np.savez_compressed(npz_path, **kw)
    print(f"Saved raw data to {npz_path}")


if __name__ == "__main__":
    main()
