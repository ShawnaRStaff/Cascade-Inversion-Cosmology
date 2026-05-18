"""Milestone 5: multi-state substrate dynamics, bounded vs unbounded sigma.

Tests whether the dynamics produces natural arrest at some p when the
substrate has higher-tier matter (sigma > 0) whose binding energy
grows with tier (z_c = 2*(sigma+1)).

Two variants run on the same seeds for direct comparison:
  - bounded: sigma_max = 4 (tier ceiling)
  - unbounded: sigma can grow without bound

Measurements per variant:
  - p(drop) trajectory: fraction of lattice that has ever toppled
  - sigma distribution at periodic snapshots
  - sigma-mean over time
  - avalanche size and duration per drop
  - terminal state of the dynamics (does p saturate?)

Run:
    .venv/bin/python -u scripts/run_milestone5_multistate.py
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

from void_cascade.sandpile_3d_multistate import (  # noqa: E402
    drive,
    initialize,
    relax,
)


L = 64
SEEDS = [1000, 1001, 1002]
N_DROPS_MAX = 600_000
SNAPSHOT_EVERY = 5_000  # take a snapshot every 5k drops


def run_variant(sigma_max: int | None, seed: int) -> dict:
    """Run one variant with one seed; snapshot periodically."""
    label = f"sigma_max={sigma_max if sigma_max is not None else 'inf'}"
    print(f"  [seed={seed}, {label}] starting...")
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
        s, T, mask = relax(state, rng, track_support=True, sigma_max=sigma_max)
        sizes[t_drop] = s
        durations[t_drop] = T
        if mask is not None:
            ever_toppled |= mask
        t_drop += 1

        if t_drop % SNAPSHOT_EVERY == 0:
            p = float(ever_toppled.mean())
            sigma_hist = np.bincount(state.sigma.ravel())
            mean_sigma_global = float(state.sigma.mean())
            mean_sigma_in_cracked = (
                float(state.sigma[ever_toppled].mean())
                if ever_toppled.any() else 0.0
            )
            snapshots.append({
                "drop": t_drop,
                "p": p,
                "sigma_hist": sigma_hist.tolist(),
                "max_sigma": int(state.sigma.max()),
                "mean_sigma_global": mean_sigma_global,
                "mean_sigma_in_cracked": mean_sigma_in_cracked,
                "grains_lost": state.grains_lost,
                "z_sum": int(state.z.sum()),
                "mean_size_window": float(sizes[max(0, t_drop - SNAPSHOT_EVERY):t_drop].mean()),
                "max_size_window": int(sizes[max(0, t_drop - SNAPSHOT_EVERY):t_drop].max()),
            })

            # Halt early if dynamics has arrested (p stops growing meaningfully)
            if len(snapshots) >= 6:
                recent_ps = [s["p"] for s in snapshots[-6:]]
                # If p has changed less than 0.001 over the last 30k drops, halt
                if max(recent_ps) - min(recent_ps) < 0.001:
                    print(f"  [seed={seed}, {label}] arrest detected at "
                          f"p={p:.4f} after {t_drop} drops")
                    break

    wall = time.time() - t0
    final_p = float(ever_toppled.mean())
    final_max_sigma = int(state.sigma.max())
    print(f"  [seed={seed}, {label}] DONE in {wall:.0f}s, drops={t_drop}, "
          f"final p={final_p:.4f}, max sigma={final_max_sigma}")
    return {
        "seed": seed,
        "sigma_max": sigma_max,
        "label": label,
        "wall_seconds": wall,
        "drops_executed": t_drop,
        "snapshots": snapshots,
        "sizes": sizes[:t_drop],
        "durations": durations[:t_drop],
        "final_p": final_p,
        "final_max_sigma": final_max_sigma,
        "final_sigma_hist": np.bincount(state.sigma.ravel()).tolist(),
        "ever_toppled": ever_toppled,
        "final_state_z": state.z,
        "final_state_sigma": state.sigma,
    }


def main() -> None:
    print(f"M5 multi-state comparison: L={L}, seeds={SEEDS}")
    print(f"  bounded variant: sigma_max=4")
    print(f"  unbounded variant: sigma_max=None")
    print(f"  N_DROPS_MAX = {N_DROPS_MAX:,}")
    print(f"  snapshot every {SNAPSHOT_EVERY:,} drops")
    print("=" * 72)

    t_start = time.time()
    results = []
    for sigma_max in (4, None):
        print(f"\nVariant: sigma_max={sigma_max}")
        for seed in SEEDS:
            r = run_variant(sigma_max, seed)
            results.append(r)

    total = time.time() - t_start
    print(f"\nAll variants done. Total wall: {total:.0f}s")

    # --- Summary table ---
    print()
    print(f"{'seed':>5}  {'variant':>14}  {'drops':>9}  {'final_p':>9}  "
          f"{'max sigma':>10}  {'mean sigma':>11}  {'wall':>7}")
    for r in results:
        snaps = r["snapshots"]
        mean_sigma = snaps[-1]["mean_sigma_global"] if snaps else 0.0
        print(f"{r['seed']:>5}  {r['label']:>14}  {r['drops_executed']:>9d}  "
              f"{r['final_p']:>9.4f}  {r['final_max_sigma']:>10d}  "
              f"{mean_sigma:>11.3f}  {r['wall_seconds']:>7.0f}")

    # --- Plot p(drop) trajectories ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    colors_bounded = ["C0", "C0", "C0"]
    colors_unbounded = ["C3", "C3", "C3"]

    # Panel 1: p(drop)
    ax = axes[0, 0]
    for r in results:
        snaps = r["snapshots"]
        d = [s["drop"] for s in snaps]
        p = [s["p"] for s in snaps]
        color = "C0" if r["sigma_max"] == 4 else "C3"
        ax.plot(d, p, "-", color=color, alpha=0.7,
                label=r["label"] if r["seed"] == SEEDS[0] else None)
    ax.set_xlabel("drops")
    ax.set_ylabel("p (fraction ever-toppled)")
    ax.set_title("p(drop): does evolution arrest?")
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Panel 2: mean sigma in cracked region over time
    ax = axes[0, 1]
    for r in results:
        snaps = r["snapshots"]
        d = [s["drop"] for s in snaps]
        ms = [s["mean_sigma_in_cracked"] for s in snaps]
        color = "C0" if r["sigma_max"] == 4 else "C3"
        ax.plot(d, ms, "-", color=color, alpha=0.7,
                label=r["label"] if r["seed"] == SEEDS[0] else None)
    ax.set_xlabel("drops")
    ax.set_ylabel(r"mean $\sigma$ in cracked region")
    ax.set_title("Tier maturation in cracked region")
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Panel 3: final sigma distribution per variant
    ax = axes[1, 0]
    bounded_dists = [r["final_sigma_hist"] for r in results if r["sigma_max"] == 4]
    unbounded_dists = [r["final_sigma_hist"] for r in results if r["sigma_max"] is None]
    # Average each variant's distributions
    def average_dist(dists):
        max_len = max(len(d) for d in dists)
        padded = np.zeros((len(dists), max_len))
        for i, d in enumerate(dists):
            padded[i, :len(d)] = d
        return padded.mean(axis=0)
    if bounded_dists:
        avg_b = average_dist(bounded_dists)
        ax.bar(np.arange(len(avg_b)) - 0.2, avg_b / avg_b.sum(),
               width=0.4, color="C0", alpha=0.7, label="bounded (sigma_max=4)")
    if unbounded_dists:
        avg_u = average_dist(unbounded_dists)
        ax.bar(np.arange(len(avg_u)) + 0.2, avg_u / avg_u.sum(),
               width=0.4, color="C3", alpha=0.7, label="unbounded")
    ax.set_xlabel("sigma")
    ax.set_ylabel("fraction of cells")
    ax.set_title("Final sigma distribution (seed-averaged)")
    ax.grid(True, alpha=0.3)
    ax.set_yscale("log")
    ax.legend()

    # Panel 4: avalanche size mean over time
    ax = axes[1, 1]
    for r in results:
        snaps = r["snapshots"]
        d = [s["drop"] for s in snaps]
        ms = [s["mean_size_window"] for s in snaps]
        color = "C0" if r["sigma_max"] == 4 else "C3"
        ax.semilogy(d, ms, "-", color=color, alpha=0.7,
                    label=r["label"] if r["seed"] == SEEDS[0] else None)
    ax.set_xlabel("drops")
    ax.set_ylabel("mean avalanche size in window")
    ax.set_title("Avalanche activity vs drops")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    fig.tight_layout()
    outdir = REPO_ROOT / "data" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_multistate_compare_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")
    plt.close(fig)

    # --- Save raw data ---
    npz_path = outdir / f"manna_3d_multistate_data_{stamp}.npz"
    kw: dict = {"L": np.int64(L)}
    for r in results:
        label_key = (
            f"bounded_s{r['seed']}" if r["sigma_max"] == 4
            else f"unbounded_s{r['seed']}"
        )
        snaps = r["snapshots"]
        kw[f"{label_key}_drops"] = np.array([s["drop"] for s in snaps])
        kw[f"{label_key}_p"] = np.array([s["p"] for s in snaps])
        kw[f"{label_key}_mean_sigma_cracked"] = np.array(
            [s["mean_sigma_in_cracked"] for s in snaps]
        )
        kw[f"{label_key}_mean_size_window"] = np.array(
            [s["mean_size_window"] for s in snaps]
        )
        kw[f"{label_key}_max_size_window"] = np.array(
            [s["max_size_window"] for s in snaps]
        )
        kw[f"{label_key}_final_p"] = np.float64(r["final_p"])
        kw[f"{label_key}_final_max_sigma"] = np.int64(r["final_max_sigma"])
        kw[f"{label_key}_final_sigma_hist"] = np.array(r["final_sigma_hist"])
        kw[f"{label_key}_ever_toppled"] = r["ever_toppled"]
        kw[f"{label_key}_final_state_sigma"] = r["final_state_sigma"]
    np.savez_compressed(npz_path, **kw)
    print(f"Saved raw data to {npz_path}")


if __name__ == "__main__":
    main()
