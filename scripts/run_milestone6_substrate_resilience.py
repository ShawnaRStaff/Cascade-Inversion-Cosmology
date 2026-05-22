"""M6 sub-track: substrate resilience experiment.

Question (from user, 2026-05-21):
    Does the substrate require 100% fracture to structurally collapse,
    or does a smaller fraction suffice? What is the role of the
    'untouched' cells in holding the substrate together — anything,
    or nothing?

Method:
    Take a real saturated substrate state from an existing run.
    Artificially 'shatter' a configurable fraction of cells (set to
    z=0, simulating a single big event having just occurred). Then
    continue normal dynamics in two regimes:
      A. WITH drops (input continues — what our normal sim does)
      B. WITHOUT drops (input stops — the user's question scenario)
    Measure recovery / collapse over time.

Fractions tested: 0% (control), 30%, 50%, 70%, 91%, 99%
Source state: L=64 final_z from new-fleet sweep (saturated, z=0.616).

For each (fraction, with/without drops):
    - carrier fraction over time
    - event size over time (with drops only — no drops = no events)
    - time to return to z_avg = 0.616 (with drops only)
    - whether dynamics 'goes silent' (without drops)

This experiment directly tests the structural-collapse intuition:
    if the substrate can survive 91% fracture (with input), what does
    it actually take to break it?
"""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.sandpile_3d import drive, initialize, MannaState3D, relax  # noqa: E402

SHATTER_FRACTIONS = [0.0, 0.30, 0.50, 0.70, 0.91, 0.99]
RECOVERY_DROPS = 5_000
WITH_DROPS_LABEL = "drops_on"
WITHOUT_DROPS_LABEL = "drops_off"


def load_saturated_state(L: int = 48, seed: int = 14800) -> MannaState3D:
    """Load the final_z from an existing L=48 run as starting point.

    Using L=48 (not L=64) keeps each scenario tractable for local
    compute — 12 scenarios x 5k drops at L=48 should finish in ~5 min.
    """
    sweep_dir = REPO_ROOT / "data" / "outputs" / "fss_sweep_20260521_031056"
    path = sweep_dir / f"L{L}_s{seed}_final.npz"
    data = np.load(path, allow_pickle=True)
    state = MannaState3D(
        z=np.asarray(data["final_z"]).astype(np.int64).copy(),
        grains_lost=int(data["grains_lost"]),
    )
    print(f"Loaded saturated state: L={L} seed={seed}, z_avg={state.z.mean():.4f}")
    return state


def shatter(state: MannaState3D, fraction: float, rng: np.random.Generator) -> int:
    """Set a random fraction of cells to z=0 (simulating a big event).

    Returns the number of grains 'lost' to the shatter (subtracted from
    state since they're discarded).
    """
    L = state.z.shape[0]
    n_total = L ** 3
    n_to_shatter = int(fraction * n_total)
    if n_to_shatter == 0:
        return 0
    flat = state.z.ravel()
    idx = rng.choice(n_total, size=n_to_shatter, replace=False)
    grains_removed = int(flat[idx].sum())
    flat[idx] = 0
    return grains_removed


def carrier_fraction(state: MannaState3D) -> float:
    return float((state.z >= 1).mean())


def run_recovery(
    state: MannaState3D,
    drops_on: bool,
    n_drops: int,
    seed: int,
) -> dict:
    """Continue dynamics from current state for n_drops. Record metrics.

    If drops_on=False, no new grains added — we just let the system
    relax and see if any spontaneous topples happen.
    """
    rng = np.random.default_rng(seed)
    z_avg_trace = []
    carrier_trace = []
    size_trace = []
    sample_every = max(1, n_drops // 200)

    if drops_on:
        for t in range(n_drops):
            drive(state, rng)
            s, _, _ = relax(state, rng, track_support=False)
            if t % sample_every == 0:
                z_avg_trace.append(float(state.z.mean()))
                carrier_trace.append(carrier_fraction(state))
                size_trace.append(int(s))
    else:
        # No drops. Just relax any still-unstable cells once.
        s, _, _ = relax(state, rng, track_support=False)
        z_avg_trace.append(float(state.z.mean()))
        carrier_trace.append(carrier_fraction(state))
        size_trace.append(int(s))
        # Should stay frozen — record 1 more measurement to confirm.
        z_avg_trace.append(float(state.z.mean()))
        carrier_trace.append(carrier_fraction(state))
        size_trace.append(0)

    return {
        "drops_on": drops_on,
        "n_drops_run": n_drops if drops_on else 0,
        "z_avg_trace": z_avg_trace,
        "carrier_trace": carrier_trace,
        "size_trace": size_trace,
        "z_avg_final": float(state.z.mean()),
        "carrier_final": carrier_fraction(state),
    }


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"substrate_resilience_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== Substrate resilience experiment ===")
    print(f"Output dir: {out_dir}")
    print(f"Test fractions: {SHATTER_FRACTIONS}")
    print(f"Recovery: {RECOVERY_DROPS} drops per condition")
    print()

    results = []
    for fraction in SHATTER_FRACTIONS:
        # Shatter scenario A: WITH drops
        state = load_saturated_state()
        rng = np.random.default_rng(42 + int(fraction * 100))
        grains_removed = shatter(state, fraction, rng)
        carrier_post_shatter = carrier_fraction(state)
        z_avg_post_shatter = float(state.z.mean())
        print(f"  [shatter={fraction:.0%}] grains_removed={grains_removed:>8d} "
              f"post-shatter: z_avg={z_avg_post_shatter:.4f} carriers={carrier_post_shatter:.3f}")
        r_on = run_recovery(state, drops_on=True, n_drops=RECOVERY_DROPS, seed=100 + int(fraction * 100))
        print(f"    [drops_on] after {RECOVERY_DROPS} drops: z_avg={r_on['z_avg_final']:.4f} carriers={r_on['carrier_final']:.3f}")

        # Shatter scenario B: WITHOUT drops
        state = load_saturated_state()
        rng = np.random.default_rng(42 + int(fraction * 100))
        shatter(state, fraction, rng)
        r_off = run_recovery(state, drops_on=False, n_drops=0, seed=100 + int(fraction * 100))
        spontaneous_topples = r_off["size_trace"][0] if r_off["size_trace"] else 0
        print(f"    [drops_off] spontaneous_topples={spontaneous_topples} "
              f"final z_avg={r_off['z_avg_final']:.4f} (started at {z_avg_post_shatter:.4f})")

        results.append({
            "fraction_shattered": fraction,
            "grains_removed": grains_removed,
            "z_avg_post_shatter": z_avg_post_shatter,
            "carrier_post_shatter": carrier_post_shatter,
            "with_drops": {
                "z_avg_trace": r_on["z_avg_trace"],
                "carrier_trace": r_on["carrier_trace"],
                "size_trace": r_on["size_trace"],
                "z_avg_final": r_on["z_avg_final"],
                "carrier_final": r_on["carrier_final"],
            },
            "without_drops": {
                "spontaneous_topples_first_relax": int(spontaneous_topples),
                "z_avg_final": r_off["z_avg_final"],
                "carrier_final": r_off["carrier_final"],
            },
        })

    # Save raw results
    with open(out_dir / "results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults: {out_dir}/results.json")

    # Plot recovery curves (carrier fraction with drops on)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    for r in results:
        trace = r["with_drops"]["carrier_trace"]
        x = np.linspace(0, RECOVERY_DROPS, len(trace))
        ax.plot(x, trace, label=f"shattered {r['fraction_shattered']:.0%}")
    ax.axhline(0.62, color="k", linestyle=":", alpha=0.5, label="equilibrium ~62%")
    ax.axhline(0.31, color="r", linestyle=":", alpha=0.5, label="percolation threshold ~31%")
    ax.set_xlabel("drops after shatter")
    ax.set_ylabel("carrier fraction (z>=1)")
    ax.set_title("Substrate recovery WITH continued input")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    fractions = [r["fraction_shattered"] for r in results]
    carrier_post = [r["carrier_post_shatter"] for r in results]
    carrier_with_drops = [r["with_drops"]["carrier_final"] for r in results]
    carrier_without_drops = [r["without_drops"]["carrier_final"] for r in results]
    ax.plot(fractions, carrier_post, "o-", label="immediately after shatter")
    ax.plot(fractions, carrier_with_drops, "o-", label=f"after {RECOVERY_DROPS} drops (with input)")
    ax.plot(fractions, carrier_without_drops, "x--", label="without input (spontaneous only)")
    ax.axhline(0.31, color="r", linestyle=":", alpha=0.5, label="percolation threshold")
    ax.set_xlabel("fraction of substrate shattered")
    ax.set_ylabel("carrier fraction")
    ax.set_title("Recovery vs damage")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_dir / "resilience.png", dpi=150)
    print(f"Plot: {out_dir}/resilience.png")


if __name__ == "__main__":
    main()
