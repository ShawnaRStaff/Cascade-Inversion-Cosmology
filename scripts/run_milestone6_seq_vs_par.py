"""M6 sub-track: sequential vs parallel Manna comparison.

Run BOTH sequential and parallel implementations at L=24 with 5 seeds
each, fit tau_a for each. The question:
  - Does sequential match literature tau_a = 1.442 (H&P 2012)?
  - Does parallel differ?
  - If both differ from literature, our fit method or finite-L effects are
    confounding the comparison.

L=24 chosen for speed (sequential is slow at larger L). 50k drops per
run is enough to reach saturation and gather statistics.
"""

from __future__ import annotations

import json
import sys
import time
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.sandpile_3d import (  # noqa: E402
    drive,
    initialize,
    relax,
    relax_sequential,
)

L = 24
N_DROPS = 50_000
N_SEEDS = 5
TAU_A_LITERATURE = 1.442


def run_one(L: int, seed: int, n_drops: int, sequential: bool) -> dict:
    rng = np.random.default_rng(seed)
    state = initialize(L)
    ever = np.zeros((L, L, L), dtype=bool)
    unique = np.zeros(n_drops, dtype=np.int64)
    sizes = np.zeros(n_drops, dtype=np.int64)
    relax_fn = relax_sequential if sequential else relax
    t0 = time.time()
    for t in range(n_drops):
        drive(state, rng)
        s, T, mask = relax_fn(state, rng, track_support=True)
        sizes[t] = s
        if mask is not None:
            unique[t] = int(mask.sum())
            ever |= mask
    return {
        "seed": seed,
        "wall_s": time.time() - t0,
        "z_avg": float(state.z.mean()),
        "p_final": float(ever.mean()),
        "unique": unique,
        "sizes": sizes,
    }


def fit_tau_a(unique_sizes: np.ndarray, sat_start: int = None) -> dict:
    if sat_start is None:
        sat_start = len(unique_sizes) // 2
    sat = unique_sizes[sat_start:]
    nonzero = sat[sat > 0]
    if len(nonzero) < 100:
        return {"error": "insufficient events"}
    a_min, a_max = max(2, int(nonzero.min())), int(nonzero.max())
    if a_max <= a_min:
        return {"error": "no spread"}
    bins = np.logspace(np.log10(a_min), np.log10(a_max + 1), 30)
    counts, edges = np.histogram(nonzero, bins=bins)
    widths = np.diff(edges)
    centers = np.sqrt(edges[:-1] * edges[1:])
    pdf = counts / (nonzero.size * widths)
    valid = (counts > 0) & np.isfinite(pdf)
    if valid.sum() < 8:
        return {"error": "too few bins"}
    fit_idx = np.where(valid)[0]
    mask = valid.copy()
    if len(fit_idx) > 5:
        mask[fit_idx[:2]] = False
        mask[fit_idx[-3:]] = False
    if mask.sum() < 4:
        return {"error": "narrow scaling regime"}
    log_x = np.log10(centers[mask])
    log_y = np.log10(pdf[mask])
    slope, _ = np.polyfit(log_x, log_y, 1)
    return {"tau_a": float(-slope), "n_events": int(nonzero.size), "n_bins": int(mask.sum())}


def main():
    out_dir = REPO_ROOT / "data" / "outputs" / f"seq_vs_par_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== Parallel vs Sequential Manna at L={L}, N_drops={N_DROPS}, n_seeds={N_SEEDS} ===")
    print(f"Literature value: tau_a = {TAU_A_LITERATURE}")
    print()

    all_results = {"parallel": [], "sequential": []}
    for mode in ["parallel", "sequential"]:
        print(f"--- {mode.upper()} ---")
        for s in range(N_SEEDS):
            seed = 1000 + s
            res = run_one(L, seed, N_DROPS, sequential=(mode == "sequential"))
            fit = fit_tau_a(res["unique"])
            max_unique_pct = res["unique"].max() / (L ** 3) * 100
            res["fit"] = fit
            res["max_unique_pct"] = max_unique_pct
            all_results[mode].append(res)
            tau_str = f"{fit['tau_a']:.3f}" if "tau_a" in fit else f"ERROR: {fit.get('error')}"
            print(f"  s={seed}: wall={res['wall_s']:>5.0f}s z_avg={res['z_avg']:.4f} "
                  f"p={res['p_final']:.3f} max_unique%={max_unique_pct:.2f} tau_a={tau_str}")

    # Aggregate
    summary = {}
    for mode, runs in all_results.items():
        taus = [r["fit"]["tau_a"] for r in runs if "tau_a" in r["fit"]]
        maxs = [r["max_unique_pct"] for r in runs]
        zs = [r["z_avg"] for r in runs]
        summary[mode] = {
            "tau_a_mean": float(np.mean(taus)),
            "tau_a_std": float(np.std(taus)),
            "max_unique_pct_mean": float(np.mean(maxs)),
            "z_avg_mean": float(np.mean(zs)),
        }

    print()
    print(f"=== Summary ===")
    print(f"  {'mode':>12}  {'tau_a':>14}  {'max_unique%':>12}  {'z_avg':>8}")
    for mode, s in summary.items():
        print(f"  {mode:>12}  {s['tau_a_mean']:.3f} +/- {s['tau_a_std']:.3f}  "
              f"{s['max_unique_pct_mean']:>11.2f}  {s['z_avg_mean']:>8.4f}")
    print(f"  literature      {TAU_A_LITERATURE:.3f} +/- 0.012      —             0.6223")

    # Verdict
    par_tau, par_std = summary["parallel"]["tau_a_mean"], summary["parallel"]["tau_a_std"]
    seq_tau, seq_std = summary["sequential"]["tau_a_mean"], summary["sequential"]["tau_a_std"]
    par_vs_lit = abs(par_tau - TAU_A_LITERATURE) / max(par_std, 0.001)
    seq_vs_lit = abs(seq_tau - TAU_A_LITERATURE) / max(seq_std, 0.001)
    print()
    print(f"  Parallel vs literature:  {par_vs_lit:.1f} sigma gap")
    print(f"  Sequential vs literature: {seq_vs_lit:.1f} sigma gap")
    if seq_vs_lit < 2:
        print("  -> Sequential MATCHES literature; parallel-vs-sequential is the source of our 1.49 vs 1.44 discrepancy")
    elif seq_vs_lit < par_vs_lit:
        print("  -> Sequential is CLOSER to literature than parallel, but neither perfect")
    else:
        print("  -> Both differ from literature similarly; fit method or finite-L effects dominate")

    # Save
    json_out = {
        "L": L,
        "n_drops": N_DROPS,
        "n_seeds": N_SEEDS,
        "literature_tau_a": TAU_A_LITERATURE,
        "summary": summary,
    }
    with open(out_dir / "results.json", "w") as f:
        json.dump(json_out, f, indent=2)
    print(f"\nSaved: {out_dir / 'results.json'}")


if __name__ == "__main__":
    main()
