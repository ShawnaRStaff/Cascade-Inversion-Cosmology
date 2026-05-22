"""M6 sub-track: aftershock / temporal-correlation analysis.

Tests whether large avalanches cluster in time (earthquake aftershock
pattern), or whether they're temporally independent (Poisson).

Mechanism in question: after a big avalanche fires, the cells along
its boundary received the released energy and may now be primed
closer to threshold. If so, the next few drops have elevated
probability of triggering another large avalanche. Result: large
events cluster in time, NOT distributed as Poisson.

Cosmologically: this would mean the "Big Bang" isn't a single
instant but a cascade of correlated catastrophic events that
appears as one phenomenon at coarse-graining.

Analyses run on every final.npz in the NEW sweep dir (L=48, 64, 96):
  1. Per-run: autocorrelation function of sizes[t] in post-saturation
     window. Positive autocorrelation at small lags = aftershock
     signature. Zero autocorrelation = Poisson (independent events).
  2. Per-run: distribution of inter-large-event times. Define
     "large" as size > 99th percentile of post-saturation sizes.
     Fit exponential (Poisson null) vs power-law (clustered).
  3. Aggregate: mean autocorrelation and inter-event distribution
     across all seeds at each L.

Output: docs/notes/milestone6_aftershock.md + plots in
data/outputs/aftershock_analysis_TIMESTAMP/
"""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
SWEEP_DIR = REPO_ROOT / "data" / "outputs" / "fss_sweep_20260521_031056"


def load_run(path: Path) -> dict:
    """Load a final.npz and return its key fields."""
    d = np.load(path, allow_pickle=True)
    return {
        "L": int(d["L"]),
        "seed": int(d["seed"]),
        "sizes": np.asarray(d["sizes"]),
        "unique_sizes": (
            np.asarray(d["unique_sizes"]) if "unique_sizes" in d.files else None
        ),
        "snapshots": list(d["snapshots"]) if "snapshots" in d.files else [],
    }


def find_saturation_start(snapshots: list[dict]) -> int:
    """Return the drop index at which p first reached 0.95.

    Post-saturation analysis uses sizes from this drop onward.
    """
    for s in snapshots:
        if s.get("p", 0) >= 0.95:
            return int(s["drop"])
    # Fallback: use the second half of the run
    return -1


def autocorrelation(x: np.ndarray, max_lag: int = 50) -> np.ndarray:
    """Compute normalized autocorrelation up to max_lag."""
    x = x - x.mean()
    var = x.var()
    if var == 0:
        return np.zeros(max_lag + 1)
    n = len(x)
    acf = np.zeros(max_lag + 1)
    for lag in range(max_lag + 1):
        if lag == 0:
            acf[lag] = 1.0
        else:
            acf[lag] = (x[:-lag] * x[lag:]).mean() / var
    return acf


def inter_large_event_times(sizes: np.ndarray, percentile: float = 99.0) -> np.ndarray:
    """Compute time between consecutive 'large' events.

    Large = size > percentile-th percentile of nonzero sizes.
    Returns array of inter-event drop counts.
    """
    nonzero = sizes[sizes > 0]
    if len(nonzero) < 10:
        return np.zeros(0, dtype=np.int64)
    threshold = np.percentile(nonzero, percentile)
    large_idx = np.where(sizes >= threshold)[0]
    if len(large_idx) < 2:
        return np.zeros(0, dtype=np.int64)
    return np.diff(large_idx)


def fit_exponential_vs_powerlaw(intertimes: np.ndarray) -> dict:
    """Compare exponential (Poisson null) vs power-law fit.

    Returns dict with both fits and a likelihood-ratio decision.
    Power law: P(t) ~ t^(-alpha)
    Exponential: P(t) ~ exp(-t/tau)
    """
    if len(intertimes) < 10:
        return {"n": len(intertimes), "verdict": "insufficient data"}

    t = intertimes[intertimes >= 1].astype(np.float64)
    if len(t) < 10:
        return {"n": len(t), "verdict": "insufficient data after filtering"}

    # Exponential fit: MLE tau = mean(t)
    tau = t.mean()
    # Log-likelihood for exponential
    log_l_exp = -len(t) * np.log(tau) - t.sum() / tau

    # Power-law fit (Hill estimator):
    # MLE for alpha given t_min=1: alpha = 1 + N / sum(log(t))
    log_t = np.log(t)
    alpha = 1.0 + len(t) / log_t.sum()
    # Log-likelihood for power law (with t_min=1, finite domain)
    log_l_pow = len(t) * np.log(alpha - 1) - alpha * log_t.sum()

    # Vuong-style comparison
    delta = log_l_pow - log_l_exp

    if delta > 5:
        verdict = "power-law strongly favored (clustering / aftershock pattern)"
    elif delta > 2:
        verdict = "power-law mildly favored"
    elif delta < -5:
        verdict = "exponential strongly favored (Poisson / independent events)"
    elif delta < -2:
        verdict = "exponential mildly favored"
    else:
        verdict = "inconclusive (both fits comparable)"

    return {
        "n": int(len(t)),
        "exp_tau": float(tau),
        "pow_alpha": float(alpha),
        "log_l_exp": float(log_l_exp),
        "log_l_pow": float(log_l_pow),
        "delta_log_l": float(delta),
        "verdict": verdict,
    }


def analyze_run(run: dict, max_lag: int = 50, percentile: float = 99.0) -> dict:
    """Run all analyses on one (L, seed) trajectory."""
    sizes = run["sizes"]
    sat_start = find_saturation_start(run["snapshots"])
    if sat_start < 0:
        sat_start = len(sizes) // 2
    sizes_sat = sizes[sat_start:]
    if len(sizes_sat) < 100:
        return {"error": "saturation window too short"}

    acf = autocorrelation(sizes_sat, max_lag=max_lag)
    inter = inter_large_event_times(sizes_sat, percentile=percentile)
    fit = fit_exponential_vs_powerlaw(inter)

    # Excess autocorrelation at small lag (lags 1-20)
    excess_acf_small_lag = float(acf[1:21].mean())

    return {
        "L": run["L"],
        "seed": run["seed"],
        "sat_start_drop": sat_start,
        "sat_window_length": len(sizes_sat),
        "n_large_events": len(inter) + 1 if len(inter) > 0 else 0,
        "acf": acf.tolist(),
        "excess_acf_lag1_to_20": excess_acf_small_lag,
        "inter_event_fit": fit,
    }


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"aftershock_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for L in [48, 64, 96, 128]:
        files = sorted(SWEEP_DIR.glob(f"L{L}_*_final.npz"))
        if not files:
            print(f"L={L}: no data yet")
            continue
        print(f"\n=== L={L} ({len(files)} seeds) ===")
        L_results = []
        for f in files:
            run = load_run(f)
            r = analyze_run(run)
            L_results.append(r)
            results.append(r)
            if "error" not in r:
                fit = r["inter_event_fit"]
                print(
                    f"  s={r['seed']}: "
                    f"sat_window={r['sat_window_length']} drops, "
                    f"n_large={r['n_large_events']}, "
                    f"acf_lag1-20={r['excess_acf_lag1_to_20']:+.4f}, "
                    f"verdict={fit.get('verdict', 'NA')}"
                )

        # Aggregate per-L
        excess_acfs = [r["excess_acf_lag1_to_20"] for r in L_results if "error" not in r]
        if excess_acfs:
            mean_excess = float(np.mean(excess_acfs))
            print(f"  L={L} mean excess ACF (lag 1-20) = {mean_excess:+.4f}")
            if mean_excess > 0.05:
                print(f"  -> SUSTAINED POSITIVE AUTOCORRELATION (aftershock signature)")
            elif mean_excess > 0.01:
                print(f"  -> mild positive autocorrelation")
            elif abs(mean_excess) < 0.01:
                print(f"  -> consistent with independent events (Poisson)")
            else:
                print(f"  -> negative autocorrelation (anti-clustering)")

    # Save JSON
    json_path = out_dir / "results.json"
    with open(json_path, "w") as f:
        # Drop the full acf arrays for json brevity; save summary
        summary = [
            {k: v for k, v in r.items() if k != "acf"}
            for r in results
        ]
        json.dump(summary, f, indent=2)
    print(f"\nResults JSON: {json_path}")

    # Save full ACF arrays as separate npz
    npz_path = out_dir / "acf_arrays.npz"
    acf_kwargs = {}
    for r in results:
        if "error" not in r:
            acf_kwargs[f"L{r['L']}_s{r['seed']}_acf"] = np.array(r["acf"])
    np.savez_compressed(npz_path, **acf_kwargs)
    print(f"Full ACF arrays: {npz_path}")


if __name__ == "__main__":
    main()
