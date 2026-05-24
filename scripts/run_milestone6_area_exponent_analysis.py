"""M6 sub-track: area exponent fit + asymptote extrapolation.

Two analyses on the existing FSS sweep data (NO new compute):

1. Fit the avalanche AREA exponent tau_a from the empirical
   distribution P(a) ~ a^(-tau_a). Compare to the published value
   for 3D Abelian Manna (sequential): tau_a = 1.442 (Huynh & Pruessner
   2012). If our parallel-update implementation gives the same
   exponent, we're in the same universality class. If different,
   that's a real finding about parallel vs sequential dynamics.

2. Extrapolate the L -> inf asymptote of unique%(L) using the
   relation area_max ~ A * L^3 (since D_a ~ d = 3 in the literature).
   Fit A from our L = 48, 64, 96, 128 data without needing L=256.
   If A is close to the carrier fraction 0.62, that supports the
   "cascades touch all carriers but no sinks" reading.

Inputs:  data/outputs/fss_sweep_20260521_031056/L*_final.npz
Outputs: data/outputs/area_exponent_analysis_<stamp>/
         - tau_a_fit.png
         - asymptote_fit.png
         - results.json
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
SWEEP_DIR = REPO_ROOT / "data" / "outputs" / "fss_sweep_20260521_031056"

# Published canonical value for 3D Abelian Manna (sequential), Huynh & Pruessner 2012
TAU_A_LITERATURE = 1.442
TAU_A_LIT_ERR = 0.012


def load_run(path: Path) -> dict:
    d = np.load(path, allow_pickle=True)
    return {
        "L": int(d["L"]),
        "seed": int(d["seed"]),
        "unique_sizes": np.asarray(d["unique_sizes"]),
        "snapshots": list(d["snapshots"]) if "snapshots" in d.files else [],
    }


def find_saturation_start(snapshots: list[dict]) -> int:
    for s in snapshots:
        if s.get("p", 0) >= 0.95:
            return int(s["drop"])
    return -1


def fit_power_law_log_binned(
    sizes: np.ndarray, n_bins: int = 30
) -> dict:
    """Fit P(a) ~ a^(-tau_a) via log-binned histogram + log-log linear fit.

    Excludes the rare-tail cutoff (top few bins) and the small-event
    rollover (bottom few bins) to fit only the scaling regime.
    """
    nonzero = sizes[sizes > 0]
    if len(nonzero) < 100:
        return {"error": "insufficient nonzero events", "n": len(nonzero)}

    a_min = max(2, int(nonzero.min()))
    a_max = int(nonzero.max())
    if a_max <= a_min:
        return {"error": "no spread in events"}

    bins = np.logspace(np.log10(a_min), np.log10(a_max + 1), n_bins)
    counts, edges = np.histogram(nonzero, bins=bins)
    widths = np.diff(edges)
    centers = np.sqrt(edges[:-1] * edges[1:])
    pdf = counts / (nonzero.size * widths)
    valid = (counts > 0) & np.isfinite(pdf)

    # Fit the scaling regime: skip bottom 2 and top 3 bins (rollover/cutoff)
    if valid.sum() < 8:
        return {"error": "too few valid bins", "n_bins_valid": int(valid.sum())}

    fit_mask = valid.copy()
    fit_idx = np.where(fit_mask)[0]
    if len(fit_idx) > 5:
        fit_mask[fit_idx[:2]] = False
        fit_mask[fit_idx[-3:]] = False

    if fit_mask.sum() < 4:
        return {"error": "scaling regime too narrow"}

    log_x = np.log10(centers[fit_mask])
    log_y = np.log10(pdf[fit_mask])
    slope, intercept = np.polyfit(log_x, log_y, 1)
    tau_fit = -slope

    # Bootstrap error
    n_boot = 100
    boot_taus = []
    rng = np.random.default_rng(42)
    for _ in range(n_boot):
        boot = rng.choice(nonzero, size=len(nonzero), replace=True)
        b_counts, _ = np.histogram(boot, bins=bins)
        b_pdf = b_counts / (boot.size * widths)
        bv = (b_counts > 0) & np.isfinite(b_pdf)
        bm = bv.copy()
        bi = np.where(bm)[0]
        if len(bi) > 5:
            bm[bi[:2]] = False
            bm[bi[-3:]] = False
        if bm.sum() < 4:
            continue
        bs, _ = np.polyfit(np.log10(centers[bm]), np.log10(b_pdf[bm]), 1)
        boot_taus.append(-bs)
    tau_err = float(np.std(boot_taus)) if boot_taus else 0.0

    return {
        "tau_a": float(tau_fit),
        "tau_a_err": tau_err,
        "n_events": int(nonzero.size),
        "n_bins_used": int(fit_mask.sum()),
        "a_min_fit": float(centers[fit_mask].min()),
        "a_max_fit": float(centers[fit_mask].max()),
        "bin_centers": centers.tolist(),
        "pdf": pdf.tolist(),
        "counts": counts.tolist(),
    }


def analyze_run(run: dict) -> dict:
    sizes = run["unique_sizes"]
    sat_start = find_saturation_start(run["snapshots"])
    if sat_start < 0:
        sat_start = len(sizes) // 2
    sizes_sat = sizes[sat_start:]
    fit = fit_power_law_log_binned(sizes_sat)
    max_area = int(sizes_sat.max()) if len(sizes_sat) > 0 else 0
    return {
        "L": run["L"],
        "seed": run["seed"],
        "sat_start_drop": int(sat_start),
        "max_area": max_area,
        "max_area_pct": float(max_area) / (run["L"] ** 3) * 100,
        "tau_a_fit": fit,
    }


def fit_asymptote(results: list[dict]) -> dict:
    """Fit area_max(L) ~ A * L^3 + corrections, extract A."""
    Ls = sorted({r["L"] for r in results})
    pairs = []
    for L in Ls:
        maxs = [r["max_area"] for r in results if r["L"] == L]
        if maxs:
            pairs.append((L, float(np.mean(maxs)), float(np.std(maxs))))
    Ls_arr = np.array([p[0] for p in pairs])
    maxs_arr = np.array([p[1] for p in pairs])
    A_per_L = maxs_arr / Ls_arr ** 3

    # Linear fit A(L) = A_inf + c/L (assume leading finite-size correction)
    x = 1.0 / Ls_arr
    if len(x) >= 2:
        slope, intercept = np.polyfit(x, A_per_L, 1)
        A_inf = float(intercept)
        c = float(slope)
    else:
        A_inf = float(A_per_L[-1])
        c = 0.0
    return {
        "L_values": Ls_arr.tolist(),
        "max_area_per_L3": A_per_L.tolist(),
        "asymptote_A_inf": A_inf,
        "finite_size_coef_c": c,
        "model": "A(L) = A_inf + c/L",
        "carrier_fraction_reference": 0.62,
    }


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"area_exponent_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)

    results = []
    print(f"=== Per-run tau_a fits ===")
    print(f"  (literature value: tau_a = {TAU_A_LITERATURE} +/- {TAU_A_LIT_ERR}, Huynh & Pruessner 2012)")
    print()
    for L in [48, 64, 96, 128]:
        files = sorted(SWEEP_DIR.glob(f"L{L}_*_final.npz"))
        if not files:
            continue
        print(f"  L={L}:")
        for f in files:
            run = load_run(f)
            r = analyze_run(run)
            results.append(r)
            fit = r["tau_a_fit"]
            if "error" in fit:
                print(f"    s={r['seed']}: ERROR ({fit['error']})")
            else:
                print(
                    f"    s={r['seed']}: tau_a = {fit['tau_a']:.3f} +/- {fit['tau_a_err']:.3f} "
                    f"(n_events={fit['n_events']}, n_bins={fit['n_bins_used']}) "
                    f"max_area_pct={r['max_area_pct']:.2f}%"
                )

    # Aggregate tau_a across all runs
    all_taus = [r["tau_a_fit"]["tau_a"] for r in results if "tau_a" in r["tau_a_fit"]]
    print()
    print(f"=== Aggregate tau_a across all 17 runs ===")
    print(f"  mean = {np.mean(all_taus):.3f}")
    print(f"  std  = {np.std(all_taus):.3f}")
    print(f"  min  = {min(all_taus):.3f}")
    print(f"  max  = {max(all_taus):.3f}")
    print()
    print(f"  Literature (Huynh & Pruessner 2012 sequential Manna 3D SC):")
    print(f"  tau_a = {TAU_A_LITERATURE} +/- {TAU_A_LIT_ERR}")
    print()
    sigma_gap = abs(np.mean(all_taus) - TAU_A_LITERATURE) / max(np.std(all_taus), 1e-6)
    if sigma_gap < 1:
        verdict = "AGREES with literature within ~1 sigma -> same universality class"
    elif sigma_gap < 2:
        verdict = "marginal agreement (1-2 sigma)"
    else:
        verdict = f"DIFFERS from literature by {sigma_gap:.1f} sigma -> parallel-update dynamics may be in different class"
    print(f"  Verdict: {verdict}")

    asymptote = fit_asymptote(results)
    print()
    print(f"=== Asymptote fit: area_max(L)/L^3 = A_inf + c/L ===")
    print(f"  L values: {asymptote['L_values']}")
    print(f"  area_max/L^3: {[f'{v:.4f}' for v in asymptote['max_area_per_L3']]}")
    print(f"  A_inf (extrapolated L -> infinity) = {asymptote['asymptote_A_inf']:.4f}")
    print(f"  finite-size coefficient c = {asymptote['finite_size_coef_c']:.4f}")
    print(f"  reference carrier fraction = 0.620")
    if 0.55 < asymptote['asymptote_A_inf'] < 0.70:
        print(f"  -> A_inf is CLOSE to carrier fraction; supports 'cascades touch all carriers' hypothesis")
    elif asymptote['asymptote_A_inf'] > 0.85:
        print(f"  -> A_inf approaches 1; cascades are essentially space-filling at large L")
    else:
        print(f"  -> A_inf is some other value; need physical interpretation")

    out = {
        "literature_tau_a": TAU_A_LITERATURE,
        "literature_tau_a_err": TAU_A_LIT_ERR,
        "our_tau_a_mean": float(np.mean(all_taus)),
        "our_tau_a_std": float(np.std(all_taus)),
        "sigma_gap_to_literature": float(sigma_gap),
        "verdict": verdict,
        "asymptote": asymptote,
        "per_run": [{k: v for k, v in r.items() if k != "tau_a_fit"} | {"tau_a": r["tau_a_fit"].get("tau_a"), "tau_a_err": r["tau_a_fit"].get("tau_a_err")} for r in results],
    }
    with open(out_dir / "results.json", "w") as f:
        json.dump(out, f, indent=2)
    print()
    print(f"Saved: {out_dir / 'results.json'}")

    # Plot 1: tau_a distributions
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax = axes[0]
    Ls = sorted({r["L"] for r in results})
    for L in Ls:
        taus_L = [r["tau_a_fit"]["tau_a"] for r in results if r["L"] == L and "tau_a" in r["tau_a_fit"]]
        if taus_L:
            ax.scatter([L] * len(taus_L), taus_L, s=80, alpha=0.7, label=f"L={L} (n={len(taus_L)})")
    ax.axhline(TAU_A_LITERATURE, color="r", linestyle="--", label=f"literature (H&P 2012): {TAU_A_LITERATURE}")
    ax.fill_between([min(Ls) - 5, max(Ls) + 5], TAU_A_LITERATURE - TAU_A_LIT_ERR,
                    TAU_A_LITERATURE + TAU_A_LIT_ERR, alpha=0.2, color="r")
    ax.set_xlabel("L")
    ax.set_ylabel(r"fitted $\tau_a$")
    ax.set_title("Area-distribution exponent per run vs literature")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Plot 2: asymptote fit
    ax = axes[1]
    x = 1.0 / np.array(asymptote["L_values"])
    y = np.array(asymptote["max_area_per_L3"])
    ax.plot(x, y, "o", markersize=10, label="data (mean per L)")
    x_fit = np.linspace(0, x.max() * 1.1, 100)
    y_fit = asymptote["asymptote_A_inf"] + asymptote["finite_size_coef_c"] * x_fit
    ax.plot(x_fit, y_fit, "-", alpha=0.6, label=f"fit: A_inf={asymptote['asymptote_A_inf']:.3f}, c={asymptote['finite_size_coef_c']:.3f}")
    ax.axhline(0.62, color="g", linestyle=":", alpha=0.6, label="carrier fraction (~62%)")
    ax.axhline(1.0, color="k", linestyle=":", alpha=0.3, label="full lattice")
    ax.set_xlabel("1/L")
    ax.set_ylabel(r"max area / $L^3$")
    ax.set_title("Asymptote extrapolation: L -> infinity")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_dir / "tau_a_and_asymptote.png", dpi=150)
    print(f"Saved: {out_dir / 'tau_a_and_asymptote.png'}")


if __name__ == "__main__":
    main()
