"""Finite-size scaling of the dynamics-driven percolation threshold p_c(L).

For each L, run many independent seeds of the 3D Manna sandpile,
each driven from an empty lattice until the cumulative ever-toppled set
first spans the box. The recorded quantity is

    p(L, seed) = fractured_fraction at the moment of first spanning

so the run is a Monte Carlo estimate of the per-realization threshold.
The mean over seeds gives p_c(L); the standard error tells us how
reliable that mean is.

The FSS extrapolation we want is

    p_c(L) = p_c_inf + a * L^(-1/nu)

with nu the correlation-length exponent for the percolation transition.
Fitting log|p_c(L) - p_c_inf| vs log(L) needs us to know p_c_inf first,
so the standard trick is: try several p_c_inf and pick the one that
linearizes the relation, or fit jointly with scipy.optimize.curve_fit.
We do the joint fit.

Reference points (NOT predictions of our model, just for orientation):
- Uncorrelated site percolation on simple cubic: p_c ~ 0.3116
- We expect a LOWER value here because toppling produces spatially
  correlated fractured sets.

Run:
    .venv/bin/python -u scripts/run_milestone3_pc_fss.py
"""

from __future__ import annotations

import sys
import time
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.percolation import check_spanning, fractured_fraction  # noqa: E402
from void_cascade.sandpile_3d import run_with_ever_toppled  # noqa: E402


# (L, n_seeds, n_drops_cap, check_every)
# n_drops_cap is generous; the run stops at first span. check_every sets
# the granularity with which we detect spanning - finer = less upward
# bias on p but more overhead.
RUN_PLAN = [
    (16, 12,   5_000,  10),
    (24, 12,  15_000,  20),
    (32, 10,  30_000,  30),
    (48,  8,  80_000,  50),
    (64,  6, 200_000,  80),
    (96,  4, 500_000, 120),
]


def find_first_spanning(
    L: int, seed: int, n_drops_cap: int, check_every: int
) -> dict:
    """Drive a fresh lattice until percolation. Return per-run metrics."""
    spanning = {"first_t": None, "first_p": None, "spans_at": None}

    def cb(t, ever, _sizes, _durations):
        r = check_spanning(ever)
        if r.percolates:
            spanning["first_t"] = t
            spanning["first_p"] = fractured_fraction(ever)
            spanning["spans_at"] = (r.spans_x, r.spans_y, r.spans_z)
            return True
        return False

    state, sizes, _durs, ever = run_with_ever_toppled(
        L=L,
        n_drops=n_drops_cap,
        seed=seed,
        check_every=check_every,
        percolation_callback=cb,
    )
    return {
        "L": L,
        "seed": seed,
        "first_t": spanning["first_t"],
        "first_p": spanning["first_p"],
        "spans_at": spanning["spans_at"],
        "n_drops_executed": int(sizes.size),
        "final_lost": int(state.grains_lost),
        "final_z_sum": int(state.z.sum()),
    }


def run_all() -> list[dict]:
    rows: list[dict] = []
    for L, n_seeds, cap, check_every in RUN_PLAN:
        print(f"L={L}  n_seeds={n_seeds}  cap={cap}  check_every={check_every}")
        t_L = time.time()
        for s_idx in range(n_seeds):
            seed = 1000 * L + s_idx
            t0 = time.time()
            row = find_first_spanning(L, seed, cap, check_every)
            row["wall"] = time.time() - t0
            rows.append(row)
            sa = row["spans_at"]
            sa_str = (
                f"{'x' if sa[0] else '-'}{'y' if sa[1] else '-'}{'z' if sa[2] else '-'}"
                if sa is not None else "----"
            )
            if row["first_p"] is None:
                print(
                    f"  seed={seed}  NO SPAN within cap={cap}  "
                    f"drops_run={row['n_drops_executed']}  "
                    f"wall={row['wall']:.1f}s"
                )
            else:
                print(
                    f"  seed={seed}  drops_to_span={row['first_t'] + 1:7d}  "
                    f"p={row['first_p']:.4f}  spans={sa_str}  "
                    f"wall={row['wall']:.1f}s"
                )
        print(f"  L={L} total wall: {time.time() - t_L:.1f}s")
    return rows


def summarize(rows: list[dict]) -> dict:
    Ls = sorted({r["L"] for r in rows})
    summary = {}
    for L in Ls:
        ps = [r["first_p"] for r in rows if r["L"] == L and r["first_p"] is not None]
        ts = [r["first_t"] + 1 for r in rows if r["L"] == L and r["first_t"] is not None]
        if not ps:
            summary[L] = None
            continue
        ps_arr = np.asarray(ps, dtype=np.float64)
        ts_arr = np.asarray(ts, dtype=np.int64)
        summary[L] = {
            "n_seeds": len(ps),
            "mean_p": float(ps_arr.mean()),
            "stderr_p": float(ps_arr.std(ddof=1) / np.sqrt(len(ps))) if len(ps) > 1 else 0.0,
            "std_p": float(ps_arr.std(ddof=1)) if len(ps) > 1 else 0.0,
            "mean_t": float(ts_arr.mean()),
            "median_t": float(np.median(ts_arr)),
        }
    return summary


def fit_constant(summary: dict) -> tuple[float, float, float, int]:
    """Weighted-mean constant model: p_c(L) = p_c for all L.

    Returns (mean_p, stderr, chi2_per_dof, dof). Weights are 1 / stderr_p^2
    per L. chi2/dof is the goodness of the flat model; values near 1
    indicate "no L-dependence detected"; values >> 1 suggest the constant
    model is too restrictive and FSS scaling should be considered.
    """
    Ls = sorted([L for L in summary if summary[L] is not None])
    p_arr = np.array([summary[L]["mean_p"] for L in Ls], dtype=np.float64)
    sigma = np.array([max(summary[L]["stderr_p"], 1e-6) for L in Ls], dtype=np.float64)
    w = 1.0 / sigma ** 2
    mean = float((p_arr * w).sum() / w.sum())
    stderr = float(1.0 / np.sqrt(w.sum()))
    dof = max(len(Ls) - 1, 1)
    chi2 = float(((p_arr - mean) ** 2 * w).sum())
    return mean, stderr, chi2 / dof, dof


def fit_pc_fss(summary: dict) -> tuple[float, float, float, float] | None:
    """Fit p_c(L) = p_c_inf + a * L^(-1/nu). Returns None if degenerate.

    Degeneracy check: if the recovered nu is implausibly large (> 20) or
    the recovered |a| is enormous compared to plausible scale (~ 1) the
    fit has collapsed to the constant model and we refuse to report it.
    """
    Ls = sorted([L for L in summary if summary[L] is not None])
    if len(Ls) < 3:
        return None
    L_arr = np.array(Ls, dtype=np.float64)
    p_arr = np.array([summary[L]["mean_p"] for L in Ls], dtype=np.float64)
    sigma = np.array([max(summary[L]["stderr_p"], 1e-6) for L in Ls], dtype=np.float64)

    def model(L, p_inf, a, nu):
        return p_inf + a * L ** (-1.0 / nu)

    p0 = [p_arr[-1], p_arr[0] - p_arr[-1], 1.0]
    try:
        popt, _pcov = curve_fit(
            model, L_arr, p_arr, p0=p0, sigma=sigma, absolute_sigma=True, maxfev=5000
        )
    except Exception as e:
        print(f"  FSS fit failed: {e}")
        return None
    p_inf, a, nu = (float(popt[0]), float(popt[1]), float(popt[2]))
    if nu > 20.0 or abs(a) > 5.0 or not (0.0 < p_inf < 1.0):
        return None
    residual = p_arr - model(L_arr, *popt)
    rms = float(np.sqrt((residual ** 2).mean()))
    return p_inf, a, nu, rms


def make_plots(rows: list[dict], summary: dict, outdir: Path, stamp: str) -> None:
    Ls = sorted([L for L in summary if summary[L] is not None])
    mean_p = np.array([summary[L]["mean_p"] for L in Ls])
    stderr_p = np.array([summary[L]["stderr_p"] for L in Ls])

    const_mean, const_stderr, chi2_per_dof, dof = fit_constant(summary)
    fss = fit_pc_fss(summary)

    # --- p vs L scatter + constant fit ---
    fig, ax = plt.subplots(figsize=(7, 5))
    for L in Ls:
        ps = [r["first_p"] for r in rows if r["L"] == L and r["first_p"] is not None]
        ax.scatter([L] * len(ps), ps, color="0.65", alpha=0.5, s=15)
    ax.errorbar(Ls, mean_p, yerr=stderr_p, fmt="o", color="C0",
                capsize=4, label="per-L mean +/- stderr")
    # Constant model
    ax.axhline(const_mean, color="C1", linestyle="-",
               label=rf"constant fit $p_c = {const_mean:.4f} \pm {const_stderr:.4f}$"
                     rf"  ($\chi^2/\nu = {chi2_per_dof:.2f}$, $\nu = {dof}$)")
    ax.fill_between([min(Ls), max(Ls)],
                    const_mean - const_stderr, const_mean + const_stderr,
                    color="C1", alpha=0.15)
    ax.axhline(0.3116, color="gray", linestyle=":",
               label=r"uncorrelated site percolation $p_c \approx 0.3116$")
    ax.set_xlabel("$L$")
    ax.set_ylabel("$p$ at first spanning")
    ax.set_title("3D Manna: dynamics-driven percolation threshold")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc="best")
    fig.tight_layout()
    path = outdir / f"manna_3d_pc_vs_L_{stamp}.png"
    fig.savefig(path, dpi=150)
    print(f"  wrote {path}")
    plt.close(fig)

    # --- 1/L plot showing flatness vs FSS extrapolation if it fits ---
    fig, ax = plt.subplots(figsize=(7, 5))
    inv_L = 1.0 / np.array(Ls, dtype=np.float64)
    ax.errorbar(inv_L, mean_p, yerr=stderr_p, fmt="o", color="C0", capsize=4,
                label="per-L mean")
    ax.axhline(const_mean, color="C1", linestyle="-",
               label=rf"constant $p_c = {const_mean:.4f}$")
    ax.fill_between([0.0, inv_L.max() * 1.05],
                    const_mean - const_stderr, const_mean + const_stderr,
                    color="C1", alpha=0.15)
    if fss is not None:
        p_inf, a, nu, rms = fss
        x_plot = np.linspace(1e-6, inv_L.max() * 1.05, 100)
        L_plot = 1.0 / x_plot
        y_plot = p_inf + a * L_plot ** (-1.0 / nu)
        ax.plot(x_plot, y_plot, "r--",
                label=rf"FSS: $p_c(\infty)={p_inf:.4f}$, $\nu={nu:.2f}$, RMS={rms:.4f}")
    else:
        ax.text(0.02, 0.05,
                "FSS fit rejected: data is consistent with no L-dependence",
                transform=ax.transAxes, fontsize=9, color="dimgray", style="italic")
    ax.axhline(0.3116, color="gray", linestyle=":",
               label=r"uncorrelated site percolation $p_c$")
    ax.set_xlabel("$1/L$")
    ax.set_ylabel("$p_c(L)$")
    ax.set_title("FSS axis: $p_c(L)$ vs $1/L$")
    ax.set_xlim(left=-0.002)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc="best")
    fig.tight_layout()
    path = outdir / f"manna_3d_pc_fss_{stamp}.png"
    fig.savefig(path, dpi=150)
    print(f"  wrote {path}")
    plt.close(fig)


def save_raw(rows: list[dict], outdir: Path, stamp: str) -> Path:
    path = outdir / f"manna_3d_pc_data_{stamp}.npz"
    Ls = np.asarray([r["L"] for r in rows])
    seeds = np.asarray([r["seed"] for r in rows])
    first_p = np.asarray([r["first_p"] if r["first_p"] is not None else np.nan
                          for r in rows])
    first_t = np.asarray([r["first_t"] if r["first_t"] is not None else -1
                          for r in rows], dtype=np.int64)
    walls = np.asarray([r["wall"] for r in rows])
    np.savez_compressed(
        path, L=Ls, seed=seeds, first_p=first_p, first_t=first_t, wall=walls
    )
    return path


def main() -> None:
    print("3D Manna percolation FSS")
    print("=" * 60)
    t_start = time.time()
    rows = run_all()
    total = time.time() - t_start

    summary = summarize(rows)
    print()
    print("Summary:")
    print(f"{'L':>4}  {'n':>3}  {'<p_c>':>8}  {'stderr':>8}  {'std':>8}  {'<drops>':>10}")
    for L in sorted(summary):
        if summary[L] is None:
            print(f"{L:>4}   --  no spans recorded")
            continue
        s = summary[L]
        print(
            f"{L:>4}  {s['n_seeds']:>3}  "
            f"{s['mean_p']:>8.4f}  {s['stderr_p']:>8.4f}  {s['std_p']:>8.4f}  "
            f"{s['mean_t']:>10.1f}"
        )

    const_mean, const_stderr, chi2_per_dof, dof = fit_constant(summary)
    print()
    print("Constant-model fit (p_c independent of L):")
    print(f"  p_c       = {const_mean:.4f} +/- {const_stderr:.4f}")
    print(f"  chi^2/dof = {chi2_per_dof:.2f}  (dof = {dof})")
    if chi2_per_dof < 3.0:
        print("  -> constant model is consistent with data; no detectable")
        print("     finite-size correction across the L values surveyed.")
    else:
        print("  -> constant model is REJECTED; real L-dependence present.")

    fit = fit_pc_fss(summary)
    print()
    if fit is None:
        print("FSS scaling fit p_c(L) = p_inf + a * L^(-1/nu):")
        print("  REJECTED as degenerate (data does not constrain nu).")
        print("  Use the constant-model p_c above as the cosmologically")
        print("  meaningful threshold.")
    else:
        p_inf, a, nu, rms = fit
        print("FSS scaling fit p_c(L) = p_inf + a * L^(-1/nu):")
        print(f"  p_c(inf) = {p_inf:.4f}")
        print(f"  a        = {a:+.4f}")
        print(f"  nu       = {nu:.3f}")
        print(f"  RMS      = {rms:.5f}")
    print()
    print("Reference: uncorrelated 3D site percolation p_c = 0.3116.")

    print()
    print(f"Total wall time: {total:.1f}s")

    outdir = REPO_ROOT / "data" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    npz = save_raw(rows, outdir, stamp)
    print(f"Saved raw data to {npz}")
    make_plots(rows, summary, outdir, stamp)


if __name__ == "__main__":
    main()
