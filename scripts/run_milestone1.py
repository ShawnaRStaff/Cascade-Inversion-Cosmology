"""Milestone 1 demo: 1D Oslo sandpile.

Runs a long simulation, extracts the avalanche size distribution from the
steady state, fits a power law to the scaling range, and writes a log-log
plot to data/outputs/. Intended to confirm tau ~ 1.55 (Christensen et al.
1996) over at least two decades of avalanche size.

Run with the repo venv:
    .venv/bin/python scripts/run_milestone1.py
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Allow running as a plain script from the repo root without installing.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from void_cascade.analysis import fit_power_law, log_binned_pdf  # noqa: E402
from void_cascade.sandpile_1d import run  # noqa: E402


def main() -> None:
    L = 128
    n_drops = 60_000
    # The pile fills at roughly L^2 / 2 grains worth of drives once you
    # account for the steady-state mean slope. Conservatively discard the
    # first ~L^2 drives so we sample only the SOC steady state.
    transient = L * L

    print(f"Running 1D Oslo: L={L}, n_drops={n_drops}, transient={transient}")
    state, sizes, durations = run(L=L, n_drops=n_drops, seed=12345)

    print(f"  mean slope at end:        {state.z.mean():.3f}")
    print(f"  grains lost off edge:     {state.grains_lost}")
    print(f"  fraction of drives s>0:   {(sizes > 0).mean():.3f}")
    print(f"  largest avalanche:        {sizes.max()}")

    s_steady = sizes[transient:]
    centers, pdf, counts = log_binned_pdf(s_steady, n_bins=40)

    # Choose the scaling range. Lower cutoff above the small-s rollover, upper
    # cutoff below the finite-size knee. The L-dependent upper cutoff is a
    # heuristic from the Oslo finite-size scaling: <s> ~ L, cutoff ~ L^D with
    # D ~ 2.25, but we leave headroom because the empirical knee is gentle.
    s_fit_min = 5.0
    s_fit_max = float(L * L) * 0.5
    tau, c = fit_power_law(centers, pdf, s_min=s_fit_min, s_max=s_fit_max)
    print(f"  fitted tau:               {tau:.3f}")
    print("  reference tau (1D Oslo):  ~1.55 (Christensen et al. 1996)")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(centers, pdf, "o", label="simulation")
    s_line = np.array([s_fit_min, s_fit_max])
    ax.loglog(
        s_line,
        10**c * s_line ** (-tau),
        "r-",
        label=rf"fit $P(s)\propto s^{{-{tau:.2f}}}$",
    )
    ax.set_xlabel("avalanche size $s$")
    ax.set_ylabel("$P(s)$")
    ax.set_title(f"1D Oslo sandpile, L={L}, steady-state avalanche distribution")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)

    outdir = Path(__file__).resolve().parents[1] / "data" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    fig_path = outdir / f"oslo_1d_avalanche_pdf_{stamp}.png"
    fig.tight_layout()
    fig.savefig(fig_path, dpi=150)
    print(f"  saved plot:               {fig_path}")


if __name__ == "__main__":
    main()
