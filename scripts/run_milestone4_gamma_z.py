"""M4 Speculation #3: compare our gamma(p) trajectory to observed gamma(z).

The two-point correlation function exponent gamma for galaxies has been
measured at many redshifts. The striking empirical fact is that gamma
stays approximately constant at ~1.7-1.9 across z = 0 to z ~ 3, while
the amplitude r_0 evolves more strongly. This is one of the well-known
near-universalities of the galaxy correlation function.

Our gamma(p) varies strongly with p. If p maps monotonically to cosmic
time / redshift, then gamma(p(z)) should also vary strongly with z.
Observations say otherwise. The decisive test: does there exist a
plausible monotonic p(z) mapping under which our gamma(p) curve
matches the observed (nearly flat) gamma(z) data?

Published galaxy correlation gamma values (compiled from the
literature, no catalog download needed):

  z ~ 0.03  gamma=1.77  Davis & Peebles 1983 (CfA)
  z ~ 0.1   gamma=1.84  Zehavi et al. 2005 (SDSS DR3)
  z ~ 0.1   gamma=1.91  Zehavi et al. 2011 (SDSS DR7, fiducial)
  z ~ 0.5   gamma=1.90  Marulli et al. 2013 (VIPERS)
  z ~ 1.0   gamma=1.66  Coil et al. 2006 (DEEP2)
  z ~ 1.0   gamma=1.86  Marulli et al. 2013 (VIPERS high-z)
  z ~ 3.0   gamma=1.80  Foucaud et al. 2003 (LBGs)

These represent different galaxy populations and selection effects,
but the slope gamma is consistently ~1.7-1.9 over a redshift range
spanning the observable universe.

Run:
    .venv/bin/python scripts/run_milestone4_gamma_z.py
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))


# Observed gamma(z) compilation
OBSERVED = [
    # (z, gamma, gamma_err, citation_short)
    (0.03, 1.77, 0.04, "Davis+Peebles 1983 CfA"),
    (0.1,  1.84, 0.05, "Zehavi+ 2005 SDSS"),
    (0.1,  1.91, 0.06, "Zehavi+ 2011 SDSS DR7"),
    (0.5,  1.90, 0.08, "Marulli+ 2013 VIPERS"),
    (1.0,  1.66, 0.12, "Coil+ 2006 DEEP2"),
    (1.0,  1.86, 0.10, "Marulli+ 2013 VIPERS z~1"),
    (3.0,  1.80, 0.15, "Foucaud+ 2003 LBGs"),
]


def latest_xi_npz() -> Path:
    outdir = REPO_ROOT / "data" / "outputs"
    candidates = sorted(outdir.glob("manna_3d_xi_data_*.npz"))
    if not candidates:
        raise FileNotFoundError("No manna_3d_xi_data_*.npz in data/outputs/")
    return candidates[-1]


def load_model_gamma(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load our gamma(p) trajectory from the latest M4 xi data.

    Reproduces the pure-power-law fit over r in [2, L/4], same as the
    fit performed in run_milestone4.py.
    """
    data = np.load(path)
    L = int(data["L"])
    p_targets = sorted([float(p) for p in data["p_targets"]])
    p_arr = []
    gamma_arr = []
    r_min_fit = 2.0
    r_max_fit = float(L) / 4.0
    for p in p_targets:
        pp = int(round(p * 100))
        r = data[f"xi_p{pp:02d}_r"]
        xi = data[f"xi_p{pp:02d}_mean"]
        np_arr = data[f"xi_p{pp:02d}_npairs"]
        sel = (r >= r_min_fit) & (r <= r_max_fit) & np.isfinite(xi) & (xi > 0) & (np_arr > 50)
        if sel.sum() < 3:
            continue
        coeffs = np.polyfit(np.log10(r[sel]), np.log10(xi[sel]), 1)
        gamma = -float(coeffs[0])
        p_arr.append(p)
        gamma_arr.append(gamma)
    return np.asarray(p_arr), np.asarray(gamma_arr)


def main() -> None:
    path = latest_xi_npz()
    print(f"Loading model gamma(p) from {path}")
    p_model, gamma_model = load_model_gamma(path)
    print(f"Model gamma(p) at {len(p_model)} p values:")
    for p, g in zip(p_model, gamma_model):
        print(f"  p={p:.2f}  gamma={g:.3f}")

    print()
    print("Observed gamma(z) from literature:")
    for z, g, ge, cite in OBSERVED:
        print(f"  z={z:5.2f}  gamma={g:.2f} +/- {ge:.2f}  {cite}")

    obs_z = np.array([o[0] for o in OBSERVED])
    obs_g = np.array([o[1] for o in OBSERVED])
    obs_e = np.array([o[2] for o in OBSERVED])

    weights = 1.0 / obs_e ** 2
    obs_weighted_mean = float((obs_g * weights).sum() / weights.sum())
    obs_weighted_se = float(1.0 / np.sqrt(weights.sum()))
    print()
    print(f"Observed weighted mean gamma = {obs_weighted_mean:.3f} +/- {obs_weighted_se:.3f}")
    print(f"Observed gamma range across z = 0.03 to 3.0: {obs_g.min():.2f} to {obs_g.max():.2f}")

    # --- The decisive question: a monotonic p(z) mapping that aligns the two ---
    # If observed gamma is nearly flat at gamma_obs, then any p(z) mapping
    # that keeps p in the narrow band where gamma_model ~ gamma_obs will
    # work. From the model data, this is p ~ 0.60 to 0.70.
    print()
    p_band_low = None
    p_band_high = None
    obs_low = obs_weighted_mean - 0.15
    obs_high = obs_weighted_mean + 0.15
    for p, g in zip(p_model, gamma_model):
        if obs_low <= g <= obs_high:
            if p_band_low is None:
                p_band_low = p
            p_band_high = p
    if p_band_low is not None:
        print(f"Model p values consistent with observed gamma "
              f"in [{obs_low:.2f}, {obs_high:.2f}]:  p in [{p_band_low:.2f}, {p_band_high:.2f}]")
    else:
        print(f"No model p value lies within +/-0.15 of observed gamma mean.")

    # --- Test 1: model variation across "observable epoch" ---
    # If the universe has spanned z=0 to z=3 and our gamma observation is
    # flat to ~0.1, then dp/dz must be small. Quantify how slow.
    print()
    print("Test: the universe has gamma roughly constant from z=0 to z=3.")
    print("If our model maps to cosmic time, dp/d(cosmic_time) must be small")
    print("in the observable epoch. From the model data, gamma changes by")
    delta_g_per_dp = []
    for i in range(1, len(p_model)):
        delta_g_per_dp.append(
            (gamma_model[i] - gamma_model[i-1]) / (p_model[i] - p_model[i-1])
        )
    delta_g_per_dp = np.array(delta_g_per_dp)
    # Around the matching p band (~0.65), find local slope dgamma/dp
    if p_band_low is not None:
        # Pick the slope in the band
        mid_p = 0.5 * (p_band_low + p_band_high)
        idx_near = int(np.argmin(np.abs(p_model[:-1] + np.diff(p_model)/2 - mid_p)))
        local_slope = float(delta_g_per_dp[idx_near])
        print(f"  d(gamma)/dp ~ {local_slope:.2f} per unit p around p ~ {mid_p:.2f}")
        # Observed: d(gamma)/d(z) approximately
        obs_slope = float(np.polyfit(obs_z, obs_g, 1)[0])
        print(f"  d(gamma)/dz observed ~ {obs_slope:+.3f} per unit z")
        # If d(gamma)/d(z)_observed = d(gamma)/dp_model * dp/dz, then
        # dp/dz = (obs slope) / (model slope)
        if local_slope != 0:
            dpdz = obs_slope / local_slope
            print(f"  Required dp/dz ~ {dpdz:+.4f}  (per unit redshift z)")
            print(f"  Over z = 0 to 3, this corresponds to Delta p ~ {dpdz * 3:+.3f}")

    # --- Plots ---
    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Left: gamma(p) model + observed band
    ax1.plot(p_model, gamma_model, "o-", color="C0", label=r"Manna $\gamma(p)$")
    ax1.axhspan(obs_g.min(), obs_g.max(), color="C1", alpha=0.15,
                label=fr"observed $\gamma$ range (0$\leq z \leq$3)")
    ax1.axhline(obs_weighted_mean, color="C1", linestyle="--",
                label=fr"observed weighted mean $\gamma = {obs_weighted_mean:.2f}$")
    if p_band_low is not None:
        ax1.axvspan(p_band_low, p_band_high, color="C2", alpha=0.15,
                    label=fr"model match band: $p \in [{p_band_low:.2f}, {p_band_high:.2f}]$")
    ax1.set_xlabel("$p$  (fractured fraction)")
    ax1.set_ylabel(r"$\gamma$  (two-point correlation slope)")
    ax1.set_title(r"Model $\gamma(p)$ vs observed range")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8)

    # Right: gamma(z) observed
    ax2.errorbar(obs_z, obs_g, yerr=obs_e, fmt="o", capsize=4, color="C1",
                 label="observed (per ref)")
    for z, g, e, cite in OBSERVED:
        ax2.annotate(cite.split()[0], (z, g),
                     xytext=(3, 3), textcoords="offset points", fontsize=7,
                     alpha=0.7)
    ax2.axhline(obs_weighted_mean, color="C1", linestyle="--",
                alpha=0.6,
                label=fr"weighted mean $\gamma = {obs_weighted_mean:.2f}$")
    ax2.set_xlabel("$z$  (cosmological redshift)")
    ax2.set_ylabel(r"$\gamma$")
    ax2.set_title(r"Observed $\gamma(z)$ for galaxy $\xi(r)$")
    ax2.set_xscale("symlog", linthresh=0.1)
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    out_path = outdir / f"manna_3d_gamma_p_vs_z_{stamp}.png"
    fig.savefig(out_path, dpi=150)
    print()
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
