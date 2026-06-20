"""Bigger-lattice front speed (FLAT-GEOMETRY BASELINE).

IMPORTANT framing: this measures how fast the catastrophe front advances in
*cells per step* on a flat square grid. That is a DYNAMICAL fact -- is the
propagation intrinsically steady, speeding up, or slowing down? It is NOT the
physical-distance / Hubble answer: straight-line distance on a flat lattice
is the wrong geometry for that (an open question -- the substrate geometry is
unjustified/probably not flat). We run it to isolate the dynamics from the
distance-metric question, with that caveat stated up front.

Bigger L gives the front room to advance many steps before hitting the
boundary (the small-L run saturated in ~20 steps). We track the leading edge
(max distance of any released cell from the centroid) per step after the tip
and fit it over the pre-saturation region.
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

from void_cascade.cascade_heat import CascadeParams  # noqa: E402
from void_cascade.cascade_heat_2d import initialize_2d, step_2d  # noqa: E402

L = 100
N_STEPS = 6000
COOLING = 0.1
HEAT = 0.12
TRACE_STEPS = 140
SEEDS = [0, 1]


def params():
    return CascadeParams(fracture_density=2.0, heat_per_crack=HEAT, diffuse=0.15,
                         cooling=COOLING, melt_heat=1.0, release_factor=0.5,
                         drive_amount=1.0, n_drive_sites=1)


def leading_edge(mask: np.ndarray) -> tuple[int, float]:
    ys, xs = np.nonzero(mask)
    n = xs.size
    if n == 0:
        return 0, 0.0
    cx, cy = xs.mean(), ys.mean()
    ext = float(np.sqrt((xs - cx) ** 2 + (ys - cy) ** 2).max())
    return n, ext


def one_run(seed):
    rng = np.random.default_rng(seed)
    st = initialize_2d(L)
    maxsw = 50 * L * L
    tip = None
    trace = []
    for t in range(N_STEPS):
        _, n_rel = step_2d(st, params(), rng, maxsw)
        if tip is None and n_rel > 0:
            tip = t
        if tip is not None:
            dt = t - tip
            if dt <= TRACE_STEPS:
                n, edge = leading_edge(st.released)
                trace.append((dt, n, edge))
            else:
                break
    return tip, trace


def fit_front(trace):
    """Fit leading edge vs time over the region before it saturates at ~L/2."""
    t = np.array([d for d, _, _ in trace], float)
    edge = np.array([e for _, _, e in trace], float)
    cap = 0.45 * L  # treat edge beyond ~0.45L as boundary-saturated
    sel = (t >= 1) & (edge < cap)
    if sel.sum() < 4:
        sel = t >= 1
    slope, intercept = np.polyfit(t[sel], edge[sel], 1)
    pred = slope * t[sel] + intercept
    r2 = 1 - ((edge[sel] - pred) ** 2).sum() / max(((edge[sel] - edge[sel].mean()) ** 2).sum(), 1e-9)
    return float(slope), float(r2), int(sel.sum())


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"cascade_heat_2d_frontspeed_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== Front speed, bigger lattice L={L} (FLAT-GEOMETRY baseline) ===\nOutput: {out_dir}\n")

    runs = []
    fig, ax = plt.subplots(figsize=(8, 6))
    for s in SEEDS:
        tip, trace = one_run(s)
        if tip is None or len(trace) < 5:
            print(f"seed {s}: did not tip / too short (tip={tip})")
            continue
        slope, r2, npts = fit_front(trace)
        runs.append({"seed": s, "tip_step": tip, "front_speed_cells_per_step": slope,
                     "fit_r2": r2, "fit_points": npts,
                     "edge_max": max(e for _, _, e in trace)})
        print(f"seed {s}: tip={tip}  front speed = {slope:.3f} cells/step  R2={r2:.3f}  "
              f"(fit over {npts} pre-saturation steps)")
        t = [d for d, _, _ in trace]; e = [ee for _, _, ee in trace]
        ax.plot(t, e, "o-", ms=3, label=f"seed {s}: {slope:.2f} cells/step (R2={r2:.2f})")

    ax.axhline(0.45 * L, color="gray", ls=":", alpha=0.6, label="boundary-saturation cap")
    ax.set_xlabel("steps after tip"); ax.set_ylabel("front leading edge (cells from centre)")
    ax.set_title(f"Flat-grid front speed, L={L} (lattice units, NOT physical distance)")
    ax.legend(); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(out_dir / "frontspeed.png", dpi=150)

    summary = {"L": L, "cooling": COOLING, "heat_per_crack": HEAT,
               "caveat": "cells/step on a FLAT grid; physical distance/Hubble needs the (open) substrate geometry",
               "runs": runs}
    with open(out_dir / "results.json", "w") as f:
        json.dump(summary, f, indent=2)
    if runs:
        sp = [r["front_speed_cells_per_step"] for r in runs]
        r2s = [r["fit_r2"] for r in runs]
        print(f"\nMean front speed {np.mean(sp):.3f} cells/step, mean R2 {np.mean(r2s):.3f}.")
        print("R2 near 1 => steady (constant) advance; low R2 => accelerating/decelerating.")
    print(f"\nResults: {out_dir}/results.json\nPlot: {out_dir}/frontspeed.png")


if __name__ == "__main__":
    main()
