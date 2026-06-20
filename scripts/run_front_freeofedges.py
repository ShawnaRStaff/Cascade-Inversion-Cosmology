"""Front, free of edges: is the spread real, or a small-box / under-driving artifact?

Runs several configs SEQUENTIALLY (no parallel heavy jobs) with FAIR per-cell
driving (n_drive_sites scaled to area, matched to the L=48 baseline rate), so
any size difference is size, not under-driving. We run them all and report
them all -- no run favored.

Honest constraints, stated up front:
- Still a FLAT square grid (distances straight-line). This removes the
  edge/boundary and under-driving confounds, NOT the flatness one.
- Substrate initialized to a constructed near-critical state (skips the slow
  fill); it self-organizes from there.

Configs:
  A  L=80  full loaded box   -- front travel vs box size (fair drive)
  B  L=160 full loaded box   -- same, bigger; does behaviour change with size?
  C  L=160 loaded DISK in empty box -- fuel test: does the front cross emptiness?

Per run we record: tip step, the front's leading edge over time, whether it
reached the box edge, and (disk) whether it crossed into the empty region.
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
from void_cascade.cascade_heat_2d import (  # noqa: E402
    avalanche_2d, diffuse_and_cool_2d, initialize_2d,
)
from void_cascade.heat_gated import release  # noqa: E402

BASELINE_CELLS = 48 * 48  # per-cell drive rate reference (L=48, 1 site/step)


def fair_sites(n_loaded_cells: int) -> int:
    return max(1, round(n_loaded_cells / BASELINE_CELLS))


def params(heat=0.10, cooling=0.10):
    return CascadeParams(fracture_density=2.0, heat_per_crack=heat, diffuse=0.15,
                         cooling=cooling, melt_heat=1.0, release_factor=0.5,
                         drive_amount=1.0, n_drive_sites=1)


def leading_edge(mask):
    ys, xs = np.nonzero(mask)
    if xs.size == 0:
        return 0, 0.0
    cx, cy = xs.mean(), ys.mean()
    return int(xs.size), float(np.sqrt((xs - cx) ** 2 + (ys - cy) ** 2).max())


def touches_box_edge(mask, pad=3):
    return bool(mask[:pad].any() or mask[-pad:].any()
               or mask[:, :pad].any() or mask[:, -pad:].any())


def run_config(label, L, region, p, n_steps, seed, disk_frac=0.30, trace=200):
    rng = np.random.default_rng(seed)
    st = initialize_2d(L)
    cx = cy = L / 2.0
    yy, xx = np.mgrid[0:L, 0:L]
    if region == "full":
        loaded = np.ones((L, L), dtype=bool)
        disk_r = None
    else:
        disk_r = disk_frac * L
        loaded = (xx - cx) ** 2 + (yy - cy) ** 2 <= disk_r ** 2
    n_loaded = int(loaded.sum())
    p = CascadeParams(**{**p.__dict__, "n_drive_sites": fair_sites(n_loaded)})

    # constructed near-critical loaded state, then settle.
    st.density[loaded] = rng.uniform(0.0, 2.0, size=n_loaded)
    avalanche_2d(st, p, rng, 50 * L * L)

    drive_flat = np.flatnonzero(loaded.ravel())
    densf = st.density.ravel()
    maxsw = 50 * L * L
    tip = None
    edge_trace = []
    edge_hit = False
    crossed_into_empty = False

    for t in range(n_steps):
        sites = rng.choice(drive_flat, size=p.n_drive_sites)
        np.add.at(densf, sites, p.drive_amount)
        avalanche_2d(st, p, rng, maxsw)
        st.heat = diffuse_and_cool_2d(st.heat, p)
        _, n_rel = release(st, p)
        if tip is None and n_rel > 0:
            tip = t
        if tip is not None:
            dt = t - tip
            n, edge = leading_edge(st.released)
            edge_trace.append((dt, n, edge))
            if touches_box_edge(st.released):
                edge_hit = True
            if disk_r is not None and edge > disk_r + 3:
                crossed_into_empty = True
            if dt >= trace:
                break

    res = {
        "label": label, "L": L, "region": region, "n_loaded_cells": n_loaded,
        "fair_drive_sites": p.n_drive_sites, "tip_step": tip,
        "ran_away": bool(st.released.any()),
        "final_released_fraction": float(st.released.mean()),
        "final_released_of_loaded": float(st.released[loaded].mean()) if n_loaded else 0.0,
        "max_front_edge": max((e for _, _, e in edge_trace), default=0.0),
        "box_half": L / 2.0, "disk_radius": disk_r,
        "reached_box_edge": edge_hit, "crossed_into_empty": crossed_into_empty,
        "edge_trace": edge_trace,
        "final_released_mask": st.released.copy(),
    }
    print(f"[{label}] L={L} {region} drive={p.n_drive_sites} tip={tip} "
          f"ran_away={res['ran_away']} max_edge={res['max_front_edge']:.1f} "
          f"(box_half={L/2:.0f}{', disk_r='+format(disk_r,'.0f') if disk_r else ''}) "
          f"reached_edge={edge_hit} crossed_empty={crossed_into_empty} "
          f"released={res['final_released_fraction']:.1%}")
    return res


def main() -> None:
    out_dir = REPO_ROOT / "data" / "outputs" / f"front_freeofedges_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"=== Front free of edges (fair drive, sequential) ===\nOutput: {out_dir}\n")

    runs = []
    runs.append(run_config("A_L80_full", 80, "full", params(), 3000, seed=10))
    runs.append(run_config("B_L160_full", 160, "full", params(), 3500, seed=11))
    runs.append(run_config("C_L160_disk", 160, "disk", params(), 3500, seed=12))

    summary = []
    for r in runs:
        s = {k: r[k] for k in ("label", "L", "region", "n_loaded_cells",
                               "fair_drive_sites", "tip_step", "ran_away",
                               "final_released_fraction", "final_released_of_loaded",
                               "max_front_edge", "box_half", "disk_radius",
                               "reached_box_edge", "crossed_into_empty")}
        summary.append(s)
    with open(out_dir / "results.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Plot: edge-vs-time for the two full boxes; final mask for the disk test.
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for r in runs:
        if r["region"] == "full" and r["edge_trace"]:
            t = [d for d, _, _ in r["edge_trace"]]; e = [ee for _, _, ee in r["edge_trace"]]
            axes[0].plot(t, e, "o-", ms=3, label=f"{r['label']} (box_half={r['box_half']:.0f})")
            axes[0].axhline(r["box_half"], ls=":", alpha=0.4)
    axes[0].set_xlabel("steps after tip"); axes[0].set_ylabel("front leading edge (cells)")
    axes[0].set_title("Front travel vs size (fair drive); dotted = box half-width")
    axes[0].legend(); axes[0].grid(True, alpha=0.3)

    disk = next((r for r in runs if r["region"] == "disk"), None)
    if disk is not None:
        axes[1].imshow(disk["final_released_mask"], cmap="hot", interpolation="nearest")
        th = np.linspace(0, 2 * np.pi, 200)
        axes[1].plot(disk["L"] / 2 + disk["disk_radius"] * np.cos(th),
                     disk["L"] / 2 + disk["disk_radius"] * np.sin(th), "c--", lw=1.5,
                     label="loaded-disk edge")
        axes[1].set_title(f"Fuel test: did the front cross into empty?\ncrossed={disk['crossed_into_empty']}")
        axes[1].legend(); axes[1].set_xticks([]); axes[1].set_yticks([])
    fig.tight_layout(); fig.savefig(out_dir / "front_freeofedges.png", dpi=150)

    print(f"\nResults: {out_dir}/results.json\nPlot: {out_dir}/front_freeofedges.png")


if __name__ == "__main__":
    main()
