"""Render a GIF of 2D Manna dynamics and a snapshot of a large avalanche.

Run with the repo venv:
    .venv/bin/python scripts/animate_2d.py

Defaults are tuned for visibility: small L for the animation so individual
cells are readable, and a larger L for the snapshot so the avalanche cluster
has interesting structure. Both files are written to data/outputs/.
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.visualize_2d import (  # noqa: E402
    record_animation,
    snapshot_large_avalanche,
)


def main() -> None:
    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    anim_path = outdir / f"manna_2d_dynamics_{stamp}.gif"
    print(f"Rendering Manna 2D dynamics animation to {anim_path}")
    record_animation(L=32, max_frames=400, seed=7, out_path=anim_path, fps=12)
    print(f"  wrote {anim_path}  ({anim_path.stat().st_size / 1024:.0f} KB)")

    snap_path = outdir / f"manna_2d_avalanche_snapshot_{stamp}.png"
    print(f"Rendering large-avalanche snapshot to {snap_path}")
    # L=96 with a tight min_area gives a recognizably extended cluster
    # without making the search too long.
    snapshot_large_avalanche(
        L=96,
        burn_in=2 * 96 * 96,
        min_area=96 * 96 // 8,  # cluster covers >=12.5% of the lattice
        max_search_drops=20_000,
        seed=11,
        out_path=snap_path,
    )
    print(f"  wrote {snap_path}  ({snap_path.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
