"""Fractal dimension of the spanning cluster from the 3D smoke run.

Loads the saved ever_toppled mask from the L=24 smoke run, isolates the
spanning cluster via 6-connectivity labeling, runs box-counting on it,
and writes a log-log plot.

Run:
    .venv/bin/python scripts/run_milestone3_dimension.py [path/to/smoke.npz]

If no path is given, picks the most recent manna_3d_smoke_*.npz in
data/outputs/.
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.cluster_geometry import box_count_dimension  # noqa: E402
from void_cascade.percolation import check_spanning  # noqa: E402
from scipy import ndimage  # noqa: E402


def latest_smoke_npz() -> Path:
    outdir = REPO_ROOT / "data" / "outputs"
    candidates = sorted(outdir.glob("manna_3d_smoke_*.npz"))
    if not candidates:
        raise FileNotFoundError("No manna_3d_smoke_*.npz in data/outputs/")
    return candidates[-1]


def largest_cluster_mask(mask: np.ndarray) -> tuple[np.ndarray, int]:
    """Return a bool mask of the single largest connected component.

    6-connectivity, matching the percolation detector.
    """
    structure = ndimage.generate_binary_structure(rank=3, connectivity=1)
    labels, _n = ndimage.label(mask, structure=structure)
    if labels.max() == 0:
        raise ValueError("empty mask")
    sizes = np.bincount(labels.ravel())
    sizes[0] = 0   # ignore background
    biggest_id = int(sizes.argmax())
    return labels == biggest_id, int(sizes[biggest_id])


def main() -> None:
    if len(sys.argv) > 1:
        path = Path(sys.argv[1])
    else:
        path = latest_smoke_npz()
    print(f"Loading {path}")
    data = np.load(path)
    ever = data["ever_toppled"]
    L = int(data["L"])
    print(f"L = {L}, ever_toppled fraction = {ever.mean():.4f}")

    span = check_spanning(ever)
    print(
        f"Spanning: x={span.spans_x} y={span.spans_y} z={span.spans_z}, "
        f"n_clusters={span.n_clusters}, largest={span.largest_cluster_size}"
    )

    largest, n_largest = largest_cluster_mask(ever)
    print(f"Isolated largest cluster: {n_largest} sites "
          f"({n_largest / ever.size:.4f} of lattice)")

    D, D_err, sizes_used, counts = box_count_dimension(largest)
    print()
    print("Box-counting on the largest cluster (spanning):")
    print(f"  box sizes used: {sizes_used.tolist()}")
    print(f"  counts:         {counts.tolist()}")
    print(f"  D_box           = {D:.3f} +/- {D_err:.3f}")
    print()
    print("Reference: a compact 3D blob has D_box ~ 3 (with finite-size")
    print("bias dragging it a few percent low). For 3D site percolation")
    print("at p_c the spanning cluster has D_f ~ 2.523.")

    outdir = REPO_ROOT / "data" / "outputs"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    fig, ax = plt.subplots(figsize=(6, 5))
    inv_b = 1.0 / sizes_used.astype(np.float64)
    ax.loglog(inv_b, counts, "o-", label="data")
    # Fitted line in log-log
    coeffs = np.polyfit(np.log10(inv_b), np.log10(counts.astype(np.float64)), 1)
    fit_y = 10 ** (coeffs[0] * np.log10(inv_b) + coeffs[1])
    ax.loglog(inv_b, fit_y, "r--", label=rf"fit slope $D_{{box}} = {D:.3f}\pm{D_err:.3f}$")
    ax.set_xlabel("$1/b$")
    ax.set_ylabel("$N(b)$  (boxes of side $b$ intersecting cluster)")
    ax.set_title(f"3D Manna spanning cluster (L={L}): box-count")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig_path = outdir / f"manna_3d_cluster_fractal_{stamp}.png"
    fig.savefig(fig_path, dpi=150)
    print(f"  wrote {fig_path}")


if __name__ == "__main__":
    main()
