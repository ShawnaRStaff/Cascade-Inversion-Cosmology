"""Milestone 3 smoke test: drive a 3D Manna lattice until percolation.

This is a small-L smoke check that the 3D dynamics + percolation detector
work end-to-end. It is not a quantitative measurement of the percolation
threshold (that requires multiple seeds, multiple L, and a careful
finite-size analysis). The point at this stage is:

  * confirm the sandpile dynamics run without errors at 3D scale,
  * confirm the cumulative ever-toppled set grows monotonically,
  * confirm the percolation detector eventually fires,
  * record what fraction of the lattice has toppled at the moment of
    first spanning, as a first rough estimate of p_c (the percolation
    occupation probability) - or rather of the dynamics-driven analog.

Outputs:
  * console summary
  * `data/outputs/manna_3d_smoke_<stamp>.npz` with the ever-toppled
    mask and a trace of (drop_idx, fractured_fraction, percolates).

Run:
    .venv/bin/python scripts/run_milestone3_smoke.py
"""

from __future__ import annotations

import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.percolation import check_spanning, fractured_fraction  # noqa: E402
from void_cascade.sandpile_3d import run_with_ever_toppled  # noqa: E402


def main() -> None:
    L = 24
    n_drops_cap = 30_000
    check_every = 200

    print(f"3D Manna smoke: L={L}, max n_drops={n_drops_cap}, "
          f"checking percolation every {check_every} drives")

    trace_t: list[int] = []
    trace_p: list[float] = []
    trace_spans: list[str] = []
    first_spanning_t: int | None = None

    def callback(t, ever, sizes, durations):
        nonlocal first_spanning_t
        p = fractured_fraction(ever)
        result = check_spanning(ever)
        trace_t.append(t)
        trace_p.append(p)
        trace_spans.append(
            f"{'x' if result.spans_x else '-'}"
            f"{'y' if result.spans_y else '-'}"
            f"{'z' if result.spans_z else '-'}"
        )
        print(
            f"  drop={t + 1:6d}  p={p:.3f}  "
            f"spans={trace_spans[-1]}  "
            f"largest={result.largest_cluster_size:6d}  "
            f"n_clusters={result.n_clusters:4d}"
        )
        if result.percolates and first_spanning_t is None:
            first_spanning_t = t
            print(f"  >>> percolation reached at drop {t + 1}, p = {p:.3f}")
            return True
        return False

    t0 = time.time()
    state, sizes, durations, ever = run_with_ever_toppled(
        L=L,
        n_drops=n_drops_cap,
        seed=42,
        check_every=check_every,
        percolation_callback=callback,
    )
    elapsed = time.time() - t0

    print()
    print(f"Wall time: {elapsed:.1f}s")
    print(f"Drops executed: {sizes.size}")
    print(f"Conservation: grains_in={sizes.size}, "
          f"sum(z)={int(state.z.sum())}, lost={state.grains_lost}, "
          f"residual={sizes.size - int(state.z.sum()) - state.grains_lost}")
    print(f"Final fractured fraction: {fractured_fraction(ever):.3f}")
    final = check_spanning(ever)
    print(
        f"Final spans: x={final.spans_x} y={final.spans_y} z={final.spans_z}, "
        f"n_clusters={final.n_clusters}, "
        f"largest={final.largest_cluster_size}"
    )

    if first_spanning_t is None:
        print()
        print("DID NOT PERCOLATE within the drop cap. Either raise the cap "
              "or lower the check_every interval.")
    else:
        print()
        print(f"First spanning at drop {first_spanning_t + 1}")

    outdir = REPO_ROOT / "data" / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = outdir / f"manna_3d_smoke_{stamp}.npz"
    np.savez_compressed(
        out_path,
        ever_toppled=ever,
        trace_t=np.asarray(trace_t, dtype=np.int64),
        trace_p=np.asarray(trace_p, dtype=np.float64),
        sizes=sizes,
        durations=durations,
        L=np.int64(L),
        first_spanning_drop=(
            np.int64(first_spanning_t)
            if first_spanning_t is not None
            else np.int64(-1)
        ),
    )
    print(f"Saved trace + final ever_toppled mask to {out_path}")


if __name__ == "__main__":
    main()
