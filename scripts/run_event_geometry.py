"""Run the per-event mask geometry experiment (M7 asymptote discriminator).

Replays the M6 FSS trajectories (same L, n_drops, seeds => bit-identical
avalanche sequences) with per-event carrier/sink classification enabled,
and saves the record-breaking event series plus the peak event's mask and
pre-event carrier field. See src/void_cascade/event_geometry.py for the
physics.

Usage:
    python scripts/run_event_geometry.py --L 48 --n-drops 250000 \
        --seeds 14800 14801 14802 14803 14804 --outdir data/outputs/event_geometry_<stamp>
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from void_cascade.event_geometry import run_peak_event_geometry


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, required=True)
    ap.add_argument("--n-drops", type=int, required=True)
    ap.add_argument("--seeds", type=int, nargs="+", required=True)
    ap.add_argument("--outdir", type=str, required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for seed in args.seeds:
        t0 = time.time()
        print(f"[{time.strftime('%F %T')}] L={args.L} seed={seed} starting", flush=True)
        out = run_peak_event_geometry(args.L, args.n_drops, seed)
        rec = out["records"]
        path = outdir / f"L{args.L}_s{seed}_eventgeom.npz"
        np.savez_compressed(
            path,
            L=args.L,
            seed=seed,
            n_drops=args.n_drops,
            sizes=out["sizes"],
            unique_sizes=out["unique_sizes"],
            rec_drop=rec["drop"],
            rec_unique=rec["unique"],
            rec_toppled_carriers=rec["toppled_carriers"],
            rec_toppled_sinks=rec["toppled_sinks"],
            rec_sink_participation=rec["sink_participation"],
            rec_carrier_coverage=rec["carrier_coverage"],
            rec_n_carriers_pre=rec["n_carriers_pre"],
            peak_mask=out["peak_mask"],
            peak_carriers_pre=out["peak_carriers_pre"],
            final_z=out["final_z"],
            grains_lost=out["grains_lost"],
        )
        peak_sp = rec["sink_participation"][-1] if len(rec["drop"]) else float("nan")
        print(
            f"[{time.strftime('%F %T')}] L={args.L} seed={seed} done in "
            f"{time.time() - t0:.0f}s  peak_unique={rec['unique'][-1]}  "
            f"peak_sink_participation={peak_sp:.3f}  -> {path}",
            flush=True,
        )


if __name__ == "__main__":
    main()
