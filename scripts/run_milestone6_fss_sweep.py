"""M6 Phase A — finite-size scaling sweep driver.

Runs the 3D Manna saturation experiment across a configurable L list
with multiple seeds per L. Each (L, seed) job runs in its own worker
process via multiprocessing.Pool, with checkpoint-resume support so
spot-instance interruption never loses progress.

Default sweep (matches the M5-data-driven plan):
    L=32   x 10 seeds
    L=48   x  5 seeds  (we already have 3 from M5 Stage 2)
    L=64   x  5 seeds  (we already have 2 from M5 Stage 5)
    L=96   x  5 seeds  (we have 1 prior run; want a proper ensemble)
    L=128  x  3 seeds  (the upper end that makes the sweep credible)

Drops-per-L scale with L^3 (more cells -> longer to saturate). Tuned
so that each run reaches p=1.0 with ~25% margin in the plateau:
    L=32   ->   100,000 drops
    L=48   ->   250,000 drops
    L=64   ->   600,000 drops
    L=96   -> 1,500,000 drops
    L=128  -> 4,000,000 drops

Outputs land under data/outputs/fss_sweep_{stamp}/ with one
{L}_s{seed}_final.npz per completed job. A summary.json at the end
captures per-job stats (peak event, final p, wall time).

Usage:
    .venv/bin/python -u scripts/run_milestone6_fss_sweep.py
      [--workers N] [--l-list 32,48,64] [--dry-run]

--dry-run prints the job list and exits. Useful for sanity check
before launching the full sweep on a paid instance.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.checkpoint import (  # noqa: E402
    CheckpointPayload,
    CheckpointPolicy,
    load_checkpoint,
    restore_rng,
    save_checkpoint,
)
from void_cascade.sandpile_3d import (  # noqa: E402
    MannaState3D,
    drive,
    initialize,
    relax,
)


DEFAULT_L_LIST = [32, 48, 64, 96, 128]
DEFAULT_SEEDS_PER_L = {32: 10, 48: 5, 64: 5, 96: 5, 128: 3}
DEFAULT_DROPS_FOR_L = {
    32: 100_000,
    48: 250_000,
    64: 600_000,
    96: 1_500_000,
    128: 4_000_000,
}
SNAPSHOT_EVERY = 5_000
CHECKPOINT_EVERY_DROPS = 50_000
CHECKPOINT_EVERY_SECONDS = 600.0


def seed_for(L: int, seed_idx: int) -> int:
    """Deterministic seed assignment so the sweep is fully reproducible."""
    return 10_000 + L * 100 + seed_idx


@dataclass
class JobResult:
    L: int
    seed: int
    drops_executed: int
    final_p: float
    peak_size: int
    peak_pct: float
    wall_seconds: float
    z_avg_final: float
    grains_lost: int


def run_one(
    L: int,
    seed: int,
    n_drops_max: int,
    ckpt_path: Path,
    out_path: Path,
    log_path: Path,
) -> JobResult:
    """Run a single (L, seed) job to completion with checkpoint-resume.

    Resumes from ckpt_path if it exists; otherwise starts fresh.
    Saves checkpoints periodically per CheckpointPolicy.
    On completion, writes the final result to out_path and returns
    a JobResult summary.
    """
    log = open(log_path, "a", buffering=1)  # line-buffered
    log.write(f"[{datetime.now().isoformat()}] L={L} seed={seed} starting\n")

    t0 = time.time()

    payload = load_checkpoint(ckpt_path)
    if payload is not None:
        log.write(
            f"[{datetime.now().isoformat()}] resuming from drop "
            f"{payload.drop}/{n_drops_max}\n"
        )
        assert payload.L == L, "checkpoint L mismatch"
        assert payload.seed == seed, "checkpoint seed mismatch"
        assert payload.n_drops_max == n_drops_max, "checkpoint n_drops_max mismatch"
        rng = restore_rng(payload.rng_state)
        state = payload.state
        ever_toppled = payload.ever_toppled
        sizes = payload.sizes
        durations = payload.durations
        snapshots = payload.snapshots
        t_drop = payload.drop
    else:
        log.write(
            f"[{datetime.now().isoformat()}] starting fresh, n_drops={n_drops_max}\n"
        )
        rng = np.random.default_rng(seed)
        state = initialize(L)
        ever_toppled = np.zeros((L, L, L), dtype=bool)
        sizes = np.zeros(n_drops_max, dtype=np.int64)
        durations = np.zeros(n_drops_max, dtype=np.int64)
        snapshots = []
        t_drop = 0

    policy = CheckpointPolicy(
        every_drops=CHECKPOINT_EVERY_DROPS,
        every_seconds=CHECKPOINT_EVERY_SECONDS,
    )

    while t_drop < n_drops_max:
        drive(state, rng)
        s, T, mask = relax(state, rng, track_support=True)
        sizes[t_drop] = s
        durations[t_drop] = T
        if mask is not None:
            ever_toppled |= mask
        t_drop += 1

        if t_drop % SNAPSHOT_EVERY == 0:
            p = float(ever_toppled.mean())
            window_sizes = sizes[max(0, t_drop - SNAPSHOT_EVERY) : t_drop]
            mean_size_w = float(window_sizes.mean())
            max_size_w = int(window_sizes.max())
            z_mean = float(state.z.mean())
            snapshots.append(
                {
                    "drop": t_drop,
                    "p": p,
                    "mean_size_window": mean_size_w,
                    "max_size_window": max_size_w,
                    "z_mean": z_mean,
                    "grains_lost": int(state.grains_lost),
                }
            )
            elapsed = time.time() - t0
            rate = t_drop / max(elapsed, 1e-9)
            log.write(
                f"[{datetime.now().isoformat()}] drop={t_drop} "
                f"p={p:.5f} <s>_win={mean_size_w:.1f} "
                f"max_s={max_size_w} z={z_mean:.3f} "
                f"rate={rate:.0f}/s elapsed={elapsed:.0f}s\n"
            )

        if policy.should_save(t_drop):
            save_checkpoint(
                ckpt_path,
                CheckpointPayload(
                    L=L,
                    seed=seed,
                    n_drops_max=n_drops_max,
                    drop=t_drop,
                    state=state,
                    rng_state=rng.bit_generator.state,
                    sizes=sizes,
                    durations=durations,
                    ever_toppled=ever_toppled,
                    snapshots=snapshots,
                ),
            )
            policy.mark_saved(t_drop)
            log.write(
                f"[{datetime.now().isoformat()}] checkpoint saved at drop={t_drop}\n"
            )

    wall = time.time() - t0

    # Final output (the artifact analysis consumes)
    np.savez_compressed(
        out_path,
        L=np.int64(L),
        seed=np.int64(seed),
        n_drops_max=np.int64(n_drops_max),
        sizes=sizes,
        durations=durations,
        ever_toppled=ever_toppled,
        final_z=state.z,
        final_p=np.float64(float(ever_toppled.mean())),
        grains_lost=np.int64(state.grains_lost),
        snapshots=np.array(snapshots, dtype=object),
    )

    result = JobResult(
        L=L,
        seed=seed,
        drops_executed=t_drop,
        final_p=float(ever_toppled.mean()),
        peak_size=int(sizes[:t_drop].max()),
        peak_pct=float(sizes[:t_drop].max()) / (L**3) * 100.0,
        wall_seconds=wall,
        z_avg_final=float(state.z.mean()),
        grains_lost=int(state.grains_lost),
    )
    log.write(
        f"[{datetime.now().isoformat()}] DONE wall={wall:.0f}s "
        f"peak={result.peak_size} ({result.peak_pct:.2f}%) "
        f"final_p={result.final_p:.5f}\n"
    )
    log.close()

    # Clean up checkpoint on successful completion (saves disk).
    if ckpt_path.exists():
        ckpt_path.unlink()

    return result


def _worker(args):
    """Pool entry point. Unpacks tuple because Pool.map prefers it."""
    return run_one(*args)


def build_job_list(
    l_list: list[int],
    seeds_per_l: dict[int, int],
    drops_for_l: dict[int, int],
    sweep_dir: Path,
) -> list[tuple]:
    jobs = []
    for L in l_list:
        for seed_idx in range(seeds_per_l[L]):
            seed = seed_for(L, seed_idx)
            jobs.append(
                (
                    L,
                    seed,
                    drops_for_l[L],
                    sweep_dir / f"L{L}_s{seed}_ckpt.npz",
                    sweep_dir / f"L{L}_s{seed}_final.npz",
                    sweep_dir / f"L{L}_s{seed}.log",
                )
            )
    return jobs


def run_single_job_mode(spec: str, sweep_dir: Path) -> None:
    """Run exactly one (L, seed) job. Used by the fleet model where each
    EC2 instance handles a single job.

    spec format: "L:seed" or "L:seed:n_drops". If n_drops omitted, uses
    DEFAULT_DROPS_FOR_L[L].
    """
    parts = spec.split(":")
    if len(parts) not in (2, 3):
        print(f"ERROR: --single-job spec must be 'L:seed' or 'L:seed:n_drops', got '{spec}'")
        sys.exit(1)
    L = int(parts[0])
    seed = int(parts[1])
    if len(parts) == 3:
        n_drops = int(parts[2])
    else:
        if L not in DEFAULT_DROPS_FOR_L:
            print(f"ERROR: L={L} has no default n_drops; pass it explicitly as 'L:seed:n_drops'")
            sys.exit(1)
        n_drops = DEFAULT_DROPS_FOR_L[L]

    sweep_dir.mkdir(parents=True, exist_ok=True)
    print(f"Single-job mode: L={L}, seed={seed}, n_drops={n_drops}")
    print(f"  output dir: {sweep_dir}")

    result = run_one(
        L=L,
        seed=seed,
        n_drops_max=n_drops,
        ckpt_path=sweep_dir / f"L{L}_s{seed}_ckpt.npz",
        out_path=sweep_dir / f"L{L}_s{seed}_final.npz",
        log_path=sweep_dir / f"L{L}_s{seed}.log",
    )

    # Write a single-job summary alongside the result.
    summary = {
        "single_job": True,
        "result": asdict(result),
    }
    with open(sweep_dir / f"L{L}_s{seed}_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(
        f"DONE L={L} seed={seed}: peak={result.peak_size} "
        f"({result.peak_pct:.2f}%) wall={result.wall_seconds:.0f}s"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--single-job",
        type=str,
        default=None,
        help="Run exactly one (L, seed) job in this process. "
             "Format: 'L:seed' or 'L:seed:n_drops'. Mutually exclusive with "
             "--l-list. Used by the fleet model.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, mp.cpu_count() - 1),
        help="parallel worker processes (default: CPU count - 1)",
    )
    parser.add_argument(
        "--l-list",
        type=str,
        default=",".join(str(L) for L in DEFAULT_L_LIST),
        help="comma-separated L values",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print the job list and exit",
    )
    parser.add_argument(
        "--sweep-dir",
        type=str,
        default=None,
        help="output dir (default: data/outputs/fss_sweep_{timestamp})",
    )
    args = parser.parse_args()

    if args.sweep_dir is None:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        sweep_dir = REPO_ROOT / "data" / "outputs" / f"fss_sweep_{stamp}"
    else:
        sweep_dir = Path(args.sweep_dir)

    if args.single_job is not None:
        run_single_job_mode(args.single_job, sweep_dir)
        return

    l_list = [int(x) for x in args.l_list.split(",")]
    for L in l_list:
        if L not in DEFAULT_SEEDS_PER_L or L not in DEFAULT_DROPS_FOR_L:
            print(
                f"ERROR: L={L} has no default seeds/drops; "
                f"edit DEFAULT_SEEDS_PER_L and DEFAULT_DROPS_FOR_L."
            )
            sys.exit(1)

    sweep_dir.mkdir(parents=True, exist_ok=True)
    jobs = build_job_list(l_list, DEFAULT_SEEDS_PER_L, DEFAULT_DROPS_FOR_L, sweep_dir)

    print(f"FSS sweep plan ({len(jobs)} jobs, {args.workers} workers)")
    print(f"  output dir: {sweep_dir}")
    print(f"{'L':>5}  {'seed':>6}  {'n_drops':>10}")
    for L, seed, n_drops, _, _, _ in jobs:
        print(f"{L:>5}  {seed:>6}  {n_drops:>10}")

    if args.dry_run:
        print("\ndry-run: not launching workers.")
        return

    print(f"\nLaunching workers at {datetime.now().isoformat()}")
    t_start = time.time()

    with mp.Pool(args.workers) as pool:
        results = pool.map(_worker, jobs)

    wall_total = time.time() - t_start
    print(f"\nAll jobs done in {wall_total:.0f}s wall.")

    summary = {
        "sweep_dir": str(sweep_dir),
        "wall_seconds_total": wall_total,
        "n_jobs": len(jobs),
        "workers": args.workers,
        "l_list": l_list,
        "results": [asdict(r) for r in results],
    }
    with open(sweep_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Summary written to {sweep_dir / 'summary.json'}")

    print(f"\n{'L':>5}  {'seed':>6}  {'peak%':>7}  {'final_p':>9}  {'wall':>7}")
    for r in results:
        print(
            f"{r.L:>5}  {r.seed:>6}  {r.peak_pct:>6.2f}%  "
            f"{r.final_p:>9.5f}  {r.wall_seconds:>6.0f}s"
        )


if __name__ == "__main__":
    main()
