"""Animation of 1D Oslo sandpile dynamics.

Captures one frame per parallel update sweep, so individual topplings are
visible as red bars in the height profile. The animation is saved as a GIF
because the environment doesn't have ffmpeg; that limits the practical frame
count to a few hundred but is enough to see the pile build up and the first
several avalanches.

Why we duplicate the relaxation logic here rather than instrumenting the
production `relax()` with a frame callback: the production path is on the hot
loop of every simulation and shouldn't carry a hook for a viz-only feature.
The duplicated logic is a few lines and is exercised by the same toppling
tests via the shared boundary rules.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from .sandpile_1d import OsloState, drive, initialize


def heights_from_slopes(z: np.ndarray) -> np.ndarray:
    """Reconstruct h_i = sum_{j>=i} z_j (open right boundary, h_L = 0)."""
    return np.cumsum(z[::-1])[::-1]


def _apply_one_sweep(state: OsloState, rng: np.random.Generator) -> np.ndarray:
    """Apply one parallel update sweep. Returns the indices that toppled."""
    z = state.z
    z_c = state.z_c
    L = z.size
    unstable = np.flatnonzero(z > z_c)
    if unstable.size == 0:
        return unstable

    if unstable.size == 1:
        i = int(unstable[0])
        if i == L - 1:
            z[i] -= 1
            z[i - 1] += 1
            state.grains_lost += 1
        elif i == 0:
            z[0] -= 2
            z[1] += 1
        else:
            z[i] -= 2
            z[i - 1] += 1
            z[i + 1] += 1
        z_c[i] = int(rng.integers(1, 3))
        return unstable

    dz = np.zeros(L, dtype=np.int64)
    right_edge_mask = unstable == (L - 1)
    interior = unstable[~right_edge_mask]
    right_edge = unstable[right_edge_mask]
    np.add.at(dz, interior, -2)
    np.add.at(dz, right_edge, -1)
    np.add.at(dz, unstable[unstable > 0] - 1, 1)
    np.add.at(dz, unstable[unstable < L - 1] + 1, 1)
    state.grains_lost += int(right_edge.size)
    z += dz
    z_c[unstable] = rng.integers(1, 3, size=unstable.size, dtype=np.int64)
    return unstable


def record_animation(
    L: int = 24,
    max_frames: int = 600,
    seed: int = 42,
    out_path: str | Path = "data/outputs/oslo_1d_animation.gif",
    fps: int = 20,
) -> Path:
    """Run an Oslo simulation and save a GIF of the height-profile dynamics.

    Frames are captured at three points: just after each grain is added (so
    the new grain is visible at site 0), and once per parallel update sweep
    during the resulting avalanche (so each cascade step is visible).

    Parameters
    ----------
    L : int
        Lattice size. Keep small (16-32) so the average steady-state avalanche
        fits in a viewable number of frames.
    max_frames : int
        Hard cap on captured frames. The simulation stops capturing once this
        many frames have been recorded, even mid-avalanche.
    seed : int
        Random seed.
    out_path : path
        Output GIF path. Parent directory is created if missing.
    fps : int
        GIF playback rate. 20 fps reads as a smooth avalanche; 10 fps reads
        as deliberate, easier to follow individual topplings.
    """
    rng = np.random.default_rng(seed)
    state = initialize(L, rng)

    frames: list[dict] = []
    drop_idx = 0
    total_topples = 0

    while len(frames) < max_frames:
        drive(state)
        drop_idx += 1
        sweep_idx = 0
        size_running = 0

        frames.append(
            {
                "h": heights_from_slopes(state.z),
                "unstable_mask": np.zeros(L, dtype=bool),
                "drop": drop_idx,
                "sweep": sweep_idx,
                "size": size_running,
                "total_topples": total_topples,
                "just_dropped": True,
            }
        )
        if len(frames) >= max_frames:
            break

        while True:
            unstable = np.flatnonzero(state.z > state.z_c)
            if unstable.size == 0:
                break
            sweep_idx += 1
            size_running += int(unstable.size)
            total_topples += int(unstable.size)

            mask_before = np.zeros(L, dtype=bool)
            mask_before[unstable] = True

            _apply_one_sweep(state, rng)

            frames.append(
                {
                    "h": heights_from_slopes(state.z),
                    "unstable_mask": mask_before,
                    "drop": drop_idx,
                    "sweep": sweep_idx,
                    "size": size_running,
                    "total_topples": total_topples,
                    "just_dropped": False,
                }
            )
            if len(frames) >= max_frames:
                break

    h_max = max(int(f["h"].max()) for f in frames)
    y_top = max(h_max + 2, 8)

    fig, ax = plt.subplots(figsize=(8, 3.5))
    x = np.arange(L)
    bars = ax.bar(x, frames[0]["h"], color="#3a6ea5", edgecolor="black", linewidth=0.4)
    ax.set_xlim(-0.5, L - 0.5)
    ax.set_ylim(0, y_top)
    ax.set_xlabel("site $i$")
    ax.set_ylabel("height $h_i$")
    ax.grid(True, axis="y", alpha=0.3)
    title = ax.set_title("")

    def update(frame_idx):
        f = frames[frame_idx]
        h = f["h"]
        for bar, hi, unst in zip(bars, h, f["unstable_mask"]):
            bar.set_height(hi)
            bar.set_color("#d64545" if unst else "#3a6ea5")
        if f["just_dropped"]:
            descr = f"drop #{f['drop']}: grain at left wall"
        elif f["sweep"] > 0:
            descr = (
                f"drop #{f['drop']}, sweep {f['sweep']}: "
                f"avalanche size so far = {f['size']}"
            )
        else:
            descr = f"drop #{f['drop']}: stable"
        title.set_text(f"{descr}   |   total topples: {f['total_topples']}")
        return list(bars) + [title]

    anim = FuncAnimation(
        fig, update, frames=len(frames), interval=int(1000 / fps), blit=False
    )

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    writer = PillowWriter(fps=fps)
    anim.save(str(out_path), writer=writer)
    plt.close(fig)
    return out_path
