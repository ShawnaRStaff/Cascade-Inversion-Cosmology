"""Visualization for the 2D Manna sandpile.

Two outputs:

1. `record_animation` — GIF showing the lattice during a small simulation.
   Each frame is one parallel update sweep; unstable cells are highlighted
   so you can see the cascade propagate. Written to disk; never displayed
   inline because the GIFs can be tens of megabytes for any meaningful
   duration.

2. `snapshot_large_avalanche` — runs until it captures a near-spanning
   avalanche, then writes a static PNG of the avalanche's support overlaid
   on the final lattice state. Useful for showing a single cluster's
   fractal structure.

All paths are returned so the caller can print them; nothing is shown
in-process.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from .sandpile_2d import MannaState, _apply_sweep, drive, initialize


def record_animation(
    L: int = 32,
    max_frames: int = 400,
    seed: int = 7,
    out_path: str | Path = "data/outputs/manna_2d_animation.gif",
    fps: int = 15,
) -> Path:
    """Run a Manna simulation and save a GIF of the lattice dynamics.

    Frame schedule: one frame immediately after each drive (so you can
    see the new grain), then one frame per parallel update sweep during
    the avalanche. The driven cell is briefly tinted in the post-drive
    frame; unstable cells in mid-avalanche frames are tinted red.

    Parameters
    ----------
    L : int
        Lattice side. Keep small (16-48) so individual cells are visible.
    max_frames : int
        Hard cap on captured frames.
    seed : int
    out_path : Path
    fps : int
    """
    rng = np.random.default_rng(seed)
    state = initialize(L)

    frames: list[dict] = []
    drop_idx = 0
    total_topples = 0

    while len(frames) < max_frames:
        i, j = drive(state, rng)
        drop_idx += 1
        sweep_idx = 0
        size_running = 0
        just_dropped_mask = np.zeros((L, L), dtype=bool)
        just_dropped_mask[i, j] = True

        frames.append(
            {
                "z": state.z.copy(),
                "highlight": just_dropped_mask,
                "highlight_kind": "drive",
                "drop": drop_idx,
                "sweep": sweep_idx,
                "size_running": size_running,
                "total_topples": total_topples,
            }
        )
        if len(frames) >= max_frames:
            break

        while True:
            unstable_mask = state.z >= 2
            if not unstable_mask.any():
                break
            sweep_idx += 1
            n_unstable = int(unstable_mask.sum())
            size_running += n_unstable
            total_topples += n_unstable

            # Snapshot the unstable set BEFORE applying the sweep so the
            # highlight matches what's about to topple.
            highlight = unstable_mask.copy()
            _apply_sweep(state, rng, toppled_mask=None)

            frames.append(
                {
                    "z": state.z.copy(),
                    "highlight": highlight,
                    "highlight_kind": "topple",
                    "drop": drop_idx,
                    "sweep": sweep_idx,
                    "size_running": size_running,
                    "total_topples": total_topples,
                }
            )
            if len(frames) >= max_frames:
                break

    z_max = max(int(f["z"].max()) for f in frames)
    z_cap = max(z_max, 3)

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    base = frames[0]["z"]
    rgb = _render_frame(base, frames[0]["highlight"], frames[0]["highlight_kind"], z_cap)
    img = ax.imshow(rgb, interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    title = ax.set_title("")

    def update(idx):
        f = frames[idx]
        img.set_data(_render_frame(f["z"], f["highlight"], f["highlight_kind"], z_cap))
        if f["highlight_kind"] == "drive":
            descr = f"drop #{f['drop']}: grain added"
        else:
            descr = (
                f"drop #{f['drop']}, sweep {f['sweep']}: "
                f"avalanche size so far = {f['size_running']}"
            )
        title.set_text(f"{descr}  |  total topples {f['total_topples']}")
        return [img, title]

    anim = FuncAnimation(
        fig, update, frames=len(frames), interval=int(1000 / fps), blit=False
    )

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    writer = PillowWriter(fps=fps)
    anim.save(str(out_path), writer=writer)
    plt.close(fig)
    return out_path


def _render_frame(
    z: np.ndarray, highlight: np.ndarray, kind: str, z_cap: int
) -> np.ndarray:
    """Compose an RGB image: grayscale by z, color overlay on highlighted cells.

    Pre-rendering avoids matplotlib re-color-mapping per frame, which is
    slow with PillowWriter; an RGB array drops straight into imshow.
    """
    L = z.shape[0]
    norm = np.clip(z.astype(np.float64) / max(z_cap, 1), 0.0, 1.0)
    gray = 0.15 + 0.65 * norm  # 0.15 (empty) to 0.80 (capped)
    rgb = np.stack([gray, gray, gray], axis=-1)
    if highlight.any():
        if kind == "drive":
            color = np.array([0.20, 0.55, 0.90])  # blue: new grain
        else:
            color = np.array([0.90, 0.25, 0.25])  # red: unstable
        rgb[highlight] = color
    return (rgb * 255).astype(np.uint8)


def snapshot_large_avalanche(
    L: int = 64,
    burn_in: int | None = None,
    min_area: int | None = None,
    max_search_drops: int = 20_000,
    seed: int = 11,
    out_path: str | Path = "data/outputs/manna_2d_avalanche_snapshot.png",
) -> Path:
    """Save a PNG of a near-spanning avalanche cluster overlaid on z.

    The driver runs through a burn-in to steady state, then keeps driving
    until an avalanche whose support area exceeds min_area is observed.
    A min_area of None defaults to L^2 / 16 (a meaningful structure but
    much smaller than the lattice itself).

    Raises
    ------
    RuntimeError
        If no qualifying avalanche occurs within max_search_drops.
    """
    from .sandpile_2d import relax  # late import to avoid cycle on Path-only use

    rng = np.random.default_rng(seed)
    state = initialize(L)
    if burn_in is None:
        burn_in = 2 * L * L
    if min_area is None:
        min_area = max(L * L // 16, 16)

    # Burn-in without tracking support.
    for _ in range(burn_in):
        drive(state, rng)
        relax(state, rng, track_support=False)

    found_mask: np.ndarray | None = None
    found_area = 0
    drops_done = 0
    while drops_done < max_search_drops:
        drive(state, rng)
        _s, _T, mask = relax(state, rng, track_support=True)
        drops_done += 1
        if mask is not None:
            a = int(mask.sum())
            if a >= min_area:
                found_mask = mask
                found_area = a
                break
    if found_mask is None:
        raise RuntimeError(
            f"No avalanche with area >= {min_area} in {max_search_drops} drops "
            f"after burn-in. Increase max_search_drops or lower min_area."
        )

    fig, ax = plt.subplots(figsize=(6, 6))
    z_view = state.z.astype(np.float64)
    ax.imshow(z_view, cmap="Greys", interpolation="nearest", vmin=0, vmax=2)
    overlay = np.zeros((L, L, 4))
    overlay[found_mask] = (0.85, 0.18, 0.18, 0.75)
    ax.imshow(overlay, interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(
        f"2D Manna L={L}: avalanche support area = {found_area} "
        f"({found_area / (L * L):.1%} of lattice)"
    )
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
