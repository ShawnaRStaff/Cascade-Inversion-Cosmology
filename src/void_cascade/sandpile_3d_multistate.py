"""Multi-state 3D Manna sandpile (Milestone 5).

The Manna model from M3 had a single state per cell (cracked / frozen).
This module extends it: each cell has an integer 'matter type' sigma in
{0, 1, 2, ...} representing the kind of substrate / matter occupying
that cell. The cosmology motivation is that real cosmic history is a
chain of new materials forming through cascading high-pressure events
(quarks -> nucleons -> nuclei -> atoms -> molecules ...), each with
greater internal stability than the previous tier.

State per cell:
- z : grain count (int)
- sigma : matter-type tier (int, in {0, 1, ...})

Toppling rule:
- Threshold scales with sigma: z_c(sigma) = 2 * (sigma + 1)
  (sigma=0 cells crack at z=2 - same as M3 Manna; higher sigma needs
   proportionally more grains, the 'binding energy' grows with tier.)
- When a cell topples it sheds exactly z_c(sigma) grains to randomly
  chosen face neighbors. Same stochastic neighbor pick as Manna.
- After toppling, sigma is updated:
    if at least one face neighbor also toppled in the same sweep:
        new_sigma = max(self_sigma, max_active_neighbor_sigma + 1)
    else:
        sigma unchanged
  This models 'fusion under heat and pressure' - matter only transforms
  to a higher tier when multiple cascading sites combine, not from
  isolated topplings.

Bounded variant:
- sigma capped at sigma_max. Analog of 'no element above iron through
  pure fusion; you need supernova conditions'. Tier ceiling.

Unbounded variant:
- sigma grows without bound. Higher tiers require more grains to crack,
  so the dynamics naturally slows as sigma accumulates.

Driving and dissipation: same as M3 Manna - bulk random drop, open
boundaries.

Open question for both variants: does the dynamics naturally arrest
at some fractured fraction p? If yes, at what p, and how does it
compare to the M4 'galaxy clustering matches' value p ~ 0.65?
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class MannaStateMS:
    """State of a 3D multi-state Manna sandpile.

    Attributes
    ----------
    z : ndarray[int], shape (L, L, L)
        Grain count per cell.
    sigma : ndarray[int], shape (L, L, L)
        Matter-type tier per cell. sigma=0 is the pristine substrate.
    grains_lost : int
        Cumulative grains lost off the boundary.
    """

    z: np.ndarray
    sigma: np.ndarray
    grains_lost: int = 0

    @property
    def L(self) -> int:
        return int(self.z.shape[0])


def initialize(L: int) -> MannaStateMS:
    """Empty L x L x L lattice with z = 0, sigma = 0 everywhere."""
    return MannaStateMS(
        z=np.zeros((L, L, L), dtype=np.int64),
        sigma=np.zeros((L, L, L), dtype=np.int64),
    )


def drive(state: MannaStateMS, rng: np.random.Generator) -> tuple[int, int, int]:
    """Add one grain at a uniformly random site."""
    L = state.L
    i = int(rng.integers(0, L))
    j = int(rng.integers(0, L))
    k = int(rng.integers(0, L))
    state.z[i, j, k] += 1
    return i, j, k


def _apply_sweep(
    state: MannaStateMS,
    rng: np.random.Generator,
    toppled_mask: np.ndarray | None = None,
    sigma_max: int | None = None,
) -> int:
    """Apply one parallel update sweep with multi-state transformation.

    Returns the number of topplings this sweep.
    """
    z_flat = state.z.ravel()
    sigma_flat = state.sigma.ravel()
    L = state.z.shape[0]
    N = z_flat.size
    LL = L * L

    # Per-cell threshold z_c = 2*(sigma+1)
    z_c_per_cell = 2 * (sigma_flat + 1)
    unstable_mask_flat = z_flat >= z_c_per_cell
    unstable_flat_idx = np.flatnonzero(unstable_mask_flat)
    n = unstable_flat_idx.size
    if n == 0:
        return 0

    if toppled_mask is not None:
        toppled_mask.ravel()[unstable_flat_idx] = True

    # Per-cell z_c for the unstable cells we are about to topple
    z_c_unstable = z_c_per_cell[unstable_flat_idx]
    total_grains = int(z_c_unstable.sum())

    # Subtract the appropriate z_c from each unstable cell
    z_flat[unstable_flat_idx] -= z_c_unstable

    # Expand to per-grain arrays: grain_origin[k] = flat index of cell
    # that emitted grain k. np.repeat handles variable z_c per cell.
    grain_origin = np.repeat(unstable_flat_idx, z_c_unstable)
    rows = grain_origin // LL
    rem = grain_origin - rows * LL
    cols = rem // L
    layers = rem - cols * L

    # Random direction per grain (6 face directions in 3D)
    dirs = rng.integers(0, 6, size=total_grains)
    ti = rows.copy()
    tj = cols.copy()
    tk = layers.copy()
    ti[dirs == 0] += 1
    ti[dirs == 1] -= 1
    tj[dirs == 2] += 1
    tj[dirs == 3] -= 1
    tk[dirs == 4] += 1
    tk[dirs == 5] -= 1

    on_lattice = (
        (ti >= 0) & (ti < L)
        & (tj >= 0) & (tj < L)
        & (tk >= 0) & (tk < L)
    )
    state.grains_lost += int(total_grains - on_lattice.sum())

    valid_flat = ti[on_lattice] * LL + tj[on_lattice] * L + tk[on_lattice]
    if valid_flat.size:
        dz = np.bincount(valid_flat, minlength=N).astype(np.int64, copy=False)
        z_flat += dz

    # --- sigma update ---
    # Cells that toppled AND have at least one face neighbor also
    # toppling this sweep get transformed:
    #   sigma_new = max(self_sigma, max_neighbor_sigma_among_active + 1)
    # Implemented via padded shifts so off-lattice neighbors contribute
    # nothing (-1 sentinel).
    unstable_3d = unstable_mask_flat.reshape(L, L, L)
    sigma_3d = state.sigma

    sig_pad = np.pad(sigma_3d, 1, constant_values=-1)
    u_pad = np.pad(unstable_3d, 1, constant_values=False)

    # has_active_neighbor[i,j,k] = OR over 6 face neighbors of u_pad
    has_active_neighbor = (
        u_pad[0:-2, 1:-1, 1:-1]
        | u_pad[2:,   1:-1, 1:-1]
        | u_pad[1:-1, 0:-2, 1:-1]
        | u_pad[1:-1, 2:,   1:-1]
        | u_pad[1:-1, 1:-1, 0:-2]
        | u_pad[1:-1, 1:-1, 2:]
    )

    # For each face direction, neighbor's sigma if neighbor is active,
    # else -1 sentinel. max over directions gives the highest active-
    # neighbor sigma (or -1 if no active neighbor).
    def directional(slc_axis: int, shift_idx: int) -> np.ndarray:
        slcs = [slice(1, -1), slice(1, -1), slice(1, -1)]
        slcs[slc_axis] = slice(0, -2) if shift_idx == 0 else slice(2, None)
        slcs_tuple = tuple(slcs)
        sig_n = sig_pad[slcs_tuple]
        u_n = u_pad[slcs_tuple]
        return np.where(u_n, sig_n, -1)

    max_neighbor_active = np.maximum.reduce([
        directional(a, s) for a in range(3) for s in (0, 1)
    ])

    transform = unstable_3d & has_active_neighbor
    new_sigma = sigma_3d.copy()
    new_sigma[transform] = np.maximum(
        sigma_3d[transform],
        max_neighbor_active[transform] + 1,
    )
    if sigma_max is not None:
        new_sigma = np.minimum(new_sigma, sigma_max)
    state.sigma = new_sigma

    return n


def relax(
    state: MannaStateMS,
    rng: np.random.Generator,
    track_support: bool = False,
    sigma_max: int | None = None,
) -> tuple[int, int, np.ndarray | None]:
    """Topple until stable. Returns (size s, duration T, support mask or None).

    Multi-state version. sigma_max=None means unbounded; otherwise cap.
    """
    L = state.L
    mask = np.zeros((L, L, L), dtype=bool) if track_support else None
    s = 0
    T = 0
    while True:
        n = _apply_sweep(state, rng, toppled_mask=mask, sigma_max=sigma_max)
        if n == 0:
            break
        T += 1
        s += n
    return s, T, mask


def run_with_ever_toppled(
    L: int,
    n_drops: int,
    seed: int | None = None,
    sigma_max: int | None = None,
    check_every: int = 1,
    callback=None,
) -> tuple[MannaStateMS, np.ndarray, np.ndarray, np.ndarray]:
    """Drive n_drops and maintain a cumulative ever_toppled mask.

    Returns (state, sizes, durations, ever_toppled).
    `callback(t, ever_toppled, state, sizes, durations)` may return True
    to halt early.
    """
    rng = np.random.default_rng(seed)
    state = initialize(L)
    sizes = np.zeros(n_drops, dtype=np.int64)
    durations = np.zeros(n_drops, dtype=np.int64)
    ever_toppled = np.zeros((L, L, L), dtype=bool)
    for t in range(n_drops):
        drive(state, rng)
        s, T, mask = relax(state, rng, track_support=True, sigma_max=sigma_max)
        sizes[t] = s
        durations[t] = T
        if mask is not None:
            ever_toppled |= mask
        if callback is not None and (t + 1) % check_every == 0:
            if callback(t, ever_toppled, state, sizes, durations):
                sizes = sizes[: t + 1]
                durations = durations[: t + 1]
                break
    return state, sizes, durations, ever_toppled
