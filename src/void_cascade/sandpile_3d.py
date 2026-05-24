"""3D Manna stochastic sandpile on a cubic lattice with open boundaries.

Direct extension of the 2D Manna model to three dimensions. A site (i, j, k)
is unstable when z >= z_c = 2; an unstable site loses 2 grains, each of
which is sent independently to one of the 6 nearest neighbors (+/-x, +/-y,
+/-z) chosen uniformly at random. Grains targeted off-lattice are lost.

Reference values for 3D Manna (Lubeck 2000, Pastor-Satorras & Vespignani
2000): tau ~ 1.35, D ~ 3.36. We do not perform a full FSS sweep in this
project; the cosmologically relevant observable is the percolation of the
ever-toppled set, not the steady-state moment exponents. The exponents
are reproducible from this code if needed but the runs are expensive.

Cosmological interpretation (per docs/model_summary.md): each toppling
event is read as a fracture in the substrate, and the union of all sites
that have toppled at least once is the cumulative fractured region. The
percolation transition - when that region first contains a connected
cluster spanning the box - is the cosmologically meaningful event.

State convention
----------------
- Site indices (i, j, k) with each in [0, L).
- z[i, j, k] is an integer >= 0; stable iff z < 2.
- All six faces are open: grains sent off-lattice are lost.
- Toppling: z[i,j,k] -= 2, then two grains each go to one of the 6
  neighbors {(i+/-1,j,k), (i,j+/-1,k), (i,j,k+/-1)} chosen independently
  and uniformly. Both grains can land on the same neighbor.

Implementation note
-------------------
We work on the flat view of z and use np.bincount to sum the +1
contributions from incoming grains. Same optimization as in 2D; it is
~50x faster than np.add.at on the same target sites.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class MannaState3D:
    """State of a 3D Manna sandpile.

    Attributes
    ----------
    z : ndarray[int], shape (L, L, L)
        Local grain count per site.
    grains_lost : int
        Cumulative grains that fell off the boundary. Conservation:
        grains_in == sum(z) + grains_lost must always hold.
    """

    z: np.ndarray
    grains_lost: int = 0

    @property
    def L(self) -> int:
        return int(self.z.shape[0])


def initialize(L: int) -> MannaState3D:
    """Empty L x L x L lattice with z = 0 everywhere."""
    return MannaState3D(z=np.zeros((L, L, L), dtype=np.int64))


def drive(state: MannaState3D, rng: np.random.Generator) -> tuple[int, int, int]:
    """Add one grain at a uniformly random site. Returns the drive location."""
    L = state.L
    i = int(rng.integers(0, L))
    j = int(rng.integers(0, L))
    k = int(rng.integers(0, L))
    state.z[i, j, k] += 1
    return i, j, k


def _apply_sweep(
    state: MannaState3D,
    rng: np.random.Generator,
    toppled_mask: np.ndarray | None,
) -> int:
    """Apply one parallel update sweep. Returns the number of topplings.

    Each currently unstable site topples once, losing 2 grains; the two
    grains are sent to independently sampled nearest neighbors among the
    6 directions {+/-x, +/-y, +/-z}. Grains targeted off-lattice are
    added to state.grains_lost.

    If toppled_mask is provided, set toppled_mask[i, j, k] = True for
    every site that toppled this sweep (used to accumulate avalanche
    support and, more importantly for the cosmology mapping, the
    ever-toppled set used by the percolation detector).
    """
    z_flat = state.z.ravel()
    L = state.z.shape[0]
    N = z_flat.size
    LL = L * L

    unstable_flat = np.flatnonzero(z_flat >= 2)
    n = unstable_flat.size
    if n == 0:
        return 0

    if toppled_mask is not None:
        toppled_mask.ravel()[unstable_flat] = True

    # Flat -> (i, j, k) decode. flat = i*LL + j*L + k.
    i_idx = unstable_flat // LL
    rem = unstable_flat - i_idx * LL
    j_idx = rem // L
    k_idx = rem - j_idx * L

    # Unique indices => direct in-place subtract is safe and fast.
    z_flat[unstable_flat] -= 2

    # Two independent grains per toppling.
    dirs = rng.integers(0, 6, size=2 * n)
    ti = np.tile(i_idx, 2)
    tj = np.tile(j_idx, 2)
    tk = np.tile(k_idx, 2)
    # 0: +x, 1: -x, 2: +y, 3: -y, 4: +z, 5: -z
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
    state.grains_lost += int(2 * n - on_lattice.sum())

    valid_flat = ti[on_lattice] * LL + tj[on_lattice] * L + tk[on_lattice]
    if valid_flat.size:
        dz = np.bincount(valid_flat, minlength=N).astype(np.int64, copy=False)
        z_flat += dz

    return n


def relax(
    state: MannaState3D,
    rng: np.random.Generator,
    track_support: bool = False,
) -> tuple[int, int, np.ndarray | None]:
    """Topple until stable. Returns (size s, duration T, support mask or None).

    `track_support=True` accumulates a (L, L, L) bool mask of every site
    that toppled at least once during this avalanche. The cost is the
    per-drive allocation of an L^3 bool, which is the dominant overhead
    at large L; only enable it when the caller needs it.
    """
    L = state.L
    mask = np.zeros((L, L, L), dtype=bool) if track_support else None
    s = 0
    T = 0
    while True:
        n = _apply_sweep(state, rng, mask)
        if n == 0:
            break
        T += 1
        s += n
    return s, T, mask


def relax_sequential(
    state: MannaState3D,
    rng: np.random.Generator,
    track_support: bool = False,
) -> tuple[int, int, np.ndarray | None]:
    """Sequential (one-at-a-time) toppling, matching Huynh & Pruessner 2012.

    Picks ONE unstable site uniformly at random, topples it (z -= 2, two
    grains distributed independently and uniformly to its 6 neighbors,
    possibly the same neighbor twice). Repeats until no unstable sites.

    This is the canonical Abelian Manna model dynamics used in the
    published literature. Our default `relax()` uses parallel updates
    (all unstable sites topple in one sweep). Abelian property
    guarantees identical FINAL state and identical unique-cells-toppled
    set, but the SEQUENCE differs — so dynamic exponents (avalanche
    size distribution P(s), area distribution P(a), duration) can
    differ between the two implementations.

    Use this function to verify our model is in the same universality
    class as the published literature.
    """
    L = state.L
    mask = np.zeros((L, L, L), dtype=bool) if track_support else None
    s = 0
    T = 0
    while True:
        # Find unstable sites.
        unstable_flat = np.flatnonzero(state.z.ravel() >= 2)
        if unstable_flat.size == 0:
            break
        # Pick ONE at random.
        idx = unstable_flat[int(rng.integers(0, unstable_flat.size))]
        LL = L * L
        i = idx // LL
        rem = idx - i * LL
        j = rem // L
        k = rem - j * L
        # Topple.
        state.z[i, j, k] -= 2
        if mask is not None:
            mask[i, j, k] = True
        # Distribute 2 grains to independently sampled neighbors.
        for _ in range(2):
            d = int(rng.integers(0, 6))
            ni, nj, nk = i, j, k
            if d == 0:
                ni += 1
            elif d == 1:
                ni -= 1
            elif d == 2:
                nj += 1
            elif d == 3:
                nj -= 1
            elif d == 4:
                nk += 1
            else:
                nk -= 1
            if 0 <= ni < L and 0 <= nj < L and 0 <= nk < L:
                state.z[ni, nj, nk] += 1
            else:
                state.grains_lost += 1
        s += 1
        T += 1
    return s, T, mask


def run(
    L: int,
    n_drops: int,
    seed: int | None = None,
) -> tuple[MannaState3D, np.ndarray, np.ndarray]:
    """Drive an empty cubic lattice for n_drops steps. Records s and T."""
    rng = np.random.default_rng(seed)
    state = initialize(L)
    sizes = np.zeros(n_drops, dtype=np.int64)
    durations = np.zeros(n_drops, dtype=np.int64)
    for t in range(n_drops):
        drive(state, rng)
        s, T, _ = relax(state, rng, track_support=False)
        sizes[t] = s
        durations[t] = T
    return state, sizes, durations


def run_with_ever_toppled(
    L: int,
    n_drops: int,
    seed: int | None = None,
    check_every: int = 1,
    percolation_callback=None,
) -> tuple[MannaState3D, np.ndarray, np.ndarray, np.ndarray]:
    """Drive the lattice and maintain a cumulative `ever_toppled` boolean field.

    This is the cosmologically relevant driver: the union of all sites
    that have toppled at least once across the entire history of drives
    is the candidate "fractured region." A separate percolation detector
    (see percolation.py) checks whether that region spans the box.

    `percolation_callback`, if provided, is called as
        percolation_callback(drop_idx, ever_toppled, sizes, durations)
    every `check_every` drives. Return True from the callback to halt
    the run early (used by the percolation-time finder).

    Returns
    -------
    state, sizes, durations, ever_toppled
        ever_toppled is a (L, L, L) bool ndarray.
    """
    rng = np.random.default_rng(seed)
    state = initialize(L)
    sizes = np.zeros(n_drops, dtype=np.int64)
    durations = np.zeros(n_drops, dtype=np.int64)
    ever_toppled = np.zeros((L, L, L), dtype=bool)
    for t in range(n_drops):
        drive(state, rng)
        s, T, mask = relax(state, rng, track_support=True)
        sizes[t] = s
        durations[t] = T
        if mask is not None:
            # Union into the cumulative set. The per-avalanche mask is
            # discarded after this OR; we only retain the cumulative one.
            ever_toppled |= mask
        if percolation_callback is not None and (t + 1) % check_every == 0:
            if percolation_callback(t, ever_toppled, sizes, durations):
                # Trim arrays to the actually-completed range.
                sizes = sizes[: t + 1]
                durations = durations[: t + 1]
                break
    return state, sizes, durations, ever_toppled
