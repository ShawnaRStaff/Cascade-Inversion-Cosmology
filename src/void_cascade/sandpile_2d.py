"""2D Manna stochastic sandpile.

Implements the Manna model (Manna 1991) on a square lattice with open
boundaries: a site is unstable when z >= z_c = 2; an unstable site loses
exactly 2 grains, each of which is sent independently to one of the 4
nearest neighbors chosen uniformly at random. Grains whose target lies
off-lattice are lost.

We use Manna rather than the deterministic Bak-Tang-Wiesenfeld model for
two reasons:

1. BTW in 2D exhibits multiscaling (Tebaldi, De Menech & Stella 1999):
   the avalanche-size distribution is not a clean simple-scaling power
   law, and a single-tau fit is misleading. Manna obeys simple scaling
   with tau ~ 1.275, D ~ 2.76 (Lubeck 2000, Chessa et al. 1999), so the
   FSS toolkit built for the 1D Oslo model drops in unchanged.

2. Manna sits in the same universality class as the SOC models used in
   the cosmological-percolation literature this project draws on
   (Carfora & Marzuoli 2023).

Driving
-------
A single grain is added at a uniformly random bulk site at each step.
Random driving in 2D Manna gives the same universality class as
center- or corner-driving but reaches steady state with fewer drops.

State convention
----------------
- Site indices (i, j) with i = row, j = column, both in [0, L).
- z[i, j] is an integer >= 0; stable iff z[i, j] < 2.
- All four edges are open: a grain sent off-lattice is lost.
- A toppling event at (i, j) is: z[i, j] -= 2, then two grains each go
  to one of (i-1, j), (i+1, j), (i, j-1), (i, j+1) chosen independently
  and uniformly. Both grains can land on the same neighbor.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class MannaState:
    """State of a 2D Manna sandpile.

    Attributes
    ----------
    z : ndarray[int], shape (L, L)
        Local grain count per site.
    grains_lost : int
        Cumulative grains that fell off the boundary. The conservation
        relation grains_in == sum(z) + grains_lost must always hold.
    """

    z: np.ndarray
    grains_lost: int = 0

    @property
    def L(self) -> int:
        return int(self.z.shape[0])


def initialize(L: int) -> MannaState:
    """Empty L x L lattice with z = 0 everywhere."""
    return MannaState(z=np.zeros((L, L), dtype=np.int64))


def drive(state: MannaState, rng: np.random.Generator) -> tuple[int, int]:
    """Add one grain at a uniformly random site. Returns the drive location.

    Random bulk driving is the standard Manna choice. The driven site is
    returned so that callers can record it (useful for animations).
    """
    L = state.L
    i = int(rng.integers(0, L))
    j = int(rng.integers(0, L))
    state.z[i, j] += 1
    return i, j


def _apply_sweep(
    state: MannaState,
    rng: np.random.Generator,
    toppled_mask: np.ndarray | None,
) -> int:
    """Apply one parallel update sweep. Returns the number of topplings.

    Each currently unstable site topples once, losing 2 grains; the two
    grains are sent to independently sampled nearest neighbors. Grains
    targeted off-lattice are added to state.grains_lost.

    If toppled_mask is provided, set toppled_mask[i, j] = True for every
    site that toppled this sweep (used to accumulate avalanche support).

    Implementation note: we work on the flat view of z and use
    np.bincount to sum the +1 contributions from incoming grains.
    np.bincount is ~50x faster than np.add.at on the same target sites,
    and dominates the hot loop in steady state where every drive triggers
    an avalanche with many sweeps.
    """
    z_flat = state.z.ravel()
    L = state.z.shape[0]
    N = z_flat.size

    unstable_flat = np.flatnonzero(z_flat >= 2)
    n = unstable_flat.size
    if n == 0:
        return 0

    if toppled_mask is not None:
        toppled_mask.ravel()[unstable_flat] = True

    rows = unstable_flat // L
    cols = unstable_flat % L

    # Unique indices => direct in-place subtract is safe and fast.
    z_flat[unstable_flat] -= 2

    # Two independent grains per toppling. Build the flat target index
    # for all 2n grains in one shot.
    dirs = rng.integers(0, 4, size=2 * n)
    target_rows = np.tile(rows, 2)
    target_cols = np.tile(cols, 2)
    target_rows[dirs == 0] -= 1   # up
    target_rows[dirs == 1] += 1   # down
    target_cols[dirs == 2] -= 1   # left
    target_cols[dirs == 3] += 1   # right

    on_lattice = (
        (target_rows >= 0)
        & (target_rows < L)
        & (target_cols >= 0)
        & (target_cols < L)
    )
    state.grains_lost += int(2 * n - on_lattice.sum())

    valid_flat = target_rows[on_lattice] * L + target_cols[on_lattice]
    if valid_flat.size:
        # bincount returns int64 on 64-bit numpy; minlength keeps shape stable.
        dz = np.bincount(valid_flat, minlength=N).astype(np.int64, copy=False)
        z_flat += dz

    return n


def relax(
    state: MannaState,
    rng: np.random.Generator,
    track_support: bool = False,
) -> tuple[int, int, np.ndarray | None]:
    """Topple until stable. Returns (size s, duration T, support mask or None).

    Parameters
    ----------
    track_support : bool
        If True, also return a boolean mask of every site that toppled at
        least once during this avalanche. The avalanche area a is then
        mask.sum(); the radius of gyration is computed from the masked
        coordinates. Off by default because the allocation of an L x L
        bool array per drive is a noticeable overhead for fast surveys.

    Notes
    -----
    The size s is the total number of toppling events (counted with
    multiplicity; a site that becomes unstable and topples three times in
    one avalanche contributes 3). The area a (track_support=True) counts
    each toppled site once.
    """
    L = state.L
    mask = np.zeros((L, L), dtype=bool) if track_support else None
    s = 0
    T = 0
    while True:
        n = _apply_sweep(state, rng, mask)
        if n == 0:
            break
        T += 1
        s += n
    return s, T, mask


def run(
    L: int,
    n_drops: int,
    seed: int | None = None,
) -> tuple[MannaState, np.ndarray, np.ndarray]:
    """Drive an empty Manna lattice for n_drops steps. Records s and T.

    Cluster geometry (area, radius of gyration, fractal dimension) is not
    tracked here because the per-drive support mask is expensive. Use
    run_with_clusters for that.
    """
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


def run_with_clusters(
    L: int,
    n_drops: int,
    seed: int | None = None,
    burn_in: int | None = None,
) -> tuple[MannaState, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Like run, but also records per-avalanche area and radius of gyration.

    The radius of gyration R_g is computed over the set of distinct sites
    that toppled during the avalanche:

        R_g = sqrt( <r^2> - <r>^2 )    with r = (i, j)

    Stable drives (no toppling) get area 0 and R_g 0.

    Tracking the avalanche support costs one (L, L) boolean per drive and
    a small per-avalanche reduction. For L <= 256 this is fine. We skip
    tracking during the optional burn_in to amortize the cost.

    Returns
    -------
    state, sizes, durations, areas, rgyr
        rgyr is float64 with NaN where area < 2 (R_g undefined for a
        single point in the conventional definition; callers should mask).
    """
    rng = np.random.default_rng(seed)
    state = initialize(L)
    sizes = np.zeros(n_drops, dtype=np.int64)
    durations = np.zeros(n_drops, dtype=np.int64)
    areas = np.zeros(n_drops, dtype=np.int64)
    rgyr = np.full(n_drops, np.nan, dtype=np.float64)

    transient = 0 if burn_in is None else int(burn_in)
    for t in range(n_drops):
        drive(state, rng)
        # During burn-in we don't need cluster geometry, save the alloc.
        track = t >= transient
        s, T, mask = relax(state, rng, track_support=track)
        sizes[t] = s
        durations[t] = T
        if mask is not None:
            a = int(mask.sum())
            areas[t] = a
            if a >= 2:
                idx = np.argwhere(mask).astype(np.float64)
                mean = idx.mean(axis=0)
                var = ((idx - mean) ** 2).sum(axis=1).mean()
                rgyr[t] = float(np.sqrt(var))
            elif a == 1:
                rgyr[t] = 0.0
    return state, sizes, durations, areas, rgyr
