"""1D Oslo rice-pile sandpile model.

Implements the Oslo model (Christensen, Corral, Frette, Feder, Jossang 1996),
a stochastic-threshold sandpile that exhibits genuine 1D self-organized
criticality with avalanche size exponent tau ~ 1.55. The deterministic 1D BTW
model is not critical (its steady state is trivial), so Oslo is the canonical
choice for 1D SOC. The only difference from BTW is a per-site threshold z_c
that is redrawn uniformly from {1, 2} after every toppling at that site; that
quenched-then-refreshed disorder breaks the deterministic propagation pattern
that kills 1D BTW.

State convention
----------------
- Sites are indexed i = 0, ..., L-1 (left-to-right).
- z[i] is the local slope z_i = h_i - h_{i+1}, with h_L = 0 (open right edge).
- Driving adds grains at i = 0 (the closed left wall).
- Toppling rule for an interior site i: z[i] -= 2, z[i-1] += 1, z[i+1] += 1.
- At i = 0 (wall): z[0] -= 2, z[1] += 1; no left neighbor exists.
- At i = L-1 (open edge): z[L-1] -= 1, z[L-2] += 1; one grain leaves the system.
- After every toppling, the toppling site draws a fresh z_c from {1, 2}.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class OsloState:
    """State of a 1D Oslo rice-pile.

    Attributes
    ----------
    z : ndarray[int], shape (L,)
        Local slopes z_i = h_i - h_{i+1}.
    z_c : ndarray[int], shape (L,)
        Per-site critical slope thresholds, each in {1, 2}.
    grains_lost : int
        Cumulative count of grains that have fallen off the right edge.
        Useful for conservation checks: grains_in - grains_lost == sum_i h_i,
        where h_i = sum_{j>=i} z_j.
    """

    z: np.ndarray
    z_c: np.ndarray
    grains_lost: int = 0

    @property
    def L(self) -> int:
        return int(self.z.size)


def initialize(L: int, rng: np.random.Generator) -> OsloState:
    """Empty pile with random thresholds in {1, 2}.

    All slopes start at 0. Thresholds are drawn uniformly from {1, 2}; this is
    the standard Oslo choice (probability 1/2 each). Changing the threshold
    probability shifts the steady-state mean slope but does not change the
    universality class.
    """
    z = np.zeros(L, dtype=np.int64)
    z_c = rng.integers(1, 3, size=L, dtype=np.int64)
    return OsloState(z=z, z_c=z_c)


def drive(state: OsloState) -> None:
    """Add one grain at the left wall: z[0] += 1.

    In rice-pile terms this is h_0 -> h_0 + 1, which raises the wall slope by
    one and leaves all other slopes unchanged.
    """
    state.z[0] += 1


def relax(state: OsloState, rng: np.random.Generator) -> tuple[int, int]:
    """Topple unstable sites until the lattice is stable.

    Uses parallel updates: at each sweep, every currently unstable site topples
    once and then draws a fresh threshold. The avalanche duration T is the
    number of parallel sweeps, the avalanche size s is the total number of
    individual topplings.

    Parallel updates are standard for the Oslo model (Christensen et al. 1996)
    and give a well-defined T. Sequential updates would give the same size s
    and the same scaling exponents but a less useful definition of T.

    Returns
    -------
    s : int
        Total number of topplings in this avalanche.
    T : int
        Number of parallel update sweeps (0 if no toppling occurred).
    """
    z = state.z
    z_c = state.z_c
    L = z.size

    s = 0
    T = 0
    # Scratch buffer for the vectorized path. Allocated lazily so that
    # avalanches that only ever have one site unstable per sweep (the common
    # case in 1D Oslo) never pay the O(L) allocation cost.
    dz: np.ndarray | None = None

    while True:
        unstable = np.flatnonzero(z > z_c)
        n_unstable = unstable.size
        if n_unstable == 0:
            break
        T += 1
        s += int(n_unstable)

        if n_unstable == 1:
            # Fast path. Almost every parallel sweep in 1D Oslo has a single
            # unstable site, and the np.add.at overhead dominates if we use
            # the vectorized path here. Doing scalar updates avoids the
            # ~microseconds-per-call cost of np.add.at on size-1 arrays.
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
            continue

        # Vectorized path for sweeps with multiple unstable sites. We use
        # np.add.at rather than `dz[idx] += 1` because two unstable sites
        # separated by one stable site both push a grain to the site
        # between them, and plain indexed += would silently drop one of
        # the contributions when the same target index appears twice.
        if dz is None:
            dz = np.zeros(L, dtype=np.int64)
        else:
            dz.fill(0)

        at_right_edge_mask = unstable == (L - 1)
        interior_topples = unstable[~at_right_edge_mask]
        right_edge_topples = unstable[at_right_edge_mask]

        # Toppling site loses grains: -2 in the bulk, -1 at the open edge.
        np.add.at(dz, interior_topples, -2)
        np.add.at(dz, right_edge_topples, -1)

        # Left neighbor gains 1, except for site 0 which has no left neighbor.
        left_targets = unstable[unstable > 0] - 1
        np.add.at(dz, left_targets, 1)

        # Right neighbor gains 1, except when the topple is itself at the edge.
        right_targets = unstable[unstable < L - 1] + 1
        np.add.at(dz, right_targets, 1)

        # A topple at site L-1 ejects exactly one grain off the right edge.
        state.grains_lost += int(right_edge_topples.size)

        z += dz

        # Redraw thresholds at every site that toppled this sweep.
        z_c[unstable] = rng.integers(1, 3, size=unstable.size, dtype=np.int64)

    return s, T


def run(
    L: int,
    n_drops: int,
    seed: int | None = None,
) -> tuple[OsloState, np.ndarray, np.ndarray]:
    """Drive an empty Oslo pile for n_drops grain additions.

    Each addition is followed by a complete relaxation. The number of returned
    avalanche records equals n_drops.

    Returns
    -------
    state : OsloState
        Final state after all drives.
    sizes : ndarray[int], shape (n_drops,)
        Avalanche sizes s in temporal order.
    durations : ndarray[int], shape (n_drops,)
        Avalanche durations T in temporal order.
    """
    rng = np.random.default_rng(seed)
    state = initialize(L, rng)
    sizes = np.zeros(n_drops, dtype=np.int64)
    durations = np.zeros(n_drops, dtype=np.int64)
    for t in range(n_drops):
        drive(state)
        s, T = relax(state, rng)
        sizes[t] = s
        durations[t] = T
    return state, sizes, durations
