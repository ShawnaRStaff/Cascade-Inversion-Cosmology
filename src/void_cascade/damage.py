"""Damage geometries for the substrate-resilience scrutiny.

The old resilience experiment removed substrate cells *uniformly at
random*. A real fracture event does not: it spreads through neighbouring
cells, like a crack running across glass. The shape of the damage changes
whether the leftover substrate stays in one connected piece, so we keep
the two geometries separate and compare them.

These are pure functions: shape + fraction + random generator in, a
boolean mask out (True = destroyed). They hold no state.

`connectivity` of the leftover (the complement of a damage mask) is
measured with `percolation.check_spanning`, which already exists; we do
not duplicate it here.
"""

from __future__ import annotations

import numpy as np


def damage_count(shape: tuple[int, int, int], fraction: float) -> int:
    """How many cells a given fraction of the lattice corresponds to.

    Truncated to a whole number of cells (floor).
    """
    n_total = int(np.prod(shape))
    return int(fraction * n_total)


def random_damage_mask(
    shape: tuple[int, int, int],
    fraction: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Destroy a fraction of cells chosen uniformly at random.

    Returns a boolean array of `shape`; True marks a destroyed cell.
    This reproduces the damage model the original resilience experiment
    used (scattered, spatially uncorrelated).
    """
    n_total = int(np.prod(shape))
    n = damage_count(shape, fraction)
    flat = np.zeros(n_total, dtype=bool)
    if n > 0:
        idx = rng.choice(n_total, size=n, replace=False)
        flat[idx] = True
    return flat.reshape(shape)


def _on_lattice_neighbors(
    cell: tuple[int, int, int],
    shape: tuple[int, int, int],
) -> list[tuple[int, int, int]]:
    """The 6 face-neighbours of `cell` that lie inside the lattice."""
    i, j, k = cell
    L0, L1, L2 = shape
    out = []
    if i + 1 < L0:
        out.append((i + 1, j, k))
    if i - 1 >= 0:
        out.append((i - 1, j, k))
    if j + 1 < L1:
        out.append((i, j + 1, k))
    if j - 1 >= 0:
        out.append((i, j - 1, k))
    if k + 1 < L2:
        out.append((i, j, k + 1))
    if k - 1 >= 0:
        out.append((i, j, k - 1))
    return out


def connected_damage_mask(
    shape: tuple[int, int, int],
    fraction: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Destroy a single connected blob of cells (a crack proxy).

    Grown Eden-style: start from one random seed cell, then repeatedly
    add a randomly chosen face-neighbour on the current frontier until
    the target number of cells is reached. The result is one 6-connected
    region (no scattered specks), which is the spatial signature a real
    spreading fracture would leave.

    Returns a boolean array of `shape`; True marks a destroyed cell.
    """
    mask = np.zeros(shape, dtype=bool)
    n = damage_count(shape, fraction)
    if n <= 0:
        return mask

    L0, L1, L2 = shape
    seed = (int(rng.integers(L0)), int(rng.integers(L1)), int(rng.integers(L2)))
    in_blob = {seed}
    mask[seed] = True

    # Frontier = candidate cells touching the blob but not yet in it.
    # Kept as a list (for uniform random pick) plus a set (for membership).
    frontier = _on_lattice_neighbors(seed, shape)
    frontier_set = set(frontier)

    while len(in_blob) < n and frontier:
        pick = int(rng.integers(len(frontier)))
        cell = frontier[pick]
        # Swap-pop to remove the chosen cell in O(1).
        frontier[pick] = frontier[-1]
        frontier.pop()
        frontier_set.discard(cell)
        if cell in in_blob:
            continue
        in_blob.add(cell)
        mask[cell] = True
        for nb in _on_lattice_neighbors(cell, shape):
            if nb not in in_blob and nb not in frontier_set:
                frontier.append(nb)
                frontier_set.add(nb)

    return mask
