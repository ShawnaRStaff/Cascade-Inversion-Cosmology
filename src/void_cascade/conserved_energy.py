"""Conserved-energy front: energy is considered, not imposed.

Instead of inventing heat (heat_per_crack) and destroying it (cooling -> an
outside sink, which is a maker's box again), energy is CONSERVED and only
changes form:

    potential  (stored darkness / fuel, at rest)
       <->  kinetic  (heat / motion)

A cell COMBUSTS when its total energy (potential + kinetic) crosses a
threshold: it converts its potential into kinetic. That kinetic then spreads
to neighbours -- a shock of energy that runs AHEAD of the front and PRIMES the
not-yet-combusted substrate toward its own threshold. The priming-ahead is the
front engine (a combustion / detonation front), not an imposed propagation
rule. Nothing is destroyed; energy that spreads is conserved.

Honest question: with energy conserved (no imposed sink), does the front
sustain by priming the substrate ahead, or does the released energy dilute
faster than it can prime the next layer (fizzle)? Built so it CAN fizzle.

Conservation law (the cornerstone test):
    energy_in (darkness driven in)  ==  sum(potential) + sum(kinetic)
                                        + boundary_lost
Boundary loss is itself a sink / box artifact -- to be removed by running this
on the edgeless growing substrate next. Here we TRACK it so conservation is
exact and provable.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class EnergyParams:
    flip_threshold: float    # combust when potential+kinetic >= this
    release_fraction: float  # fraction of potential converted to kinetic on combust
    diffuse: float           # fraction of a cell's kinetic that spreads per step (<=1)
    drive_amount: float      # darkness (potential) added per driven site per step
    n_drive_sites: int


@dataclass
class EnergyState:
    potential: np.ndarray   # stored energy (fuel), >= 0
    kinetic: np.ndarray     # motion/heat energy, >= 0
    flipped: np.ndarray     # bool, has this cell combusted
    energy_in: float = 0.0       # cumulative darkness driven in
    boundary_lost: float = 0.0   # cumulative kinetic lost off the (open) edge


def initialize_2d(L: int) -> EnergyState:
    return EnergyState(
        potential=np.zeros((L, L), dtype=float),
        kinetic=np.zeros((L, L), dtype=float),
        flipped=np.zeros((L, L), dtype=bool),
    )


def drive(state: EnergyState, p: EnergyParams, rng: np.random.Generator) -> None:
    """Absorb darkness into potential (the only energy input)."""
    flat = state.potential.ravel()
    sites = rng.integers(0, flat.size, size=p.n_drive_sites)
    np.add.at(flat, sites, p.drive_amount)
    state.energy_in += p.n_drive_sites * p.drive_amount


def ignite(state: EnergyState, p: EnergyParams, radius_frac: float = 0.04) -> None:
    """Inject enough kinetic at the centre to start combustion (no energy
    bookkeeping cheat: we add it to energy_in so conservation still holds)."""
    R, C = state.kinetic.shape
    cy, cx = R // 2, C // 2
    r = max(1, int(radius_frac * min(R, C)))
    add = p.flip_threshold * 2.0
    patch = state.kinetic[cy - r:cy + r + 1, cx - r:cx + r + 1]
    state.kinetic[cy - r:cy + r + 1, cx - r:cx + r + 1] = patch + add
    state.energy_in += add * patch.size


def combust(state: EnergyState, p: EnergyParams) -> int:
    """Cells over threshold convert potential -> kinetic (energy conserved),
    once. The kinetic then primes neighbours via diffusion."""
    total = state.potential + state.kinetic
    fire = (total >= p.flip_threshold) & (~state.flipped)
    n = int(fire.sum())
    if n:
        released = p.release_fraction * state.potential[fire]
        state.kinetic[fire] += released
        state.potential[fire] -= released
        state.flipped[fire] = True
    return n


def diffuse_kinetic(state: EnergyState, p: EnergyParams) -> None:
    """Spread kinetic to the 4 neighbours, conserving it; energy that would go
    off-lattice is counted as boundary loss (a sink to be removed by going
    edgeless). NO cooling -- nothing is destroyed in the bulk."""
    k = state.kinetic
    send = p.diffuse * k
    per = send / 4.0
    new = k - send
    new[1:, :] += per[:-1, :]
    new[:-1, :] += per[1:, :]
    new[:, 1:] += per[:, :-1]
    new[:, :-1] += per[:, 1:]
    lost = float(per[0, :].sum() + per[-1, :].sum() + per[:, 0].sum() + per[:, -1].sum())
    state.kinetic = new
    state.boundary_lost += lost


def step(state: EnergyState, p: EnergyParams, rng: np.random.Generator) -> int:
    drive(state, p, rng)
    n = combust(state, p)
    diffuse_kinetic(state, p)
    return n


def total_energy(state: EnergyState) -> float:
    return float(state.potential.sum() + state.kinetic.sum())


def conservation_residual(state: EnergyState) -> float:
    """Should be ~0: energy_in - (potential + kinetic + boundary_lost)."""
    return state.energy_in - (total_energy(state) + state.boundary_lost)
