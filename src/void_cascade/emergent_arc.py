"""Emergent arc (collapse driver = fracture-voids): nothing imposed but local rules.

One fluid+fuel grid. The whole arc emerges from a single temperature-gated
fuel<->energy exchange (no imposed implosion, no imposed ignition, no new force):

  - BUILDUP: fuel slowly accumulates everywhere (absorbing darkness).
  - FRACTURE (cold + loaded): when a cell's fuel passes a threshold and it is
    cold, it fractures -- internal energy is moved into fuel, so its PRESSURE
    DROPS, opening a low-pressure VOID. Surrounding higher-pressure material
    flows IN (collapse). This is the implosion, emergent, driven only by the
    pressure gradient the fracture creates.
  - COMBUST (hot): where inflow has compressed and heated the gas past an
    ignition point, fuel is released back into energy -- pressure spikes ->
    detonation -> expansion.

Fracture and combustion are the two directions of the same exchange, set by
temperature (cold absorbs into fuel + collapses; hot releases + expands).

Energy conserved: energy_in == E_fluid + fuel + boundary_lost. Honest scope:
2D reacting Euler (Lax-Friedrichs); abstract units; flat-local grid.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from void_cascade.material_motion_2d import (
    GAMMA, internal_energy, lax_friedrichs_step, max_wave_speed,
)


@dataclass(frozen=True)
class ArcParams:
    drive_rate: float        # fuel (darkness) absorbed per cell per step
    fracture_fuel: float     # fuel threshold to fracture (cold)
    void_frac: float         # fraction of internal energy moved to fuel on fracture
    internal_floor: float    # keep internal energy (pressure) above this
    e_ign: float             # specific internal energy to combust (hot)


@dataclass
class ArcState:
    rho: np.ndarray
    momx: np.ndarray
    momy: np.ndarray
    E: np.ndarray
    fuel: np.ndarray
    fractured: np.ndarray
    energy_in: float = 0.0
    boundary_lost: float = 0.0


def initialize(L, rho0=1.0, P0=1.0):
    E0 = P0 / (GAMMA - 1.0)
    st = ArcState(
        rho=np.full((L, L), rho0), momx=np.zeros((L, L)), momy=np.zeros((L, L)),
        E=np.full((L, L), E0), fuel=np.zeros((L, L)), fractured=np.zeros((L, L), bool),
    )
    st.energy_in = float(st.E.sum() + st.fuel.sum())
    return st


def drive(state, p, rng, n_sites):
    flat = state.fuel.ravel()
    sites = rng.integers(0, flat.size, size=n_sites)
    np.add.at(flat, sites, p.drive_rate)
    state.energy_in += n_sites * p.drive_rate


def fracture(state, p):
    """Cold, loaded cells fracture: internal energy -> fuel (pressure drops ->
    void). Conserves energy (E_fluid + fuel)."""
    e_spec = internal_energy(state.rho, state.momx, state.momy, state.E) / state.rho
    internal = internal_energy(state.rho, state.momx, state.momy, state.E)
    frac = (state.fuel >= p.fracture_fuel) & (~state.fractured) & (e_spec < p.e_ign)
    dE = np.where(frac, p.void_frac * np.maximum(internal - p.internal_floor, 0.0), 0.0)
    state.E -= dE
    state.fuel += dE
    state.fractured |= frac
    return int(frac.sum())


def combust(state, p):
    """Hot cells release fuel back into energy (pressure up -> detonation)."""
    e_spec = internal_energy(state.rho, state.momx, state.momy, state.E) / state.rho
    fire = (e_spec >= p.e_ign) & (state.fuel > 0.0)
    released = np.where(fire, state.fuel, 0.0)
    state.E += released
    state.fuel -= released
    return int(fire.sum())


def step(state, p, rng, dx, cfl, n_sites):
    drive(state, p, rng, n_sites)
    fracture(state, p)
    combust(state, p)
    s = max_wave_speed(state.rho, state.momx, state.momy, state.E)
    dt = cfl * dx / max(s, 1e-9)
    state.rho, state.momx, state.momy, state.E = lax_friedrichs_step(
        state.rho, state.momx, state.momy, state.E, dx, dt)
    return dt


def conservation_residual(state):
    return state.energy_in - (float(state.E.sum() + state.fuel.sum()) + state.boundary_lost)
