"""Buildup -> breakdown -> tip: the substrate igniting itself on its OWN stress.

Shawna's mechanism, force-free (NO gravity, NO imposed inflow). A cold substrate
slowly absorbs darkness (drive). When a cell is over-stressed it BREAKS DOWN:
it sheds stress to its neighbours (a Manna chain reaction) and converts a little
of that stress to HEAT right there. The breaking clusters, so the heat clusters.

The heat comes from the BREAKDOWN itself -- not from anything falling in. With a
slow enough drive the substrate can sit cold forever; with a fast-enough drive it
sits quiet a long time and then SUDDENLY tips to plasma. The cold control (no
stress->heat) can never tip, so a tip is never an artefact of rigging.

Energy is conserved EXACTLY: the avalanche sheds stress with periodic wrap (no
edge, nothing falls off the world), and stress->heat and combustion only move
energy between `load` and `E`. So at all times:

    E_fluid + load  ==  E_fluid_0 + (total stress driven in)

This module is the engine only: buildup -> breakdown -> tip. The fall-in (the
structure losing its cold rigidity and giving way, which then INTENSIFIES the
breakdown -- and maybe erupts back out) is the next layer, built on top of this.

Honest scope: 2D reacting fluid (Lax-Friedrichs, diffusive); abstract units;
deterministic 4-neighbour shedding (a sandpile variant -- NOT the random-2 Manna
used for the tau=1.27 measurement, but it clusters and conserves cleanly).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from void_cascade.material_motion_2d import (
    GAMMA, internal_energy, lax_friedrichs_step, max_wave_speed,
)


@dataclass(frozen=True)
class BreakdownParams:
    thr: float = 2.0           # stress threshold for a cell to break down
    hpc: float = 0.1           # stress converted to LOCAL heat per break
    e_ign: float = 2.5         # specific internal energy to combust (plasma)
    drive_sites: int = 4       # cells that absorb darkness per step
    drive_amount: float = 1.0  # stress added per site per step
    cfl: float = 0.4           # fluid timestep safety factor
    max_sweeps: int = 2000     # cap on avalanche relaxation sweeps per step


def drive(load, rng, n_sites, amount):
    """Absorb darkness: add `amount` of stress to `n_sites` random cells. Pure."""
    out = load.copy()
    flat = out.ravel()
    np.add.at(flat, rng.integers(0, flat.size, size=n_sites), amount)
    return out


def breakdown(load, rho, momx, momy, E, thr, hpc, e_ign, max_sweeps=2000):
    """Manna chain reaction with stress->heat. A COLD, over-stressed cell breaks:
    it sheds `thr` of stress equally to its 4 neighbours (periodic wrap) and turns
    up to `hpc` of stress into LOCAL heat. Relaxes until no cold over-stressed cell
    remains (or the sweep cap). Conserves load+E exactly. Returns (load, E)."""
    load = load.copy()
    E = E.copy()
    for _ in range(max_sweeps):
        e_spec = internal_energy(rho, momx, momy, E) / rho
        over = (load >= thr) & (e_spec < e_ign)
        if not over.any():
            break
        shed = np.where(over, thr, 0.0)
        nb = (np.roll(shed, 1, 0) + np.roll(shed, -1, 0)
              + np.roll(shed, 1, 1) + np.roll(shed, -1, 1)) / 4.0
        load = load - shed + nb                 # periodic shed: conserves exactly
        take = np.where(over, np.minimum(hpc, load), 0.0)
        load = load - take                      # stress ...
        E = E + take                            # ... becomes local heat (clustered)
    return load, E


def combust(rho, momx, momy, E, load, e_ign):
    """Hot cells (specific internal energy >= e_ign) release their remaining stress
    into heat. Returns (E, load, ignited)."""
    e_spec = internal_energy(rho, momx, momy, E) / rho
    fire = (e_spec >= e_ign) & (load > 0.0)
    released = np.where(fire, load, 0.0)
    return E + released, load - released, bool(fire.any())


def total_energy(E, load):
    """The conserved quantity: fluid energy + stored stress."""
    return float(E.sum() + load.sum())


def step(rho, momx, momy, E, load, p, rng):
    """One step: absorb -> break down -> combust -> let the fluid move."""
    load = drive(load, rng, p.drive_sites, p.drive_amount)
    load, E = breakdown(load, rho, momx, momy, E, p.thr, p.hpc, p.e_ign, p.max_sweeps)
    E, load, ignited = combust(rho, momx, momy, E, load, p.e_ign)
    s = max_wave_speed(rho, momx, momy, E)
    dt = p.cfl / max(s, 1e-9)
    rho, momx, momy, E = lax_friedrichs_step(rho, momx, momy, E, 1.0, dt)
    return rho, momx, momy, E, load, ignited


def run_buildup_tip(L=80, steps=3000, params=None, seed=0, sample_every=15, P0=1.0):
    """Cold, empty substrate -> slow drive -> does it tip itself to plasma, and when?
    Returns the tip step (or None), the peak temperature, the conservation residual,
    and sampled traces of temperature and mean stress over time."""
    p = params if params is not None else BreakdownParams()
    rng = np.random.default_rng(seed)
    load = np.zeros((L, L))
    rho = np.ones((L, L))
    momx = np.zeros((L, L))
    momy = np.zeros((L, L))
    E = np.full((L, L), P0 / (GAMMA - 1.0))
    e0 = total_energy(E, load)
    driven = 0.0
    tip_step = None
    max_temp = float((internal_energy(rho, momx, momy, E) / rho).max())
    t_axis, temp_trace, load_trace = [], [], []
    for t in range(steps):
        rho, momx, momy, E, load, ignited = step(rho, momx, momy, E, load, p, rng)
        driven += p.drive_sites * p.drive_amount
        if tip_step is None and ignited:
            tip_step = t
        cur = float((internal_energy(rho, momx, momy, E) / rho).max())
        max_temp = max(max_temp, cur)
        if t % sample_every == 0:
            t_axis.append(t)
            temp_trace.append(cur)
            load_trace.append(float(load.mean()))
    return {
        "L": L, "steps": steps, "tip_step": tip_step,
        "max_temp": max_temp,
        "energy_residual": total_energy(E, load) - (e0 + driven),
        "t_axis": t_axis, "temp_trace": temp_trace, "load_trace": load_trace,
    }
