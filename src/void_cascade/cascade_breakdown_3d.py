"""3D buildup -> breakdown -> tip engine.

Extends cascade_breakdown (2D) to a cubic L×L×L grid with 6 face-neighbours.
The physics is identical: the substrate absorbs darkness (drive), cells that
accumulate load above threshold break (chain reaction), converting stress to
local heat. The front propagates in all 3 directions.

Key differences from 2D:
  - 6 neighbours (±x ±y ±z) instead of 4. Shed distributes thr/6 each.
  - LF step uses material_motion_3d (5 conserved vars: rho, momx, momy, momz, E).
  - CFL conservative bound dt = cfl / max_wave_speed with cfl ≤ 0.4/sqrt(3).

Everything else (Honest Code, exact conservation, cold control) is inherited
from the 2D design.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from void_cascade.material_motion_3d import (
    GAMMA,
    internal_energy_3d,
    lax_friedrichs_step_3d,
    max_wave_speed_3d,
)


@dataclass(frozen=True)
class BreakdownParams3d:
    thr: float = 2.0
    hpc: float = 0.1
    e_ign: float = 2.5
    drive_sites: int = 4
    drive_amount: float = 1.0
    cfl: float = 0.23        # 0.4/sqrt(3) ≈ 0.23 for 3D stability
    max_sweeps: int = 2000


def drive_3d(load, rng, n_sites, amount):
    """Add stress to n_sites random cells. Pure."""
    out = load.copy()
    flat = out.ravel()
    np.add.at(flat, rng.integers(0, flat.size, size=n_sites), amount)
    return out


def breakdown_3d(load, rho, momx, momy, momz, E, thr, hpc, e_ign, max_sweeps=2000):
    """Manna chain reaction on a 3D cubic grid with 6-neighbour shedding.

    Cells above threshold that are cold shed `thr` stress evenly to 6 neighbours
    (thr/6 each, periodic wrap) and convert up to hpc of residual load to heat.
    Conserves load + E exactly.
    """
    load = load.copy()
    E    = E.copy()
    for _ in range(max_sweeps):
        e_spec = internal_energy_3d(rho, momx, momy, momz, E) / rho
        over   = (load >= thr) & (e_spec < e_ign)
        if not over.any():
            break
        shed = np.where(over, thr, 0.0)
        nb = (np.roll(shed, 1, 0) + np.roll(shed, -1, 0)
            + np.roll(shed, 1, 1) + np.roll(shed, -1, 1)
            + np.roll(shed, 1, 2) + np.roll(shed, -1, 2)) / 6.0
        load = load - shed + nb
        take = np.where(over, np.minimum(hpc, load), 0.0)
        load = load - take
        E    = E    + take
    return load, E


def combust_3d(rho, momx, momy, momz, E, load, e_ign):
    """Hot cells release remaining load into heat. Returns (E, load, ignited)."""
    e_spec   = internal_energy_3d(rho, momx, momy, momz, E) / rho
    fire     = (e_spec >= e_ign) & (load > 0.0)
    released = np.where(fire, load, 0.0)
    return E + released, load - released, bool(fire.any())


def total_energy_3d(E, load):
    return float(E.sum() + load.sum())


def step_3d(rho, momx, momy, momz, E, load, p, rng):
    """One step: drive -> breakdown -> combust -> fluid transport (3D LF)."""
    load              = drive_3d(load, rng, p.drive_sites, p.drive_amount)
    load, E           = breakdown_3d(load, rho, momx, momy, momz, E,
                                     p.thr, p.hpc, p.e_ign, p.max_sweeps)
    E, load, ignited  = combust_3d(rho, momx, momy, momz, E, load, p.e_ign)
    s  = max_wave_speed_3d(rho, momx, momy, momz, E)
    dt = p.cfl / max(s, 1e-9)
    rho, momx, momy, momz, E = lax_friedrichs_step_3d(rho, momx, momy, momz, E, 1.0, dt)
    return rho, momx, momy, momz, E, load, ignited


def correlation_length_3d(field):
    """Radial spatial correlation length of a 3D field via FFT autocorrelation.

    Returns the smallest r where the radial ACF drops below 1/e.
    Returns 0.0 for uniform fields. Pure.
    """
    f = field - field.mean()
    var = float((f ** 2).mean())
    if var < 1e-12:
        return 0.0
    L = field.shape[0]
    ft = np.fft.rfftn(f)
    acf3d = np.fft.irfftn(ft * ft.conj(), s=field.shape, axes=(0, 1, 2)).real / (L ** 3 * var)
    acf3d = np.roll(np.roll(np.roll(acf3d, L // 2, 0), L // 2, 1), L // 2, 2)
    cz, cy, cx = np.mgrid[0:L, 0:L, 0:L] - L // 2
    r = np.sqrt(cz ** 2 + cy ** 2 + cx ** 2).astype(int).ravel()
    r_max = L // 2
    acf_r = np.bincount(r, weights=acf3d.ravel(), minlength=r_max + 1)[:r_max]
    cnt   = np.bincount(r, minlength=r_max + 1)[:r_max].clip(1)
    acf_r = acf_r / cnt
    below = np.where(acf_r < 1.0 / np.e)[0]
    return float(below[0]) if len(below) > 0 else float(r_max)


def run_onset_measurement_3d(L=20, steps=3000, params=None, seed=0, sample_every=15, P0=1.0):
    """Track three onset observables during 3D buildup-to-tip.

    Same contract as run_onset_measurement (2D): clocks when load.std,
    corr_length, and n_hot first become detectable on a shared timeline.
    Returns dict with t_axis, load_std, corr_lengths, n_hot, tip_step.
    """
    p   = params if params is not None else BreakdownParams3d()
    rng = np.random.default_rng(seed)
    load = np.zeros((L, L, L))
    rho  = np.ones((L, L, L))
    momx = np.zeros((L, L, L))
    momy = np.zeros((L, L, L))
    momz = np.zeros((L, L, L))
    E    = np.full((L, L, L), P0 / (GAMMA - 1.0))
    tip_step = None
    t_axis, load_std_series, corr_lengths, n_hot_series = [], [], [], []

    for t in range(steps):
        rho, momx, momy, momz, E, load, ignited = step_3d(
            rho, momx, momy, momz, E, load, p, rng,
        )
        if tip_step is None and ignited:
            tip_step = t
        if t % sample_every == 0:
            e_spec = internal_energy_3d(rho, momx, momy, momz, E) / rho
            t_axis.append(t)
            load_std_series.append(float(load.std()))
            corr_lengths.append(correlation_length_3d(load))
            n_hot_series.append(int((e_spec >= p.e_ign).sum()))

    return {
        "L": L, "steps": steps, "tip_step": tip_step,
        "t_axis": t_axis,
        "load_std": load_std_series,
        "corr_lengths": corr_lengths,
        "n_hot": n_hot_series,
    }


def run_buildup_tip_3d(L=30, steps=3000, params=None, seed=0, sample_every=20, P0=1.0):
    """Cold 3D substrate -> slow drive -> tip to plasma (or not).

    Returns tip_step (or None), max_temp, energy_residual, and sampled traces.
    """
    p   = params if params is not None else BreakdownParams3d()
    rng = np.random.default_rng(seed)
    load = np.zeros((L, L, L))
    rho  = np.ones((L, L, L))
    momx = np.zeros((L, L, L))
    momy = np.zeros((L, L, L))
    momz = np.zeros((L, L, L))
    E    = np.full((L, L, L), P0 / (GAMMA - 1.0))
    e0      = total_energy_3d(E, load)
    driven  = 0.0
    tip_step = None
    max_temp = float((internal_energy_3d(rho, momx, momy, momz, E) / rho).max())
    t_axis, temp_trace, load_trace = [], [], []

    for t in range(steps):
        rho, momx, momy, momz, E, load, ignited = step_3d(
            rho, momx, momy, momz, E, load, p, rng
        )
        driven += p.drive_sites * p.drive_amount
        if tip_step is None and ignited:
            tip_step = t
        cur = float((internal_energy_3d(rho, momx, momy, momz, E) / rho).max())
        max_temp = max(max_temp, cur)
        if t % sample_every == 0:
            t_axis.append(t)
            temp_trace.append(cur)
            load_trace.append(float(load.mean()))

    return {
        "L": L, "steps": steps, "tip_step": tip_step,
        "max_temp": max_temp,
        "energy_residual": total_energy_3d(E, load) - (e0 + driven),
        "t_axis": t_axis, "temp_trace": temp_trace, "load_trace": load_trace,
    }
