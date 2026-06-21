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


def correlation_length(field):
    """Radial spatial correlation length of a 2D field, via FFT autocorrelation.

    Returns the smallest r where the normalized radial ACF drops below 1/e.
    Returns 0.0 for uniform fields (no spatial variation -> no structure).
    Pure function: no side effects.
    """
    f = field - field.mean()
    var = float((f ** 2).mean())
    if var < 1e-12:
        return 0.0
    L = field.shape[0]
    ft = np.fft.rfft2(f)
    acf2d = np.fft.irfft2(ft * ft.conj(), s=field.shape).real / (L * L * var)
    acf2d = np.roll(np.roll(acf2d, L // 2, 0), L // 2, 1)
    cy, cx = np.mgrid[0:L, 0:L] - L // 2
    r = np.sqrt(cy ** 2 + cx ** 2).astype(int).ravel()
    r_max = L // 2
    acf_r = np.bincount(r, weights=acf2d.ravel(), minlength=r_max + 1)[:r_max]
    cnt   = np.bincount(r, minlength=r_max + 1)[:r_max].clip(1)
    acf_r = acf_r / cnt
    below = np.where(acf_r < 1.0 / np.e)[0]
    return float(below[0]) if len(below) > 0 else float(r_max)


def run_onset_measurement(L=80, steps=3000, params=None, seed=0, sample_every=15, P0=1.0):
    """Track three onset observables during the buildup-to-tip arc.

    At each sample step records:
        n_over      : cells with load >= thr (proxy for avalanche activity)
        corr_length : radial correlation length of the load field (spatial structure)
        n_hot       : cells above ignition temperature (ignition front)

    Returns dict with t_axis, n_over, corr_lengths, n_hot, tip_step.
    """
    p = params if params is not None else BreakdownParams()
    rng = np.random.default_rng(seed)
    load = np.zeros((L, L))
    rho  = np.ones((L, L))
    momx = np.zeros((L, L))
    momy = np.zeros((L, L))
    E    = np.full((L, L), P0 / (GAMMA - 1.0))
    tip_step = None
    t_axis, load_std_series, corr_lengths, n_hot_series = [], [], [], []
    for t in range(steps):
        rho, momx, momy, E, load, ignited = step(rho, momx, momy, E, load, p, rng)
        if tip_step is None and ignited:
            tip_step = t
        if t % sample_every == 0:
            e_spec = internal_energy(rho, momx, momy, E) / rho
            t_axis.append(t)
            # breakdown always clears cells above threshold, so n_over is always 0.
            # load.std() is the honest proxy: grows as stress clusters (SOC onset).
            load_std_series.append(float(load.std()))
            corr_lengths.append(correlation_length(load))
            n_hot_series.append(int((e_spec >= p.e_ign).sum()))
    return {
        "L": L, "steps": steps, "tip_step": tip_step,
        "t_axis": t_axis,
        "load_std": load_std_series,
        "corr_lengths": corr_lengths,
        "n_hot": n_hot_series,
    }


def step_rigid(rho, momx, momy, E, load, p, rng, rigidity=0.0):
    """Like step() but cold loaded cells express stored load as elastic pressure.

    A fraction min(rigidity, 1) of each cold cell's load is temporarily lent to
    the fluid energy field before the LF step, making cold loaded cells stiffer
    (higher pressure). After LF, the lent energy returns to load in cells still
    cold; cells that crossed the ignition threshold during transport keep it as
    thermal energy (the crystal melted -> elastic energy became heat).

    Conservation: E + load is EXACTLY conserved (same as step). rigidity=0 ->
    numerically identical to step().
    """
    load = drive(load, rng, p.drive_sites, p.drive_amount)
    load, E = breakdown(load, rho, momx, momy, E, p.thr, p.hpc, p.e_ign, p.max_sweeps)
    E, load, ignited = combust(rho, momx, momy, E, load, p.e_ign)

    if rigidity > 0.0:
        e_spec = internal_energy(rho, momx, momy, E) / rho
        cold_pre = (e_spec < p.e_ign).astype(float)
        frac = min(rigidity, 1.0)
        lend = frac * load * cold_pre      # cold cells lend elastic energy to fluid
        E_eff = E + lend
        load_eff = load - lend
    else:
        E_eff, load_eff = E, load

    s = max_wave_speed(rho, momx, momy, E_eff)
    dt = p.cfl / max(s, 1e-9)
    rho, momx, momy, E_out = lax_friedrichs_step(rho, momx, momy, E_eff, 1.0, dt)

    if rigidity > 0.0:
        e_spec_out = internal_energy(rho, momx, momy, E_out) / rho
        cold_post = (e_spec_out < p.e_ign).astype(float)
        # Return lent energy to load for cells still cold; cap at E_out so E ≥ 0.
        # Global conservation holds regardless of cap:
        #   E_final + load_final = (E_out - lend_return) + (load_eff + lend_return)
        #                        = E_out + load_eff = E_eff + load_eff = E + load ✓
        lend_return = np.minimum(frac * load_eff * cold_post, E_out)
        E_final = E_out - lend_return
        load_final = load_eff + lend_return
    else:
        E_final, load_final = E_out, load_eff

    return rho, momx, momy, E_final, load_final, ignited


def step_cool(rho, momx, momy, E, load, p, rng, cooling=0.0, e_floor=1.5):
    """Like step() but hot cells bleed heat to the void each step.

    After LF transport, any cell with specific internal energy above `e_floor`
    loses `cooling` fraction of that excess:
        delta_E = cooling * max(e_spec - e_floor, 0) * rho

    Energy leaves the system -- total E+load decreases. Full accounting:
        E + load + cumulative_cooled == E0 + total_driven.

    cooling=0 -> numerically identical to step(); cooled_this_step=0.
    Returns (rho, momx, momy, E, load, ignited, cooled_this_step).
    """
    load = drive(load, rng, p.drive_sites, p.drive_amount)
    load, E = breakdown(load, rho, momx, momy, E, p.thr, p.hpc, p.e_ign, p.max_sweeps)
    E, load, ignited = combust(rho, momx, momy, E, load, p.e_ign)
    s = max_wave_speed(rho, momx, momy, E)
    dt = p.cfl / max(s, 1e-9)
    rho, momx, momy, E = lax_friedrichs_step(rho, momx, momy, E, 1.0, dt)

    cooled = 0.0
    if cooling > 0.0:
        e_spec = internal_energy(rho, momx, momy, E) / rho
        excess = np.maximum(e_spec - e_floor, 0.0)
        # cooling <= 1 guarantees e_spec stays >= e_floor; cap at 1.0 for safety.
        delta_E = min(cooling, 1.0) * excess * rho
        E = E - delta_E
        cooled = float(delta_E.sum())

    return rho, momx, momy, E, load, ignited, cooled


def run_cooling_arc(
    L=80, steps=3000, params=None, seed=0, sample_every=20, P0=1.0, cooling=0.0,
):
    """Cold substrate -> buildup -> tip (or not), with per-step heat loss.

    Uses step_cool throughout. Tracks cumulative energy lost to cooling for
    full accounting: E + load + total_cooled == E0 + total_driven.

    Returns tip_step, max_temp, energy_residual (should be ~0 with accounting),
    total_cooled, t_axis, temp_trace, load_trace.
    """
    p         = params if params is not None else BreakdownParams()
    e_floor   = P0 / (GAMMA - 1.0)
    rng       = np.random.default_rng(seed)
    load      = np.zeros((L, L))
    rho       = np.ones((L, L))
    momx      = np.zeros((L, L))
    momy      = np.zeros((L, L))
    E         = np.full((L, L), e_floor)
    e0        = total_energy(E, load)
    driven    = 0.0
    cooled    = 0.0
    tip_step  = None
    max_temp  = float((internal_energy(rho, momx, momy, E) / rho).max())
    t_axis, temp_trace, load_trace = [], [], []

    for t in range(steps):
        rho, momx, momy, E, load, ignited, dc = step_cool(
            rho, momx, momy, E, load, p, rng, cooling=cooling, e_floor=e_floor,
        )
        driven += p.drive_sites * p.drive_amount
        cooled += dc
        if tip_step is None and ignited:
            tip_step = t
        cur = float((internal_energy(rho, momx, momy, E) / rho).max())
        max_temp = max(max_temp, cur)
        if t % sample_every == 0:
            t_axis.append(t)
            temp_trace.append(cur)
            load_trace.append(float(load.mean()))

    # Full accounting: E + load + cooled should equal E0 + driven
    residual = total_energy(E, load) + cooled - (e0 + driven)
    return {
        "L": L, "steps": steps, "tip_step": tip_step,
        "max_temp": max_temp,
        "energy_residual": residual,
        "total_cooled": cooled,
        "t_axis": t_axis, "temp_trace": temp_trace, "load_trace": load_trace,
    }


def step_melt(rho, momx, momy, E, load, p, rng, melt_frac=0.0):
    """Like step() but cold cells absorb their post-LF kinetic energy as load.

    After the fluid transport step, cells that are still cold (e_spec < e_ign)
    convert a fraction `melt_frac` of their kinetic energy into mechanical stress
    (load) rather than keeping it as motion. Hot cells flow freely.

    Physical picture: cold crystalline cells are stiff walls. When hot plasma
    pushes against them, they don't flow -- they get LOADED. That load can then
    trigger more breakdown in the next step, potentially amplifying the cascade.

    Conservation: E + load is EXACTLY conserved. melt_frac=0 -> identical to step().

    The KE removed from cold cells = 0.5*(momx^2+momy^2)/rho * cold * melt_frac.
    Momentum is scaled by sqrt(1-melt_frac) for cold cells so the KE reduction
    matches exactly. E decreases by that amount; load increases by that amount.
    """
    load = drive(load, rng, p.drive_sites, p.drive_amount)
    load, E = breakdown(load, rho, momx, momy, E, p.thr, p.hpc, p.e_ign, p.max_sweeps)
    E, load, ignited = combust(rho, momx, momy, E, load, p.e_ign)

    s = max_wave_speed(rho, momx, momy, E)
    dt = p.cfl / max(s, 1e-9)
    rho, momx, momy, E = lax_friedrichs_step(rho, momx, momy, E, 1.0, dt)

    if melt_frac > 0.0:
        e_spec = internal_energy(rho, momx, momy, E) / rho
        cold = (e_spec < p.e_ign).astype(float)
        KE_cell = 0.5 * (momx ** 2 + momy ** 2) / rho
        absorbed = melt_frac * KE_cell * cold        # energy to absorb per cell
        # Scale momentum: KE ∝ mom^2, so scaling mom by sqrt(1-f) removes fraction f.
        # Capped at 1.0 to avoid domain error if melt_frac > 1.
        mom_scale = np.where(cold.astype(bool),
                             np.sqrt(max(1.0 - melt_frac, 0.0)), 1.0)
        momx  = momx  * mom_scale
        momy  = momy  * mom_scale
        # Cap absorbed at E so E never goes negative cell-by-cell.
        absorbed = np.minimum(absorbed, E)
        E    = E    - absorbed
        load = load + absorbed

    return rho, momx, momy, E, load, ignited


def run_melt_arc(
    L=80, steps=3000, params=None, seed=0, sample_every=20, P0=1.0, melt_frac=0.0,
):
    """Cold substrate -> buildup -> tip, with melt-gated rigidity.

    Identical to run_buildup_tip but uses step_melt throughout. Tracks n_hot
    at each sample so the cascade shape (slow vs fast intensification) can be
    compared across melt_frac values.

    Returns tip_step, max_temp, energy_residual, n_hot_trace, t_axis.
    """
    p   = params if params is not None else BreakdownParams()
    rng = np.random.default_rng(seed)
    load = np.zeros((L, L))
    rho  = np.ones((L, L))
    momx = np.zeros((L, L))
    momy = np.zeros((L, L))
    E    = np.full((L, L), P0 / (GAMMA - 1.0))
    e0      = total_energy(E, load)
    driven  = 0.0
    tip_step = None
    max_temp = float((internal_energy(rho, momx, momy, E) / rho).max())
    t_axis, n_hot_trace, temp_trace, load_trace = [], [], [], []

    for t in range(steps):
        rho, momx, momy, E, load, ignited = step_melt(
            rho, momx, momy, E, load, p, rng, melt_frac=melt_frac,
        )
        driven += p.drive_sites * p.drive_amount
        if tip_step is None and ignited:
            tip_step = t
        cur = float((internal_energy(rho, momx, momy, E) / rho).max())
        max_temp = max(max_temp, cur)
        if t % sample_every == 0:
            e_spec = internal_energy(rho, momx, momy, E) / rho
            t_axis.append(t)
            n_hot_trace.append(int((e_spec >= p.e_ign).sum()))
            temp_trace.append(cur)
            load_trace.append(float(load.mean()))

    return {
        "L": L, "steps": steps, "tip_step": tip_step,
        "max_temp": max_temp,
        "energy_residual": total_energy(E, load) - (e0 + driven),
        "t_axis": t_axis,
        "n_hot_trace": n_hot_trace,
        "temp_trace": temp_trace,
        "load_trace": load_trace,
    }


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


def run_full_arc(
    L=80, steps_buildup=3000, steps_after=600, params=None, seed=0,
    sample_every=25, P0=1.0, rigidity=0.0,
):
    """Cold substrate -> slow buildup -> tip -> eruption front, all one continuous run.

    Uses step_rigid throughout (reduces to step when rigidity=0). The same drive
    and params continue after the tip -- no imposed phase change. Spatial shape is
    captured via per-snapshot scalar summaries (n_hot, max_temp, max_speed, mean_load).

    Returns:
        tip_step        int or None
        snapshots       list[dict] (t, n_hot, max_temp, max_speed, mean_load, phase)
        energy_residual float   E+load - (E0 + total_driven); should be ~0
        final_fields    dict    rho, momx, momy, E, load at end of run
    """
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
    snapshots = []

    def _snap(t, phase):
        e_spec = internal_energy(rho, momx, momy, E) / rho
        speed = np.sqrt((momx / rho) ** 2 + (momy / rho) ** 2)
        return {
            "t": t,
            "n_hot": int((e_spec >= p.e_ign).sum()),
            "max_temp": float(e_spec.max()),
            "max_speed": float(speed.max()),
            "mean_load": float(load.mean()),
            "phase": phase,
        }

    total_steps = steps_buildup + steps_after
    for t in range(total_steps):
        phase = "buildup" if t < steps_buildup else "eruption"
        rho, momx, momy, E, load, ignited = step_rigid(
            rho, momx, momy, E, load, p, rng, rigidity=rigidity,
        )
        driven += p.drive_sites * p.drive_amount

        if tip_step is None and ignited:
            tip_step = t
            snapshots.append(_snap(t, "tip"))
        elif t % sample_every == 0:
            snapshots.append(_snap(t, phase))

    return {
        "L": L,
        "tip_step": tip_step,
        "snapshots": snapshots,
        "energy_residual": total_energy(E, load) - (e0 + driven),
        "final_fields": {"rho": rho, "momx": momx, "momy": momy, "E": E, "load": load},
    }
