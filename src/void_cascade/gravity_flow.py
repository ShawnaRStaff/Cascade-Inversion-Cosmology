"""Self-gravitating reacting fluid: gravity as the collapse driver.

The force-free collapse drivers failed to focus a violent enough inversion;
gravity (which exists) does it cleanly -- self-attraction focuses overdensities
to a hot core that self-ignites. This is the proper version, with real
gravitational potential energy tracked so the conservation law is honest:

    total = E_fluid (internal + kinetic) + fuel + grav_PE   (~conserved)

where grav_PE = 1/2 * sum(rho * phi), phi from Poisson's eqn (div^2 phi =
4 pi G rho), solved by FFT on the periodic grid. Gravity does work on the
fluid (rho * v . g), sourced from grav_PE -- so as the gas collapses, grav_PE
drops and fluid energy rises, total held.

Combust: hot gas (specific internal energy >= e_ign) releases fuel into energy
(detonation). 2D reacting Euler (Lax-Friedrichs) + self-gravity. Honest scope:
abstract units; LF diffusive; gravity+LF conserve to ~percent, not machine
precision (reported, not claimed exact).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from void_cascade.material_motion_2d import (
    GAMMA, internal_energy, lax_friedrichs_step, max_wave_speed,
)


@dataclass(frozen=True)
class GravParams:
    G: float          # gravitational constant (sets collapse strength)
    e_ign: float      # specific internal energy to combust
    soften: float = 3.0  # gravitational softening length (cells) -- no r->0 singularity


def gravity_potential(rho, G, soften=3.0):
    """Solve div^2 phi = 4 pi G rho on the periodic grid via FFT, with Gaussian
    softening (smooths gravity over ~soften cells, removing the r->0 singularity
    that otherwise injects spurious energy at the collapse)."""
    L = rho.shape[0]
    k = 2.0 * np.pi * np.fft.fftfreq(L)
    kx, ky = np.meshgrid(k, k, indexing="ij")
    k2 = kx * kx + ky * ky
    k2[0, 0] = 1.0
    soft = np.exp(-0.5 * k2 * soften * soften)
    rho_k = np.fft.fft2(rho - rho.mean())   # zero-mean source (periodic Poisson)
    phi_k = -4.0 * np.pi * G * rho_k / k2 * soft
    phi_k[0, 0] = 0.0
    return np.real(np.fft.ifft2(phi_k))


def gravity_accel(rho, G, soften=3.0):
    phi = gravity_potential(rho, G, soften)
    gx = -(np.roll(phi, -1, 0) - np.roll(phi, 1, 0)) / 2.0
    gy = -(np.roll(phi, -1, 1) - np.roll(phi, 1, 1)) / 2.0
    return gx, gy, phi


def grav_PE(rho, phi):
    return 0.5 * float(np.sum(rho * phi))


def combust(rho, momx, momy, E, fuel, e_ign):
    e_spec = internal_energy(rho, momx, momy, E) / rho
    fire = (e_spec >= e_ign) & (fuel > 0.0)
    rel = np.where(fire, fuel, 0.0)
    return E + rel, fuel - rel, int(fire.sum())


def step(rho, momx, momy, E, fuel, p, dx, cfl):
    # combust (hot release)
    E, fuel, _ = combust(rho, momx, momy, E, fuel, p.e_ign)
    # gravity kick. Use the MIDPOINT velocity for the work term so the energy
    # added to E exactly equals the kinetic-energy change from the kick:
    #   dKE = 1/2 rho (u_new^2 - u^2) = rho g dt * 1/2 (u + u_new) = work.
    gx, gy, _phi = gravity_accel(rho, p.G, p.soften)
    s = max_wave_speed(rho, momx, momy, E)
    g_s = float(np.max(np.sqrt(gx * gx + gy * gy)))
    dt = cfl * dx / max(s + np.sqrt(g_s * dx), 1e-9)  # CFL incl. gravity
    u, v = momx / rho, momy / rho
    momx_new = momx + rho * gx * dt
    momy_new = momy + rho * gy * dt
    u_mid = 0.5 * (u + momx_new / rho)
    v_mid = 0.5 * (v + momy_new / rho)
    E = E + rho * (u_mid * gx + v_mid * gy) * dt      # = exact KE change from kick
    momx, momy = momx_new, momy_new
    # hydro
    rho, momx, momy, E = lax_friedrichs_step(rho, momx, momy, E, dx, dt)
    return rho, momx, momy, E, fuel, dt


def total_energy(rho, momx, momy, E, fuel, G, soften=3.0):
    phi = gravity_potential(rho, G, soften)
    return float(E.sum() + fuel.sum()) + grav_PE(rho, phi)


def run_gravity_collapse(L=120, G=0.5, e_ign=2.5, fuel0=3.0, rho0=1.0, P0=1.0,
                         bump=0.5, bump_r=18, steps=600, cfl=0.3, dx=1.0, band=4,
                         seed=0, snap_steps=(0,), soften=3.0) -> dict:
    """A mild overdensity in loaded fuel -> gravity collapses it -> hot core
    self-ignites -> detonation. The inversion, emergent, driven by gravity."""
    yy, xx = np.mgrid[0:L, 0:L].astype(float)
    r = np.sqrt((xx - L / 2) ** 2 + (yy - L / 2) ** 2)
    rho = rho0 * (1.0 + bump * np.exp(-(r / bump_r) ** 2))   # central overdensity
    momx = np.zeros((L, L)); momy = np.zeros((L, L))
    E = np.full((L, L), P0 / (GAMMA - 1.0))
    fuel = np.full((L, L), fuel0)
    p = GravParams(G=G, e_ign=e_ign, soften=soften)
    c = L // 2
    sl = (slice(c - band, c + band), slice(c - band, c + band))
    total0 = total_energy(rho, momx, momy, E, fuel, G, p.soften)

    t_trace, dens_c, temp_c, burned = [], [], [], []
    ignite_step = None
    snaps = {}
    for t in range(steps):
        if t in snap_steps:
            snaps[t] = rho.copy()
        before_fuel = fuel.sum()
        rho, momx, momy, E, fuel, dt = step(rho, momx, momy, E, fuel, p, dx, cfl)
        if ignite_step is None and fuel.sum() < before_fuel - 1e-9:
            ignite_step = t
        t_trace.append(t)
        dens_c.append(float(rho[sl].mean()))
        temp_c.append(float((internal_energy(rho, momx, momy, E) / rho)[sl].mean()))
        burned.append(float((fuel < fuel0 * 0.5).mean()))
    snaps[steps - 1] = rho.copy()

    total1 = total_energy(rho, momx, momy, E, fuel, G, p.soften)
    dens_c = np.asarray(dens_c)
    return {
        "L": L, "G": G, "e_ign": e_ign, "steps": steps,
        "total_energy_residual": total1 - total0,
        "total_energy_rel_residual": (total1 - total0) / (abs(total0) + 1e-9),
        "ignite_step": ignite_step,
        "peak_central_density": float(dens_c.max()),
        "central_density_initial": float(dens_c[0]),
        "peak_central_temp": float(max(temp_c)),
        "final_burned_fraction": float(burned[-1]),
        "t_trace": t_trace, "central_density": dens_c.tolist(),
        "central_temp": temp_c, "burned": burned, "snapshots": snaps,
    }
