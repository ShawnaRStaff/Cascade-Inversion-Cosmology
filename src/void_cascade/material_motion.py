"""Material motion: the implosion as a compressible fluid (1D).

The earlier models move energy, not material. The inversion in the theory is
material *falling in on itself*, heating to plasma, and rebounding out -- that
needs MOMENTUM. The minimal honest way to get it is a compressible fluid:
mass + momentum + energy, with pressure. Then implosion -> compression-heating
-> rebound EMERGES from the physics; we don't impose it.

We solve the 1D Euler equations (the continuum law for "material that moves
and has pressure") with the Lax-Friedrichs scheme (simple, robust). Conserved
variables per cell:
    rho   (mass density)
    mom   (momentum density = rho*u)
    E     (total energy density = internal + kinetic)
Pressure (ideal gas, gamma=5/3, monatomic / plasma-like):
    P = (gamma-1) * (E - 0.5*rho*u^2)
Internal energy (the "heat") = E - 0.5*rho*u^2; compression converts kinetic
into it (material slamming together heats up), rebound converts it back --
the potential<->kinetic shift, now via real fluid motion.

Periodic boundaries -> mass, momentum, energy conserved to machine precision.
Honest scope: this is generic compressible-fluid physics (the continuum limit
of "moving material with pressure"), connected to the theory as the substrate
imploding/rebounding. 1D first.
"""

from __future__ import annotations

import numpy as np

GAMMA = 5.0 / 3.0


def primitives(rho: np.ndarray, mom: np.ndarray, E: np.ndarray):
    u = mom / rho
    P = (GAMMA - 1.0) * (E - 0.5 * mom * u)
    return u, P


def internal_energy(rho: np.ndarray, mom: np.ndarray, E: np.ndarray) -> np.ndarray:
    """Heat content: total energy minus bulk kinetic."""
    return E - 0.5 * mom * mom / rho


def _fluxes(rho, mom, E):
    u, P = primitives(rho, mom, E)
    return mom, mom * u + P, (E + P) * u


def max_wave_speed(rho, mom, E) -> float:
    u, P = primitives(rho, mom, E)
    c = np.sqrt(np.maximum(GAMMA * P / rho, 0.0))
    return float(np.max(np.abs(u) + c))


def lax_friedrichs_step(rho, mom, E, dx: float, dt: float):
    """One periodic Lax-Friedrichs update (conservative flux form)."""
    Fr, Fm, FE = _fluxes(rho, mom, E)

    def upd(U, F):
        return 0.5 * (np.roll(U, -1) + np.roll(U, 1)) - (dt / (2.0 * dx)) * (
            np.roll(F, -1) - np.roll(F, 1)
        )

    return upd(rho, Fr), upd(mom, Fm), upd(E, FE)


def converging_initial_conditions(N: int, u0: float, rho0: float = 1.0, P0: float = 1.0):
    """Material flowing toward the centre (an implosion), smooth and periodic
    (velocity zero at the centre and at the seam)."""
    x = np.arange(N)
    c0 = N / 2.0
    u = -u0 * np.sin(2.0 * np.pi * (x - c0) / N)  # converges toward centre
    rho = np.full(N, rho0, dtype=float)
    mom = rho * u
    E = P0 / (GAMMA - 1.0) + 0.5 * rho * u * u
    return rho, mom, E


def run_implosion(N: int = 200, u0: float = 1.0, steps: int = 400, cfl: float = 0.4,
                  dx: float = 1.0, band: int = 5, seed: int | None = None) -> dict:
    """Implode (converging flow) and watch the centre compress, heat, rebound."""
    rho, mom, E = converging_initial_conditions(N, u0)
    c0 = N // 2
    sl = slice(c0 - band, c0 + band)
    mass0 = float(rho.sum())
    energy0 = float(E.sum())
    mom0 = float(mom.sum())

    t_trace, dens_c, heat_c, vel_c = [], [], [], []
    min_pressure = np.inf
    for t in range(steps):
        s = max_wave_speed(rho, mom, E)
        dt = cfl * dx / max(s, 1e-9)
        rho, mom, E = lax_friedrichs_step(rho, mom, E, dx, dt)
        u, P = primitives(rho, mom, E)
        min_pressure = min(min_pressure, float(P.min()))
        t_trace.append(t)
        dens_c.append(float(rho[sl].mean()))
        heat_c.append(float(internal_energy(rho, mom, E)[sl].mean()))
        vel_c.append(float(u[sl].mean()))

    dens_c = np.asarray(dens_c)
    peak_i = int(dens_c.argmax())
    return {
        "N": N, "u0": u0, "steps": steps,
        "mass_residual": float(rho.sum() - mass0),
        "energy_residual": float(E.sum() - energy0),
        "momentum_residual": float(mom.sum() - mom0),
        "min_pressure": float(min_pressure),
        "min_density": float(rho.min()),
        "peak_central_density": float(dens_c.max()),
        "peak_density_step": peak_i,
        "central_density_initial": float(dens_c[0]),
        "central_heat_initial": float(heat_c[0]),
        "central_heat_at_peak": float(heat_c[peak_i]),
        "rebounded": bool(peak_i < steps - 1 and dens_c[-1] < dens_c[peak_i]),
        "t_trace": t_trace, "central_density": dens_c.tolist(),
        "central_heat": heat_c, "central_velocity": vel_c,
    }
