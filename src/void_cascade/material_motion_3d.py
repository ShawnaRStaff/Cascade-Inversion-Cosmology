"""Material motion in 3D: compressible Euler via Lax-Friedrichs.

Direct extension of material_motion_2d to 3D. Conserved variables per cell:
rho, momx, momy, momz, E. Pressure P=(gamma-1)*(E - 0.5*(px²+py²+pz²)/rho).
Periodic boundaries -- mass / momentum / energy are conserved exactly.

The LF average uses 6 face-neighbours (±x ±y ±z), keeping the same structure
as the 2D scheme. CFL condition: dt * max_wave_speed / dx < 1/sqrt(3).
"""

from __future__ import annotations

import numpy as np

GAMMA = 5.0 / 3.0


def internal_energy_3d(rho, momx, momy, momz, E):
    return E - 0.5 * (momx ** 2 + momy ** 2 + momz ** 2) / rho


def primitives_3d(rho, momx, momy, momz, E):
    u = momx / rho
    v = momy / rho
    w = momz / rho
    P = (GAMMA - 1.0) * (E - 0.5 * (momx * u + momy * v + momz * w))
    return u, v, w, P


def max_wave_speed_3d(rho, momx, momy, momz, E) -> float:
    u, v, w, P = primitives_3d(rho, momx, momy, momz, E)
    c = np.sqrt(np.maximum(GAMMA * P / rho, 0.0))
    return float(np.max(np.abs(u) + c))   # conservative bound (not full 3D norm)


def _flux_x_3d(rho, momx, momy, momz, E):
    u, v, w, P = primitives_3d(rho, momx, momy, momz, E)
    return momx, momx * u + P, momy * u, momz * u, (E + P) * u


def _flux_y_3d(rho, momx, momy, momz, E):
    u, v, w, P = primitives_3d(rho, momx, momy, momz, E)
    return momy, momx * v, momy * v + P, momz * v, (E + P) * v


def _flux_z_3d(rho, momx, momy, momz, E):
    u, v, w, P = primitives_3d(rho, momx, momy, momz, E)
    return momz, momx * w, momy * w, momz * w + P, (E + P) * w


def lax_friedrichs_step_3d(rho, momx, momy, momz, E, dx, dt):
    """One LF step for 3D compressible Euler with periodic boundaries.

    Each cell is updated as the average of its 6 face-neighbours minus a
    flux-divergence term. Conserves sum of each conserved variable exactly.
    """
    Fx = _flux_x_3d(rho, momx, momy, momz, E)
    Fy = _flux_y_3d(rho, momx, momy, momz, E)
    Fz = _flux_z_3d(rho, momx, momy, momz, E)
    U  = (rho, momx, momy, momz, E)

    def upd(Ui, Fxi, Fyi, Fzi):
        avg = (np.roll(Ui, -1, 0) + np.roll(Ui, 1, 0)
             + np.roll(Ui, -1, 1) + np.roll(Ui, 1, 1)
             + np.roll(Ui, -1, 2) + np.roll(Ui, 1, 2)) / 6.0
        dFx = np.roll(Fxi, -1, 0) - np.roll(Fxi, 1, 0)
        dFy = np.roll(Fyi, -1, 1) - np.roll(Fyi, 1, 1)
        dFz = np.roll(Fzi, -1, 2) - np.roll(Fzi, 1, 2)
        return avg - (dt / (2.0 * dx)) * (dFx + dFy + dFz)

    return tuple(upd(U[i], Fx[i], Fy[i], Fz[i]) for i in range(5))
