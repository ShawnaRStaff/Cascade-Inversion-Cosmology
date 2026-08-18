"""Try all three collapse drivers, honestly, favoring none.

Each must answer: does it produce a FOCUSED collapse that compresses the gas
enough to SELF-IGNITE (reach the ignition temperature with NO imposed
converging flow and NO imposed ignition) and detonate?

  (a) fracture-cascade : loaded cells fracture and the fracture CASCADES to
      loaded neighbours -> a coherent low-pressure void -> substrate collapses
      into it. No new force.
  (b) gravity          : self-gravity (FFT Poisson) pulls overdensities
      together -> focuses automatically. (Energy: gravity does work; we report
      fluid energy only, so the residual here is NOT a conservation claim.)
  (c') seeded void     : a pre-existing central void; loaded substrate collapses
      into it (focused by geometry, but the void is seeded/imposed).

Reports for each: did it self-ignite, peak compression, peak temperature
(specific internal energy) vs the ignition threshold.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from void_cascade.material_motion_2d import (  # noqa: E402
    GAMMA, internal_energy, lax_friedrichs_step, max_wave_speed,
)

E_IGN = 2.5
FUEL0 = 3.0
P0 = 1.0
E0 = P0 / (GAMMA - 1.0)        # initial internal energy (specific = 1.5 at rho=1)


def combust(rho, momx, momy, E, fuel, e_ign):
    e_spec = internal_energy(rho, momx, momy, E) / rho
    fire = (e_spec >= e_ign) & (fuel > 0.0)
    rel = np.where(fire, fuel, 0.0)
    return E + rel, fuel - rel, int(fire.sum())


def _evolve(rho, momx, momy, E, fuel, steps, e_ign, force=None):
    cb_tot = 0
    max_rho = float(rho.max())
    max_es = float((internal_energy(rho, momx, momy, E) / rho).max())
    for _ in range(steps):
        E, fuel, nc = combust(rho, momx, momy, E, fuel, e_ign)
        cb_tot += nc
        s = max_wave_speed(rho, momx, momy, E)
        dt = 0.4 / max(s, 1e-9)
        if force is not None:
            ax, ay = force(rho)
            momx = momx + rho * ax * dt
            momy = momy + rho * ay * dt
            E = E + (momx * ax + momy * ay) * dt  # work by the force
        rho, momx, momy, E = lax_friedrichs_step(rho, momx, momy, E, 1.0, dt)
        max_rho = max(max_rho, float(rho.max()))
        max_es = max(max_es, float((internal_energy(rho, momx, momy, E) / rho).max()))
    return cb_tot, max_rho, max_es


def fracture_cascade(L=120, steps=500, e_ign=E_IGN, seed=0):
    rng = np.random.default_rng(seed)
    rho = np.ones((L, L)); momx = np.zeros((L, L)); momy = np.zeros((L, L))
    E = np.full((L, L), E0); fuel = np.zeros((L, L)); fractured = np.zeros((L, L), bool)
    drive_rate, fracture_fuel, void_frac, floor = 0.04, 1.0, 0.95, 0.1
    cb_tot = 0; max_rho = 1.0; max_es = 1.5
    for _ in range(steps):
        # buildup at random sites (spatial variation)
        flat = fuel.ravel(); flat[rng.integers(0, L * L, size=L * L // 15)] += drive_rate
        # fracture with cascade: relax until no loaded cold cell triggers
        for _sweep in range(8):
            e_spec = internal_energy(rho, momx, momy, E) / rho
            frac = (fuel >= fracture_fuel) & (~fractured) & (e_spec < e_ign)
            if not frac.any():
                break
            internal = internal_energy(rho, momx, momy, E)
            dE = np.where(frac, void_frac * np.maximum(internal - floor, 0.0), 0.0)
            E -= dE; fuel += dE; fractured |= frac
            # cascade: a fractured cell pushes a little fuel to neighbours (loads them)
            push = np.zeros((L, L)); pf = 0.3 * fuel * frac
            push[1:] += pf[:-1]; push[:-1] += pf[1:]; push[:, 1:] += pf[:, :-1]; push[:, :-1] += pf[:, 1:]
            fuel = fuel - pf + push - 0  # conserve fuel in the push (4 neighbours get pf/... approx)
        E, fuel, nc = combust(rho, momx, momy, E, fuel, e_ign); cb_tot += nc
        s = max_wave_speed(rho, momx, momy, E); dt = 0.4 / max(s, 1e-9)
        rho, momx, momy, E = lax_friedrichs_step(rho, momx, momy, E, 1.0, dt)
        max_rho = max(max_rho, float(rho.max())); max_es = max(max_es, float((internal_energy(rho, momx, momy, E) / rho).max()))
    return cb_tot, max_rho, max_es


def gravity_collapse(L=120, steps=500, e_ign=E_IGN, G=0.5, seed=0):
    rng = np.random.default_rng(seed)
    rho = 1.0 + 0.2 * rng.standard_normal((L, L)); rho = np.clip(rho, 0.2, None)
    momx = np.zeros((L, L)); momy = np.zeros((L, L)); E = np.full((L, L), E0); fuel = np.full((L, L), FUEL0)
    k = 2 * np.pi * np.fft.fftfreq(L)
    kx, ky = np.meshgrid(k, k, indexing="ij"); k2 = kx * kx + ky * ky; k2[0, 0] = 1.0

    def force(r):
        rk = np.fft.fft2(r - r.mean())
        phik = -4 * np.pi * G * rk / k2; phik[0, 0] = 0
        ax = np.real(np.fft.ifft2(-1j * kx * phik))
        ay = np.real(np.fft.ifft2(-1j * ky * phik))
        return ax, ay

    return _evolve(rho, momx, momy, E, fuel, steps, e_ign, force=force)


def seeded_void(L=120, steps=500, e_ign=E_IGN, void_r=22):
    yy, xx = np.mgrid[0:L, 0:L].astype(float)
    r = np.sqrt((xx - L / 2) ** 2 + (yy - L / 2) ** 2)
    void = r < void_r
    rho = np.where(void, 0.1, 1.0); momx = np.zeros((L, L)); momy = np.zeros((L, L))
    P = np.where(void, 0.1, 2.0); E = P / (GAMMA - 1.0); fuel = np.full((L, L), FUEL0)
    return _evolve(rho, momx, momy, E, fuel, steps, e_ign)


def main() -> None:
    print(f"=== Collapse drivers (ignition threshold e_ign={E_IGN}, start temp 1.5) ===\n")
    print(f"{'driver':<20}{'self-ignited?':>14}{'peak compression':>18}{'peak temp':>12}")
    for name, fn in [("fracture-cascade", fracture_cascade),
                     ("gravity", gravity_collapse),
                     ("seeded-void", seeded_void)]:
        cb, mr, me = fn()
        print(f"{name:<20}{('YES' if cb > 0 else 'no'):>14}{mr:>17.2f}x{me:>12.2f}")
    print(f"\n(self-ignited = combustion occurred from emergent heating; peak temp must")
    print(f" exceed {E_IGN} for ignition. start temp 1.5.)")


if __name__ == "__main__":
    main()
