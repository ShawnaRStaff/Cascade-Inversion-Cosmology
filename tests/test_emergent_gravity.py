"""Emergent gravity: can attraction arise from cascade mechanics alone?

Shawna's hunch ('motion made gravity'): gravitational-like focusing might
emerge from cascade dynamics without explicitly imposing Poisson gravity.

Three experiments, letting the model decide each one:

  1. Pressure direction: hot core -> LF step -> outward blast or inward pull?
  2. Load dispersion: central load overdensity -> Manna breakdown -> focuses or spreads?
  3. Density focusing: G=0 (pure cascade) vs G=0.5 (imposed Poisson) -- which focuses?

Expected honest result: cascade = diffusive SOC spreading + outward blast pressure.
Gravitational focusing requires an explicit attractive potential; it does NOT emerge
from cascade mechanics in this model.

If any test here FAILS to show dispersal / outward blast, that would be the
emergent-gravity signal we'd want to investigate further.
"""
import numpy as np
import pytest

from void_cascade.cascade_breakdown import (
    BreakdownParams,
    breakdown,
    step,
)
from void_cascade.gravity_flow import run_gravity_collapse
from void_cascade.material_motion_2d import (
    GAMMA,
    internal_energy,
    lax_friedrichs_step,
    max_wave_speed,
)


def _cold_uniform(L, P0=1.0):
    rho  = np.ones((L, L))
    momx = np.zeros((L, L))
    momy = np.zeros((L, L))
    E    = np.full((L, L), P0 / (GAMMA - 1.0))
    return rho, momx, momy, E


def _central_load(L, load_center=5.0, load_bg=0.0, r_core=4):
    yy, xx = np.mgrid[0:L, 0:L].astype(float)
    r = np.sqrt((xx - L / 2) ** 2 + (yy - L / 2) ** 2)
    return np.where(r < r_core, load_center, load_bg)


# ------------------------------------------------------------------
# 1. Pressure direction
# ------------------------------------------------------------------

def test_hot_core_blast_is_outward():
    """Hot core creates HIGH pressure -> momentum points AWAY from core.

    Gravitational attraction would require LOW pressure at the core pulling
    material IN. The hot cascade does the opposite: blast OUT.

    LF indexing: axis 0 = rows = x-direction; axis 1 = cols = y-direction.
    momx is momentum along axis 0 (rows). The hot patch is
    E[c-2:c+2, c-2:c+2], so rows c-2..c+1 are hot and row c+2 is the first
    cold row. After one LF step, the x-flux difference between the hot row
    c+1 and the cold row c+3 drives a positive momx at row c+2 (outward
    in the +x direction).  Checking at column c (centre column) ensures both
    y-direction neighbours of (c+2, c) are also cold -> clean 1D-like signal.
    """
    L = 32
    c = L // 2
    rho, momx, momy, E = _cold_uniform(L)
    E[c - 2 : c + 2, c - 2 : c + 2] = 8.0   # hot patch

    s  = max_wave_speed(rho, momx, momy, E)
    dt = 0.4 / max(s, 1e-9)
    _, momx2, _, _ = lax_friedrichs_step(rho, momx, momy, E, 1.0, dt)

    # Row c+2 (first cold row on the + side): x-pressure wave → positive momx
    assert momx2[c + 2, c] > 0.0, (
        f"first cold row above hot patch must be pushed outward (+x); "
        f"momx at (c+2, c) = {momx2[c+2, c]:.6f}"
    )
    # Row c-3 (first cold row on the - side): pushed outward in -x
    assert momx2[c - 3, c] < 0.0, (
        f"first cold row below hot patch must be pushed outward (-x); "
        f"momx at (c-3, c) = {momx2[c-3, c]:.6f}"
    )


def test_cold_surroundings_do_not_pull_inward():
    """Even the cold shell surrounding the hot core does not create inward pull.

    After one LF step from a hot core in cold surroundings, the net momentum
    at a midpoint between the core and the boundary is outward (not inward).
    Inward would imply attractive dynamics -- that is what we are testing for.
    """
    L = 32
    c = L // 2
    rho, momx, momy, E = _cold_uniform(L, P0=1.0)
    E[c - 3 : c + 3, c - 3 : c + 3] = 6.0   # hot core

    s  = max_wave_speed(rho, momx, momy, E)
    dt = 0.4 / max(s, 1e-9)
    _, momx2, momy2, _ = lax_friedrichs_step(rho, momx, momy, E, 1.0, dt)

    # Midpoint between core edge (c+3) and boundary (c+8): should move rightward
    assert momx2[c, c + 6] >= 0.0, (
        "midpoint shell should not be pulled toward the core (no inward pull)"
    )


# ------------------------------------------------------------------
# 2. Load dispersion (SOC Manna)
# ------------------------------------------------------------------

def test_central_load_overdensity_disperses():
    """SOC Manna avalanches DIFFUSE load away from a central overdensity.

    The Manna rule sheds load to RANDOM neighbours -- no preferred direction,
    no bias toward already-loaded cells. An overdensity drains to its
    neighbours and spreads outward. This is the opposite of gravitational
    focusing (which would concentrate load further at the centre).
    """
    L = 32
    rho, momx, momy, E = _cold_uniform(L)
    load = _central_load(L, load_center=6.0, load_bg=0.0, r_core=4)

    p   = BreakdownParams(drive_sites=0, drive_amount=0.0, thr=2.0, hpc=0.1)
    rng = np.random.default_rng(0)

    c   = L // 2
    sl  = (slice(c - 3, c + 3), slice(c - 3, c + 3))
    initial_center_load = float(load[sl].mean())

    for _ in range(300):
        rho, momx, momy, E, load, _ = step(rho, momx, momy, E, load, p, rng)

    final_center_load = float(load[sl].mean())
    assert final_center_load < initial_center_load * 0.9, (
        f"SOC diffusion should drain load from centre: "
        f"initial {initial_center_load:.3f} -> final {final_center_load:.3f}"
    )


def test_load_spreads_to_periphery():
    """Load that starts at centre must appear in the immediately adjacent ring.

    Complement to the above: SOC diffusion not only drains the centre but
    carries load OUTWARD into neighbouring cells. A purely attractive dynamic
    would SUPPRESS outflow.

    Key mechanics: without driving, the cascade finishes in the FIRST step
    (all cells drop below threshold), and load is frozen thereafter. The
    check must therefore target the ring immediately outside the initial ball
    (r = r_core to r_core+4), which the first-step avalanche always reaches.
    Use hpc=0 so load is conserved (cannot disappear into heat).
    """
    L  = 32
    r_core = 3
    rho, momx, momy, E = _cold_uniform(L)
    yy, xx = np.mgrid[0:L, 0:L].astype(float)
    r = np.sqrt((xx - L / 2) ** 2 + (yy - L / 2) ** 2)
    load = np.where(r < r_core, 4.0, 0.0)

    # hpc=0: load is conserved (no conversion to heat), so it MUST spread
    p   = BreakdownParams(drive_sites=0, drive_amount=0.0, thr=2.0, hpc=0.0)
    rng = np.random.default_rng(1)

    # Ring immediately adjacent to the hot ball (always reachable in 1 cascade)
    annulus = (r >= r_core) & (r < r_core + 4)
    initial_ring_load = float(load[annulus].mean())   # = 0.0

    for _ in range(10):   # the cascade completes in step 1; a few extra steps suffice
        rho, momx, momy, E, load, _ = step(rho, momx, momy, E, load, p, rng)

    final_ring_load = float(load[annulus].mean())
    assert final_ring_load > initial_ring_load, (
        f"SOC avalanche must spread load into the adjacent ring: "
        f"initial {initial_ring_load:.4f} -> final {final_ring_load:.4f}"
    )


# ------------------------------------------------------------------
# 3. Density focusing: G=0 cascade vs G>0 gravity
# ------------------------------------------------------------------

def test_pure_cascade_does_not_focus_density():
    """Without Poisson gravity (G=0), the central overdensity does NOT grow.

    With G=0, gravity_flow.py runs only LF hydro + combust (no gravity kicks).
    The hot blast from combust pushes material OUTWARD, so peak_central_density
    should not significantly exceed the initial value.

    Contrast with test_run_gravity_collapse_central_density_rises (G=0.5),
    where gravity drives density up -- that result requires imposed gravity.
    """
    r = run_gravity_collapse(
        L=64, G=0.0, e_ign=2.5, fuel0=3.0, steps=200, cfl=0.3, seed=0
    )
    ratio = r["peak_central_density"] / r["central_density_initial"]
    # Blast may transiently compress the central band slightly, so allow 30% up.
    # Gravitational focusing at G=0.5 gives ~2x or more.
    assert ratio < 1.5, (
        f"pure cascade (G=0) should NOT significantly focus density; "
        f"got peak/initial = {ratio:.3f}. If > 1.5, investigate whether "
        f"this is a genuine emergence signal."
    )


def test_gravity_focuses_cascade_does_not():
    """Imposed gravity focuses; pure cascade disperses. Documents the contrast.

    G=0.5 central density rises (test_run_gravity_collapse_central_density_rises).
    G=0.0 central density does not rise comparably.
    The gap between them is the contribution of imposed Poisson gravity;
    nothing in the cascade mechanics closes that gap.
    """
    r_no_grav = run_gravity_collapse(
        L=64, G=0.0, e_ign=2.5, fuel0=3.0, steps=200, cfl=0.3, seed=0
    )
    r_gravity  = run_gravity_collapse(
        L=64, G=0.5, e_ign=2.5, fuel0=3.0, steps=200, cfl=0.3, seed=0
    )

    focus_no_grav = r_no_grav["peak_central_density"] / r_no_grav["central_density_initial"]
    focus_gravity  = r_gravity["peak_central_density"]  / r_gravity["central_density_initial"]

    assert focus_gravity > focus_no_grav * 1.5, (
        f"gravity (G=0.5) should focus density much more than pure cascade (G=0): "
        f"gravity ratio={focus_gravity:.3f}, cascade ratio={focus_no_grav:.3f}"
    )
