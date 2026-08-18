# Dimensional calibration: the bridge to real numbers (and why it's the rate-limiter)

**UPDATE 2026-08-17 (v2, supersedes the check below):** the 2026-06-19
"16x c" result conflated the model's TWO timescales — it converted
per-SWEEP front speeds with the per-DROP anchor (2.6 Myr). With sweep
time and drop time separated (`src/void_cascade/calibration.py`, 8
tests; `scripts/run_calibration.py` v2), the honest result is sharper:

**The four requirements (a) material causality on avalanche fronts,
(b) SOC timescale separation, (c) saturated era within 13.8 Gyr,
(d) 8.5 Mpc/cell from the xi(r) match are mutually exclusive by ~7
orders of magnitude.** Measured peak-event durations (T = 1319 ->
8039 sweeps for L=48 -> 128, T ~ L^1.84) under the causal sweep-time
floor (41.6 Myr at 8.5 Mpc cells) give single events lasting 55-334
Gyr — 4-24x the age of the universe, worsening with L. Conversely,
SOC separation plus the epoch match (~4,600 yr/drop) caps the sweep
below ~1 yr, which with causality caps the cell at ~0.2 pc, 4x10^7
below the LSS anchor.

Exactly one requirement must be dropped; the three options are worked
in the script output. The framework-preferred reading drops (a):
avalanche fronts are activity/pattern fronts (phase-velocity-like, no
material transport, superluminality permitted), while material bounds
apply only to the fluid-layer experiments. Under that reading, all
remaining constraints coexist with zero violations, and **front
superluminality becomes a REQUIRED claim of the framework, not an
embarrassment to hide** — it must be stated in any writeup.

Energy anchor (first pinned candidate): identifying z = 0.616
grains/cell with the vacuum energy density at 8.5 Mpc cells gives
1 grain = 1.6x10^61 J; the L=96 peak event moves 2.8x10^67 J (~0.2% of
the observable universe's mass-energy), and an observable-universe-
sized box needs L ~ 1650 (the L=128 box spans 1.09 Gpc).

---

**Date:** 2026-06-19. Status: framework + a check, not a result. The model is
dimensionless; reaching physical numbers needs anchors that are FREE until
pinned by data. `scripts/run_calibration.py` does the illustrative arithmetic.

## What needs anchoring
The model produces only dimensionless quantities (e.g. avalanche exponent
tau=1.27; saturation density z=0.616; join-up p_c=0.177; front speeds in
cells/step; critical loading ~0.7; critical fuel ~1.0; implosion compression
1.79x, heat spike 3.84x). To get physical numbers we need THREE anchors:

- **length:** Mpc per lattice cell
- **time:** years per step
- **energy:** Joules (or Kelvin) per model energy unit

## The honest problem: anchors are free, and the existing ones don't fit
Using the only anchors we have -- the *free* ones M4 set by matching galaxy
xi(r) (1 cell ~ 8.5 Mpc; 1 step ~ 2.6 Myr) -- the arithmetic gives:

    1 cell/step  =  3.20e6 km/s  =  10.7 x the speed of light

So **every model speed of order 1 cell/step is ~5-16x light-speed** under these
anchors:
- conserved-energy front (~0.54 cells/step): 5.8 x c
- detonation front (~1.0 cells/step): 10.7 x c
- implosion *material* velocity (u0 ~1.5 cells/step): 16 x c

Implications, stated straight:
1. An **activity/pattern front** (nothing material moving -- our "activity
   wavefront" label) *may* be superluminal, like a phase velocity. So the
   expansion front being >c is not automatically fatal.
2. **Material motion is not** -- the implosion fluid at 16x c is impossible.
   So the M4 anchors are **not mutually consistent with material velocities**:
   to make the implosion sub-light you'd need ~16x larger time/step (or
   smaller Mpc/cell). The length and time anchors are over-constrained.
3. The **energy anchor is entirely unpinned** -- nothing yet ties model energy
   to a temperature or an energy density.

## Why this is the rate-limiter
Calibration is **underdetermined** (free anchors) AND the existing free anchors
already **strain** (superluminal material). Until length, time, AND energy are
pinned from real observations -- e.g. a temperature for the plasma spike, an
energy density for z=0.616, a real length scale that's also time-consistent --
**no quantitative cosmological claim is possible**, only shape/dimensionless
comparisons. This is exactly what the project's earlier notes flagged: the
shape-only ceiling.

## What pinning each anchor would take
- **length:** a structure scale matched to a survey (M4's xi(r) r_0 is a start,
  but it's a fit, not a derivation).
- **time:** an event-rate or expansion-rate matched to an observed timescale --
  and it must be *jointly* consistent with the length anchor (the superluminal
  result shows it currently isn't).
- **energy:** the plasma spike's temperature, or z=0.616 as a physical energy
  density (the first quantitative constant the model offers as a candidate).

## Honest bottom line
The calibration *check* is itself informative: it says the current free anchors
imply superluminal material, so they must be wrong/revised, and the energy
anchor is missing. Grounding the model to reality is a real, unfinished task --
and it gates everything downstream. We did not fake numbers to dodge it.
