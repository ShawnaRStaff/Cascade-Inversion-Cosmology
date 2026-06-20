# Literature grounding: what's anchored vs. what's our own spin

**Date:** 2026-06-19. Purpose: keep the research honest. For every piece of
the model, say plainly whether it rests on published work or is our own
original synthesis. Having a spin is fine; pretending the spin is established
is not.

Three buckets: **grounded** (published work supports it), **thematic**
(a published program shares the idea but with a different mechanism), and
**original** (ours — sometimes with general-physics analogs, but not in the
cited cosmology sources).

---

## Grounded in the literature

- **SOC avalanche dynamics (the sandpile core).** Bak–Tang–Wiesenfeld 1987;
  Manna stochastic sandpile. We use the validated Manna rule, and M1–M3
  reproduced the *published exponents* (1D Oslo τ≈1.55, 2D Manna τ≈1.27, 3D
  Manna τ≈1.35, conservation D(2−τ)≈1). This is grounded with numbers.
  **Caveat found 2026-06-19:** the grounding check on the *new 1D coupled*
  model (`cascade_soc_check.py`) did NOT reproduce clean SOC — it gave
  tau~0.78 with an abnormally large mean avalanche (~5x L). That is the
  signature of 1D Manna being a degenerate sandpile case (exactly why this
  project used the Oslo rule, not Manna, in 1D — see M1). So `cascade_heat`
  in 1D is fine for asking the heat-tipping *question*, but it is NOT a clean
  SOC demonstrator, and its quantitative heat numbers (critical heat-making,
  how easily it tips) are likely exaggerated by the over-large 1D avalanches.
  Proper SOC grounding — and trustworthy heat numbers — needs 2D Manna (where
  M2 validated tau~1.27). This merges with the planned 2D heat-front build.
- **Percolation framework.** Stauffer & Aharony (percolation theory);
  cosmological phase-transition percolation Guth–Tye 1980, Turner–Weinberg–
  Widrow 1992, Gould–Tenkanen 2021. The *framework* is theirs; our measured
  cascade-driven threshold (p_c≈0.177, ~half the random 0.31) is our result.
- **SOC applied to cosmology — the precedent exists.** Moffat 1997 ("A
  Self-Organized Critical Universe") is a literal sandpile cosmology;
  Carfora–Marzuoli 2023 ("Primordial Cosmology from SOC"); Aschwanden 2011
  (astrophysical SOC). So "SOC ↔ cosmology" is a real, published research
  program — we are not inventing that link.
- **Galaxy two-point correlation comparison.** Peebles' ξ(r) slope γ≈1.8.
  M4's γ=1.80 at p≈0.65 is our most direct contact with observation, and it
  aligns with the *spectral* target the SOC-cosmology literature actually
  uses (see drift note below).
- **Mass-function comparison.** Press–Schechter (small-mass slope ~2.0);
  M4's τ_s≈2.05 near p_c lands here.

## Thematic (shared idea, different mechanism)

- **Cascade vacuum relaxation — Lukash & Mikheeva 2024 (arXiv:2506.03226).**
  Named in the README as a foundation, but read directly it describes a
  *sequential relaxation of vacuum fields through cosmological epochs* — an
  ordered chain, **not** SOC avalanches, no sandpile, no singularity. Shares
  the "vacuum energy in steps" theme; the mechanism is not ours. Honest
  status: thematic inspiration, not a mechanism we implement.
- **Non-singular bounce — Ashtekar–Pawlowski–Singh 2006 (LQC bounce).** A
  published precedent for "no singularity, a bounce instead," which matches
  the researcher's inversion picture in spirit. Different mechanism (loop
  quantum gravity vs. substrate fracture).

## Original (our spin — general-physics analogs, not in the cosmology sources)

- **Fracture → percolation → "Big Bang" mapping.** Reading each topple as a
  fracture and the percolation of the ever-toppled set as the inversion.
  This specific observable is the project's own move; the SOC-cosmology
  sources measure spectra/geometry, not fracture-percolation.
- **Heat-gated fracture / thermal-runaway catastrophe** (`cascade_heat.py`).
  Cold "freezes breaks in place," heat releases them, clustered-avalanche
  heat tips a runaway. Strong **general-physics** analogs — thermal runaway,
  fatigue-to-failure, first-order transition latent heat, detonation — but
  **no cosmology source** does this. Original.
- **Absolute-zero "looks intact" healing without bonds.** Cold = no motion =
  breaks can't separate, so it appears whole. Bond-free, consistent with our
  no-bonds finding. Our framing.
- **Expansion as a propagating front through infinite substrate.** The
  project's own "expansion-as-propagation" contribution; aligns in spirit
  with the bubble-nucleation / phase-transition-front family (Guth–Tye), but
  the substrate-front formulation is ours. (Sketched, not built.)

---

## The honest drift to watch

The published SOC-cosmology programs (Moffat, Carfora–Marzuoli) put the
cosmological *observable* in the **fluctuation spectrum / multifractal
geometry / dark-matter genesis**. Our recent heat-runaway work measures a
**catastrophe-dynamics** observable the sources do not. So while the SOC
*foundation* is well-grounded, our newest results have moved away from what
the literature actually compares against. Our strongest literature/
observational alignment remains the **spectral side (M4, γ=1.80)**, not the
heat catastrophe.

## What grounding still requires

- **Dimensional calibration.** Every cosmological comparison so far is
  shape-only (slopes, thresholds), not magnitude. Turning the model's
  numbers (z=0.616, p_c=0.177, the heat/melt scales) into physical units is
  the long-standing rate-limiter for any quantitative literature or
  observational test. Untouched.
- **Observational test of the heat catastrophe.** Internally validated by
  tests, but no comparison to any real data yet. Until then it's a
  mechanism, not a prediction.

## One-line status

Foundation (SOC, percolation, sandpile): **grounded, with reproduced
exponents.** Cosmological framing and the heat catastrophe: **original spin
with general-physics analogs, not yet anchored to the sources' observables
or to data.**
