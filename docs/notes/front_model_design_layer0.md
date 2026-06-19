# Design sketch: propagating-front model, Layer 0

**Date:** 2026-06-19
**Status:** Design only. No code yet. This is the contract we'd build to.

Plain language. A reminder of why we're here: the self-healing sandpile
could only stay active by recharging, and we proved (M5 closure) that
recharge was a fake engine. The ongoing expansion in the theory feeds
itself a different way — a wave that burns through fresh, pre-loaded
substrate. Layer 0 is the smallest honest model that can answer ONE
question:

> Does a self-feeding, one-way wave sustain itself on the substrate's own
> stored fuel — or does it fizzle out?

Built so it CAN fizzle. If we rig it to always propagate, we've learned
nothing (that was the M5 mistake).

---

## The picture in one paragraph

Every cell starts holding some stored "darkness" (fuel) and is "intact."
A wave starts at one spot. When a cell flips, it is done forever (one-way,
no healing — the honest rule we just validated) and it dumps its stored
fuel into its still-intact neighbours. If a neighbour collects enough fuel
to cross a threshold, it flips too, dumping its own fuel onward. The wave
lives only as long as each flip delivers enough to ignite the next. Like a
flame: it keeps moving if the fuel ahead is rich enough, and dies if it
isn't.

---

## Data contract (data is just data — flat, typed)

State (all arrays the same shape; start 1D, then 2D/3D):

- `load`      : array of numbers >= 0. Stored fuel per cell. **Not a free
                knob** — taken from the measured saturated substrate
                (z_avg = 0.616; cells are 0 or 1, ~62% loaded). The buildup
                phase the sandpile *legitimately* modelled (M1-M3) produces
                this loaded substrate; the front model consumes it.
- `flipped`   : bool array. One-way: False -> True only, never back.
- `received`  : array of numbers >= 0. Fuel delivered to a cell so far by
                flipped neighbours; compared against the threshold.

Parameters (declared, principled, never tuned to win):

- `theta`     : the flip threshold — how much received fuel ignites a cell.
                We **sweep** this to find the critical value, rather than
                picking one that propagates.

Bookkeeping: `front_position` (furthest flipped extent), `step`, and an
`energy_lost` counter for fuel that falls off the edge.

## The one spread rule (the engine)

One sweep:
1. Find intact cells whose `received >= theta`.
2. Flip them (mark `flipped`, permanently).
3. Each flipped cell distributes its `load` to its still-intact neighbours,
   added to their `received`. **Energy is conserved** — a cell releases
   exactly what it held, no more. (This conservation is the main anti-cheat:
   if flips could create fuel from nothing, propagation would be rigged.)
4. Repeat until no intact cell is over threshold (wave stalls) or the wave
   reaches the domain edge (wave sustains / we extend the domain).

That's it. No bonds. No external driving (the wave runs on stored fuel, not
on us poking it). No healing.

## Room to grow

Start with a long loaded region; the front advances into fresh fuel. We
approximate "infinite" by making the domain long enough that the front
never reaches the end during the run, or by following the front with a
moving window. **Watch for finite-size artifacts** — we got burned by those
before (the L-ceiling), so any "it sustains" claim must be checked against
domain size.

---

## The fizzle guardrail (built in from the start)

- **Conservation**: flips release exactly their stored fuel. No creation.
- **theta is swept, not chosen.** We locate the critical threshold
  `theta*` where behaviour flips from fizzle to sustain. Where the real
  substrate's fuel (0.616) sits relative to `theta*` is the RESULT, not an
  input.
- **Sanity test that it CAN die**: at high `theta`, the wave must stall. If
  it never dies for any `theta`, the model is rigged and we'd know the code
  is wrong, not the physics.

A likely honest twist: because the fuel is patchy (~38% of cells are empty),
whether the wave keeps finding ignitable cells is partly a **percolation**
question on the fuel — does a rich-enough connected path exist? That ties
straight back to our percolation work and is exactly the kind of
un-riggable condition we want.

---

## What we measure (so it can fail)

1. **Sustain or fizzle?** Front position vs time: advances forever, or stalls.
2. **The critical threshold `theta*`** — and which side the real substrate's
   fuel (0.616) lands on.
3. **Front speed** — position vs time: steady, speeding up, slowing? This is
   the seed of an expansion-vs-distance relation — i.e. the old Milestone-6
   "expansion-as-propagation / Hubble" question lives here.
4. (2D/3D, later) **front shape / roughness** — does it look like the cosmic
   web?

## Honest limits of Layer 0 (say them up front)

- This tests **propagation only** — NOT the falling-in / implosion. There is
  still no motion, no gravity, no material moving inward. Layer 0 is the
  *expansion-as-wavefront* part of the theory, which (per the README) is
  explicitly "the wavefront of activity, not material moving." The implosion
  is a separate, bigger, motion-requiring layer we do not touch yet.
- It is **sandpile-adjacent** and we should admit that: the new things are
  (a) irreversibility, (b) pre-loaded fuel instead of outside poking, and
  (c) self-triggering. The honest novelty is "no external driving — it runs
  on stored fuel." We should be alert to it being a relabeled sandpile and
  test for real differences.
- No heat, radiation, matter-formation, or expanding space. Later layers,
  only if a result demands them.

## Build order (mirrors M1 -> M3, cheap first)

1. **1D** with fuel drawn from the measured saturation (0/1 at ~62%): the
   cleanest, fastest test of sustain-vs-fizzle and the critical `theta*`.
2. **2D/3D** using the real saved saturated field for front shape and the
   speed/Hubble relation.
3. Tests first each step: conservation, one-way invariant, must-be-able-to-
   fizzle (high theta stalls), sustains for low theta, critical theta*
   exists.

## Not a rewrite

Reuse the irreversible-flip rule (built + tested today), the connectivity
tools, the lattice plumbing, and the saved saturated field as fuel input.
Extend-and-evolve, not start-over.
