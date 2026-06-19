# Design sketch: heat-gated substrate (cold-frozen → heat runaway)

**Date:** 2026-06-19
**Status:** BUILT — see `heat_model_findings.md` for results. The binary
version here became `heat_gated.py`; the damage-coupled (scattered →
cascading) version that actually reproduces the narrative is `cascade_heat.py`.
This note is kept as the original design contract.

Plain language. This addresses a **different beat** of the storyline than
the front sketch. The front was about the *expansion* (a wave through
fresh fuel) — that comes *after* the catastrophe. THIS model is about the
beat *before* it: the long cold buildup, and whether it quietly stays
frozen forever or eventually **tips into a runaway**.

It reconciles the two results we thought were in conflict:
- self-healing sandpile (M5) → never stops  = the **cold** regime
- irreversible fracture → one quiet pass     = the **hot** regime
They are two ends of a **temperature axis** we never had. The switch is
heat. Absolute zero is its floor; a melt/runaway point is its ceiling.

The one question it answers:

> Does the slow, cold buildup sit frozen and intact-looking forever — or
> does the heat from all that breaking eventually outrun how fast heat
> escapes, cross the melt point, and run away into catastrophe?

Built so it CAN just stay cold forever. If it always runs away no matter
what, it's rigged (the M5 lesson again).

---

## The picture in one paragraph

Every cell starts at **absolute zero** — no motion, frozen. It slowly soaks
up "darkness" and gets **denser**. When a cell gets dense enough it
**cracks** — but here's the key move: while it's cold, a crack is **frozen
in place.** The break is real, but with no motion nothing separates, so the
cell still *looks* intact. Cracks pile up invisibly. Each break gives off a
little **heat**, and heat **spreads and leaks away.** As long as heat leaks
out as fast as breaking makes it, the substrate stays cold and frozen — for
eons. But if breaking ever makes heat faster than it can leak, a spot
crosses the **melt point**: motion returns, and all the **frozen-up cracks
there suddenly let go at once** — it shatters, dumps a burst of heat, and
that heat can push neighbours over the melt point too. A **thermal
runaway** — the storm.

---

## Data contract (data is just data — flat, typed)

Per cell (start 1D, then 2D/3D):

- `density`  : number >= 0. Accumulated "darkness." Rises slowly over time.
- `cracks`   : integer >= 0. How many times this cell has fractured —
               **accumulates** (your "fracture, repeatedly, not shatter").
               This is the hidden damage that's frozen in place while cold.
- `heat`     : number >= 0, **floored at 0 (absolute zero).** Local
               temperature. Starts at 0 everywhere.

Parameters (declared, principled, never tuned to force a runaway):

- `fracture_density` : how dense before a cell cracks.
- `melt_heat`        : the temperature above which cold no longer holds —
                       frozen cracks "express" (separate / shatter).
- `heat_per_crack`   : heat released by one break.
- `dissipation`      : how fast heat spreads to neighbours and leaks off the
                       edge. **This is the guardrail knob** — enough of it
                       and the substrate can stay cold forever.

## The rules (plain)

1. **Soak up darkness.** Each step, `density` rises slowly (the eons).
2. **Fracture (cold-safe).** When `density >= fracture_density`, the cell
   cracks: `cracks += 1`, density partly relieved, a little `heat` added.
   While `heat < melt_heat`, that crack is **frozen in place** — it counts
   as hidden damage but nothing separates and nothing else happens.
3. **Heat spreads and leaks.** Each step, heat diffuses to neighbours and
   some is lost at the boundary (`dissipation`). Heat can never go below 0.
4. **Melt / release (the runaway).** If a cell's `heat >= melt_heat`, its
   accumulated `cracks` **let go at once**: it shatters and dumps a burst of
   heat (bigger if more cracks were stored). That heat spreads — and can put
   neighbours over `melt_heat` too.

No bonds. No outside pull. The "healing" is **not** re-gluing — it's simply
that **below the melt point, breaks can't express as separation** (no
motion). That's bond-free and matches the no-bonds finding.

## The self-feeding engine (NOT the recharge we threw out)

> break → heat → crosses melt point → frozen cracks release → more breaking
> → more heat → ...

A **thermal runaway**. It feeds itself through heat, exactly as the
researcher wanted ("feeds itself, just not through self-healing").

## The "must be able to stay cold" guardrail

- **Heat must dissipate** (spread + leak). The genuine, un-rigged question
  is a *race*: heat made by breaking vs heat lost. If loss wins, it stays
  cold forever; if making wins, it runs away. We do NOT pre-decide which.
- **Thresholds are declared and swept, not tuned to win.** As with the
  front's critical threshold, there should be a **critical line** between
  "stays cold" and "runs away"; finding it is the result.
- **Sanity test that it CAN stay cold:** with strong dissipation it must
  never run away. If it always runs away, the code is rigged and we'd know.
- **Absolute-zero floor:** heat >= 0; start everything at 0. The dial has a
  real bottom — it's a bounded variable, not a free knob.

## What we measure (so it can fail)

1. **Stay frozen, or run away?** (peak heat over time; fraction "released"
   vs still-frozen.)
2. **The long-quiet-then-sudden shape.** Does damage pile up invisibly for
   a long flat stretch and then tip suddenly — your narrative — or not?
3. **Hidden damage at the tip:** how many cracks were frozen in place when
   it finally let go (the "looked intact the whole time" claim).
4. **The critical line** between stay-cold and runaway (in terms of the
   make-vs-leak balance), and where principled values land on it.

## Honest limits (say them up front)

- This models the **buildup → tipping point** ONLY. It does **not** do the
  collapse/inversion (inward motion) or the expansion (the front). Those
  stay separate, later, bigger.
- Heat here is a **simple spreading number**, not real thermodynamics or
  radiation. A deliberate simplification.
- Still a **fixed lattice** (no growing domain). The runaway might be what
  *leads to* collapse, but we are not modelling the collapse itself.
- **Earned, but chosen:** the model's two regimes told us *a switch* is
  missing; WE chose temperature as that switch (good physical reason, your
  story, the cold/heat physics). Owning that is what keeps it honest.

## Build order (cheap first, like M1 → M3)

1. **1D**: a line that slowly densifies and cracks — does it stay cold, or
   run away? Find the critical make-vs-leak line. Fast, clearest.
2. **2D/3D**: does the runaway spread as a front? does the long-quiet-then-
   sudden shape hold? where do principled thresholds land?
3. Tests first each step: heat never below 0; with strong dissipation it
   stays cold (never runs away); cracks accumulate while frozen; a runaway,
   once started, releases the stored cracks; the critical line exists.

## Not a rewrite

Reuse the irreversible-fracture state machinery (cracked/fractured tracking),
the lattice plumbing, and the slow-drive idea. Add one new field — `heat`,
with its absolute-zero floor — and the melt-point rule. Extend, don't
restart.

## How this sits next to the front sketch

Two beats of the same arc:
- **This (heat-gated):** the cold buildup tipping into the catastrophe.
- **Front sketch:** the expansion that feeds on fresh substrate afterwards.
Between them sits the collapse/inversion, which needs *motion* and is the
one neither sketch covers yet. The heat model is the more central first
test — it's the one that asks whether the quiet cold even *tips* at all.
