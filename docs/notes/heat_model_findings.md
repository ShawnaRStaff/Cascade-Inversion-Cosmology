# Heat-model findings: cold buildup → sudden catastrophe (2026-06-19)

**Status:** First results. Buildup → tipping point only (no collapse, no
expansion yet). Plain language.

## Why we built this

M5's "permanent, never-ending catastrophe" was overturned — it ran on
self-healing (cells recharge and re-fire), which is the wrong engine for
fracture (see `open_holes_scrutiny.md` hole #G). The researcher's picture
needs the ongoing-ness to come from somewhere else. Working through it, the
candidate became: **heat.** A cold substrate accumulates hidden damage; if
the heat that breaking makes ever outruns how fast heat escapes, it tips
into a runaway. The temperature axis (absolute zero floor → melt point)
reconciles the two earlier results: self-healing = the cold regime,
irreversible = the hot regime.

## What we built (tests-first, all passing)

- `src/void_cascade/heat_gated.py` — first cut. Cells crack independently
  from a slow drive (scattered), each crack adds a fixed dab of heat
  ("one and done"); heat diffuses and cools; above a melt point, frozen
  cracks release. (+ `tests/test_heat_gated.py`)
- `src/void_cascade/cascade_heat.py` — the fix. Fracturing is **Manna-
  coupled**: a fracture sheds 2 grains to random neighbours, so breaks
  **cascade in avalanches** that grow from scattered (sparse early lattice)
  to clustered (filled, critical). Heat is made along the avalanche, so it
  arrives **concentrated**. (+ `tests/test_cascade_heat.py`)

A real bug was caught on the way (thanks to a "are you sure it's slow?"
challenge): the first coupled version set the shed amount equal to the
threshold, so a single grain became a lattice-spanning random walk that
fractured everything (~L² per step) — not a bounded avalanche. Fixed by
using the validated Manna rule (threshold 2, shed 2). Avalanches are now
bounded and self-organising, exactly as in M1–M3.

## The results

**1. Scattered ("one and done") heat barely tips.** `run_heat_gated_1d.py`:
the runaway only happens at essentially **zero cooling** (critical cooling
~0.003). Any real heat loss to the surrounding cold keeps it frozen. Heat
dribbled thinly across space and time is easily mopped up.

**2. Clustered (avalanche) heat tips robustly.** `run_cascade_heat_1d.py`,
heat-per-fracture 0.15: tips at **every cooling tested, 0.0 → 0.7** — a
>200× jump over the scattered model. Concentration is the whole game. But
at this strength it tips almost instantly (right when cascading begins), so
no quiet buildup.

**3. Heat-making strength is the real dial — and it gives the narrative
shape.** `run_cascade_heat_critical.py`, cooling fixed at a realistic 0.1,
sweeping heat-per-fracture:

| heat/fracture | tips at step (of 4000) |
|---|---|
| 0.000 | never (stays cold) |
| **0.005** | **3331 (83% in)** — peak heat **~12,500** |
| 0.010 | 103 |
| 0.020–0.030 | 54 |
| 0.050–0.080 | 42 |
| 0.120–0.150 | 22 |

Right at the critical edge (heat ≈ 0.005), the substrate stays **quiet for
83% of the run**, looking intact while hidden damage piles up, then
**suddenly lets go** in a huge burst. That is the cold-eons-then-sudden-
catastrophe shape, **under realistic cooling** — exactly the researcher's
story. The closer to the edge, the longer the quiet; the longer the quiet,
the more violent the release (the cell hoards tens of thousands of hidden
cracks, then dumps them at once).

## Honest caveats

- **Critical heat-making is very low (~0.005 at cooling 0.1).** So under
  realistic cooling, *almost any* nonzero heat-making eventually tips it —
  the "stays cold forever" regime is a knife-edge at heat ≈ 0. This fits
  the theory (the catastrophe is *eventually inevitable* after enough
  buildup), but stated plainly: tipping is the **generic** outcome, not a
  fine-tuned one.
- **Single seed** for the headline 3331 tip; confirm across seeds before
  leaning on the exact number.
- **Scope:** this is the **buildup → catastrophe tip ONLY.** No collapse/
  inversion (needs motion), no expansion (the front). Those remain ahead.
- **Earned but chosen:** the two earlier regimes told us *a switch* was
  missing; we chose temperature. Guardrail held throughout — with no heat
  source it stays cold (verified), so nothing is rigged to tip.

## Where it leaves the arc

We have now, with tested models: self-healing shown to be a fake engine
(M5 overturned); scattered heat shown too weak; **avalanche-clustered heat
shown to reproduce the cold-buildup → sudden-catastrophe under realistic
cooling.** The next beats — the collapse/inversion (motion) and the
expansion (front, sketched in `front_model_design_layer0.md`) — are still
unbuilt.

## Update 2026-06-19: 2D, on validated SOC — and the catastrophe is a FRONT

The 1D results above are on degenerate 1D Manna (SOC check gave tau~0.78), so
the 1D heat numbers are qualitatively right but quantitatively exaggerated.
Rebuilt in 2D (`cascade_heat_2d.py`, `run_cascade_heat_2d.py`) on the
M2-validated Manna rule:

- **Grounded:** 2D avalanche exponent **tau = 1.275** — dead-on the validated
  2D Manna value. The foundation is genuine SOC; 2D heat numbers are
  trustworthy.
- **Long quiet still holds:** at realistic cooling (0.1), it stayed quiet to
  step **1836 of 3000**, then tipped.
- **The catastrophe SPREADS AS A FRONT.** The released region grows from a
  small seed outward — 14 → 102 → 173 → 324 → 658 → 1342 cells over 16 steps
  (snapshots in `front_2d.png`: compact blobs that merge and sweep across to
  86%). It does NOT pop everywhere at once.

Why this matters: the tipping is **not simultaneous — it propagates
spatially.** That is "expansion as a propagating wavefront of activity"
(the README's own phrase), and it **emerged from the heat dynamics**, not
built in. It links the buildup→catastrophe beat to the expansion beat.

Caveats: single seed; front *speed* not characterized; it's the **activity
boundary** spreading, not material motion — so the collapse/inversion
(material falling in, needs momentum) is still unaddressed.

## CORRECTION 2026-06-19 (later): the front is a FINITE-SIZE ARTIFACT

Robustness + a bigger lattice (`run_cascade_heat_2d_robust.py`,
`run_cascade_heat_2d_frontspeed.py`) overturn the "spreads as a front" claim:

- **L=40:** every seed tips (~step 1100, tight) and the release fills ~86%.
- **L=100:** one seed **never tipped** in 6000 steps; the other tipped late
  (step 5176) and the release stayed a **tiny localized blob (edge ~3.7
  cells, speed ~0, R²≈0.04)** — it did NOT propagate.

What happened: in a *small* lattice heat diffuses across the whole thing, so
it sits near-uniformly hot — when one cell tips, the rest are already near
melt, so the release goes off almost everywhere at once and *looks* like a
fast front filling the box. At L=100 heat **dilutes into the surrounding
cold**, so a hot spot barely forms and any tip is **quenched** before it can
spread.

So: the **expansion-as-front reading is NOT supported at realistic size**, and
even the easy *tipping* was partly a small-L effect (at L=100 it barely
fires). Same pattern as M5's finite-size ceiling — dramatic small, dissolves
at scale. The buildup→tip mechanism is real but its size-robustness is now in
serious doubt; the front claim is retracted.

## RESOLUTION 2026-06-19 (edgeless test): the front IS real — the "artifact" was a confound

The "finite-size artifact" conclusion above was itself premature (caught by
Shawna): the L=100 run drove only 1 grain/step regardless of size, so the big
lattice was **under-driven ~6x per cell** — AND it had hard edges that reflect
energy. Both confound the answer.

Removing both — `cascade_heat_grow.py` / `run_cascade_heat_grow.py`: an
**edgeless, growing substrate** (fresh pre-loaded substrate padded ahead of
the front, so no edge, no reflection, no wrap), loaded fairly and ignited at
one spot — gives a clean answer:

- The front **sustains and advances at CONSTANT speed** (~0.95 cells/step),
  **linear** front-radius vs time, **identical across 3 seeds**, growing the
  domain 6x and still going when we capped it (never touched an edge).

So the catastrophe **does** propagate as a sustained front once the
edge/under-driving confounds are removed. The earlier "fizzle at scale" was
those confounds, not the physics. Constant speed = linear distance-time = the
**prerequisite for a Hubble-like relation** — *within flat-local lattice units*
(geometry still open). We ignited it (didn't wait for spontaneous tip), and
it's the *activity* wavefront, not material motion.

Net: front claim **un-retracted, with the right caveats**. The honest chain was
looked-real(small) -> artifact-worry -> caught-confounds -> real-when-done-right.
