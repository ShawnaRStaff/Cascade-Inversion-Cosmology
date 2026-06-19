# Open holes and things to examine closely

**Date:** 2026-06-19
**Status:** Live list. These are weak spots in our own tests and claims.
We write them down honestly so we don't lose them and so we don't fool
ourselves. Plain language on purpose.

A quick note on a word used a lot below:

- **"Holds together" / "connected":** imagine a cracked window. You can
  ask two different questions about it. (1) "How much glass is left?"
  That is just *counting*. (2) "Is the leftover glass still one whole
  sheet that reaches edge to edge, or has it fallen into separate
  shards?" That is *connected-ness*. These are not the same question.
  A window can still have lots of glass left (question 1 says "fine")
  while having already broken into shards (question 2 says "broken").

---

## A. Holes in the substrate-resilience experiment (the main focus)

This is the experiment that asked: "if a big event breaks 91% of the
substrate, does the leftover 9% hold it together?" It answered "no." We
now think that answer cannot be trusted, for these reasons:

1. **It counted, it did not check holding-together.** The whole
   conclusion came from counting how many cells still had something in
   them. It never once checked whether the leftover cells were still
   joined into one piece. The code to check joined-up-ness already
   exists in the project (`percolation.check_spanning`) and this
   experiment simply never used it. So it answered the *counting*
   question and called it the *holding-together* question.

2. **It compared against the wrong number.** It judged the leftovers
   against a "31%" line. That 31% is the break-point for *randomly*
   scattered dots. Our own substrate breaks at about 18%, not 31% —
   and even that 18% is about joined-up-ness, not about counting. So
   the comparison was wrong twice over.

3. **It broke the substrate in the wrong shape.** The experiment
   destroyed cells scattered *randomly* all over. But a real crack does
   not appear as random scattered specks — it spreads through a
   *connected* path, like a crack running across glass. Random damage
   and crack-shaped damage leave totally different leftovers, so a
   conclusion from random damage does not carry over to real cracks.

4. **"Nothing happens on its own" was rigged from the start.** The
   experiment "broke" cells only by *emptying* them. In this model a
   cell only acts when it is *over-full*. Emptying cells can never make
   anything over-full, so "nothing happened on its own" was guaranteed
   before the run even started — it was not a discovery. Worth saying
   plainly: an *implosion* (the substrate falling inward and cells
   getting MORE crammed) is something this model cannot even represent.

5. **"Broken = empty" is too simple.** When the substrate settles, cells
   are a mix — roughly 38% empty, 62% holding one. The experiment
   pretended a broken region is all-empty, which is not what the model
   actually leaves behind.

6. **It disagrees with an earlier result and nobody reconciled them.**
   An earlier experiment (Milestone 4, frozen-region work) checked the
   *holding-together* question the right way and found the leftover
   substrate was STILL one connected piece reaching edge-to-edge even
   when only 20% was left. The resilience experiment (counting, random
   damage) says 9% cannot hold. One says "still whole at 20% left,"
   the other says "broken by the time 9% is left." They point opposite
   ways. This contradiction is the real thing to dig into.

---

## B. Holes in the wider testing program

7. **We never tested the one-event-at-a-time assumption.** Our sim drops
   one grain, lets the whole chain reaction finish, then drops the next.
   The real universe would not wait politely. We do not know if our big
   results are real or just a side effect of this take-turns rule. A
   May-24 run looked like it tested this but actually tested a different
   thing, so the question is still open.

8. **We cannot pin the "does one event eventually cover everything"
   curve from the data we have.** Four data points were not enough — the
   fitted answer ranged anywhere from 45% to nonsense. Only bigger runs
   could settle it, and we should be honest that bigger runs are chasing
   a single all-covering event, which our own results suggest may not be
   how this substrate behaves.

9. **Two ways of measuring the same exponent disagree.** One method and
   another method give different answers, and one of them reports an
   impossible confidence. One is broken; we have not figured out which.

---

## C. Honest limits of the model (NOT things to "fix" yet)

10. **The model has no gravity, no pushing, no heat, no light, and no
    way to form matter.** It only moves grain-counts around. That means
    it cannot show the falling-inward, the heat, or the blasting-outward
    that the working mental picture is built on. This is a limit to
    state honestly in reports — NOT a thing to bolt on. We only add a
    new ingredient when the model itself clearly fails in a way that
    demands it. It has not done that.

---

## D. Mismatches with the published research we are built on

11. **The papers point at a different thing to measure.** The
    self-organized-criticality cosmology papers (Moffat; Carfora &
    Marzuoli) look at the *pattern/spectrum* of structure, not at
    "does a broken region join up edge-to-edge." Our join-up measure is
    our own invention and has the least backing in the literature.

12. **Our headline source paper is not even the same kind of model.**
    Lukash & Mikheeva (the cascade-vacuum paper) describes an orderly,
    step-by-step settling of fields — not the avalanche dynamics of a
    sandpile. They share the word "cascade"; the mechanism is different.

13. **The theory was never written down.** The two documents the project
    treats as required reading (`model_summary.md`, `glossary.md`) do
    not exist in the repo. That missing write-up is probably why "are we
    measuring the right thing" keeps coming back.

---

## E. Record-keeping holes

14. Project memory was empty; there was no write-up for the May-24
    session; the README still says we are at Milestone 2 when we are at
    Milestone 6. (Memory and this file are the start of fixing that.)

---

## Where this points (no decisions made — researcher decides)

The cleanest, cheapest thing to look at first is hole #1 + #6 together:
re-ask the resilience question the *right* way (does the leftover still
join up edge-to-edge?) on the substrate we already have on disk, and try
both random damage and crack-shaped damage. That uses tools we already
have, runs in seconds on this machine, and settles a real contradiction.
It does not add anything to the model.

---

## F. New hole found 2026-06-19: a cell can crack over and over, but the
## model has no idea of "how broken" or "collapsed"

Checked in the code (`sandpile_3d.py`). Two separate facts:

- **Yes, one cell cracks many times.** A cell that gets too full empties
  out, later refills from new input or from its neighbors, and cracks
  again — over and over, forever. Within a single chain-reaction it can
  even crack, get refilled by neighbors mid-reaction, and crack again.
  (The Milestone-6 work measured this: cells crack 2–3x more often than
  the number of distinct cells involved.) So the model definitely does
  NOT limit a cell to one crack.

- **But "fractured" is only ever an on/off stamp.** For the structure
  question, the model keeps a simple yes/no flag per cell: "has this cell
  ever cracked?" Once yes, it stays yes — cracking again changes nothing.
  There is no "slightly broken" vs "badly broken," no counter of how much
  a cell has been damaged, and crucially **no state for a cell that has
  collapsed in on itself.** A cell is either untouched or stamped-cracked,
  full stop.

Why this matters: the thing the researcher wants to study — losing enough
substrate that it shatters and *falls in* — needs a place in the model
for "collapsed." There isn't one. A cracked cell just keeps doing the
same up-empty-refill-crack loop as every other cracked cell. So the model
can tell us *when* and *where* cracking happens and how it joins up, but
it has nowhere to put "and then it caved in."

## Where the sand model ends (the boundary, in plain words)

The researcher's theory runs in stages:
1. The substrate slowly gets denser (soaking up its own "darkness").
2. The densest spot cracks first.
3. Cracking becomes common and spreads.
4. The cracks join up across the whole thing.
5. **Enough substrate is lost that the structure shatters and falls in.**
6. **Islands collapse into themselves; something becomes of what fell in.**

**The sand model faithfully covers stages 1–4.** Density piling up = grains
piling up. Densest-spot-cracks-first = fullest-cell-topples-first.
Spreading cracks = chain reactions. Cracks joining up = the join-up
("percolation") result, which happens at about 18% cracked.

**The sand model goes silent at stages 5–6** — exactly where the theory
gets interesting. It has no weight, no pull, no "falling in," and (per
section F) no "collapsed" state for anything to fall into. It only ever
spreads grains outward; it never pulls them inward.

This is NOT a failure that licenses bolting collapse onto the sandpile.
It means the collapse question lives outside what a grain-spreading model
can represent, and studying it would need a different kind of tool. That
is a researcher decision, not an automatic next step. The disciplined
last thing to squeeze from the *current* model is the cracked-window test
(section A1/A6): it can find the point where the leftover breaks into
islands — i.e. locate the edge of stage 5 — even though it cannot show
what happens past that edge.

---

## Results of the 2026-06-19 local runs (Experiments 1 and 2)

Two small local runs, no AWS, no new physics. New code: `damage.py`
(+ tests), `scripts/run_resilience_connectivity.py`,
`scripts/run_crack_dynamics.py`. Outputs under `data/outputs/`.

**Experiment 1 — the cracked-window test.** Whether the leftover holds
together depends ENTIRELY on the shape of the damage:
- Scattered-random damage: leftover stays one piece up to ~50% removed,
  then shatters into dust by ~70% removed (biggest leftover piece just
  3.6% at 70%, 0.1% at 91%).
- Connected crack-shaped damage: leftover stays ONE connected spanning
  piece even at 91% removed (biggest piece 81%); only breaks at 99%.
- The old test's number (carrier density) read 5.5% at 91% and called
  the substrate "dead." Connectivity says the opposite for realistic
  (connected) damage: the leftover 9% is still one whole piece.

Conclusion: the old "9% can't hold it together" was an artifact of using
*random* damage. For a single realistic connected event, the remnant
DOES hold together. This also reconciles hole #6 — Milestone 4 saw the
leftover span at 20% remaining because its damage was connected; the
resilience run saw collapse because its damage was random. Both correct.

**Experiment 2 — does a crack heal or feed itself?** In the current model
it HEALS. A carved connected crack refilled from 0% back to the normal
62% over 60k drops; the biggest empty region grew at the instant of
carving then shrank back to the background level. With input off, zero
topples. The model spreads and heals; it never grows a crack on its own.
(Side-finding: the settled substrate already contains one giant connected
web of empty cells, ~31% of the lattice, spanning the box — "holes" are a
permanent feature of the rest state, not only made by big events.)

Experiment 2 is also the concrete evidence behind hole #G below: the
healing is real and measured, and it is exactly the behaviour that a real
(irreversible) fracture should NOT have.

---

## G. The big one (2026-06-19): the model self-heals, but a real fracture would not

Researcher's point, and it may be the most important one on this page:
**in real life a fractured piece of substrate does not heal back up.** A
crack stays cracked. But our model heals — a cracked cell empties,
refills from neighbours and new input, and cracks again, forever,
behaving exactly like a cell that never cracked (Experiment 2 measured
this: a carved crack came all the way back to 62%).

Why this is not a small detail:

- **Self-healing is baked into what a sandpile IS.** It only keeps firing
  because cells recharge. Remove the recharge and it is a different
  machine.
- **Milestone 5's headline depends on it.** M5's "the Big Bang never
  ended — it is a permanent regime of catastrophic events" happens only
  because cells heal and re-fire over and over. If a cracked cell instead
  stayed cracked (dead, or a permanent sink), the cascade would most
  likely sweep through once, crack everything, and then GO QUIET — a
  one-time event, not an eternal plateau.
- **Irreversible fracture fits the researcher's own theory better.** The
  working picture is ONE inversion (everything fractures and falls in)
  then ongoing expansion — NOT endless repeats. Self-healing produces
  endless repeats; irreversible fracture would produce a single sweep. So
  the more faithful assumption may also be the one that matches the
  theory.

Discipline note: this is NOT bolting on complexity. It questions whether a
CORE assumption (cracked cells recharge) matches the thing we claim to
model (fracture is irreversible). Testing it means a different — arguably
simpler — base rule: a cracked cell does not recharge (becomes a
permanent sink, or is removed, or gets a changed threshold). Researcher
decision, not automatic.

Honest caveat: I am GUESSING "one sweep then quiet." I have been wrong
before (bet the event-size ceiling was intrinsic; the run proved it
finite-size). The run decides, not the guess. But if it holds, it would
overturn M5's central claim — exactly the kind of thing we chase.

### ANSWER (measured 2026-06-19): M5's permanence was a healing artifact

Ran irreversible fracture at L=48, 400k drops, both modes (sink, hole),
tests first. Code: `src/void_cascade/irreversible.py` (+ tests, 7 pass),
`scripts/run_m5_irreversible.py`. Output:
`data/outputs/m5_irreversible_20260619_135411/`.

Result, decisive:
- Biggest single event in the WHOLE run: ~33 cells = **0.03% of the
  lattice**. M5 (with healing) had events sustained at **50–71%**. About
  2000x smaller.
- By the end, events average ~0.025 cells — effectively **silent**.
- Each cell fractured ~once (107,970 topplings vs volume 110,592); fracture
  reached 97.6% then wound down to quiet. Grain conservation exact.
- Sink and hole identical (grain bookkeeping doesn't change the toppling).

Conclusion: **M5's "permanent, never-ending catastrophe regime" does not
survive once fracture cannot heal.** It was entirely a product of the
self-healing recharge. With irreversible fracture the substrate does a
single quiet sweep and stops.

This both (a) overturns M5's headline and (b) vindicates the researcher's
own correction that the ongoing-ness is NOT from self-healing — the only
thing making the model bang forever was the healing, the very engine the
theory rejects. The real ongoing storm (if any) must come from a mechanism
this model lacks: a front advancing into fresh, infinite substrate
(expansion-as-propagation), not a fixed box that recharges. That is a
different tool (bucket 2), and a researcher decision.
