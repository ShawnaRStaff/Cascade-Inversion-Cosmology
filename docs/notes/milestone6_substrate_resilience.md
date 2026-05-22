# M6 substrate resilience experiment

**Date:** 2026-05-21
**Status:** Complete

## Question

User intuition: in real materials/structures, 9% remnant doesn't hold
anything together. A human at 9% organ function dies. A building with
91% of load-bearing supports out collapses. Glass at 91% damage fails.

Does our substrate model show the same? Or does the 9% somehow hold
the lattice together?

## Method

Starting state: L=48 saturated z field from an existing FSS run
(seed=14800, z_avg=0.6159).

Artificially "shatter" a configurable fraction of cells (set to z=0,
simulating a single big avalanche having just occurred). Test six
fractions: 0%, 30%, 50%, 70%, 91%, 99%.

For each fraction, two conditions:
- **Drops ON**: continue normal dynamics (input continues, 5000 drops)
- **Drops OFF**: no new grains, just relax once and observe

Measure: carrier fraction (cells at z≥1), z_avg, spontaneous topples.

## Results

| Shatter % | z_avg post | Carrier % | Drops-on recovery (5k drops) | Drops-off spontaneous topples |
|---|---|---|---|---|
| 0%  | 0.6159 | 61.6% | 0.6154 | 0 |
| 30% | 0.4324 | 43.2% | 0.4750 | 0 |
| 50% | 0.3078 | 30.8% | 0.3518 | 0 |
| 70% | 0.1857 | 18.6% | 0.2303 | 0 |
| 91% | 0.0551 | **5.5%** | 0.1002 | 0 |
| 99% | 0.0064 | 0.6% | 0.0516 | 0 |

Plot: `data/outputs/substrate_resilience_20260521_220701/resilience.png`
Raw: `data/outputs/substrate_resilience_20260521_220701/results.json`

## Findings

### 1. Without input, the substrate is fully inert

Zero spontaneous topples in every drops-off condition. The substrate
has no internal energy to release on its own — after relaxation
completes, no cell is above threshold. Stop input, and the substrate
just sits in whatever state it was left in. **No implosion. No
collapse cascade. No aftershocks.**

### 2. A 91% event reduces carrier fraction from 62% to 5.5%

Far below the 3D percolation threshold (~31%). The substrate's
cascading capacity is **structurally destroyed** by a single event
of this magnitude. The remaining 9% intact cells (5.5% carriers +
3.5% sinks) cannot propagate cascades.

User intuition confirmed: **9% remnant does not hold the substrate
together.**

### 3. Recovery with input is slow and linear

About +0.045 carrier-fraction per 5000 drops, regardless of starting
shatter level. Full recovery from 91% shatter back to z=0.616 takes
~60,000 drops. Recovery is driven entirely by random drops gradually
re-priming sinks to carriers.

### 4. The substrate is a stress reservoir, not a load-bearing structure

In this model, the cells don't push on each other gravitationally
or mechanically. There is no "weight" being held up. Events are
triggered exclusively by external input (drops). Without input,
no events. Damage state persists indefinitely.

## What this confirms

- 9% remnant does NOT hold the substrate together (real-physics
  intuition is correct)
- 91% fracture leaves substrate carrier density well below
  percolation threshold — structurally non-functional
- Recovery requires continuous input; without it, the substrate
  stays where it was left

## What our model cannot capture (limits)

- **Implosion / falling inward**: requires gravity or spatial
  pressure, neither in this model. The substrate doesn't
  mechanically collapse — it just sits at low z.
- **Violent explosion**: the energy released when 91% fractures is
  tracked as grain count, not as heat/radiation/blast. Real physics
  would have that energy propagating outward as radiation.
- **Re-aggregation**: the model doesn't include any matter-formation
  rule. The released energy doesn't condense into anything; it just
  redistributes to neighboring cells or off the boundary.

## Implication for the cosmological reading

If the substrate had a finite energy budget (no continuous input
analog), then after a 91% event the substrate would sit at near-zero
carrier density forever. Universe goes silent after the catastrophe.

But the real universe is full of ongoing high-energy events.
Implications:
1. Continuous input must exist in the cosmological substrate
   (analogous to the user's "substrate sub-structure feeds energy
   upward" hypothesis), OR
2. The 91% catastrophe is itself the post-bang state we observe,
   and the ongoing events we measure are the slow recovery toward
   equilibrium, OR
3. The model is missing the radiation/heat propagation that would
   make the energy "go somewhere" after release rather than just
   sit as low-z cells

These are testable in extensions but not in the current model.

## Bottom line

The experiment directly answered the user's question: in our model,
the substrate does NOT have any "9% backbone" that holds it
together when 91% fractures. The carrier network is broken; the
substrate is structurally non-functional. Only continuous input
restores it.

This sharpens the cosmological reading: whatever is providing
ongoing energy input to the cosmological substrate is the thing
holding the universe in its active, eventful state. Without that
input, the substrate (and the universe) would freeze.
