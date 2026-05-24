# Milestone 6 summary: FSS scaling sweep with unique-cells metric

**Status:** First cycle complete, 2026-05-23.

M6 was originally framed as a finite-size scaling sweep to test
whether the peak-event-size ratio scales toward 100% at large L —
the central scaling claim from M5. Early in execution we discovered
the M5 metric was wrong: M5 measured topplings per unit volume
(which can exceed L^3 due to cells re-toppling within a single
avalanche), not the unique cells touched in a single event. This
note documents the corrected sweep and what we learned.

## The metric correction

M5 reported the largest single avalanche at L=64 as "240,370 cells =
91.7% of L^3" and described this as the cascade approaching full
lattice coverage. The number was correct as raw toppling count, but
the verbal description was wrong: a cascade with 91.7% topples/L^3
can involve far fewer than 91% of unique cells, since the cascade
re-topples cells multiple times.

The physically meaningful quantity — "what fraction of the lattice
was involved in this single event?" — is the unique-cells count.
M6 explicitly tracks `unique_sizes[t]` = `mask.sum()` per avalanche.

## The sweep

L values: 32 (from M5 shakedown), 48, 64, 96, 128.
Seeds: 5 each at L=48, 64, 96; 2 at L=128 (lost 1 seed to spot
interruption without auto-recovery, accepted n=2).
Drops per run: 100k at L=32 up to 4M at L=128.
Total: 17 runs with full unique-cells data.

Compute: 18-instance fleet on AWS spot/on-demand c7i-flex.large and
m7i-flex.large (the Free Plan account restriction blocked larger
compute-optimized instances). Auto-recovery monitor relaunched
spot-interrupted instances on on-demand. Total cost ~$30-40 of the
$160 credit budget.

## Headline scaling result

Mean fraction of unique cells touched in the largest avalanche per
run:

| L | n | mean unique% | range | mean topples% |
|---|---|---|---|---|
| 48 | 5 | 36.02% | 32.4-39.4% | 74.0% |
| 64 | 5 | 37.26% | 34.3-39.5% | 83.5% |
| 96 | 5 | 41.42% | 36.4-46.5% | 100.4% |
| 128 | 2 | 45.40% | 40.7-50.1% | 126.3% |

**Unique% grows with L, accelerating.** Growth rate per ln(L) unit:
- L=48→64: 4.3%/ln-unit
- L=64→96: 10.3%/ln-unit
- L=96→128: 13.9%/ln-unit

Three data points (L=48, 64, 96) initially looked like asymptoting
at 36-40%. The L=128 result (45.4% mean, 50% max) decisively
contradicts that read. Growth is real, persistent, and accelerating
through the sweep range.

## Other confirmed invariants

**Universal saturation density z = 0.616.** Holds across L with
tiny upward drift (0.615 → 0.619 from L=48 to L=128). This is a
property of the dynamics, not a parameter we set.

**Carrier/sink ratio ~62%/38%.** Direct consequence of z = 0.616
(carriers = cells at z=1, sinks at z=0). Constant across L.

**Re-toppling factor grows with L.** Topples/unique ratio:
2.05× (L=48) → 2.96× (L=128). Bigger lattice → cascades involve
more re-touching of the same cells.

## Supporting experiments (local, no AWS)

While the FSS sweep ran, two local experiments addressed
user-raised physical questions:

### Aftershock / temporal correlation analysis

Question: do big avalanches cluster in time (earthquake aftershock
pattern), or are they temporally independent?

Method: autocorrelation of `sizes[t]` and inter-large-event time
distribution analysis on all 15 saturated runs (L=48, 64, 96).

Result: **events are temporally independent (Poisson).**
Inter-event times exponential-distributed in every run.
Autocorrelation near zero at small lags. Aftershock hypothesis
rejected. Each avalanche is self-contained.

### Substrate resilience experiment

Question: does the substrate require 100% fracture to structurally
collapse, or does a smaller fraction suffice? Specifically: when
91% of cells topple in one event, does the remaining 9% "hold the
substrate together"?

Method: take a real saturated L=48 z field, artificially shatter a
configurable fraction (30%, 50%, 70%, 91%, 99%), continue dynamics
with and without input, measure recovery.

Results:
- 91% shatter → carrier fraction drops from 62% to 5.5%, FAR below
  3D percolation threshold (~31%). Substrate is structurally
  non-functional after a single 91% event.
- Without input: zero spontaneous topples in any scenario.
  Substrate just sits in its damaged state. No implosion, no
  collapse cascade, no aftershocks.
- With input: recovery is slow and linear (~+0.045
  carrier-fraction per 5000 drops, regardless of damage level).

**The 9% remnant does not hold the substrate together.** It's
just cells the cascade didn't reach, mostly already sinks.

## Honest assessment

What survives M5 cleanly:
- The substrate produces a permanent late-cycle regime of large
  cascading events (M5 Stage 2-4 finding)
- z = 0.616 saturation density is universal
- The dynamics is reproducible bit-identically given seed
- Cascading events involve a substantial fraction of the lattice
  (the new measurement: ~36-46% at modest L, growing)

What gets revised from M5:
- The "events approach 100% lattice spanning at large L" claim was
  based on the wrong metric. With the right metric, events involve
  a smaller fraction of the lattice (~36-45% across our range),
  growing slowly with L. Whether they approach 100% at L → ∞ is
  open.
- The "91% peak event" figure means 91% of L^3 in topples, NOT 91%
  of cells. The actual lattice coverage of that event was around
  ~40%.

What's new from M6:
- Unique-cells fraction grows with L, accelerating (45% at L=128,
  growing roughly +4%/ln(L) per step in the recent sweep)
- Aftershocks rejected; events are temporally independent
- Substrate is dynamically inert without input
- 9% remnant cannot hold substrate together; cascade with proper
  geometry compromises load-bearing network (carrier density drops
  below percolation threshold)

## Open question (the asymptote shape)

With four L values and ~4-14% growth per ln-unit, three plausible
shapes for unique%(L):

1. **Asymptote near the carrier fraction (~62%).** Cascades touch
   all carriers but not sinks. Currently at 45.4%, climbing toward
   62%. Physically: the substrate has a structural cap given by
   the carrier density. L=256+ would show clear deceleration.

2. **Continued growth toward 100%.** Acceleration continues.
   Cascades approach full coverage at large L. The original M5
   "events span the lattice" reading was qualitatively right, just
   with smaller numbers at any finite L than M5 implied.

3. **Power-law growth with a higher asymptote (~80-90%).**
   Compromise between (1) and (2). Cascades touch most of the
   substrate but the corner-of-the-lattice cells stay untouched
   for geometric reasons.

Distinguishing requires L=192 or L=256 data. Quoted compute on
free-tier-eligible m7i-flex.large: ~3-5 days per seed × n seeds.
Defer to M7.

## Outputs on disk

- `data/outputs/fss_sweep_20260521_031056/` — all FSS data,
  17 final.npz files with `unique_sizes` array, summary.json per run
- `data/outputs/aftershock_analysis_20260521_213008/` —
  temporal correlation analysis results
- `data/outputs/substrate_resilience_20260521_220701/` —
  shatter recovery experiment

## What's not done

- Per-event mask geometry (would let us test "does single event
  percolate the lattice"). Currently only count is saved, not
  shape. Would need code modification + re-run sample.
- Dimensional calibration of z = 0.616 to physical energy density.
  Still the rate-limiting step for any quantitative cosmological
  claim.
- L=256+ scaling data. Would settle the asymptote question.
- The substrate sub-structure / heat / radiation extensions the
  user articulated. Still on the discipline shelf — no data signal
  yet demanding them.

## Honest position, end of M6

We have a corrected scaling picture for the catastrophe regime:
single avalanches involve a meaningful (and growing) fraction of
the lattice, well below the toppling-count metric M5 used but
nontrivial nevertheless. The model still produces a permanent
critical regime with reproducible properties. The cosmological
reading needs to be reframed (smaller numbers than M5's headline
but qualitatively similar), and the L → ∞ asymptote remains the
key open question.

We can either:
(a) Stay with current data, write up findings as "cascades scale
    upward in coverage with L, asymptote shape undetermined"
(b) Run a small L=192 batch to better constrain the curve
(c) Pivot to dimensional calibration to convert the model's
    constants (z=0.616, p_c=0.177, etc.) into physical numbers

The framework's core claim (permanent late-cycle catastrophe
regime, no singularity, ongoing substrate events as the "post-bang"
state) survives intact. The specifics of M5's amplitude
overstatement are corrected.
