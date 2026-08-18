# Milestone 6 summary: FSS scaling sweep with unique-cells metric

**Status:** First cycle complete, 2026-05-23. Second cycle
(full n=5 L=128 sample + L=192 partial) analyzed 2026-08-06 —
see "Second cycle update" section; it revises the acceleration
claim and the asymptote-experiment design.

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

---

## Second cycle update (2026-08-06)

### New data since first cycle

- Three more L=128 seeds completed (22801, 22803, 22804), giving
  the full n=5 sample at every L in {48, 64, 96, 128}.
- One L=192 seed (29200, 13.5M drops) ran May 22 – Aug 4 through
  seven spot/shutdown relaunch cycles. The AWS account exhausted
  its credits and froze on Aug 4 (~21:00 UTC) with the run at
  92.9% (drop 12,547,882). The last checkpoint synced locally
  (drop 12,387,742, 91.8%) was subsequently destroyed by a
  monitor-script malfunction that launched duplicate instances and
  overwrote local + S3 checkpoints with near-scratch state. The
  surviving L=192 quantitative result is the Jul 17 analysis
  figure: peak unique fraction **>= 50.1%** (a floor, not a final
  value). Full run logs to 12.5M drops confirm z = 0.620 held
  throughout.

### Corrected headline table (n=5 everywhere)

| L | n | mean unique% (±SEM) | std | range | mean topples% | retopple |
|---|---|---|---|---|---|---|
| 48 | 5 | 36.02 ± 1.48 | 3.31 | 32.4-39.4% | 74.0% | 2.05x |
| 64 | 5 | 37.26 ± 0.84 | 1.88 | 34.3-39.5% | 83.5% | 2.24x |
| 96 | 5 | 41.42 ± 1.66 | 3.72 | 36.4-46.5% | 100.4% | 2.42x |
| 128 | 5 | 43.04 ± 2.05 | 4.58 | 38.6-50.1% | 115.9% | 2.69x |

The first-cycle L=128 row (45.40%, n=2) happened to include the
highest seed (50.1%); the full sample pulls the mean down to
43.04%.

### The acceleration claim is retracted

Growth per ln(L) unit, with propagated SEM:

- L=48→64:  4.3 ± 5.9
- L=64→96:  10.3 ± 4.6
- L=96→128: 5.6 ± 9.2

The first-cycle value of 13.9 for the last step was an n=2
artifact. With full samples, no slope differs significantly from
its neighbors; a seed-bootstrap gives P(96→128 slope exceeds
64→96 slope) = 0.34 — a coin flip. The defensible statement is:
**unique% grows with L at roughly 4-10 points per ln-unit; the
curvature (accelerating vs decelerating) is not resolved by n=5.**

### All three asymptote models fit indistinguishably

Weighted fits of unique%(L) = A·(1 − a·L^−b) to the four means:

| model | chi²/dof | predicts L192 | predicts L256 |
|---|---|---|---|
| A = 62% (carrier cap) | 0.73/2 | 45.2% | 46.8% |
| A = 100% (full coverage) | 0.59/2 | 45.9% | 47.8% |
| A free | 0.53/1 | 46.3% | 48.5% |

Two consequences, both uncomfortable and both important:

1. **The models differ by ~1 point at L=192 and ~1.7 points at
   L=256, while single-seed scatter is ~4.6 points (std at
   L=128).** Separating the means to that precision at L=256
   would need on the order of (4.6/0.5)² ≈ 80+ seeds. The
   first-cycle claim that "L=256 would settle the asymptote" is
   therefore also retracted: brute-force FSS on the mean cannot
   decide this question at feasible cost.
2. The stranded L=192 seed could never have decided it either —
   its floor of >=50.1% sits ~1σ above every model's prediction,
   consistent with all three (L=128 produced a 50.1% seed too).
   Losing the final 7% of that run cost us essentially nothing
   scientifically. (Censoring note: in 6 of 20 completed runs the
   peak event landed after the 91.8% drop mark, so the true final
   value had a ~30% chance of exceeding the floor.)

### The sharper observable: peak unique / carrier fraction

Carrier fraction is flat at 61.5-61.9% across L, so
mean-peak-unique as a fraction of carriers rises 58.5% → 60.4% →
67.0% → 69.5% (L=48→128). Under the carrier-cap hypothesis this
ratio saturates below 1; under full-coverage it crosses 1 (events
must eventually touch sink cells too). This reframing does not
escape the seed-scatter problem by itself, but it points to the
discriminator that does:

**Per-event mask geometry** (listed under "not done" since the
first cycle). Recording the actual cell set of the peak event
answers directly: does the event touch sinks at all, or is it
confined to the carrier network? Confinement is a structural
yes/no per event, not a noisy mean — a handful of instrumented
runs at L=64-96 (cheap, local or minimal cloud) likely decides
between hypothesis 1 and hypotheses 2/3 without any large-L
brute force.

### Revised position, end of second cycle

- Growth of unique% with L: confirmed, ~4-10 points/ln-unit.
- Acceleration: retracted (n=2 artifact).
- Asymptote: undetermined, and undeterminable by mean-based FSS
  at feasible seed counts. The per-event mask experiment is the
  rational next step (M7 candidate), ahead of any further large-L
  spending.
- z = 0.616-0.620 universality: reconfirmed at every L including
  the L=192 partial.
- The L=192 s29200 trajectory remains resumable in principle from
  a local same-trajectory checkpoint at drop 3.5M (deterministic
  dynamics), but the analysis above removes the scientific case
  for finishing it.
