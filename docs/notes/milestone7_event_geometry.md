# Milestone 7: per-event mask geometry — the asymptote discriminator

**Status:** Complete, 2026-08-12. All 15 M6 trajectories (L=48/64/96,
5 seeds each, identical seeds and drop counts) replayed locally with
per-event carrier/sink classification. Zero cloud cost. Code:
`src/void_cascade/event_geometry.py`; runner:
`scripts/run_event_geometry.py`; tests: `tests/test_event_geometry.py`
(8 tests, incl. bit-identical trajectory reproduction against the
uninstrumented driver). Data:
`data/outputs/event_geometry_20260808_211945/` (15 npz).

## The question

M6's second-cycle analysis showed the L→∞ coverage asymptote cannot
be decided by comparing mean peak sizes (candidate models differ by
~1.7 points at L=256 vs ~4.6-point seed scatter). The structural
alternative: classify every cell a peak avalanche topples by its
state immediately before the event — carrier (z=1) or sink (z=0).
A stable Manna configuration contains only these two states, so the
partition is exact. A sink can only topple if it receives two grains
mid-event; sink participation therefore measures genuine mid-event
recruitment, not pre-loaded structure.

## Result 1: the carrier-cap hypothesis is falsified

Peak events are NOT confined to the carrier network. In every one of
the 15 runs, 27-31% of the peak event's footprint was sink cells
before the event began. The "asymptote = 62% because only carriers
can participate" picture (hypothesis 1 of the M6 summary) is dead —
structurally, per event, with no statistical ambiguity.

## Result 2: sink recruitment grows with L, significantly

Peak-event geometry (mean ± SEM over 5 seeds):

| L | peak unique% | sink participation σ | carrier coverage |
|---|---|---|---|
| 48 | 36.02 ± 1.48 | 0.2816 ± 0.0035 | 0.419 ± 0.016 |
| 64 | 37.26 ± 0.84 | 0.2909 ± 0.0038 | 0.428 ± 0.010 |
| 96 | 41.42 ± 1.66 | 0.3055 ± 0.0027 | 0.465 ± 0.018 |

Regression over all 15 seeds: σ(L) = 0.148 + 0.0345·ln L, slope
SE 0.0066 (t = 5.2). No deceleration is visible in range (the
per-step slopes are 0.033 and 0.035 per ln-unit). Within runs, σ
also rises with event size at fixed L, and at fixed relative event
size σ still rises with L — participation depends on both.

Carrier coverage at peak grows too (0.42 → 0.46) but its linear-in-
ln L extrapolation reaches 1.0 only at L ~ 3×10^5: peak events do
not sweep the carrier network at any accessible scale. Coverage
growth is joint recruitment of carriers and sinks, not carrier-
network exhaustion.

## The exact reformulation of the asymptote question

If a peak event covered the entire lattice, its sink participation
would equal the global sink fraction: σ = 1 − c = 0.38 (c = 0.62,
constant across L). So:

- full-coverage asymptote  ⇔  σ(L) → 0.38
- capped asymptote U∞ < 1  ⇔  σ(L) → σ∞ < 0.38, giving
  U∞ ≈ c/(1−σ∞) if carrier coverage → 1 on the same trajectory

σ is bounded above by 0.38, and the measured linear-in-ln L trend
hits 0.38 at L ≈ 840 — so the trend MUST bend by there; the open
question is whether it bends below 0.38 (ceiling 89-95%) or
saturates exactly at it (100%). Current-value ceiling arithmetic:
σ∞ = 0.305 → 89%; σ∞ = 0.34 → 94%; σ∞ = 0.38 → 100%. Hypothesis 1
(62%) is excluded; hypotheses 2 and 3 both survive, now separated
by a single scalar with a hard bound.

## Why this observable makes large-L runs affordable after all

σ has ~6x smaller seed scatter than mean coverage (per-seed std
~0.008 vs ~4.6 points). At L=192, "still rising at 0.0345/ln-unit"
predicts σ ≈ 0.329 vs ≈ 0.31 if already saturating — a separation
of ~2.5 per-seed std. Three seeds at L=192 measuring σ give a
~4σ discrimination; five give ~6σ. Contrast mean-coverage FSS,
which needed ~80 seeds at L=256 for the same job. The M6-era
conclusion "brute-force large-L is futile" applies to the mean-
coverage observable, not to σ.

Cost note: L=192 at 13.5M drops is ~4-5 months single-core locally
per seed — not a local job. On cloud (spot, with the monitor bug
fixed first) it is the original ~$25/seed. Decision deferred; no
current commitment. L=128 locally (~12 days/seed, 5 parallel) would
give only ~1-2σ separation — weak; if more data is wanted, L=192
is the right buy.

## Cosmological reading

The inversion event is not limited to the pre-loaded (carrier)
fraction of the substrate: cascades recruit unloaded regions in
proportion growing with scale. At minimum ~89% of the substrate is
implicated in the largest events at large L; full-substrate
inversion (100%) remains open. Either way the M5-era qualitative
claim — the catastrophe regime engulfs the bulk of the substrate,
not a minority network — is restored, with honest numbers.
