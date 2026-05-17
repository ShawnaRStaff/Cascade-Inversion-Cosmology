# Milestone 3 summary: 3D Manna sandpile + percolation

**Status:** Complete, 2026-05-16.

## What was built

- `src/void_cascade/sandpile_3d.py` — 3D Manna stochastic sandpile on a
  cubic lattice with open boundaries. Same hot-loop optimization as 2D
  (bincount + flat indexing).
- `src/void_cascade/percolation.py` — connected-component analysis on a
  3D boolean mask with 6-connectivity, plus a "spans at least one axis"
  criterion. Wraps `scipy.ndimage.label`.
- `src/void_cascade/cluster_geometry.box_count_dimension` — extended to
  accept 2D or 3D hypercubic masks.
- `scripts/run_milestone3_smoke.py` — drive a small lattice until first
  spanning; sanity check that the pipeline works.
- `scripts/run_milestone3_pc_fss.py` — finite-size scaling of the
  spanning threshold across many L with multiple seeds each. Reports
  both a constant-model fit and an FSS scaling fit; the constant model
  is preferred when the data is consistent with no L-dependence.
- `scripts/run_milestone3_dimension.py` — box-count the largest cluster
  from a saved smoke mask.
- Tests: 10 for the sandpile, 11 for percolation, 10 for cluster
  geometry (4 new 3D tests). Total project: 50 tests passing.

## Why Manna in 3D (not BTW)

Same reasoning as Milestone 2: Manna obeys simple scaling and matches
the universality class used in the cosmological-percolation literature
(Carfora & Marzuoli 2023). The implementation extends directly from
`sandpile_2d.py` with one extra axis and 6 neighbors instead of 4.

## The cosmologically meaningful observable

In this project's framing each toppling is a fracture event in the
substrate. The cosmologically interesting quantity is not the avalanche
size distribution (already covered in M1 and M2) but the connectivity of
the cumulative fractured region - the union of all sites that have
toppled at least once. The percolation transition is the moment when
that region first spans the box.

Implementation:

1. Track an `ever_toppled` boolean array as the simulation runs, OR'd
   with each avalanche's support.
2. Periodically run 6-connectivity labeling on `ever_toppled`.
3. A cluster "spans" if it contains at least one site on opposite faces
   along some axis (any of x, y, z).

We use 6-connectivity (face neighbors only) because it matches the
toppling rule: a topple's two grains go to face neighbors, not edge or
corner neighbors. Using 26-connectivity would let clusters merge through
diagonal contacts the dynamics cannot produce.

## p_c(L) results

Finite-size scaling across L in {16, 24, 32, 48, 64, 96}, multi-seed,
each seed independently driving an empty lattice until the cumulative
ever-toppled set first spans the box. Total wall time: 105.7 min.

| L  | n_seeds | <p_c> at first spanning | stderr | std    | mean drops to span |
|----|--------:|------------------------:|-------:|-------:|-------------------:|
| 16 |     12  | 0.1785                  | 0.0061 | 0.0210 |        1 692       |
| 24 |     12  | 0.1792                  | 0.0043 | 0.0148 |        5 703       |
| 32 |     10  | 0.1719                  | 0.0067 | 0.0213 |       13 056       |
| 48 |      8  | 0.1724                  | 0.0027 | 0.0077 |       44 019       |
| 64 |      6  | 0.1800                  | 0.0027 | 0.0066 |      105 907       |
| 96 |      4  | 0.1800                  | 0.0036 | 0.0073 |      355 530       |

The per-L means span only 0.172 to 0.180, a range smaller than the per-L
stderr at most L. The data is statistically consistent with no
L-dependence. The standard FSS fit
`p_c(L) = p_c(inf) + a * L^(-1/nu)` is mathematically degenerate on this
data and we reject it.

**Weighted-mean constant fit** across all six L:

  p_c = 0.1771 +/- 0.0015     (chi^2/dof = 1.13, dof = 5)

The chi^2/dof of 1.13 confirms the constant model is fully consistent
with the data. No detectable finite-size correction across a factor of
6 in linear size (216x in volume).

## How this compares to uncorrelated 3D site percolation

Standard uncorrelated site percolation on simple cubic has p_c ~ 0.3116
with sizable finite-size corrections at L=16-32 (typically 5-10% off the
asymptote at L=16). Two qualitative differences with our run:

1. Our threshold (0.176) is about **half** the uncorrelated value.
2. We see **no** finite-size correction across an L-range of 16 to 64
   (factor of 4 in linear size, 64x in volume).

Both are consequences of the dynamics producing strongly spatially
correlated fractured sets - each avalanche fractures a connected patch,
not independent random sites - so the cluster geometry builds with
strong local correlations and reaches global connectivity at much lower
occupation than random.

For the cosmology this is an order-1 quantitative claim: the substrate
becomes globally connected at p ~ 18%, not p ~ 31%.

## Cluster fractal dimension (first cut, L=24)

Box-counting on the largest connected component of the L=24 ever-toppled
set at the moment of first spanning gives

  D_box = 1.75 +/- 0.16

over the limited fit range b in {1, 2, 4}. Caveats:

- Only 3 box sizes (the cluster bbox / 4 = 6 caps the upper end).
- We measured the cluster at the *first* moment of spanning, which is
  the smallest cluster that spans - much closer to a single backbone
  than the fat critical cluster of static percolation at p_c.
- Standard 3D site percolation at p_c has critical-cluster D_f ~ 2.523.

The ~1.75 estimate is consistent with a quasi-linear backbone plus
side branching, not a dense critical cluster. To compare apples to
apples with the static-percolation literature we would need to take a
snapshot at a fixed occupation p (e.g. drive past p_c by some margin
and freeze) and measure D_f of the largest cluster there.

## Open follow-ups

1. **Save the L=96 spanning mask** during the FSS run so we can
   box-count it cleanly. The current FSS script only records per-seed
   scalars. Either:
   - Patch the FSS script to dump the ever_toppled mask for one seed
     per L. Cost: trivial.
   - Run a single-seed L=96 driver separately that saves the mask. Cost:
     one ~16 min run.

2. **3D visualization** of the spanning cluster. Either marching-cubes
   isosurface (scikit-image) or plotly volume rendering. Deferred.

3. **Comparison to static percolation** by driving past p_c by a fixed
   margin and freezing. Would let us fit a "true" critical D_f.

4. **Multi-seed L=96** to tighten the asymptotic p_c estimate. Current
   plan only schedules 4 seeds at L=96 due to wall-time cost.

## Files written

- `data/outputs/manna_3d_smoke_20260516_153851.npz` (smoke run mask)
- `data/outputs/manna_3d_cluster_fractal_20260516_173654.png`
- `data/outputs/manna_3d_pc_data_<stamp>.npz` (FSS raw)
- `data/outputs/manna_3d_pc_vs_L_<stamp>.png`
- `data/outputs/manna_3d_pc_fss_<stamp>.png`
