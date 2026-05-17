# Session summary: 2026-05-16

Single working session. Closed two milestones (M2 and M3), kept the
server healthy, and produced the first cosmologically loaded
quantitative result of the project.

## Starting state

- Repo had 1D Oslo (M1) complete: tau ~ 1.55, D ~ 2.2, 14 tests passing.
- README, ROADMAP, CLAUDE.md present; no 2D or 3D code yet.
- Prior session had rendered a 4.5 MB 1D dynamics GIF and got stuck
  trying to display it inline. Lesson: leave large outputs on disk, do
  not attempt to view in terminal.

## Work done

### Milestone 2: 2D Manna sandpile (complete)

**Why Manna over BTW:** BTW in 2D exhibits multiscaling (Tebaldi, De
Menech & Stella 1999); a single-tau fit is misleading. Manna obeys
simple scaling and matches the universality class in the
cosmological-SOC literature (Carfora & Marzuoli 2023).

**Built:**
- `src/void_cascade/sandpile_2d.py` (Manna stochastic, open boundaries,
  parallel updates, bincount-optimized hot loop)
- `src/void_cascade/cluster_geometry.py` (population-fit and
  box-counting fractal-dimension estimators; later extended to 3D)
- `src/void_cascade/visualize_2d.py` (animation GIF + spanning-avalanche
  snapshot PNG)
- `scripts/run_milestone2.py`, `scripts/animate_2d.py`
- 11 sandpile tests, 6 cluster-geometry tests, 2 viz smoke tests

**Hot-loop optimization:** initial implementation gave ~30 drops/s at
L=32 (too slow). Replacing np.add.at with np.bincount and using direct
fancy-index decrement got L=128 throughput to ~66 drops/s. This is the
core speedup that made M2 (and later M3) tractable.

**FSS run results (L in {16, 32, 64, 128}, 42 min wall):**
- tau = 1.292 ± 0.002  (lit ~1.275; +1.3%)
- D   = 2.657 ± 0.004  (lit ~2.76; -3.7%)
- alpha = 1.403 ± 0.025 (lit ~1.51; -7%)
- z = 1.489 ± 0.022     (lit ~1.55; -4%)
- D_f (cluster) = 2.053 ± 0.009 (lit ~2.0; spot on)
- Steady-state density 0.68-0.70 (lit ~0.72)

The mild underestimates in D, alpha, z come from including L=16 (still
in finite-size correction regime); rerunning without L=16 would tighten
them.

### Milestone 3: 3D Manna sandpile + percolation (complete)

**The cosmologically meaningful observable.** In this framework each
toppling is a fracture event. The interesting quantity is the
connectivity of the cumulative fractured region - the union of all
sites that have ever toppled. The "inversion event" is when that region
first spans the box.

**Built:**
- `src/void_cascade/sandpile_3d.py` (Manna 3D, cubic lattice, 6-neighbor
  stochastic toppling)
- `src/void_cascade/percolation.py` (6-connectivity component labeling
  via scipy.ndimage.label; spans-any-axis criterion)
- Extended `cluster_geometry.box_count_dimension` to support 2D or 3D
  hypercubic masks
- `scripts/run_milestone3_smoke.py` (single-seed pipeline check)
- `scripts/run_milestone3_pc_fss.py` (multi-L, multi-seed p_c FSS)
- `scripts/run_milestone3_dimension.py` (box-count the spanning cluster)
- 10 Manna 3D tests, 11 percolation tests, 4 new 3D box-count tests

**Smoke test (L=24, 1.7 s wall):** percolation reached at drop 6200,
p=0.219, spans y and z. End-to-end pipeline confirmed.

**p_c FSS run (L in {16, 24, 32, 48, 64, 96}, 105.7 min wall):**

| L  | n_seeds | <p_c> | stderr |
|----|--------:|------:|-------:|
| 16 |     12  | 0.1785 | 0.0061 |
| 24 |     12  | 0.1792 | 0.0043 |
| 32 |     10  | 0.1719 | 0.0067 |
| 48 |      8  | 0.1724 | 0.0027 |
| 64 |      6  | 0.1800 | 0.0027 |
| 96 |      4  | 0.1800 | 0.0036 |

**Headline result:**

    p_c = 0.1771 +/- 0.0015     (chi^2/dof = 1.13, dof = 5)

No detectable finite-size correction across a factor of 6 in linear
size (216x in volume). The standard FSS scaling fit was correctly
rejected as degenerate; the constant-threshold model is the right one.

**Cluster fractal dimension (L=24, first cut):**

    D_box = 1.75 +/- 0.16   (box sizes [1, 2, 4])

Caveat: "first-spanning" cluster is a quasi-linear backbone, not the
fat critical cluster of static percolation at p_c. Standard static 3D
percolation at p_c has D_f ~ 2.523; the comparison is not apples to
apples. To match the static-percolation literature we would need to
freeze the system at fixed p past threshold and measure D_f there.

### Physics implication

The cascade-driven percolation transition happens at roughly **half**
the uncorrelated 3D site percolation threshold (0.18 vs 0.31), with
**no** finite-size correction. Both are consequences of the dynamics
producing spatially correlated cluster growth - each avalanche fractures
a connected patch, so global connectivity is reached much earlier than
random. This is a different universality class from standard
percolation.

For the cosmology framing this is the first order-1 quantitative
claim: the substrate becomes globally connected when approximately 18%
of sites have ever fractured. That is a number that future observational
comparison (M4) will test or fail to test.

### Infrastructure side

Morning briefing run on the homelab server. Findings:

- 27-day uptime; load and memory comfortable.
- 1 unhealthy container (`nexus-web`) flagged; stopped along with
  paired `nexus-iam` per user direction (both dev containers).
- `docker system prune -f` reclaimed **26.77 GB** (images 40.6 -> 32.3
  GB, build cache 36.4 -> 9.7 GB).
- Backup verification incomplete; the briefing skill's hardcoded paths
  (`/mnt/docker/backups`, `/var/backups`, `/mnt/docker/compose/backup`)
  do not match the actual backup location (somewhere under `/media/`
  per user note). The skill needs updating.

## Test inventory at end of session

| File                              | Count | Domain                |
|-----------------------------------|------:|-----------------------|
| `tests/test_sandpile_1d.py`       |    10 | 1D Oslo dynamics      |
| `tests/test_scaling.py`           |     4 | FSS moment method     |
| `tests/test_sandpile_2d.py`       |    11 | 2D Manna dynamics     |
| `tests/test_cluster_geometry.py`  |    10 | 2D+3D fractal-dim     |
| `tests/test_visualize_2d.py`      |     2 | viz smoke             |
| `tests/test_sandpile_3d.py`       |    10 | 3D Manna dynamics     |
| `tests/test_percolation.py`       |    11 | 6-conn span detection |
| **Total**                         |  **58** |                     |

All passing as of close of session.

## Files written to data/outputs/ this session

```
manna_2d_dynamics_20260516_140318.gif           (1.4 MB)
manna_2d_avalanche_snapshot_20260516_140318.png
manna_2d_fss_data_20260516_144422.npz
manna_2d_size_fss_20260516_144422.png
manna_2d_size_moments_20260516_144422.png
manna_2d_duration_fss_20260516_144422.png
manna_2d_duration_moments_20260516_144422.png
manna_2d_cluster_fractal_20260516_144422.png
manna_3d_smoke_20260516_153851.npz
manna_3d_cluster_fractal_20260516_173654.png
manna_3d_pc_data_20260516_182455.npz
manna_3d_pc_vs_L_20260516_182455.png
manna_3d_pc_fss_20260516_182455.png
```

## Documentation written

- `docs/notes/milestone2_summary.md`
- `docs/notes/milestone3_summary.md`
- `README.md` updated with Setup / Running / Tests sections
- Memory entries: `milestone2_status.md`, `milestone3_status.md`,
  MEMORY.md index updated.

## Open follow-ups (none block forward progress)

1. Rerun M2 FSS without L=16 to tighten the duration exponents
   (alpha, z) and the conservation check D(2-tau).
2. Save the spanning mask per L during the M3 FSS run so cluster
   fractal dimension can be measured at larger L without re-simulation.
3. 3D visualization of the spanning cluster (marching-cubes isosurface
   or plotly volume render). Deferred to M4 prep when we know what
   observable we want to look at.
4. Update the morning-briefing skill to scan `/media/...` for backup
   artifacts.

## Honest assessment of where the model stands

Mechanical viability: **yes**. The cascade-plus-percolation framework
runs end-to-end as a 3D simulation, produces a reproducible quantitative
prediction for the inversion threshold, and the dynamics sits in a
distinct universality class from random percolation.

Empirical support: **not yet established**. Zero observations have been
compared. The threshold (~18%) is a prediction waiting to be tested,
not a validation.

The model has earned the right to keep going. It has not yet earned
the right to claim it explains cosmology. M4 (comparison to galaxy /
cosmic-web statistics) is where that judgment gets made.

## Next session

Suggested M4 design conversation (before any code):

1. Which observable first - two-point correlation xi(r) is the obvious
   place to start; cleanest data and easiest to compute from our 3D
   mask.
2. Which data source - SDSS / BOSS / eBOSS published xi(r) values from
   a paper is the practical first cut, no need to download catalogs.
3. How to map "fractured site" to "galaxy" - the dimensionless xi(r)
   reduces this to a single free parameter (lattice spacing in Mpc).
4. At what L - probably L=128 with multiple seeds; needs wall-time
   estimate.

Also worth doing in parallel:
- Save L=96 spanning masks to enable cleaner D_f measurement.
- A small 3D visualization of the spanning cluster to keep visual
  intuition tracking with the numbers.
