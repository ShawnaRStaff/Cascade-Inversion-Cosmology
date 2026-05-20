# Session summary: 2026-05-19

Long session. Closed M5 (inversion event / catastrophe plateau) through
Stage 5 (larger-L finite-size test), set up AWS for the next-cycle
finite-size scaling sweep, queued M6 phase-A prep work.

## Starting state

- M4 cycle 1 complete (cosmological comparison: ξ(r) γ=1.80,
  τ_s=2.05 at p_c, γ(z) compatible-with-tuning).
- M5 Stages 1-4 complete from prior session: inversion-approach
  signatures, natural saturation, post-saturation state, drive past
  arrest. L=48, peak single event 71.3% of L^3.
- Open question entering today: is the 71% peak ceiling at L=48
  finite-size, or intrinsic?
- Task list had M5 Stage 4 follow-ups queued (transient comparison,
  z distribution, larger-L test).

## Work done

### M5 Stage 5: larger-L finite-size test (L=64)

Built and ran `scripts/run_milestone5_larger_L.py`: L=64, 2 seeds
(3000, 3001), no auto-arrest, 600k drops each. Same Manna rule as
Stage 4, just a bigger lattice.

**Headline result: the 71% ceiling is finite-size.**
- Seed 3000 peak event: 240,370 cells = **91.7% of L^3** at drop 255k
- Seed 3001 peak event: 211,391 cells = **80.6% of L^3** at drop 185k
- z=0.616 saturation density reproduced bit-identically to L=48
- Both seeds reached p=1.00000 at drop ~445k
- No 100% event in 1.2M total drops across both seeds at L=64

**Reproducibility (accidental but valuable):** Stage 5 inadvertently
ran twice due to a context-switch leaving an orphan process from a
prior session running when I "relaunched." The two completely
independent process invocations produced byte-identical .npz outputs
(SHA256 verified). Bit-for-bit reproducibility given seed, confirmed
across separate processes on the same hardware.

### M5 supporting analyses

- `scripts/run_milestone5_transient_comparison.py` — fits saturation
  event-size distribution and compares to observed cosmic-transient
  luminosity-function slopes (GRBs, FRBs, UHECRs, solar flares,
  magnetar giant flares). Result: Manna τ_s ≈ 1.42, observed range
  1.5-2.0. Close but not exact; suggestive only.
- `scripts/run_milestone5_z_distribution.py` — characterizes per-cell
  z statistics and spatial autocorrelation at saturation. Result: ~38%
  z=0 / ~62% z=1, autocorrelation length ~1.24 lattice units. Spatially
  unstructured saturation state.

### M5 writeup updated

`docs/notes/milestone5_summary.md` extended with Stage 5 section,
updated "Where this leaves the project" to call out the universal
invariants identified (z=0.616 across L=48 and L=64, p_c=0.177 no
L-dep), and added "Next-cycle priorities" section ranking the
data-driven follow-ups.

### Project framing conversations

User and I worked through several framing questions during the
long sim runs:
- The "wave interference" hypothesis: can multiple cascades combine
  for a final blow? Conclusion: not in the current serial-cascade
  model; would need a parallel-cascade extension (queued as #48).
- User's epistemic position: not claiming the model is right, but
  exploring outside the dominant GR/string frameworks. Both points
  about Einstein "wrong vs incomplete" and the sociology of physics
  research engaged honestly.
- Honest assessment of where the project stands: solid Tier-1
  research foundation (working code, real data, documented findings),
  potentially Tier-2 publishable observations (L-scaling, z-invariance,
  transient comparison), not yet Tier-3 (novel cosmological
  contribution — needs dimensional calibration + distinguishing
  prediction + observational test).
- Discipline statement saved to memory: no model extensions absent
  specific data signal demanding them.

### AWS setup for next-cycle FSS sweep

User has $159.96 of AWS credits, expiring 2026-10-18. Set up:
- Region chosen
- IAM user `cosmology-sim` planned (in progress)
- SSH keypair `cosmology-sim-key.pem` downloaded to project root
- Access keys CSV downloaded to project root
- `.gitignore` updated to block *.pem, *_accessKeys.csv, .env, .aws/,
  and generic credential patterns
- File permissions on credentials set to 0600
- macOS metadata files (._*) cleaned up from project root
- Verified credentials are gitignored via `git check-ignore`

## Tasks completed in session

- #44 Update M5 writeup with Stages 1-4
- #45 z distribution at saturation analysis
- #46 Larger-L test for 100% events (Stage 5)
- #47 Compare catastrophe spectrum to cosmic transients

## Tasks created for next-cycle (M6 Phase A)

- #48 Parallel cascade extension (M6 candidate, blocked by Phase A)
- #49 Add checkpoint-resume to sandpile sim
- #50 Build FSS sweep driver script
- #51 Build EC2 launch + sync automation
- #52 M6 Phase A — local prep for FSS sweep (umbrella)

## Where this leaves the project

M5 is closed. M6 is queued. The decision for next session is whether
to start with:
- **Phase A local prep** (#49 checkpoint-resume → #50 sweep driver →
  #51 EC2 automation → smoke-test locally) — about 1 day of dev work,
  zero AWS spend, prepares for Phase B/C
- Or some other priority the user identifies after a night's sleep

The agreed sequence:
- Phase A: local prep (1 day, $0)
- Phase B: AWS shakedown on L ∈ {32, 48, 64} (~3 hr wall, ~$3)
- Phase C: real FSS sweep at L=96 and L=128 (~5 days wall, ~$90)
- Phase D: back to local for analysis, fit scaling law, writeup

Estimated total compute cost: ~$95 of $160 credits.

## Files added or modified in session

Added:
- `scripts/run_milestone5_larger_L.py`
- `scripts/run_milestone5_transient_comparison.py`
- `scripts/run_milestone5_z_distribution.py`
- `docs/notes/milestone5_summary.md`
- `docs/notes/session_summary_2026-05-19.md` (this file)
- Data outputs in `data/outputs/` (two Stage 5 result sets, SHA256
  identical; transient and z-distribution plots)

Modified:
- `.gitignore` (credential protection)
- `scripts/run_milestone5_multistate.py` (per-snapshot progress prints)

## Commits

- `2e4f3da` M5 Stage 5: L=64 finite-size test + writeup + supporting analyses
- `5861568` gitignore: protect AWS credentials and macOS metadata

Both pushed to `origin/main`.

## How to pick up cold tomorrow

1. Read `docs/notes/milestone5_summary.md` (Stages 1-5 of M5, full).
2. Read this file (where we left off + next-cycle plan).
3. Check task list (`TaskList`) — pending tasks #48-#52 are the
   queue for M6 Phase A.
4. Verify AWS credits / IAM / keypair still ready (`aws sts
   get-caller-identity` with the IAM creds will confirm access).
5. Start with task #49 (checkpoint-resume) — it's the foundation for
   everything else and isolated enough to start fresh on.

That should be enough to be productive within 15 minutes of sitting
down.
