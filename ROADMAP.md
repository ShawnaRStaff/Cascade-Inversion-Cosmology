# Roadmap
 
## Implementation plan for the void cascade simulation
 
This document tracks the concrete simulation milestones. The goal is to convert the qualitative model into a working 3D sandpile that produces structures comparable to the observed cosmic web. Each milestone has deliverables, key questions, and success criteria.
 
---
 
## Milestone 1: 1D sandpile prototype
 
### Purpose
Verify that the basic SOC dynamics work and that the code produces the expected power-law avalanche distribution. This is preparation, not the cosmological model yet.
 
### Deliverables
- Python module implementing a 1D sandpile with configurable parameters
- Tracking of avalanche events: size, duration, time of occurrence
- Plot showing avalanche size distribution on log-log axes
- Verification that the slope matches the known 1D SOC exponent
### Key questions
- Does the system actually self-organize to a critical state without fine-tuning?
- How long does it take to reach criticality?
- What's the right way to drive the system (uniform addition, random addition, edge addition)?
### Success criteria
- Avalanche size distribution shows a clear power law over at least two decades
- The exponent is consistent with published 1D SOC values
- The code is clean enough to extend to 2D without rewriting
### Estimated effort
Few days. This is mostly to build intuition and infrastructure.
 
---
 
## Milestone 2: 2D sandpile with visualization
 
### Purpose
Move to a dimension where the avalanche geometry is visible and interesting. Verify that fractal avalanche clusters emerge.
 
### Deliverables
- Python module implementing a 2D sandpile on a square lattice
- Animation showing avalanche propagation in real time
- Static visualizations of representative avalanche clusters
- Measurement of the fractal dimension of avalanche clusters
- Comparison of the avalanche size distribution to the 1D case
### Key questions
- Do the avalanche clusters look fractal?
- How does the dimension of the lattice affect the universal exponents?
- What does the steady-state look like visually?
### Success criteria
- Visual output shows clear avalanche dynamics with clusters of all sizes
- Fractal dimension of clusters is measured and reported
- The avalanche size distribution exponent is measured and matches the published 2D value
### Estimated effort
One to two weeks. The visualization is the main work.
 
---
 
## Milestone 3: 3D sandpile with cosmological reinterpretation
 
### Purpose
This is the cosmological version. Reinterpret cells as substrate units, "grains" as accumulated darkness, topplings as fracture events. Track the percolation behavior of the fractured cluster.
 
### Deliverables
- Python module implementing a 3D sandpile on a cubic lattice (or alternative geometry)
- Tracking of the largest connected fractured cluster over time
- Detection of the percolation transition
- 3D visualization of the fractured region (using plotly or similar)
- Measurement of the percolation threshold and the cluster's fractal dimension
### Key questions
- At what fraction of fractured cells does percolation occur?
- What does the connected fractured region look like at percolation?
- How does the cluster's fractal dimension compare to the observed cosmic web?
- What's the lifetime of the system from initialization to percolation?
### Success criteria
- The percolation transition is clearly identified
- The 3D fractal dimension of the fractured cluster is measured
- The cosmic web fractal dimension (about 2.0 to 2.5 from observations) is in the same range as the simulation produces, or the difference is quantified
### Estimated effort
Two to four weeks, depending on how much visualization work is done.
 
---
 
## Milestone 4: Comparison to observations
 
### Purpose
Test whether the simulated cascade reproduces observed cosmic structure statistics. This is the first real empirical test.
 
### Deliverables
- Power spectrum of the simulated fractured region
- Cluster mass function
- Two-point correlation function
- Comparison plots against observational data from SDSS and similar surveys
### Key questions
- Does the simulated structure match the cosmic web on any scales?
- Where does it match and where does it diverge?
- Are the divergences telling us something about the substrate geometry?
### Success criteria
- Quantitative comparison is produced for at least one statistical measure
- Any matching or divergence is documented and explained
### Estimated effort
Four to eight weeks. This is the hardest milestone because it requires understanding observational data well enough to compare against.
 
---
 
## Milestone 5: Address the inversion event
 
### Purpose
Currently the simulation tracks the slow build phase but not the catastrophic inversion. Add a representation of the percolation transition as a global event that "inverts" the substrate.
 
### Deliverables
- Implementation of an inversion rule when the connected fractured region exceeds threshold
- Tracking of what happens to the lattice state at and after inversion
- Comparison of pre-inversion and post-inversion structure
### Key questions
- What's the right rule for what "inversion" does to the lattice?
- Does the post-inversion structure look like the observed universe?
- Can the inversion produce the CMB-like uniformity?
### Success criteria
- The inversion is well-defined operationally
- Post-inversion statistics are computed and reported
### Estimated effort
Open-ended. This is where the model is least specified.
 
---
 
## Milestone 6: Address the redshift problem
 
### Purpose
Derive a quantitative version of expansion-as-propagation that reproduces Hubble's law without invoking metric expansion.
 
### Deliverables
- A geometric or kinematic derivation of an effective Hubble-like relation
- A prediction for how the redshift-distance relation would differ from the standard cosmology, if at all
- A test plan for distinguishing the two interpretations observationally
### Key questions
- Can a propagating wavefront in a static substrate reproduce H_0 ≈ 70 km/s/Mpc?
- Does this interpretation predict any observable deviations from the standard FLRW redshift-distance relation?
- How does the model handle the Type Ia supernova observations that established the accelerating expansion?
### Success criteria
- A quantitative model exists that produces Hubble-like behavior
- The model is either consistent with observations or has an identified failure mode
### Estimated effort
Open-ended. This is the make-or-break empirical question.
 
---
 
## Working principles
 
This is a personal research project. The following principles apply:
 
1. **Working code beats elegant theory.** A simulation that runs and produces output is more useful than an equation that hasn't been computed.
2. **Be honest about uncertainty.** When something is hand-wavy, label it as such in comments and notes.
3. **Compare to observations as early as possible.** It's easy to fall in love with a simulation that doesn't match reality. Cross-reference often.
4. **Document the decisions.** When choosing a lattice type or threshold rule, write down why. Future-you will not remember.
5. **Don't oversell.** The goal is to learn whether the model has merit, not to prove it does.
## Open theoretical questions
 
These are unresolved and might force changes to the implementation plan:
 
- What physically is being accumulated as "darkness" in the substrate?
- What triggers the first fracture in a perfectly still substrate?
- Is the lattice geometry physically motivated or arbitrary?
- What is the geometry of "inside out" at the inversion event?
- How does the model handle the observed homogeneity of the CMB?
Each of these is a candidate for a research note in `docs/notes/`.