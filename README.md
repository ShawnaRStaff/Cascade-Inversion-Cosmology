# Void Cascade Cosmology Model
 
A computational research project exploring whether the observed structure of the universe can be reproduced by self-organized criticality dynamics in a structured infinite substrate, without invoking general relativity as the foundation or a singular Big Bang event.
 
## Core proposition
 
The universe emerges from a cascade of microscopic fracture events in a pre-existing structured void, not from a singularity. What we observe as cosmic expansion is the propagating wavefront of activity through this substrate, not the metric stretching of space itself. The "Big Bang" is the catastrophic percolation transition when the cascade becomes globally connected, not an act of creation.
 
This model is qualitative in its current form. The goal of the repository is to develop a quantitative version through simulation and to test its predictions against observations of the cosmic web, the CMB, and large-scale structure.
 
## Physics frameworks being synthesized
 
The model combines three published research programs:
 
1. **Cascade vacuum relaxation** (Lukash & Mikheeva, 2024, arXiv:2506.03226) for the energy-step structure of vacuum evolution
2. **Self-organized criticality** (Bak, Tang & Wiesenfeld, 1987; Moffat, 1997; Carfora & Marzuoli, 2023) for microscopic avalanche dynamics
3. **Percolation theory of cosmological phase transitions** (Guth & Tye, 1980; Turner, Weinberg & Widrow, 1992; Gould & Tenkanen, 2021) for the catastrophic transition event
Original contributions of this model: the microscopic-fracture-to-effective-field bridge, the expansion-as-propagation interpretation, the inversion event as structural failure rather than singularity, and the substrate-accumulated "darkness" as a candidate component of dark matter.
 
## What this repository will contain
 
Code:
 
- A 3D sandpile simulation with cosmological reinterpretation (in progress)
- Analysis tools for avalanche size distributions, fracture cluster geometry, and percolation behavior
- Comparison utilities against observational cosmic web statistics
Documents:
 
- The model summary and glossary (see `docs/`)
- Research notes as the model develops (see `docs/notes/`)
- The technical roadmap for the simulation work (`ROADMAP.md`)
## Status

Milestone 1 (1D sandpile prototype) is complete: an Oslo rice-pile model with finite-size scaling across L = 32, 64, 128, 256 reproduces the published exponents (tau ~ 1.55, D ~ 2.2, conservation D(2-tau) ~ 1). Milestone 2 (2D sandpile with visualization) is in progress: the Manna stochastic model is implemented and the FSS run is being collected. Milestone 3 (3D + cosmological reinterpretation) has not started.

## Setup

This project targets Python 3.11+. From the repo root:

```
python -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Dependencies are numpy, scipy, matplotlib, and pytest. plotly is reserved for 3D work in Milestone 3.

## Running the simulations

All scripts write their output to `data/outputs/` with timestamped filenames. Nothing opens an interactive window; results are plain PNG / GIF / NPZ files you open yourself.

1D Oslo sandpile (Milestone 1):

```
.venv/bin/python scripts/run_milestone1.py       # single L, avalanche PDF + fit
.venv/bin/python scripts/run_milestone1_fss.py   # full FSS across four L
.venv/bin/python scripts/animate_1d.py           # GIF of pile dynamics
```

2D Manna sandpile (Milestone 2):

```
.venv/bin/python scripts/run_milestone2.py       # FSS + cluster fractal dimension
.venv/bin/python scripts/animate_2d.py           # dynamics GIF + large-avalanche snapshot
```

Tests (under a minute for the geometry tests; the sandpile tests take a couple minutes because they include short steady-state runs):

```
.venv/bin/python -m pytest
```
 
## How to read this repository
 
Start with the documents in `docs/` for the theoretical framework, then `ROADMAP.md` for the implementation plan. Code in `src/` will be added as milestones are reached.
 
## Approach to development
 
This is a personal research project, not a polished software product. The repository emphasizes clarity of thinking over engineering perfection. Code will be readable, documents will state what is solid versus what is hand-wavy, and the model will be revised as evidence accumulates or contradicts.
 
## License
 
To be decided. 
 
## References
 
The three foundational papers are the priority reading list:
 
- Lukash, V. N. & Mikheeva, E. V. (2024). "Cascade Relaxation of the Gravitating Vacuum as a Generator of the Evolving Universe." arXiv:2506.03226
- Bak, P., Tang, C. & Wiesenfeld, K. (1987). "Self-organized criticality: An explanation of 1/f noise." Physical Review Letters 59, 381
- Carfora, M. & Marzuoli, A. (2023). "Primordial Cosmology from Self-Organized Criticality."
Secondary reading:
 
- Moffat, J. W. (1997). "A Self-Organized Critical Universe." arXiv:gr-qc/9702014
- Ashtekar, A., Pawlowski, T. & Singh, P. (2006). On the Big Bounce in loop quantum cosmology
- Aschwanden, M. J. (2011). "Self-Organized Criticality in Astrophysics: The Statistics of Nonlinear Processes in the Universe"
- Hagedorn, R. (1965). On the Hagedorn temperature and hadronic matter