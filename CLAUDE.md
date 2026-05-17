Claude · MD
Copy

# CLAUDE.md
# This file is read by Claude Code when working in this repository.
# It provides context, conventions, and current state for the project.
 
## Project overview
 
This repository develops a computational model for an alternative cosmology framework called Cascade Inversion Cosmology (working name; may change). The proposition: the universe emerged from a cascade of microscopic fracture events in a structured infinite substrate, with the Big Bang reinterpreted as a percolation-driven inversion event rather than a singularity. The model synthesizes three published physics frameworks (cascade vacuum relaxation, self-organized criticality, percolation theory of cosmological phase transitions) with original contributions on expansion-as-propagation and inversion as structural failure.
 
The full theoretical framework is in `docs/model_summary.md`. The vocabulary mapping to standard physics is in `docs/glossary.md`. The roadmap is in `ROADMAP.md`.
 
## Current focus
 
Milestone 1 of the roadmap: build a minimal 1D sandpile simulation in Python that reproduces the known SOC dynamics and power-law avalanche size distribution. This is preparation work to verify the code infrastructure before scaling to 2D and 3D.
 
The goal at this stage is functional code with clear outputs, not engineering perfection.
 
## Tech stack
 
- Language: Python 3.11 or newer
- Core libraries: numpy, scipy, matplotlib for prototypes; plotly for 3D visualization later
- Testing: pytest
- Notebook environment: Jupyter for exploration; production code lives in `src/`
- Type hints: yes, use them
- Formatting: black, line length 100
## Code style
 
- Write code as if a physicist with moderate Python skills will read it, not just other software engineers
- Comments should explain why, not what
- Variable names should match the physics where possible (n0 for substrate density, p for percolation probability, tau for the avalanche exponent, etc.)
- Functions should do one thing and be small
- Prefer pure functions; avoid hidden state unless there's a real reason
- Avoid premature optimization; readable beats clever
## Project structure
 
```
void-cascade-cosmology/
├── README.md
├── ROADMAP.md
├── CLAUDE.md (this file)
├── docs/
│   ├── model_summary.md
│   ├── glossary.md
│   └── notes/
├── src/
│   └── void_cascade/
│       ├── __init__.py
│       ├── sandpile_1d.py
│       ├── sandpile_2d.py
│       └── sandpile_3d.py
├── tests/
├── notebooks/
└── data/
```
 
## How to engage with me
 
I am the researcher developing this model. I am also a software engineer; I know how to code. What I need from you is help thinking through the physics, suggesting clean implementations, and catching errors in either domain.
 
Push back when something doesn't make sense. If I propose a simulation approach that won't produce the output I'm describing, tell me. If I confuse two physics concepts, correct me. If a piece of code I write is wrong, fix it.
 
Don't soften assessments. I want to know if an approach is bad before I spend a day implementing it.
 
## Defaults when I haven't specified
 
- Use numpy arrays over Python lists for anything numerical
- Save outputs to `data/outputs/` with timestamped filenames
- Plot results with matplotlib; use plotly only for 3D
- Add a docstring to every function explaining what physics it implements
- When in doubt about a physics parameter, use the values from Bak, Tang & Wiesenfeld 1987 or Aschwanden 2011 as defaults
- Write tests for any function that does math; not for plotting or I/O
## Things to avoid
 
- Adding dependencies beyond the core stack unless I ask
- Overengineering with classes and abstractions before the code needs them
- Generating boilerplate that I haven't asked for (no setup.py until I want to package, no extensive logging until I have a reason)
- Treating this as crank physics. The model is grounded in serious recent literature. If you're skeptical of an aspect, say so specifically; don't dismiss generally.
- Suggesting I consult "a physicist" as if I'm not capable of evaluating ideas myself.
## Conventions
 
- Avalanche size: `s`
- Avalanche duration: `T`
- Lattice site values: array `z`
- Threshold: `z_c`
- Power law exponent: `tau`
- Percolation threshold: `p_c`
- Fractal dimension: `D`
- Time step: `dt`
When in doubt, ask. When obvious, just do it.
 
## Reading priority
 
If you have time to read background, prioritize in this order:
 
1. `docs/model_summary.md` (the theoretical framework)
2. `ROADMAP.md` (what we're building and why)
3. `docs/glossary.md` (vocabulary mapping)
4. The Bak, Tang & Wiesenfeld 1987 paper (foundational SOC)
5. Lukash & Mikheeva 2024 arXiv:2506.03226 (cascade vacuum relaxation)
6. Carfora & Marzuoli 2023 (SOC applied to primordial cosmology)