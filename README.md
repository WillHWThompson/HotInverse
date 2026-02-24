# HotInverse

![Status](https://img.shields.io/badge/status-research%20artifact-blue)
![Language](https://img.shields.io/badge/language-Julia-informational)
![Focus](https://img.shields.io/badge/focus-inverse%20optimization-success)

Inverse optimization framework for HOT-style facility placement and spatial robustness experiments.

## Project Highlights
- Julia package-style source organization with modular algorithms and data structures.
- Scripted experiment drivers for optimization sweeps and constraint evaluation.
- Integrated data + plotting directories for reproducible analysis outputs.
- Lightweight overview notebook for reviewer onboarding.

## What This Project Is
`HotInverse` is a Julia research codebase for optimization and simulation of spatial facility layouts under constraints and perturbations. It includes:
- core optimization/data structures in `src/`,
- experiment drivers in `scripts/`,
- data assets and intermediate outputs,
- notebook-style exploratory analysis.

## Repository Layout
- `src/` — package source (`HotInverse.jl`) and algorithm modules
- `scripts/` — experiment and optimization entry points
- `data/` — experiment inputs and generated outputs
- `notebooks/` — analysis notebooks
- `plots/` — generated visual outputs
- `test/` — Julia tests

## Quickstart
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Typical run pattern:
```bash
julia --project=. scripts/run_num_fac_optimization.jl
```

## Example Notebook
For a high-level usage walkthrough:
- `notebooks/quickstart_overview.ipynb`

## Notes for Reviewers
This repository is maintained as a research artifact and is intended to present complete project structure, methodology, and workflow. Some historical scripts and datasets reflect ongoing experimentation rather than polished production APIs.

## Project Status
- Maturity: advanced research prototype
- Language: Julia (with DrWatson-style project organization)
- Focus: methodological experimentation and inverse spatial optimization

## Author
Will Thompson
