# Project progress — Dopant Diffusion and Thermal-Budget Simulation

## Original objective

Build a scientifically transparent Python package that evolves one-dimensional boron chemical-concentration profiles in silicon through ramp-and-hold thermal schedules. The project must compare finite-source and fixed-surface-concentration boundary cases, quantify process sensitivity, and validate numerical results against analytical diffusion solutions and conservation relationships.

## Agreed V1 scope

One configurable dilute-Fickian boron-in-silicon diffusivity model; piecewise-linear thermal schedules; a conservative one-dimensional finite-volume solver with implicit integration; two boundary-condition cases; dose, threshold-depth and thermal-budget diagnostics; validated figures and process reports; package, CLI, tests and CI.

Excluded from V1: electrical activation; concentration-dependent or defect-mediated diffusion; transient-enhanced diffusion; oxidation effects; two-dimensional geometry; calibrated industrial process prediction.

## Architecture

| Component | Role |
| --- | --- |
| `materials.py` | Boron-in-silicon Arrhenius diffusivity and parameter checks. |
| `thermal.py` | Ramp/hold schedule representation and integrated diffusivity. |
| `profiles.py` | Initial profiles and analytical reference solutions. |
| `solver.py` | Conservative finite-volume implicit PDE solver. |
| `diagnostics.py` | Dose, flux balance, threshold-depth and convergence diagnostics. |
| `workflow.py` / `cli.py` | Reproducible config-driven simulations and reports. |
| `plotting.py` | Scientific figures. |

## Milestones

| Milestone | Status | Evidence |
| --- | --- | --- |
| Recover original roadmap and inspect Projects 1–3 | VERIFIED, PUBLISHED | Original portfolio record and GitHub repositories checked. |
| Create repository, scope, model documentation and recovery record | VERIFIED, PUBLISHED | Initial GitHub commit `d2f3a9a`. |
| Implement material and thermal-schedule core | IN PROGRESS | Local implementation awaiting tests and publication. |
| Implement conservative diffusion solver and analytical references | NOT STARTED | — |
| Validate conservation, analytical cases and convergence | NOT STARTED | — |
| Add process comparison, figures and reproducible reports | NOT STARTED | — |
| Package, CLI, documentation and CI | NOT STARTED | — |
| Final QA and v0.1.0 release | NOT STARTED | — |

## Scientific validation

No numerical PDE result has yet been validated. The initial documentation defines the planned analytical reference cases: a reflected Gaussian under no-flux boundaries and the complementary-error-function solution for constant surface concentration.

## Test status

No automated tests have been run for Project 4 yet.

## Known issues

- The selected diffusivity is a representative dilute-Fickian parameterisation, not a calibration of a particular furnace or implant process.
- Thermal schedule and material modules are untested locally at this checkpoint.

## Remaining work

1. Add tests for Arrhenius and thermal-budget behavior.
2. Implement and validate the conservative solver.
3. Build diagnostics, figures, reports and CLI.
4. Perform package/CI/release QA.

## Latest verified GitHub checkpoint

`d2f3a9a` — initial GitHub-created scaffold. The local additions in the working tree are newer and must be tested and published before the next milestone is marked complete.
