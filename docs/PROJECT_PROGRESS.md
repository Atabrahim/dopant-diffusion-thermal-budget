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
| Create repository, scope, model documentation and recovery record | VERIFIED, PUBLISHED | Repository and baseline documentation are published. |
| Implement material and thermal-schedule core | VERIFIED, PUBLISHED | Arrhenius boron model, schedule integration and six focused tests are published. |
| Implement conservative diffusion solver and analytical references | VERIFIED, PUBLISHED | Sparse finite-volume backward-Euler solver and Gaussian/erfc references are published. |
| Validate conservation, analytical cases and convergence | VERIFIED, PUBLISHED | 15 tests pass locally; analytical comparison and mesh/time refinement are included. |
| Add process comparison, figures and reproducible reports | VERIFIED, PUBLISHED | Reproduction script generates three inspected figures, two profile reports and six sensitivity rows; remote tree matches local files. |
| Package, CLI, documentation and CI | IN PROGRESS | JSON-driven CLI and documentation are published; CI and clean installation remain. |
| Final QA and v0.1.0 release | NOT STARTED | — |

## Scientific validation

The sparse finite-volume solver was tested against independent analytical solutions at constant diffusivity. On a 4 µm domain with 800 cells and 1,000 implicit substeps, the reflected-Gaussian relative L2 error is `1.3647e-4`; the finite-source relative dose-balance error is `6.4e-14`. For a constant surface concentration on a 4 µm domain with 1,200 cells, the erfc relative L2 error below 1.5 µm is `1.8656e-4`; the flux-balance error relative to final dose is `6.02e-13`.

The suite also checks that refining both the mesh and the backward-Euler substeps reduces analytical error. The Gaussian validation domain was enlarged from 2 µm to 4 µm after diagnosis showed that the far no-flux boundary was contaminating a semi-infinite analytical comparison. This is recorded as a validation correction, not hidden by loosening the test.

## Test status

`19 passed` locally with the process workflow and CLI tests. Ruff lint and formatting checks pass; the documented CLI and reproduction script run end-to-end.

## Known issues

- The selected diffusivity is a representative dilute-Fickian parameterisation, not a calibration of a particular furnace or implant process.
- The supplied Arrhenius parameters are illustrative and configurable; they do not represent a calibrated furnace/implant flow.
- The current solver assumes spatially uniform diffusivity, so it cannot model concentration- or defect-dependent diffusion.

## Remaining work

1. Publish this progress update and the figure/report provenance notes.
2. Configure GitHub Actions and verify a successful run.
3. Build sdist/wheel, install in a clean environment, verify README commands and release QA.

## Latest verified GitHub checkpoint

`1d7ea5f` — process reports, figures, CLI, tests and README published to `main`; local/remote file trees compared. This progress update is the next change to publish.
