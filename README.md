# Dopant Diffusion and Thermal-Budget Simulation

**A conservative one-dimensional Python simulator for boron diffusion in silicon during ramp-and-hold thermal processes.**

> Development status: the repository scaffold and fixed Version 1 scientific scope are published. The numerical solver, validation suite, process reports and final release follow as separately verified milestones.

## Why this project

Dopant diffusion converts a spatial dopant distribution into a device-relevant concentration profile. In semiconductor processing, the temperature-time history of an anneal controls how far an implanted or pre-deposited dopant spreads. That directly affects shallow junctions, sheet resistance and process integration windows.

This project models the **chemical concentration** of one dopant (boron) in silicon. It deliberately does not convert concentration to electrically active carrier density: activation, clustering, transient-enhanced diffusion, point defects and concentration-dependent diffusivity require a more detailed model.

## Fixed V1 scope

- A temperature-dependent Arrhenius diffusivity for boron in silicon.
- Piecewise-linear ramp-and-hold temperature schedules.
- A one-dimensional, conservative finite-volume discretisation of Fick's second law.
- Implicit time integration for finite-source and fixed-surface-concentration boundary cases.
- Dose, concentration-threshold depth and thermal-budget diagnostics.
- Analytical validation against reflected-Gaussian and complementary-error-function diffusion solutions.
- Mesh/time-step convergence, positivity and dose/flux-balance checks.
- A documented command-line workflow that generates reproducible figures and process-comparison reports.

Future work is intentionally excluded from V1: concentration-dependent diffusivity, electrical activation, oxidation-enhanced diffusion, transient-enhanced diffusion, two-dimensional effects and a full TCAD process flow.

## Scientific model

For a spatially uniform but time-dependent diffusivity, the chemical concentration \(C(x,t)\) obeys Fick's second law:

\[
\frac{\partial C}{\partial t}=\frac{\partial}{\partial x}\left[D(T)\frac{\partial C}{\partial x}\right].
\]

Here \(x\) is depth in metres, \(t\) is time in seconds, \(C\) is concentration in m\(^{-3}\), and \(D\) is diffusivity in m\(^2\) s\(^{-1}\). The simplified Arrhenius model is

\[
D(T) = D_0\exp\left(-\frac{E_a}{k_B T}\right),
\]

where \(D_0\) is the prefactor, \(E_a\) the activation energy, \(k_B\) Boltzmann's constant and \(T\) absolute temperature. See [the model documentation](docs/MODEL.md) for assumptions, boundary conditions, units and references.

## Planned results

The completed release will include concentration-versus-depth profiles, temperature schedules, thermal-budget comparisons, dose-balance diagnostics and analytical/numerical agreement figures.

## Installation and usage

Installation and command-line examples will be added only after they are verified against the packaged implementation.

## License

The original software is released under the [MIT License](LICENSE).
