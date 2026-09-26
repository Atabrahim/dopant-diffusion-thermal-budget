# Numerical validation

## Purpose

The package separates validation of the numerical diffusion solver from interpretation of the synthetic thermal case study. The analytical references below assume a semi-infinite, constant-diffusivity medium. They are not experimental boron profiles.

## Independent reference cases

### Finite source

A dose-normalised reflected Gaussian is evolved with zero flux at both computational boundaries. At constant diffusivity, its analytical standard deviation follows

$$
\sigma^2(t) = \sigma_0^2 + 2Dt.
$$

The validation compares the numerical profile with that reflected-Gaussian solution on a 4 µm domain.

### Constant surface concentration

For a semi-infinite material initially at zero concentration with fixed surface concentration $C_s$,

$$
C(x,t) = C_s\,\operatorname{erfc}\left(\frac{x}{2\sqrt{Dt}}\right).
$$

The numerical solver uses a fixed concentration at the surface and a zero deep boundary. The comparison is restricted to depth below 1.5 µm, well inside the validated 4 µm domain.

## Recorded results

| Test quantity | Verified result |
| --- | ---: |
| Reflected-Gaussian relative L2 error, 800 cells and 1,000 implicit steps | $1.3647\times10^{-4}$ |
| Finite-source relative dose-balance error | $6.4\times10^{-14}$ |
| Constant-surface erfc relative L2 error, 1,200 cells | $1.8656\times10^{-4}$ |
| Fixed-surface relative flux-balance error | $6.02\times10^{-13}$ |

The automated suite also checks that jointly refining the spatial mesh and backward-Euler substeps reduces Gaussian-reference error, and that time-step refinement reduces the same error at fixed mesh.

## Boundary-domain diagnosis

A first Gaussian comparison on a 2 µm no-flux domain showed approximately 1.05% relative L2 error that did not improve with mesh refinement. A parameter sweep showed the cause was physical boundary contamination: a zero-flux wall is not the semi-infinite boundary assumed by the Gaussian formula. Enlarging the comparison domain to 4 µm reduced the error to $1.3647\times10^{-4}$ at the documented discretisation. The test was corrected by changing the domain, not by loosening tolerance.

## Conservation diagnostics

For no-flux finite-source cases, the solver reports

$$
\Delta Q = Q_{\mathrm{final}} - Q_{\mathrm{initial}} -
\int_0^{t}(J_{\mathrm{left}}-J_{\mathrm{right}})\,dt,
$$

where $Q=\int C\,dx$ is dose per unit area and $J=-D\,\partial C/\partial x$. The finite-volume scheme reports this balance rather than silently assuming conservation. Fixed-concentration boundaries are expected to exchange dose with the domain, so their flux balance—not equality of initial and final dose—is tested.

## What validation does not establish

Analytical agreement and balance checks establish implementation correctness only in the stated dilute-Fickian regime. They do not validate the default Arrhenius parameters for a specific furnace, implant energy, dose, or wafer. They also do not include activation, defects, high-concentration effects, or process metrology.
