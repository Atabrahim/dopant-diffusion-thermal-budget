# Physical model and assumptions

## Chemical concentration transport

The package models the chemical boron concentration \(C(x,t)\) in silicon using

\[
\frac{\partial C}{\partial t}=\frac{\partial}{\partial x}\left[D(T)\frac{\partial C}{\partial x}\right].
\]

| Symbol | Meaning | SI unit |
| --- | --- | --- |
| \(C\) | chemical boron concentration | m\(^{-3}\) |
| \(x\) | depth below the surface | m |
| \(t\) | time | s |
| \(D\) | boron chemical diffusivity | m\(^2\) s\(^{-1}\) |
| \(T\) | absolute temperature | K |

For V1, \(D\) is spatially uniform and depends only on the supplied thermal history:

\[
D(T)=D_0\exp\left[-\frac{E_a}{k_B T}\right].
\]

The default boron-in-silicon parameterisation is \(D_0=0.76\ \mathrm{cm^2\,s^{-1}}\) and \(E_a=3.46\ \mathrm{eV}\), exposed as inputs rather than universal material constants. Reported values depend on concentration, interstitial/vacancy conditions, measurement technique and the applicable temperature range.

## Boundary-condition cases

### Finite source

A buried, dose-normalised reflected Gaussian is evolved with zero flux at both boundaries. The surface condition is \(\partial C/\partial x=0\), so total simulated dose changes only through explicitly reported numerical flux balance. This is a useful model for limited-dose drive-in, not an exact implantation simulator.

### Constant surface concentration

The surface concentration is fixed to \(C_s\) and the deep boundary is held at zero when the domain is chosen sufficiently long. On a semi-infinite domain starting from zero concentration, the analytical result is

\[
C(x,\Theta)=C_s\operatorname{erfc}\left(\frac{x}{2\sqrt{\Theta}}\right),
\qquad
\Theta=\int_0^t D[T(\tau)]\,d\tau.
\]

The finite computational domain is therefore selected wide enough that the far-boundary approximation is quantitatively checked.

## Numerical method

Cell-centred finite volumes discretise diffusive face fluxes. Backward Euler advances the solution in diffusion time, using a sparse linear solve. The method is unconditionally stable for this linear problem; that does **not** remove the need for mesh and time-step convergence checks.

## What this model does not claim

Chemical concentration is not the same as activated dopant density or carrier concentration. The model omits concentration-dependent diffusion, transient-enhanced diffusion, clustering, segregation, oxidation effects, implant damage and two-dimensional geometry. It is an educational, traceable process-model foundation rather than an industrial TCAD replacement.

## References

1. S. M. Sze and K. K. Ng, *Physics of Semiconductor Devices*, 3rd ed., Wiley, 2007, diffusion fundamentals and process context.
2. S. Wolf and R. N. Tauber, *Silicon Processing for the VLSI Era, Vol. 1*, 2nd ed., Lattice Press, 2000, dopant diffusion data and process assumptions.
3. J. Crank, *The Mathematics of Diffusion*, 2nd ed., Oxford University Press, 1975, Gaussian and complementary-error-function solutions.
