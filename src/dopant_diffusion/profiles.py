"""Initial concentration profiles and analytical diffusion reference cases."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.special import erfc


def _validated_depth_m(depth_m: ArrayLike) -> NDArray[np.float64]:
    depth = np.asarray(depth_m, dtype=float)
    if depth.ndim != 1 or len(depth) < 2:
        raise ValueError("depth_m must be a one-dimensional array with at least two values.")
    if np.any(~np.isfinite(depth)) or np.any(depth < 0.0):
        raise ValueError("depth_m must contain finite, non-negative values.")
    if np.any(np.diff(depth) <= 0.0):
        raise ValueError("depth_m must be strictly increasing.")
    return depth


def reflected_gaussian_m3(
    depth_m: ArrayLike,
    dose_m2: float,
    projected_range_m: float,
    standard_deviation_m: float,
) -> NDArray[np.float64]:
    """Return a dose-normalised Gaussian reflected at a no-flux surface.

    The mirror term makes the profile integrate to ``dose_m2`` on the ideal
    half-line and imposes zero concentration gradient at the surface. ``dose_m2``
    is an areal chemical dose (m⁻²), so the returned concentration is m⁻³.
    """

    depth = _validated_depth_m(depth_m)
    parameters = (dose_m2, projected_range_m, standard_deviation_m)
    if not all(np.isfinite(value) for value in parameters):
        raise ValueError("Gaussian-profile parameters must be finite.")
    if dose_m2 <= 0.0 or projected_range_m < 0.0 or standard_deviation_m <= 0.0:
        raise ValueError(
            "Dose and standard deviation must be positive; projected range is non-negative."
        )

    denominator = np.sqrt(2.0 * np.pi) * standard_deviation_m
    forward = np.exp(-0.5 * ((depth - projected_range_m) / standard_deviation_m) ** 2)
    mirror = np.exp(-0.5 * ((depth + projected_range_m) / standard_deviation_m) ** 2)
    return np.asarray(dose_m2 * (forward + mirror) / denominator, dtype=float)


def reflected_gaussian_solution_m3(
    depth_m: ArrayLike,
    dose_m2: float,
    projected_range_m: float,
    initial_standard_deviation_m: float,
    thermal_budget_m2: float,
) -> NDArray[np.float64]:
    """Analytical no-flux solution for a reflected Gaussian initial profile.

    For uniform-in-space diffusivity, the thermal history enters through
    ``thermal_budget_m2 = ∫D(T)dt`` and the variance evolves as
    ``σ² = σ₀² + 2 Θ``.
    """

    if not np.isfinite(thermal_budget_m2) or thermal_budget_m2 < 0.0:
        raise ValueError("thermal_budget_m2 must be finite and non-negative.")
    standard_deviation_m = np.sqrt(initial_standard_deviation_m**2 + 2.0 * thermal_budget_m2)
    return reflected_gaussian_m3(depth_m, dose_m2, projected_range_m, standard_deviation_m)


def constant_surface_erfc_solution_m3(
    depth_m: ArrayLike, surface_concentration_m3: float, thermal_budget_m2: float
) -> NDArray[np.float64]:
    """Return the semi-infinite constant-surface-concentration solution.

    The result is ``C_s erfc[x / (2 sqrt(Θ))]`` for initially dopant-free
    material, where ``Θ = ∫D(T)dt``. It is an analytical *reference* for the
    numerical finite-domain model, not experimental data.
    """

    depth = _validated_depth_m(depth_m)
    if not np.isfinite(surface_concentration_m3) or surface_concentration_m3 <= 0.0:
        raise ValueError("surface_concentration_m3 must be finite and positive.")
    if not np.isfinite(thermal_budget_m2) or thermal_budget_m2 < 0.0:
        raise ValueError("thermal_budget_m2 must be finite and non-negative.")
    if thermal_budget_m2 == 0.0:
        return np.zeros_like(depth)
    argument = depth / (2.0 * np.sqrt(thermal_budget_m2))
    return np.asarray(surface_concentration_m3 * erfc(argument), dtype=float)
