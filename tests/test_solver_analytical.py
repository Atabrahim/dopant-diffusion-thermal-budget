"""Scientific validation tests for the diffusion solver."""

from __future__ import annotations

import numpy as np
import pytest

from dopant_diffusion.materials import BoronInSilicon
from dopant_diffusion.profiles import (
    constant_surface_erfc_solution_m3,
    reflected_gaussian_m3,
    reflected_gaussian_solution_m3,
)
from dopant_diffusion.solver import (
    BoundaryCondition,
    SpatialGrid,
    integrated_dose_m2,
    solve_diffusion,
)
from dopant_diffusion.thermal import ThermalSchedule, ThermalSegment


class ConstantDiffusivity:
    """Test-only diffusivity model independent of temperature."""

    def __init__(self, value_m2_per_s: float) -> None:
        self.value_m2_per_s = value_m2_per_s

    def diffusivity_m2_per_s(self, temperature_K: float | np.ndarray) -> np.ndarray:
        return np.asarray(temperature_K, dtype=float) * 0.0 + self.value_m2_per_s


def _isothermal_schedule(duration_s: float) -> ThermalSchedule:
    return ThermalSchedule((ThermalSegment(1_273.15, 1_273.15, duration_s),))


def test_no_flux_solver_conserves_discrete_dose() -> None:
    grid = SpatialGrid(length_m=2.0e-6, cells=400)
    initial = reflected_gaussian_m3(grid.centres_m, 2.0e16, 4.0e-7, 4.0e-8)
    result = solve_diffusion(
        initial,
        grid,
        _isothermal_schedule(100.0),
        ConstantDiffusivity(2.0e-15),
        substeps_per_segment=100,
    )

    relative_dose_change = (
        abs(result.final_dose_m2 - result.initial_dose_m2) / result.initial_dose_m2
    )
    assert relative_dose_change < 1.0e-12
    assert result.cumulative_boundary_dose_change_m2 == pytest.approx(0.0, abs=1.0e-18)
    assert abs(result.dose_balance_error_m2) / result.initial_dose_m2 < 1.0e-12
    assert np.min(result.concentration_m3) >= 0.0


def test_reflected_gaussian_agrees_with_analytical_solution() -> None:
    grid = SpatialGrid(length_m=4.0e-6, cells=800)
    dose_m2 = 2.0e16
    projected_range_m = 4.0e-7
    initial_sigma_m = 4.0e-8
    duration_s = 100.0
    diffusivity_m2_per_s = 2.0e-15
    initial = reflected_gaussian_m3(grid.centres_m, dose_m2, projected_range_m, initial_sigma_m)
    result = solve_diffusion(
        initial,
        grid,
        _isothermal_schedule(duration_s),
        ConstantDiffusivity(diffusivity_m2_per_s),
        substeps_per_segment=1_000,
    )
    reference = reflected_gaussian_solution_m3(
        grid.centres_m,
        dose_m2,
        projected_range_m,
        initial_sigma_m,
        diffusivity_m2_per_s * duration_s,
    )

    relative_l2_error = np.linalg.norm(result.concentration_m3 - reference) / np.linalg.norm(
        reference
    )
    assert relative_l2_error < 0.001


def test_constant_surface_case_agrees_with_erfc_far_from_domain_boundary() -> None:
    grid = SpatialGrid(length_m=4.0e-6, cells=1_200)
    surface_concentration_m3 = 1.0e25
    duration_s = 100.0
    diffusivity_m2_per_s = 2.0e-15
    result = solve_diffusion(
        np.zeros(grid.cells),
        grid,
        _isothermal_schedule(duration_s),
        ConstantDiffusivity(diffusivity_m2_per_s),
        left_boundary=BoundaryCondition("fixed_concentration", surface_concentration_m3),
        right_boundary=BoundaryCondition("fixed_concentration", 1.0e-30),
        substeps_per_segment=1_000,
    )
    reference = constant_surface_erfc_solution_m3(
        grid.centres_m, surface_concentration_m3, diffusivity_m2_per_s * duration_s
    )
    comparison = grid.centres_m < 1.5e-6
    relative_l2_error = np.linalg.norm(
        result.concentration_m3[comparison] - reference[comparison]
    ) / np.linalg.norm(reference[comparison])

    assert relative_l2_error < 0.015
    assert abs(result.dose_balance_error_m2) / result.final_dose_m2 < 1.0e-10


def test_refining_the_mesh_reduces_analytical_error() -> None:
    errors: list[float] = []
    for cells in (200, 400):
        grid = SpatialGrid(length_m=4.0e-6, cells=cells)
        initial = reflected_gaussian_m3(grid.centres_m, 2.0e16, 4.0e-7, 4.0e-8)
        result = solve_diffusion(
            initial,
            grid,
            _isothermal_schedule(100.0),
            ConstantDiffusivity(2.0e-15),
            substeps_per_segment=4_000,
        )
        reference = reflected_gaussian_solution_m3(grid.centres_m, 2.0e16, 4.0e-7, 4.0e-8, 2.0e-13)
        errors.append(
            float(np.linalg.norm(result.concentration_m3 - reference) / np.linalg.norm(reference))
        )

    assert errors[1] < errors[0]


def test_refining_backward_euler_steps_reduces_analytical_error() -> None:
    grid = SpatialGrid(length_m=4.0e-6, cells=800)
    initial = reflected_gaussian_m3(grid.centres_m, 2.0e16, 4.0e-7, 4.0e-8)
    reference = reflected_gaussian_solution_m3(grid.centres_m, 2.0e16, 4.0e-7, 4.0e-8, 2.0e-13)
    errors: list[float] = []

    for substeps in (1_000, 4_000):
        result = solve_diffusion(
            initial,
            grid,
            _isothermal_schedule(100.0),
            ConstantDiffusivity(2.0e-15),
            substeps_per_segment=substeps,
        )
        errors.append(
            float(np.linalg.norm(result.concentration_m3 - reference) / np.linalg.norm(reference))
        )

    assert errors[1] < errors[0]


def test_solver_rejects_invalid_input_shape() -> None:
    grid = SpatialGrid(length_m=1.0e-6, cells=20)
    with pytest.raises(ValueError, match="shape"):
        solve_diffusion(
            np.zeros(19),
            grid,
            _isothermal_schedule(1.0),
            BoronInSilicon(),
        )


def test_integrated_dose_requires_non_negative_concentration() -> None:
    grid = SpatialGrid(length_m=1.0e-6, cells=20)
    with pytest.raises(ValueError, match="non-negative"):
        integrated_dose_m2(-np.ones(grid.cells), grid)
