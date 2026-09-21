"""Conservative finite-volume solver for one-dimensional Fickian diffusion."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.sparse import csc_matrix, diags
from scipy.sparse.linalg import spsolve

from .thermal import DiffusivityModel, ThermalSchedule

BoundaryKind = Literal["no_flux", "fixed_concentration"]


@dataclass(frozen=True)
class SpatialGrid:
    """Uniform cell-centred one-dimensional grid.

    ``length_m`` is the simulated depth and ``cells`` is the number of finite
    volumes. Concentrations live at cell centres, not at boundary faces.
    """

    length_m: float
    cells: int

    def __post_init__(self) -> None:
        if not np.isfinite(self.length_m) or self.length_m <= 0.0:
            raise ValueError("length_m must be finite and positive.")
        if self.cells < 3:
            raise ValueError("cells must be at least 3.")

    @property
    def spacing_m(self) -> float:
        """Finite-volume cell width in m."""

        return self.length_m / self.cells

    @property
    def centres_m(self) -> NDArray[np.float64]:
        """Cell-centre depths in m."""

        return (np.arange(self.cells, dtype=float) + 0.5) * self.spacing_m


@dataclass(frozen=True)
class BoundaryCondition:
    """One boundary condition for the finite-volume domain.

    A fixed concentration represents a Dirichlet boundary at a cell face. A
    no-flux boundary represents a zero concentration gradient. Concentrations
    are chemical concentrations in m⁻³.
    """

    kind: BoundaryKind = "no_flux"
    concentration_m3: float = 0.0

    def __post_init__(self) -> None:
        if self.kind not in {"no_flux", "fixed_concentration"}:
            raise ValueError("kind must be 'no_flux' or 'fixed_concentration'.")
        if not np.isfinite(self.concentration_m3) or self.concentration_m3 < 0.0:
            raise ValueError("concentration_m3 must be finite and non-negative.")
        if self.kind == "fixed_concentration" and self.concentration_m3 <= 0.0:
            raise ValueError("A fixed concentration must be positive.")


@dataclass(frozen=True)
class DiffusionResult:
    """Final concentration profile and numerical conservation diagnostics."""

    grid: SpatialGrid
    concentration_m3: NDArray[np.float64]
    thermal_budget_m2: float
    time_steps: int
    initial_dose_m2: float
    final_dose_m2: float
    cumulative_boundary_dose_change_m2: float

    @property
    def dose_balance_error_m2(self) -> float:
        """Final-minus-initial dose minus discrete boundary contribution (m⁻²)."""

        return self.final_dose_m2 - self.initial_dose_m2 - self.cumulative_boundary_dose_change_m2


def _validate_concentration_m3(
    concentration_m3: ArrayLike, grid: SpatialGrid
) -> NDArray[np.float64]:
    concentration = np.asarray(concentration_m3, dtype=float)
    if concentration.shape != (grid.cells,):
        raise ValueError(f"concentration_m3 must have shape ({grid.cells},).")
    if np.any(~np.isfinite(concentration)) or np.any(concentration < 0.0):
        raise ValueError("concentration_m3 must contain finite, non-negative values.")
    return concentration


def integrated_dose_m2(concentration_m3: ArrayLike, grid: SpatialGrid) -> float:
    """Return cell-integrated areal chemical dose in m⁻²."""

    concentration = _validate_concentration_m3(concentration_m3, grid)
    return float(np.sum(concentration) * grid.spacing_m)


def _operator_and_source(
    grid: SpatialGrid, left: BoundaryCondition, right: BoundaryCondition
) -> tuple[csc_matrix, NDArray[np.float64]]:
    """Build the finite-volume Laplacian and constant boundary source.

    ``operator @ concentration + source`` has units m⁻⁵ when concentration is
    m⁻³. The finite-volume formulation gives exact discrete conservation when
    both boundaries have zero flux.
    """

    cells = grid.cells
    inverse_spacing_squared = 1.0 / grid.spacing_m**2
    lower = np.ones(cells - 1, dtype=float)
    diagonal = -2.0 * np.ones(cells, dtype=float)
    upper = np.ones(cells - 1, dtype=float)
    source = np.zeros(cells, dtype=float)

    if left.kind == "no_flux":
        diagonal[0] = -1.0
    else:
        diagonal[0] = -3.0
        source[0] = 2.0 * left.concentration_m3 * inverse_spacing_squared

    if right.kind == "no_flux":
        diagonal[-1] = -1.0
    else:
        diagonal[-1] = -3.0
        source[-1] = 2.0 * right.concentration_m3 * inverse_spacing_squared

    operator = diags((lower, diagonal, upper), offsets=(-1, 0, 1), format="csc")
    return csc_matrix(operator * inverse_spacing_squared), source


def _boundary_dose_increment_m2(
    concentration_m3: NDArray[np.float64],
    grid: SpatialGrid,
    left: BoundaryCondition,
    right: BoundaryCondition,
    thermal_budget_increment_m2: float,
) -> float:
    """Return net dose injected through boundaries during one implicit step."""

    inverse_half_spacing = 2.0 / grid.spacing_m
    left_influx = 0.0
    right_influx = 0.0
    if left.kind == "fixed_concentration":
        left_influx = inverse_half_spacing * (left.concentration_m3 - concentration_m3[0])
    if right.kind == "fixed_concentration":
        right_influx = inverse_half_spacing * (right.concentration_m3 - concentration_m3[-1])
    return float(thermal_budget_increment_m2 * (left_influx + right_influx))


def solve_diffusion(
    initial_concentration_m3: ArrayLike,
    grid: SpatialGrid,
    schedule: ThermalSchedule,
    diffusivity: DiffusivityModel,
    *,
    left_boundary: BoundaryCondition = BoundaryCondition(),
    right_boundary: BoundaryCondition = BoundaryCondition(),
    substeps_per_segment: int = 50,
) -> DiffusionResult:
    """Advance a chemical concentration profile through a thermal schedule.

    Backward Euler is applied in integrated diffusivity, ``Θ = ∫D[T(t)]dt``, not by treating a
    temperature ramp as a single arbitrary endpoint temperature. It is stable
    for this linear model, but temporal and spatial convergence remain required
    validation checks.
    """

    if substeps_per_segment < 1:
        raise ValueError("substeps_per_segment must be at least 1.")
    concentration = _validate_concentration_m3(initial_concentration_m3, grid).copy()
    operator, source = _operator_and_source(grid, left_boundary, right_boundary)
    increments_m2 = schedule.segment_thermal_budgets_m2(diffusivity, substeps_per_segment)
    identity = diags(np.ones(grid.cells, dtype=float), offsets=0, format="csc")
    cumulative_boundary_dose_change_m2 = 0.0
    initial_dose_m2 = integrated_dose_m2(concentration, grid)

    for increment_m2 in increments_m2:
        system = identity - increment_m2 * operator
        right_hand_side = concentration + increment_m2 * source
        concentration = np.asarray(spsolve(system, right_hand_side), dtype=float)
        concentration_scale = max(float(np.max(concentration)), 1.0)
        if np.min(concentration) < -1.0e-12 * concentration_scale:
            raise RuntimeError(
                "Implicit diffusion step produced a non-physical negative concentration."
            )
        concentration = np.maximum(concentration, 0.0)
        cumulative_boundary_dose_change_m2 += _boundary_dose_increment_m2(
            concentration, grid, left_boundary, right_boundary, increment_m2
        )

    final_dose_m2 = integrated_dose_m2(concentration, grid)
    return DiffusionResult(
        grid=grid,
        concentration_m3=concentration,
        thermal_budget_m2=float(sum(increments_m2)),
        time_steps=len(increments_m2),
        initial_dose_m2=initial_dose_m2,
        final_dose_m2=final_dose_m2,
        cumulative_boundary_dose_change_m2=cumulative_boundary_dose_change_m2,
    )
