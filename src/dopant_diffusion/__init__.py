"""Conservative one-dimensional silicon dopant-diffusion simulation."""

from .materials import BoronInSilicon
from .solver import BoundaryCondition, SpatialGrid, solve_diffusion
from .thermal import ThermalSchedule, ThermalSegment
from .workflow import FiniteSourceCase, SimulationOutcome, run_finite_source_case

__all__ = [
    "BoronInSilicon",
    "BoundaryCondition",
    "FiniteSourceCase",
    "SimulationOutcome",
    "SpatialGrid",
    "ThermalSchedule",
    "ThermalSegment",
    "run_finite_source_case",
    "solve_diffusion",
]
