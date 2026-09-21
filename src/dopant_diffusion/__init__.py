"""Conservative one-dimensional silicon dopant-diffusion simulation."""

from .materials import BoronInSilicon
from .thermal import ThermalSchedule, ThermalSegment

__all__ = ["BoronInSilicon", "ThermalSchedule", "ThermalSegment"]
