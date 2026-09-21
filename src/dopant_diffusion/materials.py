"""Material parameters and thermally activated diffusivity models."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

BOLTZMANN_CONSTANT_EV_PER_K = 8.617_333_262_145e-5
"""Boltzmann constant in electronvolts per kelvin."""


@dataclass(frozen=True)
class BoronInSilicon:
    """Dilute, Fickian boron diffusivity parameterisation for silicon.

    Parameters are deliberately configurable because reported boron diffusivities
    depend on concentration, point-defect conditions and the experimental method.
    The V1 default is a representative intrinsic-diffusion Arrhenius law:
    ``D = 0.76 cm²/s * exp(-3.46 eV / (k_B T))``.

    This is a chemical-diffusion model only. It must not be interpreted as an
    electrical-activation model.
    """

    prefactor_m2_per_s: float = 0.76e-4
    activation_energy_eV: float = 3.46

    def __post_init__(self) -> None:
        if not np.isfinite(self.prefactor_m2_per_s) or self.prefactor_m2_per_s <= 0.0:
            raise ValueError("prefactor_m2_per_s must be finite and positive.")
        if not np.isfinite(self.activation_energy_eV) or self.activation_energy_eV <= 0.0:
            raise ValueError("activation_energy_eV must be finite and positive.")

    def diffusivity_m2_per_s(self, temperature_K: ArrayLike) -> NDArray[np.float64]:
        """Return diffusivity in m²/s at absolute temperature ``temperature_K``."""

        temperature = np.asarray(temperature_K, dtype=float)
        if np.any(~np.isfinite(temperature)) or np.any(temperature <= 0.0):
            raise ValueError("temperature_K must contain finite values greater than zero kelvin.")
        exponent = -self.activation_energy_eV / (BOLTZMANN_CONSTANT_EV_PER_K * temperature)
        return np.asarray(self.prefactor_m2_per_s * np.exp(exponent), dtype=float)
