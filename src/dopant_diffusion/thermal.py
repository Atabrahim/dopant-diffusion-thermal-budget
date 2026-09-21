"""Piecewise-linear thermal schedules and integrated thermal budget."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
from scipy.integrate import quad


class DiffusivityModel(Protocol):
    """Protocol for an object that evaluates diffusivity in m²/s."""

    def diffusivity_m2_per_s(self, temperature_K: float | np.ndarray) -> np.ndarray: ...


@dataclass(frozen=True)
class ThermalSegment:
    """One linear temperature ramp or isothermal hold.

    Temperatures are in K and duration is in s.
    """

    start_temperature_K: float
    end_temperature_K: float
    duration_s: float

    def __post_init__(self) -> None:
        values = (self.start_temperature_K, self.end_temperature_K, self.duration_s)
        if not all(np.isfinite(value) for value in values):
            raise ValueError("Thermal-segment values must be finite.")
        if self.start_temperature_K <= 0.0 or self.end_temperature_K <= 0.0:
            raise ValueError("Temperatures must be greater than zero kelvin.")
        if self.duration_s <= 0.0:
            raise ValueError("duration_s must be positive.")

    def temperature_at_fraction(self, fraction: float) -> float:
        """Return linearly interpolated temperature for fraction in [0, 1]."""

        if not 0.0 <= fraction <= 1.0:
            raise ValueError("fraction must lie in [0, 1].")
        return self.start_temperature_K + fraction * (
            self.end_temperature_K - self.start_temperature_K
        )


@dataclass(frozen=True)
class ThermalSchedule:
    """Continuous sequence of linear ramps and holds."""

    segments: tuple[ThermalSegment, ...]

    def __post_init__(self) -> None:
        if not self.segments:
            raise ValueError("A thermal schedule requires at least one segment.")
        for previous, current in zip(self.segments[:-1], self.segments[1:], strict=True):
            if not np.isclose(previous.end_temperature_K, current.start_temperature_K):
                raise ValueError(
                    "Adjacent thermal segments must have matching endpoint temperatures."
                )

    @property
    def duration_s(self) -> float:
        """Total schedule duration in s."""

        return float(sum(segment.duration_s for segment in self.segments))

    def time_temperature_arrays(
        self, points_per_segment: int = 100
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return sampled time (s) and temperature (K) arrays for plotting."""

        if points_per_segment < 2:
            raise ValueError("points_per_segment must be at least 2.")
        times: list[np.ndarray] = []
        temperatures: list[np.ndarray] = []
        offset_s = 0.0
        for index, segment in enumerate(self.segments):
            local_time = np.linspace(0.0, segment.duration_s, points_per_segment)
            if index:
                local_time = local_time[1:]
            fraction = local_time / segment.duration_s
            times.append(offset_s + local_time)
            temperatures.append(
                segment.start_temperature_K
                + fraction * (segment.end_temperature_K - segment.start_temperature_K)
            )
            offset_s += segment.duration_s
        return np.concatenate(times), np.concatenate(temperatures)

    def segment_thermal_budgets_m2(
        self, diffusivity: DiffusivityModel, substeps_per_segment: int
    ) -> list[float]:
        """Integrate ``∫D[T(t)]dt`` for equal-time substeps of each schedule segment.

        The returned diffusion-time increments have units m². Numerical quadrature
        prevents a high-temperature ramp from being represented by an arbitrary
        single endpoint temperature.
        """

        if substeps_per_segment < 1:
            raise ValueError("substeps_per_segment must be at least 1.")
        increments: list[float] = []
        for segment in self.segments:
            local_edges = np.linspace(0.0, segment.duration_s, substeps_per_segment + 1)
            for left_s, right_s in zip(local_edges[:-1], local_edges[1:], strict=True):

                def integrand(local_time_s: float) -> float:
                    fraction = local_time_s / segment.duration_s
                    temperature_K = segment.temperature_at_fraction(fraction)
                    return float(diffusivity.diffusivity_m2_per_s(temperature_K))

                budget, _ = quad(integrand, float(left_s), float(right_s), epsabs=1e-20)
                increments.append(float(budget))
        return increments

    def thermal_budget_m2(self, diffusivity: DiffusivityModel) -> float:
        """Return the integrated diffusivity ``∫ D[T(t)]dt`` in m²."""

        return float(sum(self.segment_thermal_budgets_m2(diffusivity, substeps_per_segment=1)))
