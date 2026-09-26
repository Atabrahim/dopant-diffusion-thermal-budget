"""Reproducible finite-source process cases, reports and sensitivity studies."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .diagnostics import concentration_crossings_m
from .profiles import reflected_gaussian_m3
from .solver import DiffusionResult, SpatialGrid, solve_diffusion
from .thermal import DiffusivityModel, ThermalSchedule, ThermalSegment


@dataclass(frozen=True)
class FiniteSourceCase:
    """Inputs for a finite-source chemical-diffusion simulation.

    All concentration and geometry values are SI. The initial profile is a
    reflected Gaussian, which represents a simplified, dose-limited profile; it
    is not an ion-implantation damage model.
    """

    name: str
    grid: SpatialGrid
    schedule: ThermalSchedule
    dose_m2: float
    projected_range_m: float
    standard_deviation_m: float
    threshold_m3: float
    substeps_per_segment: int = 200

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("A process case requires a non-empty name.")
        parameters = (
            self.dose_m2,
            self.projected_range_m,
            self.standard_deviation_m,
            self.threshold_m3,
        )
        if not all(np.isfinite(value) for value in parameters):
            raise ValueError("Finite-source case parameters must be finite.")
        if self.dose_m2 <= 0.0 or self.standard_deviation_m <= 0.0 or self.threshold_m3 <= 0.0:
            raise ValueError("Dose, standard deviation and threshold must be positive.")
        if self.projected_range_m < 0.0:
            raise ValueError("projected_range_m must be non-negative.")
        if self.substeps_per_segment < 1:
            raise ValueError("substeps_per_segment must be at least 1.")


@dataclass(frozen=True)
class SimulationOutcome:
    """Numerical result plus process diagnostics for one finite-source case."""

    case: FiniteSourceCase
    initial_concentration_m3: np.ndarray
    result: DiffusionResult
    crossing_depths_m: np.ndarray

    @property
    def max_concentration_m3(self) -> float:
        """Final maximum chemical concentration in m⁻³."""

        return float(np.max(self.result.concentration_m3))

    def report(self) -> dict[str, Any]:
        """Return a JSON-serialisable, unit-labelled process report."""

        return {
            "case": self.case.name,
            "model": {
                "quantity": "chemical boron concentration",
                "initial_profile": "reflected Gaussian",
                "transport": "one-dimensional Fickian diffusion with spatially uniform D(T)",
            },
            "units": {
                "depth": "nm",
                "concentration": "cm^-3",
                "dose": "cm^-2",
                "thermal_budget": "m^2",
            },
            "thermal_schedule": [asdict(segment) for segment in self.case.schedule.segments],
            "grid": {"length_m": self.case.grid.length_m, "cells": self.case.grid.cells},
            "input": {
                "dose_cm2": self.case.dose_m2 / 1.0e4,
                "projected_range_nm": self.case.projected_range_m * 1.0e9,
                "standard_deviation_nm": self.case.standard_deviation_m * 1.0e9,
                "threshold_cm3": self.case.threshold_m3 / 1.0e6,
            },
            "output": {
                "thermal_budget_m2": self.result.thermal_budget_m2,
                "initial_dose_cm2": self.result.initial_dose_m2 / 1.0e4,
                "final_dose_cm2": self.result.final_dose_m2 / 1.0e4,
                "dose_balance_error_cm2": self.result.dose_balance_error_m2 / 1.0e4,
                "max_concentration_cm3": self.max_concentration_m3 / 1.0e6,
                "threshold_crossings_nm": (self.crossing_depths_m * 1.0e9).tolist(),
                "implicit_steps": self.result.time_steps,
            },
        }


@dataclass(frozen=True)
class ScaledDiffusivity:
    """Wrap a diffusivity model to make a declared parameter-sensitivity case."""

    model: DiffusivityModel
    multiplier: float

    def __post_init__(self) -> None:
        if not np.isfinite(self.multiplier) or self.multiplier <= 0.0:
            raise ValueError("multiplier must be finite and positive.")

    def diffusivity_m2_per_s(self, temperature_K: float | np.ndarray) -> np.ndarray:
        """Return scaled diffusivity in m²/s."""

        return self.multiplier * np.asarray(self.model.diffusivity_m2_per_s(temperature_K))


def run_finite_source_case(
    case: FiniteSourceCase, diffusivity: DiffusivityModel
) -> SimulationOutcome:
    """Run one finite-source process case and calculate transparent diagnostics."""

    initial = reflected_gaussian_m3(
        case.grid.centres_m,
        case.dose_m2,
        case.projected_range_m,
        case.standard_deviation_m,
    )
    result = solve_diffusion(
        initial,
        case.grid,
        case.schedule,
        diffusivity,
        substeps_per_segment=case.substeps_per_segment,
    )
    crossings = concentration_crossings_m(
        case.grid.centres_m, result.concentration_m3, case.threshold_m3
    )
    return SimulationOutcome(case, initial, result, crossings)


def schedule_with_temperature_offset(
    schedule: ThermalSchedule, temperature_offset_K: float
) -> ThermalSchedule:
    """Return a schedule shifted by a declared temperature assumption in K."""

    if not np.isfinite(temperature_offset_K):
        raise ValueError("temperature_offset_K must be finite.")
    return ThermalSchedule(
        tuple(
            ThermalSegment(
                segment.start_temperature_K + temperature_offset_K,
                segment.end_temperature_K + temperature_offset_K,
                segment.duration_s,
            )
            for segment in schedule.segments
        )
    )


def sensitivity_rows(
    case: FiniteSourceCase,
    diffusivity: DiffusivityModel,
    *,
    diffusivity_multipliers: tuple[float, ...] = (0.5, 1.0, 2.0),
    temperature_offsets_K: tuple[float, ...] = (-10.0, 0.0, 10.0),
) -> list[dict[str, float | str]]:
    """Quantify declared model-input sensitivity without calling it uncertainty.

    Each row changes either an Arrhenius-diffusivity multiplier or every schedule
    temperature. The reported deepest crossing is an explicitly chosen summary;
    the full set of crossings is retained in the generated report.
    """

    rows: list[dict[str, float | str]] = []
    for multiplier in diffusivity_multipliers:
        outcome = run_finite_source_case(case, ScaledDiffusivity(diffusivity, multiplier))
        rows.append(_sensitivity_row("diffusivity multiplier", multiplier, outcome))
    for offset_K in temperature_offsets_K:
        shifted_case = FiniteSourceCase(
            name=f"{case.name}; temperature offset {offset_K:+.0f} K",
            grid=case.grid,
            schedule=schedule_with_temperature_offset(case.schedule, offset_K),
            dose_m2=case.dose_m2,
            projected_range_m=case.projected_range_m,
            standard_deviation_m=case.standard_deviation_m,
            threshold_m3=case.threshold_m3,
            substeps_per_segment=case.substeps_per_segment,
        )
        rows.append(
            _sensitivity_row(
                "temperature offset", offset_K, run_finite_source_case(shifted_case, diffusivity)
            )
        )
    return rows


def _sensitivity_row(
    parameter: str, value: float, outcome: SimulationOutcome
) -> dict[str, float | str]:
    crossings_nm = outcome.crossing_depths_m * 1.0e9
    deepest_crossing_nm = float(crossings_nm[-1]) if len(crossings_nm) else float("nan")
    return {
        "parameter": parameter,
        "value": value,
        "thermal_budget_m2": outcome.result.thermal_budget_m2,
        "max_concentration_cm3": outcome.max_concentration_m3 / 1.0e6,
        "crossing_count": float(len(crossings_nm)),
        "deepest_threshold_crossing_nm": deepest_crossing_nm,
    }


def write_outcome(output_dir: Path, outcome: SimulationOutcome) -> tuple[Path, Path]:
    """Write a unit-labelled JSON report and CSV profile for one outcome."""

    output_dir.mkdir(parents=True, exist_ok=True)
    stem = _safe_stem(outcome.case.name)
    report_path = output_dir / f"{stem}_report.json"
    profile_path = output_dir / f"{stem}_profile.csv"
    report_path.write_text(json.dumps(outcome.report(), indent=2) + "\n", encoding="utf-8")
    with profile_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=("depth_nm", "initial_concentration_cm3", "final_concentration_cm3"),
        )
        writer.writeheader()
        for depth_m, initial_m3, final_m3 in zip(
            outcome.case.grid.centres_m,
            outcome.initial_concentration_m3,
            outcome.result.concentration_m3,
            strict=True,
        ):
            writer.writerow(
                {
                    "depth_nm": depth_m * 1.0e9,
                    "initial_concentration_cm3": initial_m3 / 1.0e6,
                    "final_concentration_cm3": final_m3 / 1.0e6,
                }
            )
    return report_path, profile_path


def write_sensitivity_report(output_path: Path, rows: list[dict[str, float | str]]) -> None:
    """Write sensitivity rows to a stable CSV schema."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = (
        "parameter",
        "value",
        "thermal_budget_m2",
        "max_concentration_cm3",
        "crossing_count",
        "deepest_threshold_crossing_nm",
    )
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _safe_stem(name: str) -> str:
    """Make a predictable local report filename without changing the case name."""

    return "".join(character.lower() if character.isalnum() else "_" for character in name).strip(
        "_"
    )
