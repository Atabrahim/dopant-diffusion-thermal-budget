"""Tests for reproducible process workflows and the command-line interface."""

from __future__ import annotations

import json

import pytest

from dopant_diffusion.cli import main, run_from_config
from dopant_diffusion.materials import BoronInSilicon
from dopant_diffusion.solver import SpatialGrid
from dopant_diffusion.thermal import ThermalSchedule, ThermalSegment
from dopant_diffusion.workflow import FiniteSourceCase, run_finite_source_case, sensitivity_rows


def _case() -> FiniteSourceCase:
    return FiniteSourceCase(
        name="workflow test",
        grid=SpatialGrid(2.0e-6, 240),
        schedule=ThermalSchedule((ThermalSegment(1_273.15, 1_273.15, 60.0),)),
        dose_m2=1.0e17,
        projected_range_m=2.0e-7,
        standard_deviation_m=5.0e-8,
        threshold_m3=1.0e23,
        substeps_per_segment=40,
    )


def test_workflow_reports_all_threshold_crossings_and_dose_balance() -> None:
    outcome = run_finite_source_case(_case(), BoronInSilicon())
    report = outcome.report()

    assert len(outcome.crossing_depths_m) == 2
    assert abs(outcome.result.dose_balance_error_m2) / outcome.result.initial_dose_m2 < 1.0e-12
    assert report["units"]["concentration"] == "cm^-3"
    assert report["model"]["quantity"] == "chemical boron concentration"


def test_sensitivity_rows_are_explicit_about_changed_input() -> None:
    rows = sensitivity_rows(_case(), BoronInSilicon())

    assert len(rows) == 6
    assert {row["parameter"] for row in rows} == {"diffusivity multiplier", "temperature offset"}
    assert all(row["crossing_count"] == 2.0 for row in rows)


def test_cli_writes_reproducible_json_and_csv(tmp_path) -> None:
    configuration = {
        "case_name": "CLI test",
        "material": {"prefactor_m2_per_s": 7.6e-5, "activation_energy_eV": 3.46},
        "grid": {"length_m": 2.0e-6, "cells": 120},
        "finite_source": {
            "dose_m2": 1.0e17,
            "projected_range_m": 2.0e-7,
            "standard_deviation_m": 5.0e-8,
        },
        "threshold_m3": 1.0e23,
        "substeps_per_segment": 20,
        "thermal_schedule": [
            {"start_temperature_K": 1_273.15, "end_temperature_K": 1_273.15, "duration_s": 30.0}
        ],
    }
    config_path = tmp_path / "case.json"
    output_dir = tmp_path / "outputs"
    config_path.write_text(json.dumps(configuration), encoding="utf-8")

    assert main([str(config_path), "--output-dir", str(output_dir)]) == 0

    report = json.loads(next(output_dir.glob("*_report.json")).read_text(encoding="utf-8"))
    assert next(output_dir.glob("*_profile.csv")).is_file()
    assert report["case"] == "CLI test"


def test_invalid_configuration_is_explained() -> None:
    with pytest.raises(ValueError, match="Invalid process configuration"):
        run_from_config({"material": {}})
