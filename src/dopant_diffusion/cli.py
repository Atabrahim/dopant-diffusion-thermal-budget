"""Command-line interface for reproducible finite-source process simulations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from .materials import BoronInSilicon
from .plotting import plot_concentration_profiles, plot_temperature_schedules
from .solver import SpatialGrid
from .thermal import ThermalSchedule, ThermalSegment
from .workflow import FiniteSourceCase, run_finite_source_case, write_outcome


def build_parser() -> argparse.ArgumentParser:
    """Build the public command-line parser."""

    parser = argparse.ArgumentParser(
        prog="dopant-diffusion",
        description="Run a reproducible one-dimensional finite-source boron diffusion case.",
    )
    parser.add_argument(
        "config", type=Path, help="Path to a JSON finite-source process configuration."
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True, help="Directory for CSV and JSON outputs."
    )
    parser.add_argument(
        "--plots", action="store_true", help="Also save profile and thermal-schedule PNG figures."
    )
    return parser


def run_from_config(config: dict[str, Any]):
    """Create a material model and process case from the documented JSON schema."""

    try:
        material_config = config["material"]
        grid_config = config["grid"]
        profile_config = config["finite_source"]
        schedule_config = config["thermal_schedule"]
        material = BoronInSilicon(**material_config)
        schedule = ThermalSchedule(tuple(ThermalSegment(**segment) for segment in schedule_config))
        case = FiniteSourceCase(
            name=str(config["case_name"]),
            grid=SpatialGrid(**grid_config),
            schedule=schedule,
            threshold_m3=float(config["threshold_m3"]),
            substeps_per_segment=int(config.get("substeps_per_segment", 200)),
            **profile_config,
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"Invalid process configuration: {error}") from error
    return material, case


def main(arguments: list[str] | None = None) -> int:
    """Run the CLI and return a conventional process exit code."""

    parser = build_parser()
    namespace = parser.parse_args(arguments)
    try:
        config = json.loads(namespace.config.read_text(encoding="utf-8"))
        material, case = run_from_config(config)
        outcome = run_finite_source_case(case, material)
        report_path, profile_path = write_outcome(namespace.output_dir, outcome)
        if namespace.plots:
            plot_temperature_schedules(
                {case.name: case.schedule}, namespace.output_dir / "schedule.png"
            )
            plot_concentration_profiles([outcome], namespace.output_dir / "profile.png")
    except (OSError, ValueError, json.JSONDecodeError, RuntimeError) as error:
        parser.error(str(error))
    print(f"Wrote {report_path}")
    print(f"Wrote {profile_path}")
    return 0
