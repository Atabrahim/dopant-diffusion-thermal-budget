"""Regenerate the published V1 process-comparison reports and figures."""

from __future__ import annotations

import json
from pathlib import Path

from dopant_diffusion.materials import BoronInSilicon
from dopant_diffusion.plotting import (
    plot_concentration_profiles,
    plot_sensitivity,
    plot_temperature_schedules,
)
from dopant_diffusion.solver import SpatialGrid
from dopant_diffusion.thermal import ThermalSchedule, ThermalSegment
from dopant_diffusion.workflow import (
    FiniteSourceCase,
    run_finite_source_case,
    sensitivity_rows,
    write_outcome,
    write_sensitivity_report,
)


def main() -> None:
    """Generate the tracked V1 case study from declared model inputs."""

    repository_root = Path(__file__).resolve().parents[1]
    figure_dir = repository_root / "figures"
    report_dir = repository_root / "reports"
    material = BoronInSilicon()
    grid = SpatialGrid(length_m=2.0e-6, cells=800)
    baseline_schedule = ThermalSchedule(
        (
            ThermalSegment(298.15, 1_273.15, 90.0),
            ThermalSegment(1_273.15, 1_273.15, 1_200.0),
            ThermalSegment(1_273.15, 298.15, 90.0),
        )
    )
    cooler_schedule = ThermalSchedule(
        (
            ThermalSegment(298.15, 1_223.15, 90.0),
            ThermalSegment(1_223.15, 1_223.15, 1_200.0),
            ThermalSegment(1_223.15, 298.15, 90.0),
        )
    )
    shared = {
        "grid": grid,
        "dose_m2": 1.0e17,
        "projected_range_m": 2.0e-7,
        "standard_deviation_m": 5.0e-8,
        "threshold_m3": 1.0e23,
        "substeps_per_segment": 200,
    }
    baseline = FiniteSourceCase(
        name="1000 °C, 20 min drive-in", schedule=baseline_schedule, **shared
    )
    cooler = FiniteSourceCase(name="950 °C, 20 min drive-in", schedule=cooler_schedule, **shared)
    outcomes = [
        run_finite_source_case(baseline, material),
        run_finite_source_case(cooler, material),
    ]
    for outcome in outcomes:
        write_outcome(report_dir, outcome)
    rows = sensitivity_rows(baseline, material)
    write_sensitivity_report(report_dir / "sensitivity.csv", rows)
    comparison = {
        "kind": "synthetic model comparison; not measured process data",
        "baseline": outcomes[0].report(),
        "cooler_schedule": outcomes[1].report(),
        "sensitivity_definition": (
            "declared input sensitivity, not a probabilistic uncertainty interval"
        ),
        "sensitivity_rows": rows,
    }
    (report_dir / "process_comparison.json").write_text(
        json.dumps(comparison, indent=2) + "\n", encoding="utf-8"
    )
    plot_temperature_schedules(
        {baseline.name: baseline_schedule, cooler.name: cooler_schedule},
        figure_dir / "thermal_schedules.png",
    )
    plot_concentration_profiles(outcomes, figure_dir / "boron_profiles.png")
    plot_sensitivity(rows, figure_dir / "sensitivity.png")


if __name__ == "__main__":
    main()
