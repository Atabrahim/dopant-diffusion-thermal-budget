"""Publication-oriented plots for the documented process workflow."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np

from .thermal import ThermalSchedule
from .workflow import SimulationOutcome


def _pyplot():
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:  # pragma: no cover - exercised by installation guidance
        raise RuntimeError("Plotting requires matplotlib. Install with '.[plots]'.") from error
    return plt


def plot_temperature_schedules(schedules: dict[str, ThermalSchedule], output_path: Path) -> None:
    """Plot temperature schedules against time in minutes."""

    plt = _pyplot()
    figure, axis = plt.subplots(figsize=(7.4, 4.2), constrained_layout=True)
    for label, schedule in schedules.items():
        time_s, temperature_K = schedule.time_temperature_arrays()
        axis.plot(time_s / 60.0, temperature_K - 273.15, linewidth=2.0, label=label)
    axis.set_xlabel("Time (min)")
    axis.set_ylabel("Temperature (°C)")
    axis.set_title("Ramp-and-hold thermal schedules")
    axis.legend(frameon=False)
    axis.grid(alpha=0.25)
    _save(figure, output_path)


def plot_concentration_profiles(outcomes: Iterable[SimulationOutcome], output_path: Path) -> None:
    """Plot initial/final chemical-concentration profiles with explicit units."""

    outcomes = list(outcomes)
    if not outcomes:
        raise ValueError("At least one outcome is required for a profile plot.")
    plt = _pyplot()
    figure, axis = plt.subplots(figsize=(7.4, 4.6), constrained_layout=True)
    first = outcomes[0]
    axis.semilogy(
        first.case.grid.centres_m * 1.0e9,
        first.initial_concentration_m3 / 1.0e6,
        color="0.35",
        linestyle="--",
        linewidth=1.5,
        label="Initial reflected-Gaussian chemical profile",
    )
    for outcome in outcomes:
        axis.semilogy(
            outcome.case.grid.centres_m * 1.0e9,
            outcome.result.concentration_m3 / 1.0e6,
            linewidth=2.0,
            label=outcome.case.name,
        )
    threshold_cm3 = first.case.threshold_m3 / 1.0e6
    plotted_maximum_cm3 = max(
        float(first.initial_concentration_m3.max() / 1.0e6),
        *(float(outcome.result.concentration_m3.max() / 1.0e6) for outcome in outcomes),
        threshold_cm3,
    )
    axis.axhline(
        threshold_cm3, color="0.25", linewidth=1.0, linestyle=":", label="Reported threshold"
    )
    axis.set_xlabel("Depth below surface (nm)")
    axis.set_ylabel("Chemical boron concentration (cm$^{-3}$)")
    axis.set_title("Finite-source boron diffusion in silicon")
    # The reflected-Gaussian far-tail can be many decades below a useful plot.
    # Base the visible range on reported quantities, not the numerical tail.
    display_floor_cm3 = 1.0e13
    axis.set_ylim(display_floor_cm3, 2.0 * plotted_maximum_cm3)
    visible_depth_m = 0.0
    initial_visible = first.initial_concentration_m3 >= display_floor_cm3 * 1.0e6
    if np.any(initial_visible):
        visible_depth_m = float(first.case.grid.centres_m[initial_visible][-1])
    for outcome in outcomes:
        visible = outcome.result.concentration_m3 >= display_floor_cm3 * 1.0e6
        if np.any(visible):
            visible_depth_m = max(visible_depth_m, float(outcome.case.grid.centres_m[visible][-1]))
    if visible_depth_m:
        axis.set_xlim(
            0.0,
            min(max(case.case.grid.length_m for case in outcomes), 1.15 * visible_depth_m) * 1.0e9,
        )
    axis.legend(frameon=False, fontsize=8)
    axis.grid(alpha=0.2, which="both")
    _save(figure, output_path)


def plot_sensitivity(rows: list[dict[str, float | str]], output_path: Path) -> None:
    """Plot declared input sensitivity of the deepest threshold crossing."""

    plt = _pyplot()
    groups = {
        parameter: [row for row in rows if row["parameter"] == parameter]
        for parameter in {str(row["parameter"]) for row in rows}
    }
    figure, axes = plt.subplots(1, len(groups), figsize=(7.4, 3.8), constrained_layout=True)
    if len(groups) == 1:
        axes = [axes]
    for axis, (parameter, group) in zip(axes, sorted(groups.items()), strict=True):
        sorted_group = sorted(group, key=lambda row: float(row["value"]))
        x = [float(row["value"]) for row in sorted_group]
        y = [float(row["deepest_threshold_crossing_nm"]) for row in sorted_group]
        axis.plot(x, y, marker="o", linewidth=2.0)
        axis.set_xlabel("Multiplier" if parameter == "diffusivity multiplier" else "Offset (K)")
        axis.set_ylabel("Deepest threshold crossing (nm)")
        axis.set_title(parameter.capitalize())
        axis.grid(alpha=0.25)
    _save(figure, output_path)


def _save(figure, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=180, bbox_inches="tight")
    figure.clf()
