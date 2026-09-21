"""Tests for material parameters and thermal-history calculations."""

from __future__ import annotations

import numpy as np
import pytest

from dopant_diffusion.materials import BoronInSilicon
from dopant_diffusion.thermal import ThermalSchedule, ThermalSegment


def test_arrhenius_diffusivity_is_positive_and_increases_with_temperature() -> None:
    material = BoronInSilicon()
    diffusivity = material.diffusivity_m2_per_s(np.array([1_073.15, 1_273.15, 1_473.15]))

    assert np.all(diffusivity > 0.0)
    assert np.all(np.diff(diffusivity) > 0.0)


def test_invalid_material_parameters_are_rejected() -> None:
    with pytest.raises(ValueError, match="positive"):
        BoronInSilicon(prefactor_m2_per_s=0.0)
    with pytest.raises(ValueError, match="greater than zero"):
        BoronInSilicon().diffusivity_m2_per_s(0.0)


def test_isothermal_thermal_budget_equals_diffusivity_times_duration() -> None:
    material = BoronInSilicon()
    schedule = ThermalSchedule((ThermalSegment(1_273.15, 1_273.15, 180.0),))

    expected = float(material.diffusivity_m2_per_s(1_273.15)) * 180.0

    assert schedule.thermal_budget_m2(material) == pytest.approx(expected, rel=1e-12)


def test_ramp_thermal_budget_is_bounded_by_endpoint_diffusivities() -> None:
    material = BoronInSilicon()
    duration_s = 90.0
    schedule = ThermalSchedule((ThermalSegment(1_073.15, 1_273.15, duration_s),))
    budget = schedule.thermal_budget_m2(material)
    endpoint_diffusivities = material.diffusivity_m2_per_s(np.array([1_073.15, 1_273.15]))

    assert (
        endpoint_diffusivities.min() * duration_s
        < budget
        < endpoint_diffusivities.max() * duration_s
    )


def test_schedule_sampling_preserves_endpoints_and_continuity() -> None:
    schedule = ThermalSchedule(
        (
            ThermalSegment(298.15, 1_273.15, 90.0),
            ThermalSegment(1_273.15, 1_273.15, 600.0),
        )
    )

    time_s, temperature_K = schedule.time_temperature_arrays(points_per_segment=3)

    assert time_s[0] == pytest.approx(0.0)
    assert time_s[-1] == pytest.approx(690.0)
    assert temperature_K[0] == pytest.approx(298.15)
    assert temperature_K[-1] == pytest.approx(1_273.15)
    assert np.all(np.diff(time_s) > 0.0)


def test_discontinuous_schedule_is_rejected() -> None:
    with pytest.raises(ValueError, match="matching endpoint"):
        ThermalSchedule(
            (
                ThermalSegment(298.15, 1_273.15, 90.0),
                ThermalSegment(1_200.0, 1_200.0, 60.0),
            )
        )
