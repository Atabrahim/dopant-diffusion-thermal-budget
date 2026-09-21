"""Tests for process diagnostics."""

from __future__ import annotations

import numpy as np
import pytest

from dopant_diffusion.diagnostics import concentration_crossings_m


def test_threshold_crossings_report_both_sides_of_a_finite_source_profile() -> None:
    depth_m = np.array([0.0, 1.0, 2.0, 3.0]) * 1.0e-6
    concentration_m3 = np.array([0.0, 2.0, 2.0, 0.0]) * 1.0e20

    crossings = concentration_crossings_m(depth_m, concentration_m3, 1.0e20)

    assert crossings == pytest.approx(np.array([0.5e-6, 2.5e-6]))


def test_threshold_crossings_validate_inputs() -> None:
    with pytest.raises(ValueError, match="strictly increasing"):
        concentration_crossings_m([0.0, 1.0, 1.0], [1.0, 2.0, 3.0], 2.0)
    with pytest.raises(ValueError, match="positive"):
        concentration_crossings_m([0.0, 1.0], [1.0, 2.0], 0.0)
