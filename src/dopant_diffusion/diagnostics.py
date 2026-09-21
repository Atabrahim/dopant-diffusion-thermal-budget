"""Physically interpretable diffusion-process diagnostics."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def concentration_crossings_m(
    depth_m: ArrayLike, concentration_m3: ArrayLike, threshold_m3: float
) -> NDArray[np.float64]:
    """Return linearly interpolated depths where concentration crosses a threshold.

    A finite-source profile can have two crossings, whereas a monotonic
    constant-source profile normally has one. The function therefore reports all
    crossings instead of silently selecting a device junction definition.
    """

    depth = np.asarray(depth_m, dtype=float)
    concentration = np.asarray(concentration_m3, dtype=float)
    if depth.ndim != 1 or concentration.shape != depth.shape or len(depth) < 2:
        raise ValueError("depth_m and concentration_m3 must be matching one-dimensional arrays.")
    if np.any(~np.isfinite(depth)) or np.any(~np.isfinite(concentration)):
        raise ValueError("Depth and concentration values must be finite.")
    if np.any(np.diff(depth) <= 0.0):
        raise ValueError("depth_m must be strictly increasing.")
    if not np.isfinite(threshold_m3) or threshold_m3 <= 0.0:
        raise ValueError("threshold_m3 must be finite and positive.")

    shifted = concentration - threshold_m3
    crossings: list[float] = []
    exact = np.flatnonzero(shifted == 0.0)
    crossings.extend(float(depth[index]) for index in exact)
    sign_change_indices = np.flatnonzero(shifted[:-1] * shifted[1:] < 0.0)
    for index in sign_change_indices:
        left_fraction = -shifted[index] / (shifted[index + 1] - shifted[index])
        crossings.append(float(depth[index] + left_fraction * (depth[index + 1] - depth[index])))
    return np.asarray(sorted(set(crossings)), dtype=float)
