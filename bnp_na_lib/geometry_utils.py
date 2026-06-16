"""Shared geometry helpers for bnp_na.

This module intentionally has no DSSR or Phenix dependency.
"""
from __future__ import annotations

import math
from typing import Iterable

import numpy as np


class GeometryError(Exception):
    """Raised when a geometric operation cannot be completed safely."""


def as_vector(values: Iterable[float], name: str = "vector") -> np.ndarray:
    """Return a finite 3-vector as a NumPy array."""
    arr = np.asarray(list(values), dtype=float)
    if arr.shape != (3,) or not np.all(np.isfinite(arr)):
        raise GeometryError(f"{name} must be a finite 3-vector.")
    return arr


def normalize_vector(values: Iterable[float], name: str = "vector") -> np.ndarray:
    """Return a normalized finite 3-vector."""
    arr = as_vector(values, name=name)
    norm = float(np.linalg.norm(arr))
    if norm <= 1e-12:
        raise GeometryError(f"{name} has near-zero length.")
    return arr / norm


def axis_angle_to_matrix(axis: Iterable[float], theta_rad: float) -> np.ndarray:
    """Return a 3x3 rotation matrix from axis-angle input."""
    axis_v = normalize_vector(axis, name="rotation axis")
    theta = float(theta_rad)
    a = math.cos(theta / 2.0)
    b, c, d = -axis_v * math.sin(theta / 2.0)
    return np.array(
        [
            [a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
            [2 * (b * c + a * d), a * a + c * c - b * b - d * d, 2 * (c * d - a * b)],
            [2 * (b * d - a * c), 2 * (c * d + a * b), a * a + d * d - b * b - c * c],
        ],
        dtype=float,
    )


def rotation_matrix_from_to(v_from: Iterable[float], v_to: Iterable[float]) -> np.ndarray:
    """Return the shortest-arc rotation matrix that maps v_from onto v_to."""
    a = normalize_vector(v_from, name="source vector")
    b = normalize_vector(v_to, name="target vector")
    cross = np.cross(a, b)
    dot = float(np.dot(a, b))

    if np.isclose(dot, 1.0):
        return np.eye(3, dtype=float)

    if np.isclose(dot, -1.0):
        helper = np.array([1.0, 0.0, 0.0], dtype=float)
        if np.allclose(np.abs(a), helper):
            helper = np.array([0.0, 1.0, 0.0], dtype=float)
        axis = helper - float(np.dot(helper, a)) * a
        axis = normalize_vector(axis, name="180-degree rotation axis")
        return axis_angle_to_matrix(axis, math.pi)

    cross_norm = float(np.linalg.norm(cross))
    if cross_norm <= 1e-12:
        return np.eye(3, dtype=float)

    k = np.array(
        [
            [0.0, -cross[2], cross[1]],
            [cross[2], 0.0, -cross[0]],
            [-cross[1], cross[0], 0.0],
        ],
        dtype=float,
    )
    return np.eye(3, dtype=float) + k + k @ k * ((1.0 - dot) / (cross_norm ** 2))


def rotation_matrix_z(deg: float) -> np.ndarray:
    """Return a right-handed rotation matrix about +Z."""
    angle = math.radians(float(deg))
    c = math.cos(angle)
    s = math.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=float)


def rotation_matrix_y(deg: float) -> np.ndarray:
    """Return a right-handed rotation matrix about +Y."""
    angle = math.radians(float(deg))
    c = math.cos(angle)
    s = math.sin(angle)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]], dtype=float)
