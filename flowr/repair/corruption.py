"""Deterministic coordinate corruptions used by repair experiments."""

import math
from collections.abc import Sequence

import torch


def build_torsion_corruption(
    coords_good: torch.Tensor,
    axis_origin: int,
    axis_target: int,
    editable_atom_indices: Sequence[int],
    angle_degrees: float,
) -> torch.Tensor:
    """Rotate selected atoms around an ordered bond axis using the right-hand rule."""
    if coords_good.ndim != 2 or coords_good.shape[-1] != 3:
        raise ValueError("coords_good must have shape [n_atoms, 3]")
    n_atoms = coords_good.shape[0]
    indices = list(editable_atom_indices)
    if axis_origin == axis_target:
        raise ValueError("axis_origin and axis_target must differ")
    if not 0 <= axis_origin < n_atoms or not 0 <= axis_target < n_atoms:
        raise IndexError("axis atom index is out of range")
    if len(indices) != len(set(indices)):
        raise ValueError("editable_atom_indices must be unique")
    if any(index < 0 or index >= n_atoms for index in indices):
        raise IndexError("editable atom index is out of range")
    if axis_origin in indices or axis_target in indices:
        raise ValueError("axis atoms must remain fixed")

    origin = coords_good[axis_origin]
    axis = coords_good[axis_target] - origin
    axis_norm = torch.linalg.vector_norm(axis)
    if float(axis_norm.item()) < 1e-8:
        raise ValueError("rotation axis has zero length")
    axis = axis / axis_norm

    coords_bad = coords_good.clone()
    if not indices:
        return coords_bad

    selected = torch.tensor(indices, dtype=torch.long, device=coords_good.device)
    vectors = coords_good[selected] - origin
    angle = torch.tensor(
        math.radians(angle_degrees),
        dtype=coords_good.dtype,
        device=coords_good.device,
    )
    cos_angle = torch.cos(angle)
    sin_angle = torch.sin(angle)
    axis_rows = axis.expand_as(vectors)
    rotated = (
        vectors * cos_angle
        + torch.cross(axis_rows, vectors, dim=-1) * sin_angle
        + axis_rows
        * (vectors * axis_rows).sum(dim=-1, keepdim=True)
        * (1.0 - cos_angle)
    )
    coords_bad[selected] = origin + rotated
    return coords_bad
