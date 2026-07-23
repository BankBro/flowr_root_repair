"""Deterministic coordinate corruptions used by repair experiments."""

import math
from collections import deque
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import Lipinski


DEFAULT_TORSION_ANGLES = (
    30.0,
    -30.0,
    60.0,
    -60.0,
    90.0,
    -90.0,
    120.0,
    -120.0,
    150.0,
    -150.0,
    180.0,
)


@dataclass(frozen=True)
class RotatableBranch:
    """A deterministic smaller branch around one rotatable bond."""

    axis_origin: int
    axis_target: int
    editable_atom_indices: tuple[int, ...]


@dataclass(frozen=True)
class TorsionCandidate:
    """Coordinates and provenance for one enumerated torsion corruption."""

    axis_origin: int
    axis_target: int
    editable_atom_indices: tuple[int, ...]
    angle_degrees: float
    coords: np.ndarray


def _component_without_bond(
    mol: Chem.Mol, start: int, blocked: frozenset[int]
) -> set[int]:
    visited = {start}
    queue = deque([start])
    while queue:
        atom_index = queue.popleft()
        for neighbor in mol.GetAtomWithIdx(atom_index).GetNeighbors():
            neighbor_index = neighbor.GetIdx()
            if frozenset((atom_index, neighbor_index)) == blocked:
                continue
            if neighbor_index not in visited:
                visited.add(neighbor_index)
                queue.append(neighbor_index)
    return visited


def smaller_rotatable_branches(mol: Chem.Mol) -> list[RotatableBranch]:
    """Return the smaller heavy-atom branch for every RDKit rotatable bond."""
    if mol.GetNumConformers() == 0:
        raise ValueError("mol must contain a conformer")
    matches = mol.GetSubstructMatches(Lipinski.RotatableBondSmarts, uniquify=True)
    branches: list[RotatableBranch] = []
    for left, right in matches:
        blocked = frozenset((left, right))
        options = []
        for origin, target in ((left, right), (right, left)):
            component = _component_without_bond(mol, target, blocked)
            editable = tuple(sorted(component - {target}))
            if editable:
                options.append((len(editable), editable, origin, target))
        if not options:
            continue
        _, editable, origin, target = min(options)
        branches.append(
            RotatableBranch(
                axis_origin=origin,
                axis_target=target,
                editable_atom_indices=editable,
            )
        )
    return sorted(
        branches,
        key=lambda branch: (
            min(branch.axis_origin, branch.axis_target),
            max(branch.axis_origin, branch.axis_target),
            branch.axis_origin,
            branch.editable_atom_indices,
        ),
    )


def enumerate_torsion_candidates(
    mol: Chem.Mol,
    angles: Sequence[float] = DEFAULT_TORSION_ANGLES,
) -> list[TorsionCandidate]:
    """Enumerate deterministic rigid torsions of each smaller rotatable branch."""
    coords = np.asarray(mol.GetConformer().GetPositions(), dtype=np.float64)
    coords_tensor = torch.from_numpy(coords)
    candidates = []
    for branch in smaller_rotatable_branches(mol):
        for angle in angles:
            corrupted = build_torsion_corruption(
                coords_tensor,
                axis_origin=branch.axis_origin,
                axis_target=branch.axis_target,
                editable_atom_indices=branch.editable_atom_indices,
                angle_degrees=float(angle),
            )
            candidates.append(
                TorsionCandidate(
                    axis_origin=branch.axis_origin,
                    axis_target=branch.axis_target,
                    editable_atom_indices=branch.editable_atom_indices,
                    angle_degrees=float(angle),
                    coords=corrupted.numpy(),
                )
            )
    return candidates


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
