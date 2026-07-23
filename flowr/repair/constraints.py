"""Hard constraints for coordinate-only ligand repair."""

from typing import Any

import torch


DISCRETE_STATE_KEYS = (
    "atomics",
    "charges",
    "hybridization",
    "is_aromatic",
    "bonds",
    "interactions",
)


def _copy_value(value: Any) -> Any:
    if torch.is_tensor(value):
        return value.clone()
    if isinstance(value, list):
        return value.copy()
    if isinstance(value, dict):
        return value.copy()
    return value


def project_coordinate_only_state(
    state: dict[str, Any],
    reference: dict[str, Any],
    fixed_mask: torch.Tensor,
) -> dict[str, Any]:
    """Restore fixed coordinates and every discrete ligand attribute."""
    if "coords" not in state or "coords" not in reference:
        raise KeyError("state and reference must contain coords")
    if state["coords"].shape != reference["coords"].shape:
        raise ValueError("state and reference coordinates must have the same shape")
    if fixed_mask.shape != reference["coords"].shape[:2]:
        raise ValueError("fixed_mask must have shape [batch_size, n_atoms]")

    projected = {key: _copy_value(value) for key, value in state.items()}
    fixed_mask = fixed_mask.to(device=state["coords"].device, dtype=torch.bool)
    coords = torch.where(
        fixed_mask.unsqueeze(-1), reference["coords"], state["coords"]
    )
    if "mask" in reference:
        coords = coords * reference["mask"].to(coords.dtype).unsqueeze(-1)
    projected["coords"] = coords

    for key in DISCRETE_STATE_KEYS:
        if key in reference and torch.is_tensor(reference[key]):
            projected[key] = reference[key].clone()

    for key in ("mask", "fragment_mask", "fragment_mode"):
        if key in reference:
            projected[key] = _copy_value(reference[key])
    return projected


def states_have_same_discrete_graph(
    state: dict[str, Any], reference: dict[str, Any]
) -> bool:
    """Return whether all available discrete tensors match exactly."""
    for key in DISCRETE_STATE_KEYS:
        if key in reference and torch.is_tensor(reference[key]):
            if key not in state or not torch.equal(state[key], reference[key]):
                return False
    return True


def max_fixed_coordinate_drift(
    coords: torch.Tensor, reference_coords: torch.Tensor, fixed_mask: torch.Tensor
) -> float:
    """Return the maximum Euclidean displacement over fixed atoms."""
    if coords.shape != reference_coords.shape:
        raise ValueError("coords and reference_coords must have the same shape")
    if fixed_mask.shape != coords.shape[:-1]:
        raise ValueError("fixed_mask shape does not match coordinates")
    fixed_mask = fixed_mask.to(device=coords.device, dtype=torch.bool)
    if not fixed_mask.any():
        return 0.0
    drift = torch.linalg.vector_norm(coords - reference_coords, dim=-1)
    return float(drift[fixed_mask].max().item())
