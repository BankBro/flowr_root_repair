"""Coordinate priors for local ligand repair."""

from typing import Optional

import torch


def build_local_coordinate_prior(
    coords_bad: torch.Tensor,
    editable_mask: torch.Tensor,
    generator: Optional[torch.Generator] = None,
) -> torch.Tensor:
    """Sample an isotropic prior and align it to the editable-region centroid."""
    if coords_bad.ndim != 2 or coords_bad.shape[-1] != 3:
        raise ValueError("coords_bad must have shape [n_atoms, 3]")
    if editable_mask.shape != coords_bad.shape[:1]:
        raise ValueError("editable_mask must have shape [n_atoms]")

    editable_mask = editable_mask.to(device=coords_bad.device, dtype=torch.bool)
    coords_prior = coords_bad.clone()
    n_editable = int(editable_mask.sum().item())
    if n_editable == 0:
        return coords_prior

    noise = torch.randn(
        (n_editable, 3),
        dtype=coords_bad.dtype,
        device=coords_bad.device,
        generator=generator,
    )
    target_centroid = coords_bad[editable_mask].mean(dim=0, keepdim=True)
    coords_prior[editable_mask] = noise - noise.mean(dim=0, keepdim=True) + target_centroid
    return coords_prior
