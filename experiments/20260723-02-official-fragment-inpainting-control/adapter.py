"""Exact-mask adapter for the official FLOWR.ROOT fragment inpainting path."""

from typing import Any

import numpy as np
import torch


def fixed_first_index_map(fixed_mask: np.ndarray | torch.Tensor) -> np.ndarray:
    """Map official fixed-first output slots back to original ligand indices."""
    mask = np.asarray(torch.as_tensor(fixed_mask).cpu(), dtype=bool)
    if mask.ndim != 1:
        raise ValueError("fixed_mask must be one dimensional")
    return np.concatenate([np.flatnonzero(mask), np.flatnonzero(~mask)])


def build_exact_fragment_prior(
    interpolant: Any,
    target_ligand: Any,
    rdkit_mol: Any,
    fixed_mask: np.ndarray | torch.Tensor,
) -> tuple[Any, np.ndarray]:
    """Call the official private prior builder with the frozen true mask."""
    mask = torch.as_tensor(
        fixed_mask,
        dtype=torch.bool,
        device=target_ligand.coords.device,
    )
    if mask.ndim != 1 or mask.numel() != target_ligand.seq_length:
        raise ValueError("fixed_mask must match target_ligand length")
    if not mask.any() or mask.all():
        raise ValueError("fragment inpainting requires fixed and editable atoms")

    prior = interpolant._build_inference_prior(
        target_ligand,
        rdkit_mol,
        mode="fragment_inpainting",
        mask=mask,
        is_local=True,
    )
    output_to_reference = fixed_first_index_map(mask)
    n_fixed = int(mask.sum().item())
    expected_prior_mask = torch.zeros_like(prior.fragment_mask, dtype=torch.bool)
    expected_prior_mask[:n_fixed] = True
    if not torch.equal(prior.fragment_mask.bool(), expected_prior_mask):
        raise RuntimeError("official prior did not preserve fixed-first fragment layout")

    fixed_indices = torch.where(mask)[0]
    for field in ("coords", "atomics", "charges", "hybridization"):
        prior_value = getattr(prior, field, None)
        target_value = getattr(target_ligand, field, None)
        if prior_value is None and target_value is None:
            continue
        if prior_value is None or target_value is None:
            raise RuntimeError(f"official prior changed fixed field availability: {field}")
        if not torch.equal(prior_value[:n_fixed], target_value[fixed_indices]):
            raise RuntimeError(f"official prior changed fixed field: {field}")
    return prior, output_to_reference


def restore_quantized_fixed_coordinates(
    coords_output: np.ndarray,
    coords_bad: np.ndarray,
    fixed_mask: np.ndarray,
    output_to_reference: np.ndarray,
    tolerance: float = 1e-5,
) -> tuple[np.ndarray, float, bool]:
    """Measure physical-frame drift and snap only float32-scale fixed drift."""
    coords_output = np.asarray(coords_output, dtype=float)
    coords_bad = np.asarray(coords_bad, dtype=float)
    fixed_mask = np.asarray(fixed_mask, dtype=bool)
    output_to_reference = np.asarray(output_to_reference, dtype=int)
    if coords_output.shape != (output_to_reference.size, 3):
        raise ValueError("coords_output does not match output index map")

    fixed_output = np.flatnonzero(fixed_mask[output_to_reference])
    fixed_reference = output_to_reference[fixed_output]
    raw_drift = (
        float(
            np.linalg.norm(
                coords_output[fixed_output] - coords_bad[fixed_reference], axis=-1
            ).max()
        )
        if fixed_output.size
        else 0.0
    )
    restored = coords_output.copy()
    within_tolerance = bool(np.isfinite(raw_drift) and raw_drift <= tolerance)
    if within_tolerance:
        restored[fixed_output] = coords_bad[fixed_reference]
    return restored, raw_drift, within_tolerance
