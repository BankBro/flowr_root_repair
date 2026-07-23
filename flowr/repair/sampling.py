"""FLOWR.ROOT sampling restricted to ligand coordinates."""

from typing import Any

import torch

from flowr.repair.constraints import project_coordinate_only_state


def build_repair_times(
    batch_size: int,
    device: torch.device | str,
    coord_time: float = 0.0,
    pocket_time: float = 0.0,
) -> list[torch.Tensor]:
    """Build scalar model times with discrete ligand state fixed at time one."""
    return [
        torch.full((batch_size,), coord_time, device=device),
        torch.ones(batch_size, device=device),
        torch.full((batch_size,), pocket_time, device=device),
    ]


def _initial_self_condition(
    model: Any, prior: dict[str, Any], reference: dict[str, Any], fixed_mask: torch.Tensor
) -> dict[str, Any]:
    cond_batch: dict[str, Any] = {
        "coords": prior["coords"].clone(),
        "atomics": torch.zeros_like(prior["atomics"]),
        "bonds": torch.zeros_like(prior["bonds"]),
    }
    if getattr(model, "sc_charges", False):
        cond_batch["charges"] = torch.zeros_like(prior["charges"])
    if getattr(model, "add_feats", False) and "hybridization" in prior:
        cond_batch["hybridization"] = torch.zeros_like(prior["hybridization"])
    return project_coordinate_only_state(cond_batch, reference, fixed_mask)


def sample_oracle_repair(
    model: Any,
    ligand_bad: dict[str, Any],
    pocket_data: dict[str, Any],
    prior_coords: torch.Tensor,
    fixed_mask: torch.Tensor,
    steps: int = 100,
) -> dict[str, Any]:
    """Run deterministic Euler sampling while keeping graph and fixed atoms exact."""
    if steps <= 0:
        raise ValueError("steps must be positive")
    if prior_coords.shape != ligand_bad["coords"].shape:
        raise ValueError("prior_coords must match ligand coordinates")
    if fixed_mask.shape != ligand_bad["mask"].shape:
        raise ValueError("fixed_mask must match ligand mask")

    valid_mask = ligand_bad["mask"].bool()
    fixed_mask = fixed_mask.bool() & valid_mask
    reference = {
        key: value.clone() if torch.is_tensor(value) else value.copy() if isinstance(value, list) else value
        for key, value in ligand_bad.items()
    }
    reference["fragment_mask"] = fixed_mask.long()
    reference["fragment_mode"] = ["fragment_inpainting"] * prior_coords.shape[0]

    if torch.equal(fixed_mask, valid_mask):
        return project_coordinate_only_state(reference, reference, fixed_mask)

    if getattr(model.integrator, "use_sde_simulation", False):
        raise ValueError("coordinate-only repair requires deterministic ODE sampling")
    if float(getattr(model.integrator, "coord_noise_level", 0.0)) != 0.0:
        raise ValueError("coordinate-only repair requires zero integration noise")
    if not getattr(model, "inpainting_mode", False):
        raise ValueError("model must be loaded with fragment inpainting enabled")

    prior = {
        key: value.clone() if torch.is_tensor(value) else value.copy() if isinstance(value, list) else value
        for key, value in reference.items()
    }
    prior["coords"] = prior_coords.clone()
    prior = project_coordinate_only_state(prior, reference, fixed_mask)
    curr = {
        key: value.clone() if torch.is_tensor(value) else value.copy() if isinstance(value, list) else value
        for key, value in prior.items()
    }
    cond_batch = _initial_self_condition(model, prior, reference, fixed_mask)

    batch_size = prior_coords.shape[0]
    device = prior_coords.device
    times = build_repair_times(batch_size, device)
    time_points = torch.linspace(0.0, 1.0, steps + 1, device=device)
    step_sizes = time_points[1:] - time_points[:-1]

    with torch.inference_mode():
        pocket_equis, pocket_invs = model.gen.get_pocket_encoding(
            pocket_data["coords"],
            pocket_data["atom_names"],
            pocket_atom_charges=torch.argmax(pocket_data["charges"], dim=-1),
            pocket_bond_types=torch.argmax(pocket_data["bonds"], dim=-1),
            pocket_res_types=pocket_data["res_names"],
            pocket_atom_mask=pocket_data["mask"],
        )

        for step_size in step_sizes:
            cond = cond_batch if model.self_condition else None
            out = model(
                curr,
                pocket_data,
                times,
                cond_batch=cond,
                pocket_equis=pocket_equis,
                pocket_invs=pocket_invs,
                training=False,
            )
            predicted, cond_batch = model._get_predictions(out)
            cond_batch = project_coordinate_only_state(
                cond_batch, reference, fixed_mask
            )

            for _ in range(model.sc_recycle_steps - 1):
                out = model(
                    curr,
                    pocket_data,
                    times,
                    cond_batch=cond_batch if model.self_condition else None,
                    pocket_equis=pocket_equis,
                    pocket_invs=pocket_invs,
                    training=False,
                )
                predicted, cond_batch = model._get_predictions(out)
                cond_batch = project_coordinate_only_state(
                    cond_batch, reference, fixed_mask
                )

            next_coords = model.integrator.coord_step(
                curr["coords"],
                predicted["coords"],
                mask=reference["mask"],
                times=times[0],
                step_size=step_size,
            )
            curr = {**curr, "coords": next_coords}
            curr = project_coordinate_only_state(curr, reference, fixed_mask)
            times[0] = times[0] + step_size
            times[2] = times[2] + step_size

        final_times = [times[0] - 1e-4, times[1], times[2] - 1e-4]
        out = model(
            curr,
            pocket_data,
            final_times,
            cond_batch=cond_batch if model.self_condition else None,
            pocket_equis=pocket_equis,
            pocket_invs=pocket_invs,
            training=False,
        )
        predicted, _ = model._get_predictions(out)

    return project_coordinate_only_state(predicted, reference, fixed_mask)
