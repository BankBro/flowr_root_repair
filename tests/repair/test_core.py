import numpy as np
import pytest
import torch
from rdkit import Chem

from flowr.repair import (
    build_local_coordinate_prior,
    build_repair_times,
    build_torsion_corruption,
    classify_experiment,
    copy_mol_with_coords,
    max_fixed_coordinate_drift,
    molecular_graph_signature,
    project_coordinate_only_state,
    sample_oracle_repair,
    states_have_same_discrete_graph,
)


def test_local_prior_is_seeded_centered_and_keeps_fixed_atoms():
    coords = torch.tensor(
        [[0.0, 0.0, 0.0], [2.0, 1.0, -1.0], [4.0, 3.0, 1.0]]
    )
    editable = torch.tensor([False, True, True])
    first = build_local_coordinate_prior(
        coords, editable, torch.Generator().manual_seed(7)
    )
    second = build_local_coordinate_prior(
        coords, editable, torch.Generator().manual_seed(7)
    )

    assert torch.equal(first, second)
    assert torch.equal(first[~editable], coords[~editable])
    assert torch.allclose(first[editable].mean(0), coords[editable].mean(0))


def test_local_prior_returns_exact_copy_for_empty_editable_region():
    coords = torch.randn(4, 3)
    prior = build_local_coordinate_prior(coords, torch.zeros(4, dtype=torch.bool))
    assert torch.equal(prior, coords)
    assert prior.data_ptr() != coords.data_ptr()


def test_torsion_corruption_uses_right_hand_rule_and_keeps_axis_fixed():
    coords = torch.tensor(
        [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
        dtype=torch.float64,
    )
    corrupted = build_torsion_corruption(coords, 0, 1, [2], 90.0)

    assert torch.equal(corrupted[:2], coords[:2])
    assert torch.allclose(
        corrupted[2], torch.tensor([0.0, 1.0, 0.0], dtype=torch.float64), atol=1e-12
    )


def test_torsion_corruption_rejects_axis_atoms_as_editable():
    with pytest.raises(ValueError, match="axis atoms"):
        build_torsion_corruption(torch.zeros(3, 3), 0, 1, [1, 2], 30.0)


def _state(coords):
    return {
        "coords": coords.clone(),
        "atomics": torch.eye(3).unsqueeze(0),
        "charges": torch.ones(1, 3, 2),
        "hybridization": torch.ones(1, 3, 2),
        "bonds": torch.ones(1, 3, 3, 2),
        "mask": torch.tensor([[1, 1, 0]]),
        "fragment_mask": torch.tensor([[1, 0, 0]]),
        "fragment_mode": ["fragment_inpainting"],
    }


def test_projection_restores_fixed_coordinates_graph_and_padding():
    reference = _state(torch.tensor([[[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.0, 0.0]]]))
    state = _state(reference["coords"] + 9.0)
    state["atomics"] = torch.zeros_like(state["atomics"])
    state["bonds"] = torch.zeros_like(state["bonds"])
    fixed = torch.tensor([[True, False, False]])

    projected = project_coordinate_only_state(state, reference, fixed)

    assert torch.equal(projected["coords"][0, 0], reference["coords"][0, 0])
    assert torch.equal(projected["coords"][0, 1], state["coords"][0, 1])
    assert torch.equal(projected["coords"][0, 2], torch.zeros(3))
    assert states_have_same_discrete_graph(projected, reference)
    assert max_fixed_coordinate_drift(projected["coords"], reference["coords"], fixed) == 0.0


def test_repair_times_keep_discrete_time_at_one():
    times = build_repair_times(2, "cpu", coord_time=0.25, pocket_time=0.5)
    assert torch.equal(times[0], torch.tensor([0.25, 0.25]))
    assert torch.equal(times[1], torch.ones(2))
    assert torch.equal(times[2], torch.tensor([0.5, 0.5]))


def test_all_fixed_repair_bypasses_model():
    reference = _state(torch.randn(1, 3, 3))
    fixed = reference["mask"].bool()
    output = sample_oracle_repair(
        model=None,
        ligand_bad=reference,
        pocket_data={},
        prior_coords=reference["coords"] + 5.0,
        fixed_mask=fixed,
        steps=5,
    )
    valid = reference["mask"].bool()
    assert torch.equal(output["coords"][valid], reference["coords"][valid])
    assert torch.equal(output["coords"][~valid], torch.zeros_like(output["coords"][~valid]))
    assert states_have_same_discrete_graph(output, reference)


def _records(successes_per_case, runs_per_case=10):
    records = []
    for case_index, successes in enumerate(successes_per_case):
        for rollout in range(runs_per_case):
            records.append(
                {
                    "system_id": f"case-{case_index}",
                    "success": rollout < successes,
                }
            )
    return records


@pytest.mark.parametrize(
    ("successes", "outcome"),
    [
        ([2, 2, 2, 2, 1], "NO-GO"),
        ([2, 2, 2, 2, 2], "CONDITIONAL"),
        ([5, 5, 5, 5, 5], "GO"),
        ([10, 10, 4, 1, 0], "CONDITIONAL"),
    ],
)
def test_experiment_gate_boundaries(successes, outcome):
    case_ids = [f"case-{index}" for index in range(5)]
    result = classify_experiment(_records(successes), case_ids)
    assert result["outcome"] == outcome


def test_experiment_gate_marks_missing_rollout_incomplete():
    case_ids = [f"case-{index}" for index in range(5)]
    records = _records([5, 5, 5, 5, 5])[:-1]
    assert classify_experiment(records, case_ids)["outcome"] == "INCOMPLETE"


def test_rdkit_coordinate_copy_preserves_complete_graph():
    mol = Chem.MolFromSmiles("CCO")
    mol.AddConformer(Chem.Conformer(mol.GetNumAtoms()))
    signature = molecular_graph_signature(mol)
    copied = copy_mol_with_coords(mol, np.arange(9).reshape(3, 3))

    assert molecular_graph_signature(copied) == signature
    assert np.array_equal(copied.GetConformer().GetPositions(), np.arange(9).reshape(3, 3))
