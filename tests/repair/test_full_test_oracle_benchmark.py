import importlib.util
import sys
from pathlib import Path

import numpy as np
import torch
from rdkit import Chem

from flowr.repair import (
    fixed_first_index_map,
    restore_quantized_fixed_coordinates,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
EXPERIMENT_DIR = (
    REPO_ROOT / "experiments" / "20260723-04-full-test-oracle-repair-benchmark"
)
sys.path.insert(0, str(EXPERIMENT_DIR))


def _load_module(name: str):
    path = EXPERIMENT_DIR / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"full_benchmark_{name}", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


protocol = _load_module("protocol")
runner = _load_module("run")


def _mol(smiles: str, coords: list[tuple[float, float, float]]) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for index, position in enumerate(coords):
        conformer.SetAtomPosition(index, position)
    mol.AddConformer(conformer)
    return mol


def test_shared_official_adapter_preserves_historical_mapping_and_drift():
    fixed = np.array([True, False, True, False])
    mapping = fixed_first_index_map(fixed)
    coords_bad = np.arange(12, dtype=float).reshape(4, 3)
    coords_output = coords_bad[mapping].copy()
    coords_output[:2] += 1e-6

    restored, raw_drift, ok = restore_quantized_fixed_coordinates(
        coords_output, coords_bad, fixed, mapping
    )

    assert np.array_equal(mapping, np.array([0, 2, 1, 3]))
    assert ok and raw_drift > 0.0
    assert np.array_equal(restored[:2], coords_bad[mapping[:2]])


def test_official_adapter_normalizes_lmdb_numpy_interactions():
    system = type("System", (), {})()
    system.interactions = np.zeros((2, 3, 4), dtype=np.int64)

    normalized = runner._tensorize_interactions(system)

    assert normalized is system
    assert torch.is_tensor(normalized.interactions)
    assert normalized.interactions.shape == (2, 3, 4)


def test_official_pb_valid_requires_all_configured_fields():
    values = {field: True for field in protocol.PB_FIELDS}
    assert protocol.pb_valid(values)
    values["internal_energy"] = False
    assert not protocol.pb_valid(values)


def test_official_quality_can_pass_while_strict_identity_fails(monkeypatch):
    good = _mol("CCO", [(0, 0, 0), (1, 0, 0), (2, 0, 0)])
    generated = _mol("CCN", [(0, 0, 0), (1, 0, 0), (2, 1, 0)])
    material = {
        "good_mol": good,
        "bad_mol": good,
        "protein_mol": Chem.MolFromSmiles("C"),
        "coords_good": np.asarray(good.GetConformer().GetPositions()),
        "coords_bad": np.asarray(good.GetConformer().GetPositions()),
        "fixed_mask": np.array([True, True, False]),
        "editable_atom_indices": [2],
    }
    material["protein_mol"].AddConformer(Chem.Conformer(1))
    monkeypatch.setattr(protocol.smolRD, "mol_is_valid", lambda mol, connected: True)
    monkeypatch.setattr(protocol, "official_condition_match", lambda *args: True)
    monkeypatch.setattr(
        protocol,
        "pb_results",
        lambda *args: {field: True for field in protocol.PB_FIELDS},
    )
    monkeypatch.setattr(protocol, "evaluate_strain", lambda *args, **kwargs: [1.0])
    monkeypatch.setattr(protocol, "fixed_fragment_matches", lambda *args: (True, True))

    result = protocol.evaluate_output(
        mol_pred=generated,
        material=material,
        output_to_reference=np.arange(3),
        raw_fixed_drift=0.0,
        run_completed=True,
        buster=object(),
    )

    assert result["official_quality_success"]
    assert not result["same_molecule"]
    assert not result["strict_success"]


def test_cluster_bootstrap_keeps_complex_as_sampling_unit():
    records = []
    for system_id in ("a", "b", "c"):
        for seed in range(10):
            records.extend(
                [
                    {
                        "method": "coordinate_only",
                        "system_id": system_id,
                        "seed": seed,
                        "strict_success": True,
                    },
                    {
                        "method": "official_inpainting",
                        "system_id": system_id,
                        "seed": seed,
                        "strict_success": False,
                    },
                ]
            )

    result = runner._paired_bootstrap(records, ["a", "b", "c"], "strict_success")

    assert result["unit"] == "complex"
    assert result["coordinate_minus_official"] == 1.0
    assert result["ci95"] == [1.0, 1.0]
    assert result["conclusion"] == "coordinate_only_better"


def test_formal_seed_inventory_is_disjoint_by_test_index():
    first = {runner.FORMAL_SEED_BASE + rollout for rollout in range(10)}
    second = {
        runner.FORMAL_SEED_BASE + runner.ROLLOUTS_PER_CASE + rollout
        for rollout in range(10)
    }
    assert len(first) == len(second) == 10
    assert first.isdisjoint(second)


def test_sampling_inventory_requires_exact_paired_keys():
    cases = [{"system_id": "case-a", "test_index": 2}]
    records = []
    for method in ("coordinate_only", "official_inpainting"):
        for rollout in range(runner.ROLLOUTS_PER_CASE):
            records.append(
                {
                    "method": method,
                    "system_id": "case-a",
                    "seed": str(
                        runner.FORMAL_SEED_BASE
                        + 2 * runner.ROLLOUTS_PER_CASE
                        + rollout
                    ),
                }
            )

    runner._validate_sampling_inventory(records, cases)

    records[-1] = {**records[-1], "seed": "999"}
    try:
        runner._validate_sampling_inventory(records, cases)
    except RuntimeError as error:
        assert "frozen paired key set" in str(error)
    else:
        raise AssertionError("Mismatched paired inventory was accepted")


def test_endpoint_funnel_is_cumulative():
    passed = {
        "run_status": "completed",
        "valid": True,
        "fully_connected_valid": True,
        "condition_match": True,
        "pb_valid": True,
        "same_molecule": True,
        "fixed_atoms_ok": True,
        "fixed_bonds_ok": True,
        "fixed_coords_ok": True,
    }
    disconnected = {**passed, "fully_connected_valid": False}
    disconnected["condition_match"] = True

    funnel = runner._endpoint_funnel([passed, disconnected])
    counts = {stage["stage"]: stage["remaining"] for stage in funnel}

    assert counts["valid"] == 2
    assert counts["fully_connected_valid"] == 1
    assert counts["condition_match"] == 1
    assert counts["fixed_coords_ok"] == 1
