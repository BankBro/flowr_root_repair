import importlib.util
from pathlib import Path

import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import AllChem

from flowr.repair import (
    canonical_isomeric_smiles,
    fixed_fragment_matches,
    summarize_repair_funnel,
)


ADAPTER_PATH = (
    Path(__file__).resolve().parents[2]
    / "experiments"
    / "20260723-02-official-fragment-inpainting-control"
    / "adapter.py"
)
SPEC = importlib.util.spec_from_file_location("official_control_adapter", ADAPTER_PATH)
adapter = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(adapter)


def _embedded(smiles: str) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=17) == 0
    return Chem.RemoveHs(mol)


def test_canonical_identity_is_order_independent_and_stereo_aware():
    reference = _embedded("F[C@H](Cl)Br")
    reordered = Chem.RenumberAtoms(reference, list(reversed(range(reference.GetNumAtoms()))))
    opposite = _embedded("F[C@@H](Cl)Br")

    assert canonical_isomeric_smiles(reordered) == canonical_isomeric_smiles(reference)
    assert canonical_isomeric_smiles(opposite) != canonical_isomeric_smiles(reference)


def test_fixed_fragment_matching_uses_output_to_reference_map():
    reference = _embedded("CCCO")
    output_to_reference = np.array([0, 2, 1, 3])
    prediction = Chem.RenumberAtoms(reference, output_to_reference.tolist())

    atoms_ok, bonds_ok = fixed_fragment_matches(
        prediction,
        reference,
        output_to_reference,
        np.array([0, 2]),
    )
    assert atoms_ok
    assert bonds_ok


def test_fixed_first_index_map_and_small_drift_restoration():
    fixed = np.array([True, False, True, False])
    mapping = adapter.fixed_first_index_map(fixed)
    assert np.array_equal(mapping, np.array([0, 2, 1, 3]))

    coords_bad = np.arange(12, dtype=float).reshape(4, 3)
    coords_output = coords_bad[mapping].copy()
    coords_output[:2] += 2e-6
    restored, raw_drift, ok = adapter.restore_quantized_fixed_coordinates(
        coords_output,
        coords_bad,
        fixed,
        mapping,
    )
    assert ok
    assert raw_drift > 0.0
    assert np.array_equal(restored[:2], coords_bad[mapping[:2]])


def test_large_fixed_drift_is_not_hidden_by_restoration():
    fixed = np.array([True, False])
    mapping = adapter.fixed_first_index_map(fixed)
    coords_bad = np.zeros((2, 3))
    coords_output = coords_bad[mapping].copy()
    coords_output[0, 0] = 1e-3

    restored, raw_drift, ok = adapter.restore_quantized_fixed_coordinates(
        coords_output,
        coords_bad,
        fixed,
        mapping,
    )
    assert not ok
    assert raw_drift == 1e-3
    assert restored[0, 0] == 1e-3


def test_exact_adapter_passes_true_mask_to_official_private_builder():
    class Ligand:
        seq_length = 3
        coords = torch.arange(9).reshape(3, 3)
        atomics = torch.eye(3)
        charges = torch.eye(3)
        hybridization = torch.eye(3)

    class Interpolant:
        def _build_inference_prior(self, target, rdkit_mol, mode, mask, is_local):
            self.call = (mode, mask.clone(), is_local)
            order = torch.cat([torch.where(mask)[0], torch.where(~mask)[0]])
            prior = type("Prior", (), {})()
            prior.coords = target.coords[order].clone()
            prior.atomics = target.atomics[order].clone()
            prior.charges = target.charges[order].clone()
            prior.hybridization = target.hybridization[order].clone()
            prior.fragment_mask = torch.tensor([True, True, False])
            return prior

    interpolant = Interpolant()
    fixed = torch.tensor([True, False, True])
    _, mapping = adapter.build_exact_fragment_prior(
        interpolant, Ligand(), rdkit_mol=None, fixed_mask=fixed
    )

    assert interpolant.call[0] == "fragment_inpainting"
    assert torch.equal(interpolant.call[1], fixed)
    assert interpolant.call[2] is True
    assert np.array_equal(mapping, np.array([0, 2, 1]))


def test_funnel_is_cumulative():
    records = [
        {
            "run_status": "completed",
            "usable_output": True,
            "no_protein_clash": True,
            "geometry_ok": True,
            "same_molecule": True,
        },
        {
            "run_status": "completed",
            "usable_output": True,
            "no_protein_clash": True,
            "geometry_ok": False,
            "same_molecule": True,
        },
        {
            "run_status": "model_error",
            "usable_output": False,
            "no_protein_clash": False,
            "geometry_ok": False,
            "same_molecule": False,
        },
    ]
    funnel = summarize_repair_funnel(records)
    assert [stage["remaining"] for stage in funnel] == [3, 2, 2, 1, 1]
