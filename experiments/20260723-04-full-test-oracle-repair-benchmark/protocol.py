"""Frozen construction and evaluation protocol for the full-test benchmark."""

import hashlib
import math
import string
from dataclasses import asdict
from pathlib import Path
from typing import Any

import gemmi
import numpy as np
import pandas as pd
from posebusters import PoseBusters
from posebusters.tools.protein import get_atom_type_mask
from rdkit import Chem

from flowr.gen.utils import check_substructure_match
from flowr.repair import (
    canonical_isomeric_smiles,
    copy_mol_with_coords,
    enumerate_torsion_candidates,
    fixed_fragment_matches,
    molecular_graph_signature,
    smaller_rotatable_branches,
)
from flowr.util import rdkit as smolRD
from flowr.util.metrics import evaluate_strain


EXPERIMENT_ID = "20260723-04-full-test-oracle-repair-benchmark"
MIN_RELATIVE_DISTANCE = 0.50
MAX_RELATIVE_DISTANCE = 0.75
TARGET_RELATIVE_DISTANCE = 0.625
MIN_CLASHES = 1
MAX_CLASHES = 4
FIXED_DRIFT_TOLERANCE = 1e-5
LMDB_COORD_TOLERANCE = 5.1e-5

PB_FIELDS = [
    "mol_pred_loaded",
    "mol_cond_loaded",
    "sanitization",
    "inchi_convertible",
    "all_atoms_connected",
    "bond_lengths",
    "bond_angles",
    "internal_steric_clash",
    "aromatic_ring_flatness",
    "double_bond_flatness",
    "internal_energy",
    "protein-ligand_maximum_distance",
    "minimum_distance_to_protein",
    "minimum_distance_to_organic_cofactors",
    "minimum_distance_to_inorganic_cofactors",
    "minimum_distance_to_waters",
    "volume_overlap_with_protein",
    "volume_overlap_with_organic_cofactors",
    "volume_overlap_with_inorganic_cofactors",
    "volume_overlap_with_waters",
]

ALLOWED_BAD_PROTEIN_FAILURES = {
    "minimum_distance_to_protein",
    "volume_overlap_with_protein",
}

CONSTRUCTION_FIELDS = [
    "experiment_id",
    "test_index",
    "system_id",
    "status",
    "reason",
    "ligand_atoms",
    "baseline_pb_valid",
    "baseline_failed_fields",
    "rotatable_branch_count",
    "moderate_candidate_count",
    "pb_candidate_attempts",
    "axis_origin",
    "axis_target",
    "angle_degrees",
    "editable_atom_count",
    "editable_atom_indices",
    "fixed_atom_count",
    "num_pairwise_clashes",
    "min_relative_distance",
    "bad_failed_pb_fields",
    "max_fixed_drift",
    "good_sdf_path",
    "bad_sdf_path",
    "source_sdf_path",
    "source_cif_path",
    "source_sdf_sha256",
    "source_cif_sha256",
    "good_sdf_sha256",
    "bad_sdf_sha256",
]

EVALUATION_FIELDS = [
    "valid",
    "fully_connected_valid",
    "condition_evaluable",
    "condition_match",
    "pb_evaluable",
    "pb_valid",
    "official_quality_success",
    "same_molecule",
    "fixed_atoms_ok",
    "fixed_bonds_ok",
    "fixed_coords_ok",
    "max_fixed_drift",
    "strict_success",
    "strain_evaluable",
    "strain_energy_kcal_mol",
    "editable_rmsd_to_good",
    "all_atom_rmsd_to_good",
    *PB_FIELDS,
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_sdf(path: Path, *, sanitize: bool = True) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(
        str(path),
        removeHs=True,
        sanitize=sanitize,
        strictParsing=True,
    )
    mol = next((candidate for candidate in supplier if candidate is not None), None)
    if mol is None:
        raise ValueError(f"Could not load ligand from {path}")
    return mol


def pocket_mol_from_cif(path: Path) -> Chem.Mol:
    """Convert one centered pocket CIF to an RDKit protein molecule."""
    structure = gemmi.read_structure(str(path))
    chain_names = string.ascii_uppercase + string.ascii_lowercase + string.digits
    for model in structure:
        if len(model) > len(chain_names):
            raise ValueError(f"Too many chains for deterministic PDB conversion: {path}")
        for index, chain in enumerate(model):
            chain.name = chain_names[index]
    mol = Chem.MolFromPDBBlock(
        structure.make_pdb_string(),
        sanitize=False,
        removeHs=False,
        proximityBonding=False,
    )
    if mol is None:
        raise ValueError(f"Could not convert pocket CIF to RDKit Mol: {path}")
    return mol


def pb_results(
    buster: PoseBusters, mol: Chem.Mol, protein_mol: Chem.Mol
) -> dict[str, bool]:
    frame = buster.bust([mol], None, protein_mol)
    if len(frame) != 1:
        raise RuntimeError(f"PoseBusters returned {len(frame)} rows")
    missing = [field for field in PB_FIELDS if field not in frame.columns]
    if missing:
        raise RuntimeError(f"PoseBusters output is missing columns: {missing}")
    row = frame.iloc[0]
    return {
        field: bool(row[field]) if not pd.isna(row[field]) else False
        for field in PB_FIELDS
    }


def pb_valid(values: dict[str, Any]) -> bool:
    return bool(all(bool(values.get(field, False)) for field in PB_FIELDS))


def _relative_protein_distances(
    mol: Chem.Mol, protein_mol: Chem.Mol
) -> tuple[int, float]:
    ligand_coords = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    ligand_symbols = np.asarray([atom.GetSymbol() for atom in mol.GetAtoms()])
    ligand_heavy = ligand_symbols != "H"
    ligand_coords = ligand_coords[ligand_heavy]
    ligand_symbols = ligand_symbols[ligand_heavy]

    ignore_types = {
        "hydrogens",
        "organic_cofactors",
        "inorganic_cofactors",
        "waters",
    }
    protein_mask = np.asarray(
        get_atom_type_mask(protein_mol, ignore_types), dtype=bool
    )
    protein_coords = np.asarray(
        protein_mol.GetConformer().GetPositions(), dtype=float
    )[protein_mask]
    protein_symbols = np.asarray(
        [atom.GetSymbol() for atom in protein_mol.GetAtoms()]
    )[protein_mask]
    if not ligand_coords.size or not protein_coords.size:
        raise ValueError("Ligand or protein has no atoms for clash calculation")

    periodic_table = Chem.GetPeriodicTable()
    ligand_radii = np.asarray(
        [periodic_table.GetRvdw(str(symbol)) for symbol in ligand_symbols]
    )
    protein_radii = np.asarray(
        [periodic_table.GetRvdw(str(symbol)) for symbol in protein_symbols]
    )
    distances = np.linalg.norm(
        ligand_coords[:, None, :] - protein_coords[None, :, :], axis=-1
    )
    relative = distances / (ligand_radii[:, None] + protein_radii[None, :])
    return int((relative < MAX_RELATIVE_DISTANCE).sum()), float(relative.min())


def _base_construction_record(
    test_index: int,
    system_id: str,
    source_sdf: Path,
    source_cif: Path,
) -> dict[str, Any]:
    return {
        "experiment_id": EXPERIMENT_ID,
        "test_index": test_index,
        "system_id": system_id,
        "status": "baseline_error",
        "reason": "",
        "ligand_atoms": 0,
        "baseline_pb_valid": False,
        "baseline_failed_fields": "",
        "rotatable_branch_count": 0,
        "moderate_candidate_count": 0,
        "pb_candidate_attempts": 0,
        "axis_origin": "",
        "axis_target": "",
        "angle_degrees": "",
        "editable_atom_count": 0,
        "editable_atom_indices": "",
        "fixed_atom_count": 0,
        "num_pairwise_clashes": 0,
        "min_relative_distance": float("nan"),
        "bad_failed_pb_fields": "",
        "max_fixed_drift": float("nan"),
        "good_sdf_path": "",
        "bad_sdf_path": "",
        "source_sdf_path": str(source_sdf),
        "source_cif_path": str(source_cif),
        "source_sdf_sha256": sha256(source_sdf),
        "source_cif_sha256": sha256(source_cif),
        "good_sdf_sha256": "",
        "bad_sdf_sha256": "",
    }


def construct_clash_case(
    *,
    test_index: int,
    system: Any,
    source_sdf: Path,
    source_cif: Path,
    buster: PoseBusters,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    """Construct one deterministic isolated torsion clash or an attrition record."""
    system_id = str(system.metadata["system_id"])
    record = _base_construction_record(
        test_index, system_id, source_sdf, source_cif
    )
    try:
        good_mol = load_sdf(source_sdf, sanitize=True)
        protein_mol = pocket_mol_from_cif(source_cif)
        record["ligand_atoms"] = good_mol.GetNumAtoms()
        raw_atomics = [atom.GetAtomicNum() for atom in good_mol.GetAtoms()]
        if raw_atomics != system.ligand.atomics.tolist():
            raise ValueError("SDF and LMDB ligand atom order differ")
        coords_good = np.asarray(
            good_mol.GetConformer().GetPositions(), dtype=np.float64
        )
        coord_delta = float(
            np.max(np.abs(coords_good - system.ligand.coords.numpy()))
        )
        if coord_delta > LMDB_COORD_TOLERANCE:
            raise ValueError(f"SDF and LMDB coordinates differ by {coord_delta}")
        baseline_pb = pb_results(buster, good_mol, protein_mol)
    except Exception as error:
        record["reason"] = f"{type(error).__name__}: {error}"
        return record, None

    record["baseline_pb_valid"] = pb_valid(baseline_pb)
    record["baseline_failed_fields"] = ",".join(
        field for field in PB_FIELDS if not baseline_pb[field]
    )
    if not record["baseline_pb_valid"]:
        record["status"] = "baseline_pb_failed"
        record["reason"] = record["baseline_failed_fields"]
        return record, None

    branches = smaller_rotatable_branches(good_mol)
    record["rotatable_branch_count"] = len(branches)
    if not branches:
        record["status"] = "no_rotatable_bond"
        record["reason"] = "RDKit found no non-empty smaller rotatable branch"
        return record, None

    moderate = []
    for candidate in enumerate_torsion_candidates(good_mol):
        bad_mol = copy_mol_with_coords(good_mol, candidate.coords)
        clashes, min_relative = _relative_protein_distances(bad_mol, protein_mol)
        if (
            MIN_CLASHES <= clashes <= MAX_CLASHES
            and MIN_RELATIVE_DISTANCE <= min_relative < MAX_RELATIVE_DISTANCE
        ):
            moderate.append(
                (
                    abs(min_relative - TARGET_RELATIVE_DISTANCE),
                    abs(candidate.angle_degrees),
                    len(candidate.editable_atom_indices),
                    candidate.axis_origin,
                    candidate.axis_target,
                    candidate.angle_degrees,
                    clashes,
                    min_relative,
                    candidate,
                )
            )
    record["moderate_candidate_count"] = len(moderate)
    if not moderate:
        record["status"] = "no_moderate_candidate"
        record["reason"] = "No candidate met the frozen clash severity interval"
        return record, None

    for ranked in sorted(moderate, key=lambda item: item[:-1]):
        *_, clashes, min_relative, candidate = ranked
        record["pb_candidate_attempts"] += 1
        bad_mol = copy_mol_with_coords(good_mol, candidate.coords)
        if molecular_graph_signature(bad_mol) != molecular_graph_signature(good_mol):
            continue
        fixed_mask = np.ones(good_mol.GetNumAtoms(), dtype=bool)
        fixed_mask[np.asarray(candidate.editable_atom_indices, dtype=int)] = False
        fixed_drift = float(
            np.linalg.norm(
                candidate.coords[fixed_mask] - coords_good[fixed_mask], axis=-1
            ).max()
        )
        if fixed_drift > 1e-10:
            continue
        bad_pb = pb_results(buster, bad_mol, protein_mol)
        nonprotein_ok = all(
            bad_pb[field]
            for field in PB_FIELDS
            if field not in ALLOWED_BAD_PROTEIN_FAILURES
        )
        if nonprotein_ok and not bad_pb["minimum_distance_to_protein"]:
            record.update(
                {
                    "status": "eligible",
                    "reason": "",
                    "axis_origin": candidate.axis_origin,
                    "axis_target": candidate.axis_target,
                    "angle_degrees": candidate.angle_degrees,
                    "editable_atom_count": len(candidate.editable_atom_indices),
                    "editable_atom_indices": ",".join(
                        str(index) for index in candidate.editable_atom_indices
                    ),
                    "fixed_atom_count": int(fixed_mask.sum()),
                    "num_pairwise_clashes": clashes,
                    "min_relative_distance": min_relative,
                    "bad_failed_pb_fields": ",".join(
                        field for field in PB_FIELDS if not bad_pb[field]
                    ),
                    "max_fixed_drift": fixed_drift,
                }
            )
            material = {
                "test_index": test_index,
                "system_id": system_id,
                "good_mol": good_mol,
                "bad_mol": bad_mol,
                "protein_mol": protein_mol,
                "coords_good": coords_good,
                "coords_bad": candidate.coords,
                "fixed_mask": fixed_mask,
                "candidate": asdict(candidate) | {"coords": candidate.coords},
            }
            return record, material

    record["status"] = "no_isolated_pb_candidate"
    record["reason"] = "Moderate candidates failed non-protein PoseBusters checks"
    return record, None


def official_condition_match(
    generated: Chem.Mol,
    reference_bad: Chem.Mol,
    editable_atom_indices: list[int],
) -> bool:
    generated_copy = Chem.Mol(generated)
    reference_copy = Chem.Mol(reference_bad)
    Chem.SanitizeMol(generated_copy)
    Chem.SanitizeMol(reference_copy)
    return bool(
        check_substructure_match(
            generated_copy,
            reference_copy,
            inpainting_mode="substructure_inpainting",
            substructure_query=editable_atom_indices,
        )
    )


def _rmsd(coords: np.ndarray, reference: np.ndarray, mask: np.ndarray) -> float:
    if not np.asarray(mask, dtype=bool).any():
        return 0.0
    delta = coords[mask] - reference[mask]
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=-1))))


def evaluate_output(
    *,
    mol_pred: Chem.Mol | None,
    material: dict[str, Any],
    output_to_reference: np.ndarray,
    raw_fixed_drift: float,
    run_completed: bool,
    buster: PoseBusters,
) -> dict[str, Any]:
    """Evaluate one output under official quality and strict repair endpoints."""
    result: dict[str, Any] = {
        "valid": False,
        "fully_connected_valid": False,
        "condition_evaluable": False,
        "condition_match": False,
        "pb_evaluable": False,
        "pb_valid": False,
        "official_quality_success": False,
        "same_molecule": False,
        "fixed_atoms_ok": False,
        "fixed_bonds_ok": False,
        "fixed_coords_ok": False,
        "max_fixed_drift": raw_fixed_drift,
        "strict_success": False,
        "strain_evaluable": False,
        "strain_energy_kcal_mol": float("nan"),
        "editable_rmsd_to_good": float("nan"),
        "all_atom_rmsd_to_good": float("nan"),
        **{field: False for field in PB_FIELDS},
    }
    if not run_completed or mol_pred is None:
        return result

    result["valid"] = bool(smolRD.mol_is_valid(mol_pred, connected=False))
    result["fully_connected_valid"] = bool(
        smolRD.mol_is_valid(mol_pred, connected=True)
    )
    if result["fully_connected_valid"]:
        result["condition_evaluable"] = True
        try:
            result["condition_match"] = official_condition_match(
                mol_pred,
                material["bad_mol"],
                list(material["editable_atom_indices"]),
            )
        except Exception:
            result["condition_evaluable"] = False
        try:
            strain_mol = Chem.Mol(mol_pred)
            Chem.SanitizeMol(strain_mol)
            strain_mol = Chem.AddHs(strain_mol, addCoords=True)
            strain = float(
                evaluate_strain(
                    [strain_mol],
                    n_steps=500,
                    add_hs=False,
                    force_field_name="MMFF94s",
                    return_list=True,
                )[0]
            )
            result["strain_energy_kcal_mol"] = strain
            result["strain_evaluable"] = math.isfinite(strain)
        except Exception:
            pass

    try:
        values = pb_results(buster, mol_pred, material["protein_mol"])
        result.update(values)
        result["pb_evaluable"] = True
        result["pb_valid"] = pb_valid(values)
    except Exception:
        pass

    result["official_quality_success"] = bool(
        result["fully_connected_valid"]
        and result["condition_match"]
        and result["pb_valid"]
    )
    result["same_molecule"] = bool(
        canonical_isomeric_smiles(mol_pred)
        == canonical_isomeric_smiles(material["good_mol"])
    )

    output_to_reference = np.asarray(output_to_reference, dtype=int)
    coords_output = np.asarray(
        mol_pred.GetConformer().GetPositions(), dtype=float
    )
    n_reference = material["coords_good"].shape[0]
    mapping_ok = bool(
        output_to_reference.size == coords_output.shape[0]
        and output_to_reference.size == n_reference
        and np.array_equal(np.sort(output_to_reference), np.arange(n_reference))
    )
    if mapping_ok:
        fixed_indices = np.flatnonzero(material["fixed_mask"])
        fixed_atoms_ok, fixed_bonds_ok = fixed_fragment_matches(
            mol_pred,
            material["bad_mol"],
            output_to_reference,
            fixed_indices,
        )
        result["fixed_atoms_ok"] = fixed_atoms_ok
        result["fixed_bonds_ok"] = fixed_bonds_ok
        coords_reference_order = np.empty_like(coords_output)
        coords_reference_order[output_to_reference] = coords_output
        if fixed_indices.size:
            measured_drift = float(
                np.linalg.norm(
                    coords_reference_order[fixed_indices]
                    - material["coords_bad"][fixed_indices],
                    axis=-1,
                ).max()
            )
        else:
            measured_drift = 0.0
        result["max_fixed_drift"] = max(raw_fixed_drift, measured_drift)
        result["fixed_coords_ok"] = bool(
            np.isfinite(result["max_fixed_drift"])
            and result["max_fixed_drift"] <= FIXED_DRIFT_TOLERANCE
        )
        editable_mask = ~np.asarray(material["fixed_mask"], dtype=bool)
        result["editable_rmsd_to_good"] = _rmsd(
            coords_reference_order, material["coords_good"], editable_mask
        )
        result["all_atom_rmsd_to_good"] = _rmsd(
            coords_reference_order,
            material["coords_good"],
            np.ones(n_reference, dtype=bool),
        )

    result["strict_success"] = bool(
        result["official_quality_success"]
        and result["same_molecule"]
        and result["fixed_atoms_ok"]
        and result["fixed_bonds_ok"]
        and result["fixed_coords_ok"]
    )
    return result
