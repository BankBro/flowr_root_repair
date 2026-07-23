"""Evaluation and experiment gates for coordinate-only repair."""

from collections import Counter
from collections.abc import Iterable
from typing import Any

import numpy as np
from rdkit import Chem


IGNORED_PROTEIN_TYPES = {
    "hydrogens",
    "organic_cofactors",
    "inorganic_cofactors",
    "waters",
}


def copy_mol_with_coords(mol: Chem.Mol, coords: np.ndarray) -> Chem.Mol:
    """Copy an RDKit molecule and replace its conformer coordinates."""
    coords = np.asarray(coords, dtype=float)
    if coords.shape != (mol.GetNumAtoms(), 3):
        raise ValueError("coords must have shape [mol.GetNumAtoms(), 3]")
    copied = Chem.Mol(mol)
    if copied.GetNumConformers() == 0:
        copied.AddConformer(Chem.Conformer(copied.GetNumAtoms()))
    conformer = copied.GetConformer()
    for index, position in enumerate(coords):
        conformer.SetAtomPosition(index, tuple(float(value) for value in position))
    return copied


def molecular_graph_signature(mol: Chem.Mol) -> tuple[Any, ...]:
    """Return an atom-order-sensitive signature of the complete discrete graph."""
    atoms = tuple(
        (
            atom.GetAtomicNum(),
            atom.GetFormalCharge(),
            int(atom.GetHybridization()),
            atom.GetIsAromatic(),
            int(atom.GetChiralTag()),
        )
        for atom in mol.GetAtoms()
    )
    bonds = tuple(
        sorted(
            (
                min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
                max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
                str(bond.GetBondType()),
                bond.GetIsAromatic(),
                int(bond.GetStereo()),
            )
            for bond in mol.GetBonds()
        )
    )
    return atoms, bonds


def canonical_isomeric_smiles(mol: Chem.Mol | None) -> str | None:
    """Return a canonical stereo-aware identity after RDKit sanitization."""
    if mol is None:
        return None
    try:
        sanitized = Chem.Mol(mol)
        Chem.SanitizeMol(sanitized)
        return Chem.MolToSmiles(
            sanitized,
            canonical=True,
            isomericSmiles=True,
        )
    except Exception:
        return None


def _atom_identity(atom: Chem.Atom) -> tuple[Any, ...]:
    return (
        atom.GetAtomicNum(),
        atom.GetFormalCharge(),
        atom.GetIsotope(),
        atom.GetIsAromatic(),
    )


def _bond_identity(bond: Chem.Bond | None) -> tuple[Any, ...] | None:
    if bond is None:
        return None
    return str(bond.GetBondType()), bond.GetIsAromatic()


def fixed_fragment_matches(
    mol_pred: Chem.Mol | None,
    mol_reference: Chem.Mol,
    output_to_reference: np.ndarray,
    fixed_reference_indices: np.ndarray,
) -> tuple[bool, bool]:
    """Check atom identities and the induced fixed-fragment graph under reordering."""
    if mol_pred is None:
        return False, False

    try:
        sanitized_pred = Chem.Mol(mol_pred)
        sanitized_reference = Chem.Mol(mol_reference)
        Chem.SanitizeMol(sanitized_pred)
        Chem.SanitizeMol(sanitized_reference)
    except Exception:
        return False, False

    output_to_reference = np.asarray(output_to_reference, dtype=int)
    fixed_reference_indices = np.asarray(fixed_reference_indices, dtype=int)
    if sanitized_pred.GetNumAtoms() != output_to_reference.size:
        return False, False

    fixed_set = set(fixed_reference_indices.tolist())
    output_indices = [
        output_index
        for output_index, reference_index in enumerate(output_to_reference.tolist())
        if reference_index in fixed_set
    ]
    if len(output_indices) != fixed_reference_indices.size:
        return False, False

    atoms_ok = all(
        _atom_identity(sanitized_pred.GetAtomWithIdx(output_index))
        == _atom_identity(
            sanitized_reference.GetAtomWithIdx(int(output_to_reference[output_index]))
        )
        for output_index in output_indices
    )
    bonds_ok = True
    for left_pos, left_output in enumerate(output_indices):
        left_reference = int(output_to_reference[left_output])
        for right_output in output_indices[left_pos + 1 :]:
            right_reference = int(output_to_reference[right_output])
            pred_bond = sanitized_pred.GetBondBetweenAtoms(left_output, right_output)
            reference_bond = sanitized_reference.GetBondBetweenAtoms(
                left_reference, right_reference
            )
            if _bond_identity(pred_bond) != _bond_identity(reference_bond):
                bonds_ok = False
                break
        if not bonds_ok:
            break
    return bool(atoms_ok), bool(bonds_ok)


def _empty_pose_metrics() -> dict[str, Any]:
    return {
        "finite_coords": False,
        "sanitized": False,
        "no_protein_clash": False,
        "num_pairwise_clashes": -1,
        "min_relative_distance": float("nan"),
        "bond_lengths_ok": False,
        "bond_angles_ok": False,
        "no_internal_clash": False,
        "geometry_ok": False,
    }


def evaluate_inpainting_candidate(
    mol_pred: Chem.Mol | None,
    mol_bad: Chem.Mol,
    protein_mol: Chem.Mol,
    coords_output: np.ndarray,
    coords_bad: np.ndarray,
    fixed_reference_indices: np.ndarray,
    output_to_reference: np.ndarray,
    *,
    coords_good: np.ndarray | None = None,
    raw_fixed_drift: float | None = None,
    fixed_drift_tolerance: float = 1e-5,
    run_completed: bool = True,
) -> dict[str, Any]:
    """Evaluate native local redesign and strict same-molecule repair endpoints."""
    coords_output = np.asarray(coords_output, dtype=float)
    coords_bad = np.asarray(coords_bad, dtype=float)
    output_to_reference = np.asarray(output_to_reference, dtype=int)
    fixed_reference_indices = np.asarray(fixed_reference_indices, dtype=int)

    shape_ok = coords_output.shape == (output_to_reference.size, 3)
    mapping_ok = bool(
        output_to_reference.size == coords_bad.shape[0]
        and np.array_equal(np.sort(output_to_reference), np.arange(coords_bad.shape[0]))
    )
    if shape_ok and mapping_ok:
        coords_reference_order = np.empty_like(coords_output)
        coords_reference_order[output_to_reference] = coords_output
        calculated_fixed_drift = (
            float(
                np.linalg.norm(
                    coords_reference_order[fixed_reference_indices]
                    - coords_bad[fixed_reference_indices],
                    axis=-1,
                ).max()
            )
            if fixed_reference_indices.size
            else 0.0
        )
    else:
        coords_reference_order = np.full_like(coords_bad, np.nan, dtype=float)
        calculated_fixed_drift = float("inf")

    maximum_fixed_drift = (
        calculated_fixed_drift if raw_fixed_drift is None else float(raw_fixed_drift)
    )
    fixed_coords_ok = bool(
        np.isfinite(maximum_fixed_drift)
        and maximum_fixed_drift <= fixed_drift_tolerance
    )

    pose = evaluate_pose(mol_pred, protein_mol) if mol_pred is not None else _empty_pose_metrics()
    canonical_pred = canonical_isomeric_smiles(mol_pred)
    canonical_bad = canonical_isomeric_smiles(mol_bad)
    single_component = bool(
        canonical_pred is not None
        and mol_pred is not None
        and len(Chem.GetMolFrags(mol_pred)) == 1
    )
    fixed_atoms_ok, fixed_bonds_ok = fixed_fragment_matches(
        mol_pred,
        mol_bad,
        output_to_reference,
        fixed_reference_indices,
    )
    fixed_fragment_retained = bool(
        fixed_coords_ok and fixed_atoms_ok and fixed_bonds_ok
    )
    same_molecule = bool(
        canonical_pred is not None
        and canonical_bad is not None
        and canonical_pred == canonical_bad
    )

    editable_mask = np.ones(coords_bad.shape[0], dtype=bool)
    editable_mask[fixed_reference_indices] = False
    if coords_good is None or not shape_ok or not mapping_ok:
        editable_rmsd = float("nan")
        all_atom_rmsd = float("nan")
    else:
        coords_good = np.asarray(coords_good, dtype=float)
        editable_rmsd = _direct_rmsd(
            coords_reference_order, coords_good, editable_mask
        )
        all_atom_rmsd = _direct_rmsd(
            coords_reference_order,
            coords_good,
            np.ones(coords_bad.shape[0], dtype=bool),
        )

    usable_output = bool(
        run_completed
        and pose["finite_coords"]
        and pose["sanitized"]
        and single_component
        and fixed_fragment_retained
    )
    native_success = bool(
        usable_output and pose["no_protein_clash"] and pose["geometry_ok"]
    )
    strict_success = bool(native_success and same_molecule)
    return {
        **pose,
        "single_component": single_component,
        "fixed_atoms_ok": fixed_atoms_ok,
        "fixed_bonds_ok": fixed_bonds_ok,
        "fixed_coords_ok": fixed_coords_ok,
        "fixed_fragment_retained": fixed_fragment_retained,
        "raw_max_fixed_drift": maximum_fixed_drift,
        "same_molecule": same_molecule,
        "canonical_isomeric_smiles": canonical_pred or "",
        "editable_rmsd_to_good": editable_rmsd,
        "all_atom_rmsd_to_good": all_atom_rmsd,
        "usable_output": usable_output,
        "native_success": native_success,
        "strict_success": strict_success,
    }


def summarize_repair_funnel(records: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    """Build cumulative native and strict endpoint attrition counts."""
    records = list(records)
    predicates = [
        ("attempted", lambda record: True),
        (
            "completed_usable_fixed",
            lambda record: record.get("run_status") == "completed"
            and bool(record.get("usable_output", False)),
        ),
        ("no_protein_clash", lambda record: bool(record.get("no_protein_clash", False))),
        ("internal_geometry", lambda record: bool(record.get("geometry_ok", False))),
        ("same_molecule", lambda record: bool(record.get("same_molecule", False))),
    ]
    active = list(records)
    funnel = []
    previous = len(active)
    for stage, predicate in predicates:
        if stage != "attempted":
            active = [record for record in active if predicate(record)]
        remaining = len(active)
        funnel.append(
            {
                "stage": stage,
                "eliminated_at_stage": 0 if stage == "attempted" else previous - remaining,
                "remaining": remaining,
                "fraction_of_attempts": remaining / len(records) if records else 0.0,
            }
        )
        previous = remaining
    return funnel


def evaluate_pose(mol: Chem.Mol, protein_mol: Chem.Mol) -> dict[str, Any]:
    """Run the frozen intermolecular and internal geometry checks."""
    from posebusters.modules.distance_geometry import check_geometry
    from posebusters.modules.intermolecular_distance import (
        check_intermolecular_distance,
    )

    finite_coords = bool(np.isfinite(mol.GetConformer().GetPositions()).all())
    sanitized = False
    if finite_coords:
        try:
            sanitized_mol = Chem.Mol(mol)
            sanitized = Chem.SanitizeMol(sanitized_mol) == Chem.SanitizeFlags.SANITIZE_NONE
        except Exception:
            sanitized = False

    try:
        intermolecular = check_intermolecular_distance(
            mol,
            protein_mol,
            radius_type="vdw",
            radius_scale=1.0,
            clash_cutoff=0.75,
            ignore_types=IGNORED_PROTEIN_TYPES,
        )["results"]
    except Exception:
        intermolecular = {}

    try:
        geometry = check_geometry(mol)["results"]
    except Exception:
        geometry = {}

    no_clashes = bool(intermolecular.get("no_clashes", False))
    bond_lengths_ok = bool(geometry.get("bond_lengths_within_bounds", False))
    bond_angles_ok = bool(geometry.get("bond_angles_within_bounds", False))
    no_internal_clash = bool(geometry.get("no_internal_clash", False))
    return {
        "finite_coords": finite_coords,
        "sanitized": sanitized,
        "no_protein_clash": no_clashes,
        "num_pairwise_clashes": int(intermolecular.get("num_pairwise_clashes", -1)),
        "min_relative_distance": float(
            intermolecular.get("most_extreme_relative_distance", float("nan"))
        ),
        "bond_lengths_ok": bond_lengths_ok,
        "bond_angles_ok": bond_angles_ok,
        "no_internal_clash": no_internal_clash,
        "geometry_ok": bool(
            sanitized and bond_lengths_ok and bond_angles_ok and no_internal_clash
        ),
    }


def _direct_rmsd(coords: np.ndarray, reference: np.ndarray, mask: np.ndarray) -> float:
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return 0.0
    squared = np.sum((coords[mask] - reference[mask]) ** 2, axis=-1)
    return float(np.sqrt(np.mean(squared)))


def evaluate_repair(
    mol_pred: Chem.Mol,
    mol_bad: Chem.Mol,
    protein_mol: Chem.Mol,
    coords_pred: np.ndarray,
    coords_bad: np.ndarray,
    coords_good: np.ndarray,
    fixed_mask: np.ndarray,
    discrete_state_equal: bool,
    fixed_drift_tolerance: float = 1e-6,
) -> dict[str, Any]:
    """Evaluate one repair rollout using the frozen four-part success rule."""
    coords_pred = np.asarray(coords_pred, dtype=float)
    coords_bad = np.asarray(coords_bad, dtype=float)
    coords_good = np.asarray(coords_good, dtype=float)
    fixed_mask = np.asarray(fixed_mask, dtype=bool)
    editable_mask = ~fixed_mask

    pose = evaluate_pose(mol_pred, protein_mol)
    graph_unchanged = bool(
        discrete_state_equal
        and molecular_graph_signature(mol_pred) == molecular_graph_signature(mol_bad)
    )
    if fixed_mask.any():
        fixed_drift = float(
            np.linalg.norm(coords_pred[fixed_mask] - coords_bad[fixed_mask], axis=-1).max()
        )
    else:
        fixed_drift = 0.0
    fixed_coords_ok = bool(fixed_drift <= fixed_drift_tolerance)

    result = {
        **pose,
        "graph_unchanged": graph_unchanged,
        "fixed_coords_ok": fixed_coords_ok,
        "max_fixed_drift": fixed_drift,
        "editable_rmsd_to_good": _direct_rmsd(
            coords_pred, coords_good, editable_mask
        ),
        "all_atom_rmsd_to_good": _direct_rmsd(
            coords_pred, coords_good, np.ones_like(fixed_mask, dtype=bool)
        ),
    }
    result["success"] = bool(
        result["finite_coords"]
        and result["no_protein_clash"]
        and result["geometry_ok"]
        and result["fixed_coords_ok"]
        and result["graph_unchanged"]
    )
    return result


def classify_experiment(
    records: Iterable[dict[str, Any]],
    expected_system_ids: Iterable[str],
    rollouts_per_system: int = 10,
) -> dict[str, Any]:
    """Apply the frozen GO, NO-GO, CONDITIONAL, or INCOMPLETE gate."""
    records = list(records)
    expected_system_ids = list(expected_system_ids)
    counts = Counter(record["system_id"] for record in records)
    successes = Counter(
        record["system_id"] for record in records if bool(record.get("success", False))
    )
    complete = len(records) == len(expected_system_ids) * rollouts_per_system and all(
        counts[system_id] == rollouts_per_system for system_id in expected_system_ids
    )
    total_successes = sum(successes.values())

    if not complete:
        outcome = "INCOMPLETE"
    elif total_successes <= 9:
        outcome = "NO-GO"
    elif total_successes >= 25 and sum(
        successes[system_id] >= 3 for system_id in expected_system_ids
    ) >= 4:
        outcome = "GO"
    else:
        outcome = "CONDITIONAL"

    return {
        "outcome": outcome,
        "complete": complete,
        "total_runs": len(records),
        "total_successes": total_successes,
        "overall_success_rate": (
            total_successes / len(records) if records else 0.0
        ),
        "per_case": {
            system_id: {
                "runs": counts[system_id],
                "successes": successes[system_id],
                "success_rate": (
                    successes[system_id] / counts[system_id]
                    if counts[system_id]
                    else 0.0
                ),
            }
            for system_id in expected_system_ids
        },
    }
