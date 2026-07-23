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
