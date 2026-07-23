#!/usr/bin/env python
"""Run the frozen oracle-mask coordinate repair pilot."""

import argparse
import csv
import hashlib
import json
import subprocess
import time
from argparse import Namespace
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import torch
from rdkit import Chem

from flowr.data.dataset import PocketComplexLMDBDataset
from flowr.gen.utils import load_data_from_pdb
from flowr.repair import (
    build_local_coordinate_prior,
    build_torsion_corruption,
    classify_experiment,
    copy_mol_with_coords,
    evaluate_pose,
    evaluate_repair,
    molecular_graph_signature,
    sample_oracle_repair,
    states_have_same_discrete_graph,
)
from flowr.scriptutil import complex_transform, get_n_bond_types, load_model
from flowr.util.pocket import PocketComplexBatch


EXPERIMENT_ID = "20260723-01-oracle-mask-repair-pilot"
EXPECTED_CHECKPOINT_SHA256 = (
    "b818f41dc12ffb6bc558bb0ad997055581e07cd9e49dcac1b794ed9993c46e4c"
)
FORMAL_STEPS = 100
FORMAL_BATCH_SIZE = 2
FORMAL_SEED_BASE = 2026072300
ROLLOUTS_PER_CASE = 10

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
CASES_PATH = SCRIPT_DIR / "cases.json"
OUTPUT_DIR = REPO_ROOT / "outputs" / EXPERIMENT_ID
DATASET_ROOT = REPO_ROOT / ".data" / "datasets" / "spindr"
LMDB_DIR = DATASET_ROOT / "final_from_smol"
RAW_VAL_DIR = DATASET_ROOT / "data_prepared" / "val"
CHECKPOINT_PATH = REPO_ROOT / ".data" / "checkpoints" / "flowr_root" / "flowr_root_v2.2.ckpt"

RUN_FIELDS = [
    "experiment_id",
    "stage",
    "case_index",
    "system_id",
    "rollout_index",
    "seed",
    "steps",
    "batch_size",
    "run_status",
    "error",
    "finite_coords",
    "sanitized",
    "no_protein_clash",
    "num_pairwise_clashes",
    "min_relative_distance",
    "bond_lengths_ok",
    "bond_angles_ok",
    "no_internal_clash",
    "geometry_ok",
    "graph_unchanged",
    "fixed_coords_ok",
    "max_fixed_drift",
    "editable_rmsd_to_good",
    "all_atom_rmsd_to_good",
    "success",
    "runtime_seconds",
    "peak_gpu_memory_mb",
    "sdf_path",
]
BOOL_FIELDS = {
    "finite_coords",
    "sanitized",
    "no_protein_clash",
    "bond_lengths_ok",
    "bond_angles_ok",
    "no_internal_clash",
    "geometry_ok",
    "graph_unchanged",
    "fixed_coords_ok",
    "success",
}
INT_FIELDS = {
    "case_index",
    "rollout_index",
    "seed",
    "steps",
    "batch_size",
    "num_pairwise_clashes",
}
FLOAT_FIELDS = {
    "min_relative_distance",
    "max_fixed_drift",
    "editable_rmsd_to_good",
    "all_atom_rmsd_to_good",
    "runtime_seconds",
    "peak_gpu_memory_mb",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_default(value: Any) -> Any:
    if isinstance(value, (np.bool_, np.integer, np.floating)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value)}")


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _write_csv(path: Path, records: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow({field: record.get(field, "") for field in fields})
    temporary.replace(path)


def _read_run_records(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    records = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            parsed: dict[str, Any] = dict(row)
            for field in BOOL_FIELDS:
                parsed[field] = row[field].lower() == "true" if row[field] else False
            for field in INT_FIELDS:
                parsed[field] = int(row[field]) if row[field] else -1
            for field in FLOAT_FIELDS:
                parsed[field] = float(row[field]) if row[field] else float("nan")
            records.append(parsed)
    return records


def _load_cases() -> list[dict[str, Any]]:
    payload = json.loads(CASES_PATH.read_text(encoding="utf-8"))
    if payload["experiment_id"] != EXPERIMENT_ID:
        raise ValueError("cases.json experiment_id does not match")
    cases = sorted(payload["cases"], key=lambda item: item["case_index"])
    if [case["case_index"] for case in cases] != list(range(len(cases))):
        raise ValueError("case_index values must be contiguous and zero based")
    return cases


def _load_sdf(path: Path) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
    mol = next((candidate for candidate in supplier if candidate is not None), None)
    if mol is None:
        raise ValueError(f"Could not load ligand from {path}")
    return Chem.RemoveHs(mol)


def _load_protein(path: Path) -> Chem.Mol:
    mol = Chem.MolFromPDBFile(str(path), sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not load protein from {path}")
    return mol


def _processing_args(case: dict[str, Any]) -> Namespace:
    system_id = case["system_id"]
    return Namespace(
        add_hs=False,
        add_hs_and_optimize=False,
        kekulize=False,
        use_pdbfixer=False,
        protonate_pocket=False,
        pocket_cutoff=6.0,
        cut_pocket=True,
        max_pocket_size=1000,
        min_pocket_size=10,
        compute_interactions=False,
        pocket_type="holo",
        pdb_id=system_id,
        ligand_id=system_id,
        pdb_file=str(RAW_VAL_DIR / f"{system_id}.pdb"),
        ligand_file=str(RAW_VAL_DIR / f"{system_id}.sdf"),
    )


def _prepare_case(case: dict[str, Any]) -> dict[str, Any]:
    system_id = case["system_id"]
    good_mol = _load_sdf(RAW_VAL_DIR / f"{system_id}.sdf")
    protein_mol = _load_protein(RAW_VAL_DIR / f"{system_id}.pdb")
    system = load_data_from_pdb(
        _processing_args(case), remove_hs=True, remove_aromaticity=True
    )
    if system is None:
        raise RuntimeError(f"Official processing failed for {system_id}")

    coords_good = np.asarray(good_mol.GetConformer().GetPositions(), dtype=np.float64)
    if good_mol.GetNumAtoms() != system.ligand.seq_length:
        raise ValueError(f"Raw/model atom count mismatch for {system_id}")
    raw_atomics = [atom.GetAtomicNum() for atom in good_mol.GetAtoms()]
    if raw_atomics != system.ligand.atomics.tolist():
        raise ValueError(f"Raw/model atom order mismatch for {system_id}")
    if np.max(np.abs(coords_good - system.ligand.coords.numpy())) > 3e-5:
        raise ValueError(f"Official processing changed ligand coordinates for {system_id}")

    coords_bad_tensor = build_torsion_corruption(
        torch.from_numpy(coords_good),
        axis_origin=case["axis_origin"],
        axis_target=case["axis_target"],
        editable_atom_indices=case["editable_atom_indices"],
        angle_degrees=case["angle_degrees"],
    )
    coords_bad = coords_bad_tensor.numpy()
    bad_mol = copy_mol_with_coords(good_mol, coords_bad)
    bad_ligand = system.ligand._copy_with(
        coords=torch.as_tensor(coords_bad, dtype=system.ligand.coords.dtype)
    )
    bad_system = system._copy_with(ligand=bad_ligand)
    fixed_mask = np.ones(good_mol.GetNumAtoms(), dtype=bool)
    fixed_mask[np.asarray(case["editable_atom_indices"], dtype=int)] = False
    return {
        "case": case,
        "good_mol": good_mol,
        "bad_mol": bad_mol,
        "protein_mol": protein_mol,
        "bad_system": bad_system,
        "coords_good": coords_good,
        "coords_bad": coords_bad,
        "fixed_mask": fixed_mask,
    }


def _write_mol(path: Path, mol: Chem.Mol, properties: dict[str, Any] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    mol = Chem.Mol(mol)
    for key, value in (properties or {}).items():
        mol.SetProp(str(key), str(value))
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()


def _load_lmdb_validation_systems() -> dict[str, Any]:
    dataset = PocketComplexLMDBDataset(
        root=str(LMDB_DIR),
        remove_hs=True,
        remove_aromaticity=True,
        skip_non_valid=False,
    )
    splits = np.load(LMDB_DIR / "splits.npz")
    systems = {}
    try:
        for index in splits["idx_val"]:
            system = dataset[int(index)]
            systems[system.metadata["system_id"]] = system
    finally:
        dataset._close_env()
    return systems


def run_preflight() -> dict[str, Any]:
    cases = _load_cases()
    checkpoint_sha = _sha256(CHECKPOINT_PATH)
    if checkpoint_sha != EXPECTED_CHECKPOINT_SHA256:
        raise ValueError("Checkpoint SHA256 does not match the frozen protocol")
    lmdb_systems = _load_lmdb_validation_systems()
    baseline_records = []
    inputs_dir = OUTPUT_DIR / "inputs"

    for case in cases:
        system_id = case["system_id"]
        if system_id not in lmdb_systems:
            raise ValueError(f"{system_id} is not in the official validation split")
        material = _prepare_case(case)
        raw_atomics = [atom.GetAtomicNum() for atom in material["good_mol"].GetAtoms()]
        if raw_atomics != lmdb_systems[system_id].ligand.atomics.tolist():
            raise ValueError(f"LMDB/raw atom order mismatch for {system_id}")

        axis_bond = material["good_mol"].GetBondBetweenAtoms(
            case["axis_origin"], case["axis_target"]
        )
        if axis_bond is None:
            raise ValueError(f"Frozen axis is not a bond for {system_id}")
        if molecular_graph_signature(material["good_mol"]) != molecular_graph_signature(
            material["bad_mol"]
        ):
            raise ValueError(f"Corruption changed the graph for {system_id}")

        fixed_drift = np.linalg.norm(
            material["coords_bad"][material["fixed_mask"]]
            - material["coords_good"][material["fixed_mask"]],
            axis=-1,
        ).max()
        if fixed_drift > 1e-10:
            raise ValueError(f"Corruption moved fixed atoms for {system_id}")

        good_pose = evaluate_pose(material["good_mol"], material["protein_mol"])
        bad_pose = evaluate_pose(material["bad_mol"], material["protein_mol"])
        if not good_pose["no_protein_clash"] or not good_pose["geometry_ok"]:
            raise ValueError(f"G_good baseline failed for {system_id}: {good_pose}")
        if not bad_pose["geometry_ok"]:
            raise ValueError(f"G_bad geometry failed for {system_id}: {bad_pose}")
        if bad_pose["num_pairwise_clashes"] != case["expected_pairwise_clashes"]:
            raise ValueError(f"G_bad clash count changed for {system_id}: {bad_pose}")
        if abs(
            bad_pose["min_relative_distance"]
            - case["expected_min_relative_distance"]
        ) > 0.002:
            raise ValueError(f"G_bad minimum relative distance changed for {system_id}")

        _write_mol(
            inputs_dir / f"{system_id}__good.sdf",
            material["good_mol"],
            {"system_id": system_id, "pose": "G_good"},
        )
        _write_mol(
            inputs_dir / f"{system_id}__bad.sdf",
            material["bad_mol"],
            {
                "system_id": system_id,
                "pose": "G_bad",
                "axis": f"{case['axis_origin']}->{case['axis_target']}",
                "angle_degrees": case["angle_degrees"],
            },
        )
        for pose_name, pose in (("G_good", good_pose), ("G_bad", bad_pose)):
            baseline_records.append(
                {
                    "system_id": system_id,
                    "pose": pose_name,
                    **pose,
                }
            )

    baseline_fields = ["system_id", "pose"] + list(
        key for key in baseline_records[0] if key not in {"system_id", "pose"}
    )
    _write_csv(OUTPUT_DIR / "baselines.csv", baseline_records, baseline_fields)
    result = {
        "experiment_id": EXPERIMENT_ID,
        "status": "passed",
        "checkpoint_path": str(CHECKPOINT_PATH.relative_to(REPO_ROOT)),
        "checkpoint_sha256": checkpoint_sha,
        "lmdb_validation_count": len(lmdb_systems),
        "case_count": len(cases),
        "cases": [case["system_id"] for case in cases],
    }
    _write_json(OUTPUT_DIR / "preflight.json", result)
    print(json.dumps(result, indent=2))
    return result


def _model_args(steps: int) -> Namespace:
    return Namespace(
        ckpt_path=str(CHECKPOINT_PATH),
        data_path=str(LMDB_DIR),
        save_dir=str(OUTPUT_DIR),
        integration_steps=steps,
        ode_sampling_strategy="linear",
        corrector_iters=0,
        categorical_strategy="uniform-sample",
        cat_sampling_noise_level=1,
        coord_noise_scale=0.0,
        use_sde_simulation=False,
        use_cosine_scheduler=False,
        arch="pocket",
        pocket_type="holo",
        pocket_noise="fix",
        interaction_conditional=False,
        scaffold_hopping=False,
        scaffold_elaboration=False,
        linker_inpainting=False,
        core_growing=False,
        fragment_inpainting=True,
        fragment_growing=False,
        substructure_inpainting=False,
        substructure=None,
        graph_inpainting=None,
        lora_finetuned=False,
        lora_finetuning=False,
        freeze_layers=False,
        affinity_finetuning=False,
    )


def _load_checkpoint_model(steps: int, device: str):
    args = _model_args(steps)
    result = load_model(args)
    model, hparams = result[0], result[1]
    vocabs = {
        "vocab": result[2],
        "vocab_charges": result[3],
        "vocab_hybridization": result[4],
        "vocab_aromatic": result[5],
    }
    if hparams["remove_hs"] is not True or hparams["remove_aromaticity"] is not True:
        raise ValueError("Checkpoint representation differs from the frozen protocol")
    model = model.to(device)
    model.eval()
    return model, hparams, vocabs


def _transform_bad_system(material: dict[str, Any], vocabs: dict[str, Any]):
    return complex_transform(
        material["bad_system"],
        vocab=vocabs["vocab"],
        vocab_charges=vocabs["vocab_charges"],
        vocab_hybridization=vocabs["vocab_hybridization"],
        vocab_aromatic=vocabs["vocab_aromatic"],
        n_bonds=get_n_bond_types("uniform-sample"),
        coord_std=1.0,
        pocket_noise="fix",
        pocket_noise_std=0.0,
        use_interactions=False,
        rotate_complex=False,
    )


def _complex_batch_to_dict(batch: PocketComplexBatch) -> dict[str, Any]:
    data = {
        "coords": batch.coords(state="holo").float(),
        "atomics": batch.atomics(state="holo").float(),
        "bonds": batch.adjacency(state="holo").float(),
        "charges": batch.charges(state="holo").float(),
        "atom_names": batch.atom_names(state="holo").long(),
        "res_names": batch.res_names(state="holo").long(),
        "lig_mask": batch.lig_mask(state="holo").long(),
        "pocket_mask": batch.pocket_mask(state="holo").long(),
        "mask": batch.mask.long(),
    }
    hybridization = batch.hybridization(state="holo")
    if hybridization is not None:
        data["hybridization"] = hybridization.float()
    return data


def _to_device(data: dict[str, Any], device: str) -> dict[str, Any]:
    return {
        key: value.to(device) if torch.is_tensor(value) else value
        for key, value in data.items()
    }


def _run_batch(
    model: Any,
    transformed_system: Any,
    material: dict[str, Any],
    seeds: list[int],
    steps: int,
    device: str,
    stage: str,
    batch_size_used: int,
    output_root: Path,
) -> list[dict[str, Any]]:
    batch = PocketComplexBatch([transformed_system] * len(seeds))
    combined = _to_device(_complex_batch_to_dict(batch), device)
    ligand_bad = model.builder.extract_ligand_from_complex(combined)
    pocket_data = model.builder.extract_pocket_from_complex(combined)
    fixed_single = torch.as_tensor(material["fixed_mask"], dtype=torch.bool)
    fixed_mask = fixed_single.unsqueeze(0).expand(len(seeds), -1).to(device)
    editable_single = ~fixed_single

    prior_rows = []
    for row, seed in enumerate(seeds):
        generator = torch.Generator(device="cpu").manual_seed(seed)
        prior_rows.append(
            build_local_coordinate_prior(
                ligand_bad["coords"][row].detach().cpu(),
                editable_single,
                generator,
            )
        )
    prior_coords = torch.stack(prior_rows).to(device)

    if device.startswith("cuda"):
        torch.cuda.reset_peak_memory_stats(device)
    start = time.perf_counter()
    output = sample_oracle_repair(
        model,
        ligand_bad,
        pocket_data,
        prior_coords,
        fixed_mask,
        steps=steps,
    )
    if device.startswith("cuda"):
        torch.cuda.synchronize(device)
    elapsed = time.perf_counter() - start
    peak_memory_mb = (
        torch.cuda.max_memory_allocated(device) / (1024**2)
        if device.startswith("cuda")
        else 0.0
    )
    graph_equal = states_have_same_discrete_graph(output, ligand_bad)
    center = transformed_system.com.squeeze(0).detach().cpu().numpy()
    case = material["case"]
    records = []

    for row, seed in enumerate(seeds):
        coords_centered = output["coords"][row].detach().cpu().numpy()
        coords_physical = coords_centered + center
        mol_pred = copy_mol_with_coords(material["bad_mol"], coords_physical)
        metrics = evaluate_repair(
            mol_pred=mol_pred,
            mol_bad=material["bad_mol"],
            protein_mol=material["protein_mol"],
            coords_pred=coords_physical,
            coords_bad=material["coords_bad"],
            coords_good=material["coords_good"],
            fixed_mask=material["fixed_mask"],
            discrete_state_equal=graph_equal,
        )
        rollout_index = (
            seed - FORMAL_SEED_BASE - case["case_index"] * ROLLOUTS_PER_CASE
            if stage == "formal"
            else row
        )
        sdf_path = output_root / case["system_id"] / f"seed_{seed}.sdf"
        _write_mol(
            sdf_path,
            mol_pred,
            {
                "experiment_id": EXPERIMENT_ID,
                "stage": stage,
                "system_id": case["system_id"],
                "seed": seed,
                "success": metrics["success"],
            },
        )
        records.append(
            {
                "experiment_id": EXPERIMENT_ID,
                "stage": stage,
                "case_index": case["case_index"],
                "system_id": case["system_id"],
                "rollout_index": rollout_index,
                "seed": seed,
                "steps": steps,
                "batch_size": batch_size_used,
                "run_status": "completed",
                "error": "",
                **metrics,
                "runtime_seconds": elapsed / len(seeds),
                "peak_gpu_memory_mb": peak_memory_mb,
                "sdf_path": str(sdf_path.relative_to(REPO_ROOT)),
            }
        )
    return records


def _is_oom(error: Exception) -> bool:
    return isinstance(error, torch.OutOfMemoryError) or "out of memory" in str(error).lower()


def _failure_record(
    case: dict[str, Any], seed: int, error: Exception, batch_size: int
) -> dict[str, Any]:
    rollout_index = seed - FORMAL_SEED_BASE - case["case_index"] * ROLLOUTS_PER_CASE
    return {
        "experiment_id": EXPERIMENT_ID,
        "stage": "formal",
        "case_index": case["case_index"],
        "system_id": case["system_id"],
        "rollout_index": rollout_index,
        "seed": seed,
        "steps": FORMAL_STEPS,
        "batch_size": batch_size,
        "run_status": "model_error",
        "error": f"{type(error).__name__}: {error}",
        "success": False,
        "runtime_seconds": 0.0,
        "peak_gpu_memory_mb": 0.0,
        "sdf_path": "",
    }


def run_smoke(steps: int, device: str) -> dict[str, Any]:
    if steps not in {5, 100}:
        raise ValueError("Frozen smoke supports only 5 or 100 steps")
    run_preflight()
    case = _load_cases()[0]
    material = _prepare_case(case)
    model, _, vocabs = _load_checkpoint_model(steps, device)
    transformed = _transform_bad_system(material, vocabs)
    smoke_root = OUTPUT_DIR / "smoke" / f"steps_{steps}"
    records = _run_batch(
        model,
        transformed,
        material,
        [424242, 424243],
        steps,
        device,
        "smoke",
        2,
        smoke_root,
    )
    technical_ok = all(
        record["finite_coords"]
        and record["fixed_coords_ok"]
        and record["graph_unchanged"]
        for record in records
    )
    result = {
        "experiment_id": EXPERIMENT_ID,
        "steps": steps,
        "technical_ok": technical_ok,
        "records": records,
    }
    _write_json(smoke_root / "smoke.json", result)
    if not technical_ok:
        raise RuntimeError(f"{steps}-step smoke failed technical constraints")
    print(json.dumps(result, indent=2, default=_json_default))
    return result


def run_formal(device: str) -> dict[str, Any]:
    run_preflight()
    cases = _load_cases()
    records_path = OUTPUT_DIR / "runs.csv"
    records = _read_run_records(records_path)
    completed_keys = {(record["system_id"], record["seed"]) for record in records}
    model, _, vocabs = _load_checkpoint_model(FORMAL_STEPS, device)

    for case in cases:
        material = _prepare_case(case)
        transformed = _transform_bad_system(material, vocabs)
        all_seeds = [
            FORMAL_SEED_BASE
            + case["case_index"] * ROLLOUTS_PER_CASE
            + rollout_index
            for rollout_index in range(ROLLOUTS_PER_CASE)
        ]
        pending = [
            seed for seed in all_seeds if (case["system_id"], seed) not in completed_keys
        ]
        for offset in range(0, len(pending), FORMAL_BATCH_SIZE):
            seeds = pending[offset : offset + FORMAL_BATCH_SIZE]
            try:
                new_records = _run_batch(
                    model,
                    transformed,
                    material,
                    seeds,
                    FORMAL_STEPS,
                    device,
                    "formal",
                    FORMAL_BATCH_SIZE,
                    OUTPUT_DIR / "repairs",
                )
            except Exception as error:
                if _is_oom(error):
                    if device.startswith("cuda"):
                        torch.cuda.empty_cache()
                    new_records = []
                    for seed in seeds:
                        try:
                            new_records.extend(
                                _run_batch(
                                    model,
                                    transformed,
                                    material,
                                    [seed],
                                    FORMAL_STEPS,
                                    device,
                                    "formal",
                                    1,
                                    OUTPUT_DIR / "repairs",
                                )
                            )
                        except Exception as retry_error:
                            if _is_oom(retry_error):
                                raise RuntimeError(
                                    f"Infrastructure OOM repeated for {case['system_id']} seed {seed}"
                                ) from retry_error
                            new_records.append(
                                _failure_record(case, seed, retry_error, 1)
                            )
                else:
                    new_records = [
                        _failure_record(case, seed, error, len(seeds)) for seed in seeds
                    ]

            records.extend(new_records)
            records.sort(key=lambda item: (item["case_index"], item["seed"]))
            _write_csv(records_path, records, RUN_FIELDS)
            completed_keys.update(
                (record["system_id"], record["seed"]) for record in new_records
            )
            print(
                f"Completed {len(records)}/50 formal runs: {case['system_id']} {seeds}"
            )

    return summarize_results(device=device)


def _git_commit() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
    ).strip()


def _numeric_summary(values: list[float]) -> dict[str, float] | None:
    finite = np.asarray([value for value in values if np.isfinite(value)], dtype=float)
    if finite.size == 0:
        return None
    return {
        "min": float(np.min(finite)),
        "median": float(np.median(finite)),
        "mean": float(np.mean(finite)),
        "max": float(np.max(finite)),
    }


def summarize_results(device: str = "cuda") -> dict[str, Any]:
    cases = _load_cases()
    records = _read_run_records(OUTPUT_DIR / "runs.csv")
    expected_ids = [case["system_id"] for case in cases]
    gate = classify_experiment(records, expected_ids, ROLLOUTS_PER_CASE)
    condition_fields = [
        "no_protein_clash",
        "geometry_ok",
        "fixed_coords_ok",
        "graph_unchanged",
    ]
    failure_counts = {
        field: sum(not bool(record.get(field, False)) for record in records)
        for field in condition_fields
    }
    per_case_diagnostics = {}
    for system_id in expected_ids:
        case_records = [record for record in records if record["system_id"] == system_id]
        per_case_diagnostics[system_id] = {
            "editable_rmsd_to_good": _numeric_summary(
                [record["editable_rmsd_to_good"] for record in case_records]
            ),
            "all_atom_rmsd_to_good": _numeric_summary(
                [record["all_atom_rmsd_to_good"] for record in case_records]
            ),
            "remaining_clashes": _numeric_summary(
                [float(record["num_pairwise_clashes"]) for record in case_records]
            ),
        }

    summary = {
        "experiment_id": EXPERIMENT_ID,
        "implementation_commit": _git_commit(),
        "checkpoint_sha256": EXPECTED_CHECKPOINT_SHA256,
        "protocol": {
            "steps": FORMAL_STEPS,
            "solver": "euler",
            "schedule": "linear",
            "corrector_iters": 0,
            "rollouts_per_case": ROLLOUTS_PER_CASE,
            "formal_seed_base": FORMAL_SEED_BASE,
        },
        "gate": gate,
        "run_status_counts": dict(Counter(record["run_status"] for record in records)),
        "failure_counts": failure_counts,
        "diagnostics": {
            "editable_rmsd_to_good": _numeric_summary(
                [record["editable_rmsd_to_good"] for record in records]
            ),
            "all_atom_rmsd_to_good": _numeric_summary(
                [record["all_atom_rmsd_to_good"] for record in records]
            ),
            "runtime_seconds": _numeric_summary(
                [record["runtime_seconds"] for record in records]
            ),
            "peak_gpu_memory_mb": _numeric_summary(
                [record["peak_gpu_memory_mb"] for record in records]
            ),
            "per_case": per_case_diagnostics,
        },
        "device": (
            torch.cuda.get_device_name(device)
            if device.startswith("cuda") and torch.cuda.is_available()
            else device
        ),
    }
    _write_json(OUTPUT_DIR / "summary.json", summary)
    print(json.dumps(summary, indent=2, default=_json_default))
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["preflight", "smoke", "formal", "summarize"])
    parser.add_argument("--steps", type=int, default=5, help="Smoke steps: 5 or 100")
    parser.add_argument("--device", default="cuda")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.stage == "preflight":
        run_preflight()
    elif args.stage == "smoke":
        run_smoke(args.steps, args.device)
    elif args.stage == "formal":
        run_formal(args.device)
    else:
        summarize_results(args.device)


if __name__ == "__main__":
    main()
