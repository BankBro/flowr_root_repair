#!/usr/bin/env python
"""Run the frozen official fragment-inpainting control experiment."""

import argparse
import csv
import hashlib
import importlib.util
import json
import random
import subprocess
import sys
import time
from argparse import Namespace
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import torch
from rdkit import Chem

from adapter import (
    build_exact_fragment_prior,
    fixed_first_index_map,
    restore_quantized_fixed_coordinates,
)
from flowr.gen.utils import load_util
from flowr.repair import (
    classify_experiment,
    copy_mol_with_coords,
    evaluate_inpainting_candidate,
    summarize_repair_funnel,
)
from flowr.scriptutil import load_model
from flowr.util.pocket import PocketComplex, PocketComplexBatch


EXPERIMENT_ID = "20260723-02-official-fragment-inpainting-control"
SOURCE_EXPERIMENT_ID = "20260723-01-oracle-mask-repair-pilot"
EXPECTED_CHECKPOINT_SHA256 = (
    "b818f41dc12ffb6bc558bb0ad997055581e07cd9e49dcac1b794ed9993c46e4c"
)
EXPECTED_SOURCE_HASHES = {
    "cases.json": "bd88d952b964fb3926c459307974cdc85a3f1b5d5090a22d3082111d96b0d527",
    "runs.csv": "7c8c2ffae13cbff371e06282cb2f1590350e6d04740f848ed3abdb9e3e72cd6e",
    "summary.json": "07fbdcc015730419e81c97a63d7fe7cd00cbeb445df6bf750f9b76fbeaa4948c",
}
EXPECTED_BAD_SDF_HASHES = {
    "3rog__1__1.A__1.B": "3c3abb753877de3f42ea4738db2a31dffe591e845044ea28dc81af5bcdbeb82d",
    "4bv5__2__1.B__1.D": "cdfe075981c09c63b4553fc117393aa94e68ef18573a9fd4dda329d1e3f5fcb6",
    "4f0s__1__1.A__1.B": "5e971a4f30e7b736a3e2556345542cc82ffedcda753b11ac2ae278bcc62ce523",
    "6pvz__1__1.A__1.C": "b734cfc27047edb65eced68b95ac5484440eef8f4c42898b16b5bf79a3a69cce",
    "7ddl__1__1.A__1.T": "a78c6729f2cf23663973add269dcbe2d401a7f2c2c55e6f2abdb626aff3ad9aa",
}
FORMAL_STEPS = 100
FORMAL_SEED_BASE = 2026072300
ROLLOUTS_PER_CASE = 10
FIXED_DRIFT_TOLERANCE = 1e-5
SAMPLING_IMPLEMENTATION_COMMIT = "52e114f6fcf48f2797c55279efd2391ea3b5a424"

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
SOURCE_EXPERIMENT_DIR = REPO_ROOT / "experiments" / SOURCE_EXPERIMENT_ID
SOURCE_OUTPUT_DIR = REPO_ROOT / "outputs" / SOURCE_EXPERIMENT_ID
OUTPUT_DIR = REPO_ROOT / "outputs" / EXPERIMENT_ID
CHECKPOINT_PATH = (
    REPO_ROOT / ".data" / "checkpoints" / "flowr_root" / "flowr_root_v2.2.ckpt"
)

RUN_FIELDS = [
    "experiment_id",
    "source_experiment_id",
    "method",
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
    "single_component",
    "fixed_atoms_ok",
    "fixed_bonds_ok",
    "fixed_coords_ok",
    "fixed_fragment_retained",
    "raw_max_fixed_drift",
    "no_protein_clash",
    "num_pairwise_clashes",
    "min_relative_distance",
    "bond_lengths_ok",
    "bond_angles_ok",
    "no_internal_clash",
    "geometry_ok",
    "usable_output",
    "same_molecule",
    "canonical_isomeric_smiles",
    "editable_rmsd_to_good",
    "all_atom_rmsd_to_good",
    "native_success",
    "strict_success",
    "runtime_seconds",
    "peak_gpu_memory_mb",
    "output_sha256",
    "output_to_reference",
    "sdf_path",
]
BOOL_FIELDS = {
    "finite_coords",
    "sanitized",
    "single_component",
    "fixed_atoms_ok",
    "fixed_bonds_ok",
    "fixed_coords_ok",
    "fixed_fragment_retained",
    "no_protein_clash",
    "bond_lengths_ok",
    "bond_angles_ok",
    "no_internal_clash",
    "geometry_ok",
    "usable_output",
    "same_molecule",
    "native_success",
    "strict_success",
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
    "raw_max_fixed_drift",
    "min_relative_distance",
    "editable_rmsd_to_good",
    "all_atom_rmsd_to_good",
    "runtime_seconds",
    "peak_gpu_memory_mb",
}


def _load_source_module():
    path = SOURCE_EXPERIMENT_DIR / "run.py"
    spec = importlib.util.spec_from_file_location("frozen_oracle_repair_pilot", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load source experiment runner from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


SOURCE = _load_source_module()


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


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RUN_FIELDS, lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow({field: record.get(field, "") for field in RUN_FIELDS})
    temporary.replace(path)


def _read_records(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    records = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            parsed: dict[str, Any] = dict(row)
            for field in BOOL_FIELDS:
                parsed[field] = row.get(field, "").lower() == "true"
            for field in INT_FIELDS:
                parsed[field] = int(row[field]) if row.get(field) else -1
            for field in FLOAT_FIELDS:
                parsed[field] = float(row[field]) if row.get(field) else float("nan")
            records.append(parsed)
    return records


def _load_sdf(path: Path, sanitize: bool = False) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=sanitize)
    mol = next((candidate for candidate in supplier if candidate is not None), None)
    if mol is None:
        raise ValueError(f"Could not load molecule from {path}")
    return Chem.RemoveHs(mol, sanitize=sanitize)


def _write_mol(path: Path, mol: Chem.Mol, properties: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    copied = Chem.Mol(mol)
    for key, value in properties.items():
        copied.SetProp(str(key), str(value))
    writer = Chem.SDWriter(str(path))
    writer.write(copied)
    writer.close()


def _source_file_hashes() -> dict[str, str]:
    paths = {
        "cases.json": SOURCE_EXPERIMENT_DIR / "cases.json",
        "runs.csv": SOURCE_OUTPUT_DIR / "runs.csv",
        "summary.json": SOURCE_OUTPUT_DIR / "summary.json",
    }
    observed = {name: _sha256(path) for name, path in paths.items()}
    if observed != EXPECTED_SOURCE_HASHES:
        raise ValueError(f"Frozen source hashes changed: {observed}")
    for system_id, expected in EXPECTED_BAD_SDF_HASHES.items():
        path = SOURCE_OUTPUT_DIR / "inputs" / f"{system_id}__bad.sdf"
        observed_hash = _sha256(path)
        if observed_hash != expected:
            raise ValueError(f"Frozen G_bad hash changed for {system_id}")
    return observed


def _reevaluate_coordinate_only() -> list[dict[str, Any]]:
    source_records = SOURCE._read_run_records(SOURCE_OUTPUT_DIR / "runs.csv")
    materials = {
        case["system_id"]: SOURCE._prepare_case(case) for case in SOURCE._load_cases()
    }
    records = []
    for source_record in source_records:
        material = materials[source_record["system_id"]]
        mol_pred = _load_sdf(REPO_ROOT / source_record["sdf_path"], sanitize=False)
        coords_pred = np.asarray(mol_pred.GetConformer().GetPositions(), dtype=float)
        fixed_indices = np.flatnonzero(material["fixed_mask"])
        metrics = evaluate_inpainting_candidate(
            mol_pred=mol_pred,
            mol_bad=material["bad_mol"],
            protein_mol=material["protein_mol"],
            coords_output=coords_pred,
            coords_bad=material["coords_bad"],
            coords_good=material["coords_good"],
            fixed_reference_indices=fixed_indices,
            output_to_reference=np.arange(coords_pred.shape[0]),
            raw_fixed_drift=source_record["max_fixed_drift"],
            fixed_drift_tolerance=FIXED_DRIFT_TOLERANCE,
            run_completed=source_record["run_status"] == "completed",
        )
        if metrics["native_success"] != source_record["success"]:
            raise ValueError(
                f"Common evaluator changed native result for "
                f"{source_record['system_id']} seed {source_record['seed']}"
            )
        records.append(
            {
                "experiment_id": EXPERIMENT_ID,
                "source_experiment_id": SOURCE_EXPERIMENT_ID,
                "method": "coordinate_only",
                "stage": "formal",
                "case_index": source_record["case_index"],
                "system_id": source_record["system_id"],
                "rollout_index": source_record["rollout_index"],
                "seed": source_record["seed"],
                "steps": source_record["steps"],
                "batch_size": source_record["batch_size"],
                "run_status": source_record["run_status"],
                "error": source_record["error"],
                **metrics,
                "runtime_seconds": source_record["runtime_seconds"],
                "peak_gpu_memory_mb": source_record["peak_gpu_memory_mb"],
                "output_sha256": _sha256(REPO_ROOT / source_record["sdf_path"]),
                "output_to_reference": json.dumps(
                    list(range(coords_pred.shape[0])), separators=(",", ":")
                ),
                "sdf_path": source_record["sdf_path"],
            }
        )

    native_count = sum(record["native_success"] for record in records)
    strict_count = sum(record["strict_success"] for record in records)
    if (native_count, strict_count) != (14, 13):
        raise ValueError(
            f"Common evaluator expected coordinate-only 14/13, got "
            f"{native_count}/{strict_count}"
        )
    _write_csv(OUTPUT_DIR / "coordinate_only_runs.csv", records)
    return records


def run_preflight() -> dict[str, Any]:
    checkpoint_hash = _sha256(CHECKPOINT_PATH)
    if checkpoint_hash != EXPECTED_CHECKPOINT_SHA256:
        raise ValueError("Checkpoint SHA256 does not match the frozen protocol")
    source_hashes = _source_file_hashes()
    cases = SOURCE._load_cases()
    for case in cases:
        material = SOURCE._prepare_case(case)
        mapping = fixed_first_index_map(material["fixed_mask"])
        expected = np.concatenate(
            [
                np.flatnonzero(material["fixed_mask"]),
                np.asarray(case["editable_atom_indices"], dtype=int),
            ]
        )
        if not np.array_equal(mapping, expected):
            raise ValueError(f"Exact mask index map changed for {case['system_id']}")

    coordinate_records = _reevaluate_coordinate_only()
    result = {
        "experiment_id": EXPERIMENT_ID,
        "status": "passed",
        "checkpoint_sha256": checkpoint_hash,
        "source_hashes": source_hashes,
        "bad_sdf_hashes": EXPECTED_BAD_SDF_HASHES,
        "case_count": len(cases),
        "coordinate_only_reproduced": {
            "runs": len(coordinate_records),
            "native_successes": sum(
                record["native_success"] for record in coordinate_records
            ),
            "strict_successes": sum(
                record["strict_success"] for record in coordinate_records
            ),
        },
    }
    _write_json(OUTPUT_DIR / "preflight.json", result)
    print(json.dumps(result, indent=2, default=_json_default))
    return result


def _model_args(steps: int) -> Namespace:
    args = SOURCE._model_args(steps)
    additions = {
        "solver": "euler",
        "pocket_coord_noise_std": 0.0,
        "pocket_time": None,
        "interaction_time": None,
        "separate_pocket_interpolation": False,
        "separate_interaction_interpolation": False,
        "prior_center_file": None,
        "max_fragment_cuts": 3,
        "dataset": "spindr",
        "sample_mol_sizes": False,
        "rotation_alignment": False,
        "permutation_alignment": False,
        "anisotropic_prior": False,
        "ref_ligand_com_prior": False,
        "ref_ligand_com_noise_std": 0.05,
        "virtual_atom_p": 0.0,
        "final_inpaint": True,
    }
    for key, value in additions.items():
        setattr(args, key, value)
    return args


def _load_official_components(steps: int, device: str):
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
        raise ValueError("Checkpoint representation differs from frozen source inputs")
    transform, interpolant = load_util(
        args,
        hparams,
        vocabs["vocab"],
        vocabs["vocab_charges"],
        vocabs["vocab_hybridization"],
        vocabs["vocab_aromatic"],
    )
    model = model.to(device)
    model.eval()
    return args, model, transform, interpolant


def _complex_batch_to_dict(batch: PocketComplexBatch, state: str) -> dict[str, Any]:
    data = {
        "coords": batch.coords(state=state).float(),
        "atomics": batch.atomics(state=state).float(),
        "bonds": batch.adjacency(state=state).float(),
        "interactions": batch.interactions(state=state),
        "charges": batch.charges(state=state).float(),
        "atom_names": batch.atom_names(state=state).long(),
        "res_names": batch.res_names(state=state).long(),
        "lig_mask": batch.lig_mask(state=state).long(),
        "pocket_mask": batch.pocket_mask(state=state).long(),
        "fragment_mask": batch.fragment_mask(),
        "fragment_mode": batch.fragment_mode(),
        "mask": batch.mask.long(),
        "complex": batch._systems,
    }
    hybridization = batch.hybridization(state=state)
    if hybridization is not None:
        data["hybridization"] = hybridization.float()
    return data


def _to_device(data: dict[str, Any], device: str) -> dict[str, Any]:
    return {
        key: value.to(device) if torch.is_tensor(value) else value
        for key, value in data.items()
    }


def _reset_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed % (2**32))
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def _output_hash(output: dict[str, Any]) -> str:
    digest = hashlib.sha256()
    for key in ("coords", "atomics", "charges", "hybridization", "bonds", "mask"):
        value = output.get(key)
        if not torch.is_tensor(value):
            continue
        tensor = value.detach().cpu().contiguous()
        digest.update(key.encode("utf-8"))
        digest.update(str(tensor.dtype).encode("ascii"))
        digest.update(np.asarray(tensor.shape, dtype=np.int64).tobytes())
        digest.update(tensor.numpy().tobytes())
    return digest.hexdigest()


def _run_one(
    args: Namespace,
    model: Any,
    transform: Any,
    interpolant: Any,
    material: dict[str, Any],
    seed: int,
    steps: int,
    device: str,
    stage: str,
    output_root: Path,
    artifact_label: str | None = None,
) -> dict[str, Any]:
    _reset_seed(seed)
    if device.startswith("cuda"):
        torch.cuda.reset_peak_memory_stats(device)
    start = time.perf_counter()

    transformed = transform(material["bad_system"])
    prior_ligand, output_to_reference = build_exact_fragment_prior(
        interpolant,
        transformed.ligand,
        material["bad_mol"],
        material["fixed_mask"],
    )
    prior_system = PocketComplex(
        ligand=prior_ligand,
        apo=transformed.holo,
        interactions=transformed.interactions,
        metadata=transformed.metadata,
        fragment_mask=prior_ligand.fragment_mask,
        fragment_mode=prior_ligand.fragment_mode,
        com=transformed.com,
    )
    prior_data = _to_device(
        _complex_batch_to_dict(PocketComplexBatch([prior_system]), state="apo"),
        device,
    )
    posterior_data = _to_device(
        _complex_batch_to_dict(PocketComplexBatch([transformed]), state="holo"),
        device,
    )

    ligand_prior = model.builder.extract_ligand_from_complex(prior_data)
    ligand_prior["interactions"] = prior_data["interactions"]
    ligand_prior["fragment_mask"] = prior_data["fragment_mask"]
    ligand_prior["fragment_mode"] = prior_data["fragment_mode"]
    pocket_data = model.builder.extract_pocket_from_complex(posterior_data)
    pocket_data["interactions"] = posterior_data["interactions"]
    pocket_data["complex"] = posterior_data["complex"]
    times = [
        torch.zeros(1, device=device),
        torch.zeros(1, device=device),
        torch.zeros(1, device=device),
    ]
    output = model._generate(
        ligand_prior,
        pocket_data,
        steps=steps,
        times=times,
        strategy=args.ode_sampling_strategy,
        solver=args.solver,
        corr_iters=args.corrector_iters,
        save_traj=False,
        final_inpaint=True,
        apply_guidance=False,
        coord_noise_level=0.0,
    )
    if device.startswith("cuda"):
        torch.cuda.synchronize(device)
    runtime = time.perf_counter() - start
    peak_memory = (
        torch.cuda.max_memory_allocated(device) / (1024**2)
        if device.startswith("cuda")
        else 0.0
    )

    coords_raw = output["coords"][0].detach().cpu().numpy().astype(np.float64)
    coords_restored, raw_drift, _ = restore_quantized_fixed_coordinates(
        coords_raw,
        material["coords_bad"],
        material["fixed_mask"],
        output_to_reference,
        tolerance=FIXED_DRIFT_TOLERANCE,
    )
    raw_mol = model._generate_mols(output, sanitise=False)[0]
    mol_pred = None
    if raw_mol is not None and raw_mol.GetNumAtoms() == coords_restored.shape[0]:
        mol_pred = copy_mol_with_coords(raw_mol, coords_restored)

    fixed_indices = np.flatnonzero(material["fixed_mask"])
    metrics = evaluate_inpainting_candidate(
        mol_pred=mol_pred,
        mol_bad=material["bad_mol"],
        protein_mol=material["protein_mol"],
        coords_output=coords_restored,
        coords_bad=material["coords_bad"],
        coords_good=material["coords_good"],
        fixed_reference_indices=fixed_indices,
        output_to_reference=output_to_reference,
        raw_fixed_drift=raw_drift,
        fixed_drift_tolerance=FIXED_DRIFT_TOLERANCE,
        run_completed=True,
    )

    case = material["case"]
    rollout_index = (
        seed - FORMAL_SEED_BASE - case["case_index"] * ROLLOUTS_PER_CASE
        if stage == "formal"
        else -1
    )
    filename = artifact_label or f"seed_{seed}"
    sdf_path = output_root / case["system_id"] / f"{filename}.sdf"
    artifact_error = ""
    if mol_pred is not None:
        try:
            _write_mol(
                sdf_path,
                mol_pred,
                {
                    "experiment_id": EXPERIMENT_ID,
                    "method": "official_fragment_inpainting",
                    "stage": stage,
                    "system_id": case["system_id"],
                    "seed": seed,
                    "native_success": metrics["native_success"],
                    "strict_success": metrics["strict_success"],
                    "raw_max_fixed_drift": raw_drift,
                    "output_to_reference": json.dumps(
                        output_to_reference.tolist(), separators=(",", ":")
                    ),
                },
            )
        except Exception as error:
            artifact_error = f"SDF write failed: {type(error).__name__}: {error}"
            sdf_path = Path()
    else:
        sdf_path = Path()

    return {
        "experiment_id": EXPERIMENT_ID,
        "source_experiment_id": SOURCE_EXPERIMENT_ID,
        "method": "official_fragment_inpainting",
        "stage": stage,
        "case_index": case["case_index"],
        "system_id": case["system_id"],
        "rollout_index": rollout_index,
        "seed": seed,
        "steps": steps,
        "batch_size": 1,
        "run_status": "completed",
        "error": artifact_error,
        **metrics,
        "runtime_seconds": runtime,
        "peak_gpu_memory_mb": peak_memory,
        "output_sha256": _output_hash(output),
        "output_to_reference": json.dumps(
            output_to_reference.tolist(), separators=(",", ":")
        ),
        "sdf_path": str(sdf_path.relative_to(REPO_ROOT)) if sdf_path.parts else "",
    }


def _failure_record(
    case: dict[str, Any], seed: int, error: Exception, steps: int, stage: str
) -> dict[str, Any]:
    rollout_index = (
        seed - FORMAL_SEED_BASE - case["case_index"] * ROLLOUTS_PER_CASE
        if stage == "formal"
        else -1
    )
    return {
        "experiment_id": EXPERIMENT_ID,
        "source_experiment_id": SOURCE_EXPERIMENT_ID,
        "method": "official_fragment_inpainting",
        "stage": stage,
        "case_index": case["case_index"],
        "system_id": case["system_id"],
        "rollout_index": rollout_index,
        "seed": seed,
        "steps": steps,
        "batch_size": 1,
        "run_status": "model_error",
        "error": f"{type(error).__name__}: {error}",
        "native_success": False,
        "strict_success": False,
        "runtime_seconds": 0.0,
        "peak_gpu_memory_mb": 0.0,
        "output_to_reference": "",
        "sdf_path": "",
    }


def _is_oom(error: Exception) -> bool:
    return isinstance(error, torch.OutOfMemoryError) or "out of memory" in str(error).lower()


def run_smoke(steps: int, device: str) -> dict[str, Any]:
    if steps not in {5, 100}:
        raise ValueError("Frozen smoke supports only 5 or 100 steps")
    run_preflight()
    case = SOURCE._load_cases()[0]
    material = SOURCE._prepare_case(case)
    args, model, transform, interpolant = _load_official_components(steps, device)
    smoke_root = OUTPUT_DIR / "smoke" / f"steps_{steps}"
    records = [
        _run_one(
            args,
            model,
            transform,
            interpolant,
            material,
            seed,
            steps,
            device,
            "smoke",
            smoke_root,
            artifact_label=label,
        )
        for seed, label in (
            (424242, "seed_424242_first"),
            (424243, "seed_424243"),
            (424242, "seed_424242_repeat"),
        )
    ]
    reproducible = records[0]["output_sha256"] == records[2]["output_sha256"]
    technical_ok = bool(
        reproducible
        and all(
            record["run_status"] == "completed"
            and record["finite_coords"]
            and record["fixed_coords_ok"]
            for record in records
        )
    )
    result = {
        "experiment_id": EXPERIMENT_ID,
        "steps": steps,
        "technical_ok": technical_ok,
        "reproducible_same_seed": reproducible,
        "records": records,
    }
    _write_json(smoke_root / "smoke.json", result)
    if not technical_ok:
        raise RuntimeError(f"{steps}-step official smoke failed")
    print(json.dumps(result, indent=2, default=_json_default))
    return result


def run_formal(device: str) -> dict[str, Any]:
    run_preflight()
    cases = SOURCE._load_cases()
    records_path = OUTPUT_DIR / "official_runs.csv"
    records = _read_records(records_path)
    completed_keys = {(record["system_id"], record["seed"]) for record in records}
    args, model, transform, interpolant = _load_official_components(FORMAL_STEPS, device)

    for case in cases:
        material = SOURCE._prepare_case(case)
        for rollout_index in range(ROLLOUTS_PER_CASE):
            seed = FORMAL_SEED_BASE + case["case_index"] * 10 + rollout_index
            if (case["system_id"], seed) in completed_keys:
                continue
            try:
                record = _run_one(
                    args,
                    model,
                    transform,
                    interpolant,
                    material,
                    seed,
                    FORMAL_STEPS,
                    device,
                    "formal",
                    OUTPUT_DIR / "official_repairs",
                )
            except Exception as error:
                if _is_oom(error):
                    if device.startswith("cuda"):
                        torch.cuda.empty_cache()
                    try:
                        record = _run_one(
                            args,
                            model,
                            transform,
                            interpolant,
                            material,
                            seed,
                            FORMAL_STEPS,
                            device,
                            "formal",
                            OUTPUT_DIR / "official_repairs",
                        )
                    except Exception as retry_error:
                        if _is_oom(retry_error):
                            raise RuntimeError(
                                f"Repeated batch-size-1 OOM for {case['system_id']} seed {seed}"
                            ) from retry_error
                        record = _failure_record(
                            case, seed, retry_error, FORMAL_STEPS, "formal"
                        )
                else:
                    record = _failure_record(case, seed, error, FORMAL_STEPS, "formal")

            records.append(record)
            records.sort(key=lambda item: (item["case_index"], item["seed"]))
            _write_csv(records_path, records)
            completed_keys.add((case["system_id"], seed))
            print(
                f"Completed {len(records)}/50 official runs: "
                f"{case['system_id']} seed {seed}"
            )
    return summarize_results(device)


def _gate(records: list[dict[str, Any]], field: str) -> dict[str, Any]:
    case_ids = [case["system_id"] for case in SOURCE._load_cases()]
    routed = [{**record, "success": bool(record.get(field, False))} for record in records]
    return classify_experiment(routed, case_ids, ROLLOUTS_PER_CASE)


def _per_case(records: list[dict[str, Any]]) -> dict[str, dict[str, int]]:
    result = {}
    for case in SOURCE._load_cases():
        system_id = case["system_id"]
        selected = [record for record in records if record["system_id"] == system_id]
        result[system_id] = {
            "runs": len(selected),
            "native_successes": sum(record["native_success"] for record in selected),
            "strict_successes": sum(record["strict_success"] for record in selected),
        }
    return result


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


def _reevaluate_official_artifacts(
    records: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    materials = {
        case["system_id"]: SOURCE._prepare_case(case) for case in SOURCE._load_cases()
    }
    boolean_mismatches = []
    audited = []
    comparison_fields = sorted(BOOL_FIELDS)
    for record in records:
        if record["run_status"] != "completed" or not record.get("sdf_path"):
            audited.append(record)
            continue
        material = materials[record["system_id"]]
        sdf_path = REPO_ROOT / record["sdf_path"]
        mol_pred = _load_sdf(sdf_path, sanitize=False)
        coords_pred = np.asarray(mol_pred.GetConformer().GetPositions(), dtype=float)
        output_to_reference = np.asarray(
            json.loads(record["output_to_reference"]), dtype=int
        )
        metrics = evaluate_inpainting_candidate(
            mol_pred=mol_pred,
            mol_bad=material["bad_mol"],
            protein_mol=material["protein_mol"],
            coords_output=coords_pred,
            coords_bad=material["coords_bad"],
            coords_good=material["coords_good"],
            fixed_reference_indices=np.flatnonzero(material["fixed_mask"]),
            output_to_reference=output_to_reference,
            raw_fixed_drift=record["raw_max_fixed_drift"],
            fixed_drift_tolerance=FIXED_DRIFT_TOLERANCE,
            run_completed=True,
        )
        changed = [
            field
            for field in comparison_fields
            if bool(record.get(field, False)) != bool(metrics.get(field, False))
        ]
        if changed:
            boolean_mismatches.append(
                {
                    "system_id": record["system_id"],
                    "seed": record["seed"],
                    "fields": changed,
                }
            )
        audited.append({**record, **metrics})

    _write_csv(OUTPUT_DIR / "official_runs.csv", audited)
    audit = {
        "experiment_id": EXPERIMENT_ID,
        "status": "passed",
        "audited_completed_sdfs": sum(
            record["run_status"] == "completed" and bool(record.get("sdf_path"))
            for record in audited
        ),
        "boolean_mismatches_from_prewrite_evaluation": boolean_mismatches,
        "native_successes_after_artifact_evaluation": sum(
            record["native_success"] for record in audited
        ),
        "strict_successes_after_artifact_evaluation": sum(
            record["strict_success"] for record in audited
        ),
    }
    _write_json(OUTPUT_DIR / "artifact_audit.json", audit)
    return audited


def summarize_results(device: str = "cuda") -> dict[str, Any]:
    coordinate_records = _read_records(OUTPUT_DIR / "coordinate_only_runs.csv")
    if not coordinate_records:
        run_preflight()
        coordinate_records = _read_records(OUTPUT_DIR / "coordinate_only_runs.csv")
    official_records = _read_records(OUTPUT_DIR / "official_runs.csv")
    if official_records:
        official_records = _reevaluate_official_artifacts(official_records)
    coordinate_native = sum(record["native_success"] for record in coordinate_records)
    coordinate_strict = sum(record["strict_success"] for record in coordinate_records)
    official_native = sum(record["native_success"] for record in official_records)
    official_strict = sum(record["strict_success"] for record in official_records)

    comparison = {
        "experiment_id": EXPERIMENT_ID,
        "denominator": 50,
        "coordinate_only": {
            "native_successes": coordinate_native,
            "strict_successes": coordinate_strict,
            "per_case": _per_case(coordinate_records),
            "funnel": summarize_repair_funnel(coordinate_records),
        },
        "official_fragment_inpainting": {
            "native_successes": official_native,
            "strict_successes": official_strict,
            "per_case": _per_case(official_records),
            "funnel": summarize_repair_funnel(official_records),
        },
        "official_minus_coordinate_only_percentage_points": {
            "native": (official_native - coordinate_native) * 2.0,
            "strict": (official_strict - coordinate_strict) * 2.0,
        },
    }
    summary = {
        "experiment_id": EXPERIMENT_ID,
        "sampling_implementation_commit": SAMPLING_IMPLEMENTATION_COMMIT,
        "evaluation_implementation_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
        ).strip(),
        "checkpoint_sha256": EXPECTED_CHECKPOINT_SHA256,
        "protocol": {
            "method": "official_fragment_inpainting",
            "exact_oracle_mask": True,
            "steps": FORMAL_STEPS,
            "solver": "euler",
            "schedule": "linear",
            "corrector_iters": 0,
            "final_inpaint": True,
            "batch_size": 1,
            "rollouts_per_case": ROLLOUTS_PER_CASE,
            "formal_seed_base": FORMAL_SEED_BASE,
        },
        "native_gate_descriptive": _gate(official_records, "native_success"),
        "strict_gate": _gate(official_records, "strict_success"),
        "run_status_counts": dict(
            Counter(record["run_status"] for record in official_records)
        ),
        "diagnostics": {
            "runtime_seconds": _numeric_summary(
                [record["runtime_seconds"] for record in official_records]
            ),
            "peak_gpu_memory_mb": _numeric_summary(
                [record["peak_gpu_memory_mb"] for record in official_records]
            ),
            "raw_max_fixed_drift": _numeric_summary(
                [record["raw_max_fixed_drift"] for record in official_records]
            ),
            "editable_rmsd_to_good": _numeric_summary(
                [record["editable_rmsd_to_good"] for record in official_records]
            ),
        },
        "device": (
            torch.cuda.get_device_name(device)
            if device.startswith("cuda") and torch.cuda.is_available()
            else device
        ),
    }
    _write_json(OUTPUT_DIR / "comparison.json", comparison)
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
