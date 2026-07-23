#!/usr/bin/env python
"""Run the frozen SPINDR full-test oracle repair benchmark."""

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import random
import subprocess
import sys
import time
from argparse import Namespace
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import posebusters
import torch
from posebusters import PoseBusters
from rdkit import Chem, rdBase

from flowr.data.dataset import PocketComplexLMDBDataset
from flowr.gen.utils import load_util
from flowr.repair import (
    build_exact_fragment_prior,
    build_local_coordinate_prior,
    copy_mol_with_coords,
    restore_quantized_fixed_coordinates,
    sample_oracle_repair,
)
from flowr.scriptutil import complex_transform, get_n_bond_types, load_model
from flowr.util.pocket import PocketComplex, PocketComplexBatch

from protocol import (
    CONSTRUCTION_FIELDS,
    EVALUATION_FIELDS,
    EXPERIMENT_ID,
    FIXED_DRIFT_TOLERANCE,
    PB_FIELDS,
    construct_clash_case,
    evaluate_output,
    load_sdf,
    pb_results,
    pb_valid,
    pocket_mol_from_cif,
    sha256,
)


EXPECTED_CHECKPOINT_SHA256 = (
    "b818f41dc12ffb6bc558bb0ad997055581e07cd9e49dcac1b794ed9993c46e4c"
)
FORMAL_STEPS = 100
ROLLOUTS_PER_CASE = 10
FORMAL_SEED_BASE = 2026072300
BOOTSTRAP_SEED = 2026072304
BOOTSTRAP_REPLICATES = 10_000

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATASET_ROOT = REPO_ROOT / ".data" / "datasets" / "spindr"
LMDB_DIR = DATASET_ROOT / "final_from_smol"
PREPARED_ROOT = DATASET_ROOT / "data_prepared"
OUTPUT_DIR = REPO_ROOT / "outputs" / EXPERIMENT_ID
CASES_PATH = OUTPUT_DIR / "cases.json"
CHECKPOINT_PATH = (
    REPO_ROOT / ".data" / "checkpoints" / "flowr_root" / "flowr_root_v2.2.ckpt"
)
POSEBUSTERS_CONFIG = REPO_ROOT / "posebusters" / "config" / "dock.yml"

SAMPLING_FIELDS = [
    "experiment_id",
    "method",
    "stage",
    "test_index",
    "system_id",
    "rollout_index",
    "seed",
    "steps",
    "batch_size",
    "run_status",
    "error",
    "runtime_seconds",
    "peak_gpu_memory_mb",
    "raw_max_fixed_drift",
    "output_to_reference",
    "output_sha256",
    "sdf_path",
]
RUN_FIELDS = [*SAMPLING_FIELDS, *EVALUATION_FIELDS]
BASELINE_FIELDS = [
    "experiment_id",
    "method",
    "test_index",
    "system_id",
    *EVALUATION_FIELDS,
]
CASE_RATE_FIELDS = [
    "experiment_id",
    "test_index",
    "system_id",
    "editable_atom_count",
    "num_pairwise_clashes",
    "min_relative_distance",
    "coordinate_official_successes",
    "coordinate_official_rate",
    "inpainting_official_successes",
    "inpainting_official_rate",
    "coordinate_strict_successes",
    "coordinate_strict_rate",
    "inpainting_strict_successes",
    "inpainting_strict_rate",
]


def _relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT.resolve()))


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.bool_, np.integer, np.floating)):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(
            _json_safe(payload),
            indent=2,
            sort_keys=True,
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _csv_value(value: Any) -> Any:
    if isinstance(value, (np.bool_, np.integer, np.floating)):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    return value


def _write_csv(path: Path, records: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow(
                {field: _csv_value(record.get(field, "")) for field in fields}
            )
    temporary.replace(path)


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_mol(path: Path, mol: Chem.Mol, properties: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    copied = Chem.Mol(mol)
    for key, value in properties.items():
        copied.SetProp(str(key), str(value))
    writer = Chem.SDWriter(str(path))
    writer.write(copied)
    writer.close()


def _load_test_systems() -> list[tuple[int, Any]]:
    dataset = PocketComplexLMDBDataset(
        root=str(LMDB_DIR),
        remove_hs=True,
        remove_aromaticity=True,
        skip_non_valid=False,
    )
    splits = np.load(LMDB_DIR / "splits.npz")
    systems = []
    try:
        for test_index, dataset_index in enumerate(splits["idx_test"]):
            systems.append((test_index, dataset[int(dataset_index)]))
    finally:
        dataset._close_env()
    if len(systems) != 225:
        raise RuntimeError(f"Expected 225 test systems, found {len(systems)}")
    ids = [str(system.metadata["system_id"]) for _, system in systems]
    if len(ids) != len(set(ids)):
        raise RuntimeError("Official test split contains duplicate system IDs")
    return systems


def _system_map() -> dict[str, Any]:
    return {
        str(system.metadata["system_id"]): system
        for _, system in _load_test_systems()
    }


def _checkpoint_hash() -> str:
    observed = sha256(CHECKPOINT_PATH)
    if observed != EXPECTED_CHECKPOINT_SHA256:
        raise ValueError("Checkpoint SHA256 does not match the frozen protocol")
    return observed


def _version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        if name == "posebusters":
            return str(getattr(posebusters, "__version__", "unknown"))
        return "unknown"


def run_preflight() -> dict[str, Any]:
    systems = _load_test_systems()
    missing = []
    input_entries = []
    for test_index, system in systems:
        system_id = str(system.metadata["system_id"])
        for kind, path in (
            ("ligand_sdf", PREPARED_ROOT / f"{system_id}.sdf"),
            ("pocket_cif", PREPARED_ROOT / f"{system_id}.cif"),
        ):
            if not path.is_file():
                missing.append(_relative(path))
            else:
                input_entries.append(
                    {
                        "test_index": test_index,
                        "system_id": system_id,
                        "kind": kind,
                        "path": _relative(path),
                        "bytes": path.stat().st_size,
                        "sha256": sha256(path),
                    }
                )
    if missing:
        raise FileNotFoundError(f"Missing official test inputs: {missing}")

    manifest = {"experiment_id": EXPERIMENT_ID, "entries": input_entries}
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    _write_json(OUTPUT_DIR / "input_manifest.json", manifest)
    result = {
        "experiment_id": EXPERIMENT_ID,
        "status": "passed",
        "official_test_entries": len(systems),
        "ligand_sdfs": sum(e["kind"] == "ligand_sdf" for e in input_entries),
        "pocket_cifs": sum(e["kind"] == "pocket_cif" for e in input_entries),
        "checkpoint_sha256": _checkpoint_hash(),
        "posebusters_config_sha256": sha256(POSEBUSTERS_CONFIG),
        "input_manifest_sha256": sha256(OUTPUT_DIR / "input_manifest.json"),
        "posebusters_version": _version("posebusters"),
        "rdkit_version": rdBase.rdkitVersion,
    }
    _write_json(OUTPUT_DIR / "preflight.json", result)
    print(json.dumps(result, indent=2))
    return result


def run_construct() -> dict[str, Any]:
    run_preflight()
    systems = _load_test_systems()
    buster = PoseBusters(config="dock")
    records: list[dict[str, Any]] = []
    cases = []
    inputs_dir = OUTPUT_DIR / "inputs"
    for position, (test_index, system) in enumerate(systems, start=1):
        system_id = str(system.metadata["system_id"])
        source_sdf = PREPARED_ROOT / f"{system_id}.sdf"
        source_cif = PREPARED_ROOT / f"{system_id}.cif"
        record, material = construct_clash_case(
            test_index=test_index,
            system=system,
            source_sdf=source_sdf,
            source_cif=source_cif,
            buster=buster,
        )
        record["source_sdf_path"] = _relative(source_sdf)
        record["source_cif_path"] = _relative(source_cif)
        if material is not None:
            good_path = inputs_dir / f"{system_id}__good.sdf"
            bad_path = inputs_dir / f"{system_id}__bad.sdf"
            properties = {
                "experiment_id": EXPERIMENT_ID,
                "test_index": test_index,
                "system_id": system_id,
            }
            _write_mol(good_path, material["good_mol"], properties | {"pose": "G_good"})
            _write_mol(
                bad_path,
                material["bad_mol"],
                properties
                | {
                    "pose": "G_bad",
                    "axis": f"{record['axis_origin']}->{record['axis_target']}",
                    "angle_degrees": record["angle_degrees"],
                    "editable_atom_indices": record["editable_atom_indices"],
                },
            )
            record.update(
                {
                    "good_sdf_path": _relative(good_path),
                    "bad_sdf_path": _relative(bad_path),
                    "good_sdf_sha256": sha256(good_path),
                    "bad_sdf_sha256": sha256(bad_path),
                }
            )
            cases.append(
                {
                    key: record[key]
                    for key in (
                        "test_index",
                        "system_id",
                        "axis_origin",
                        "axis_target",
                        "angle_degrees",
                        "editable_atom_count",
                        "editable_atom_indices",
                        "fixed_atom_count",
                        "num_pairwise_clashes",
                        "min_relative_distance",
                        "good_sdf_path",
                        "bad_sdf_path",
                        "good_sdf_sha256",
                        "bad_sdf_sha256",
                        "source_cif_path",
                        "source_cif_sha256",
                    )
                }
            )
        records.append(record)
        _write_csv(OUTPUT_DIR / "construction.csv", records, CONSTRUCTION_FIELDS)
        if position % 10 == 0 or position == len(systems):
            print(
                f"Constructed {position}/{len(systems)}; eligible={len(cases)}",
                flush=True,
            )

    payload = {
        "experiment_id": EXPERIMENT_ID,
        "protocol": {
            "test_entries": 225,
            "angle_grid_degrees": [30, -30, 60, -60, 90, -90, 120, -120, 150, -150, 180],
            "clash_count_interval": [1, 4],
            "relative_distance_interval": [0.50, 0.75],
            "relative_distance_target": 0.625,
        },
        "eligible_cases": cases,
    }
    _write_json(CASES_PATH, payload)
    status_counts = dict(Counter(record["status"] for record in records))
    summary = {
        "experiment_id": EXPERIMENT_ID,
        "status": "completed",
        "test_entries": len(records),
        "eligible_cases": len(cases),
        "status_counts": status_counts,
        "cases_sha256": sha256(CASES_PATH),
        "construction_sha256": sha256(OUTPUT_DIR / "construction.csv"),
    }
    _write_json(OUTPUT_DIR / "construction_summary.json", summary)
    print(json.dumps(summary, indent=2))
    return summary


def _load_cases() -> list[dict[str, Any]]:
    if not CASES_PATH.is_file():
        raise FileNotFoundError("Run construct before loading eligible cases")
    payload = json.loads(CASES_PATH.read_text(encoding="utf-8"))
    if payload.get("experiment_id") != EXPERIMENT_ID:
        raise ValueError("cases.json experiment ID does not match")
    cases = sorted(payload["eligible_cases"], key=lambda case: case["test_index"])
    if not cases:
        raise RuntimeError("Construction produced no eligible cases")
    return cases


def _load_material(case: dict[str, Any], system: Any) -> dict[str, Any]:
    good_path = REPO_ROOT / case["good_sdf_path"]
    bad_path = REPO_ROOT / case["bad_sdf_path"]
    if sha256(good_path) != case["good_sdf_sha256"]:
        raise ValueError(f"G_good hash changed for {case['system_id']}")
    if sha256(bad_path) != case["bad_sdf_sha256"]:
        raise ValueError(f"G_bad hash changed for {case['system_id']}")
    good_mol = load_sdf(good_path, sanitize=True)
    bad_mol = load_sdf(bad_path, sanitize=True)
    protein_mol = pocket_mol_from_cif(REPO_ROOT / case["source_cif_path"])
    coords_good = np.asarray(good_mol.GetConformer().GetPositions(), dtype=np.float64)
    coords_bad = np.asarray(bad_mol.GetConformer().GetPositions(), dtype=np.float64)
    editable = np.asarray(
        [int(value) for value in case["editable_atom_indices"].split(",")],
        dtype=int,
    )
    fixed_mask = np.ones(good_mol.GetNumAtoms(), dtype=bool)
    fixed_mask[editable] = False
    bad_ligand = system.ligand._copy_with(
        coords=torch.as_tensor(coords_bad, dtype=system.ligand.coords.dtype)
    )
    bad_system = system._copy_with(ligand=bad_ligand)
    return {
        "case": case,
        "test_index": int(case["test_index"]),
        "system_id": str(case["system_id"]),
        "good_mol": good_mol,
        "bad_mol": bad_mol,
        "protein_mol": protein_mol,
        "coords_good": coords_good,
        "coords_bad": coords_bad,
        "editable_atom_indices": editable.tolist(),
        "fixed_mask": fixed_mask,
        "bad_system": bad_system,
    }


def _model_args(steps: int, *, official: bool) -> Namespace:
    args = Namespace(
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
    if official:
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


def _load_coordinate_components(steps: int, device: str):
    result = load_model(_model_args(steps, official=False))
    model, hparams = result[0], result[1]
    if hparams["remove_hs"] is not True or hparams["remove_aromaticity"] is not True:
        raise ValueError("Checkpoint representation differs from frozen protocol")
    vocabs = {
        "vocab": result[2],
        "vocab_charges": result[3],
        "vocab_hybridization": result[4],
        "vocab_aromatic": result[5],
    }
    model = model.to(device)
    model.eval()
    return model, vocabs


def _load_official_components(steps: int, device: str):
    args = _model_args(steps, official=True)
    result = load_model(args)
    model, hparams = result[0], result[1]
    if hparams["remove_hs"] is not True or hparams["remove_aromaticity"] is not True:
        raise ValueError("Checkpoint representation differs from frozen protocol")
    transform, interpolant = load_util(
        args,
        hparams,
        result[2],
        result[3],
        result[4],
        result[5],
    )
    model = model.to(device)
    model.eval()
    return args, model, transform, interpolant


def _coordinate_complex_dict(batch: PocketComplexBatch) -> dict[str, Any]:
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


def _official_complex_dict(batch: PocketComplexBatch, state: str) -> dict[str, Any]:
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


def _tensorize_interactions(system: PocketComplex) -> PocketComplex:
    """Normalize LMDB NumPy interactions for the official batch adapter."""
    if system.interactions is not None and not torch.is_tensor(system.interactions):
        system.interactions = torch.as_tensor(system.interactions)
    return system


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


def _tensor_output_hash(output: dict[str, Any]) -> str:
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


def _transform_coordinate_system(material: dict[str, Any], vocabs: dict[str, Any]):
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


def _run_coordinate_one(
    *,
    model: Any,
    vocabs: dict[str, Any],
    material: dict[str, Any],
    seed: int,
    steps: int,
    device: str,
    stage: str,
    output_root: Path,
) -> dict[str, Any]:
    _reset_seed(seed)
    transformed = _transform_coordinate_system(material, vocabs)
    combined = _to_device(
        _coordinate_complex_dict(PocketComplexBatch([transformed])), device
    )
    ligand_bad = model.builder.extract_ligand_from_complex(combined)
    pocket_data = model.builder.extract_pocket_from_complex(combined)
    fixed_single = torch.as_tensor(material["fixed_mask"], dtype=torch.bool)
    fixed_mask = fixed_single.unsqueeze(0).to(device)
    generator = torch.Generator(device="cpu").manual_seed(seed)
    prior_coords = build_local_coordinate_prior(
        ligand_bad["coords"][0].detach().cpu(), ~fixed_single, generator
    ).unsqueeze(0).to(device)

    if device.startswith("cuda"):
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
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
    runtime = time.perf_counter() - started
    peak_memory = (
        torch.cuda.max_memory_allocated(device) / (1024**2)
        if device.startswith("cuda")
        else 0.0
    )

    center = transformed.com.squeeze(0).detach().cpu().numpy()
    coords_raw = output["coords"][0].detach().cpu().numpy().astype(np.float64) + center
    fixed_indices = np.flatnonzero(material["fixed_mask"])
    raw_drift = (
        float(
            np.linalg.norm(
                coords_raw[fixed_indices] - material["coords_bad"][fixed_indices],
                axis=-1,
            ).max()
        )
        if fixed_indices.size
        else 0.0
    )
    coords_final = coords_raw.copy()
    coords_final[material["fixed_mask"]] = material["coords_bad"][material["fixed_mask"]]
    mol_pred = copy_mol_with_coords(material["bad_mol"], coords_final)
    sdf_path = output_root / material["system_id"] / f"seed_{seed}.sdf"
    _write_mol(
        sdf_path,
        mol_pred,
        {
            "experiment_id": EXPERIMENT_ID,
            "method": "coordinate_only",
            "stage": stage,
            "system_id": material["system_id"],
            "seed": seed,
            "output_to_reference": json.dumps(list(range(mol_pred.GetNumAtoms()))),
        },
    )
    rollout_index = (
        seed - FORMAL_SEED_BASE - material["test_index"] * ROLLOUTS_PER_CASE
        if stage == "formal"
        else -1
    )
    return {
        "experiment_id": EXPERIMENT_ID,
        "method": "coordinate_only",
        "stage": stage,
        "test_index": material["test_index"],
        "system_id": material["system_id"],
        "rollout_index": rollout_index,
        "seed": seed,
        "steps": steps,
        "batch_size": 1,
        "run_status": "completed",
        "error": "",
        "runtime_seconds": runtime,
        "peak_gpu_memory_mb": peak_memory,
        "raw_max_fixed_drift": raw_drift,
        "output_to_reference": json.dumps(
            list(range(mol_pred.GetNumAtoms())), separators=(",", ":")
        ),
        "output_sha256": _tensor_output_hash(output),
        "sdf_path": _relative(sdf_path),
    }


def _run_official_one(
    *,
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
) -> dict[str, Any]:
    _reset_seed(seed)
    if device.startswith("cuda"):
        torch.cuda.reset_peak_memory_stats(device)
    started = time.perf_counter()
    transformed = _tensorize_interactions(transform(material["bad_system"]))
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
        _official_complex_dict(PocketComplexBatch([prior_system]), state="apo"),
        device,
    )
    posterior_data = _to_device(
        _official_complex_dict(PocketComplexBatch([transformed]), state="holo"),
        device,
    )
    ligand_prior = model.builder.extract_ligand_from_complex(prior_data)
    ligand_prior["interactions"] = prior_data["interactions"]
    ligand_prior["fragment_mask"] = prior_data["fragment_mask"]
    ligand_prior["fragment_mode"] = prior_data["fragment_mode"]
    pocket_data = model.builder.extract_pocket_from_complex(posterior_data)
    pocket_data["interactions"] = posterior_data["interactions"]
    pocket_data["complex"] = posterior_data["complex"]
    times = [torch.zeros(1, device=device) for _ in range(3)]
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
    runtime = time.perf_counter() - started
    peak_memory = (
        torch.cuda.max_memory_allocated(device) / (1024**2)
        if device.startswith("cuda")
        else 0.0
    )

    coords_raw = output["coords"][0].detach().cpu().numpy().astype(np.float64)
    coords_final, raw_drift, _ = restore_quantized_fixed_coordinates(
        coords_raw,
        material["coords_bad"],
        material["fixed_mask"],
        output_to_reference,
        tolerance=FIXED_DRIFT_TOLERANCE,
    )
    raw_mol = model._generate_mols(output, sanitise=False)[0]
    mol_pred = None
    if raw_mol is not None and raw_mol.GetNumAtoms() == coords_final.shape[0]:
        mol_pred = copy_mol_with_coords(raw_mol, coords_final)
    if mol_pred is None:
        raise RuntimeError("Official generation did not produce a writable molecule")
    sdf_path = output_root / material["system_id"] / f"seed_{seed}.sdf"
    _write_mol(
        sdf_path,
        mol_pred,
        {
            "experiment_id": EXPERIMENT_ID,
            "method": "official_inpainting",
            "stage": stage,
            "system_id": material["system_id"],
            "seed": seed,
            "output_to_reference": json.dumps(output_to_reference.tolist()),
        },
    )
    rollout_index = (
        seed - FORMAL_SEED_BASE - material["test_index"] * ROLLOUTS_PER_CASE
        if stage == "formal"
        else -1
    )
    return {
        "experiment_id": EXPERIMENT_ID,
        "method": "official_inpainting",
        "stage": stage,
        "test_index": material["test_index"],
        "system_id": material["system_id"],
        "rollout_index": rollout_index,
        "seed": seed,
        "steps": steps,
        "batch_size": 1,
        "run_status": "completed",
        "error": "",
        "runtime_seconds": runtime,
        "peak_gpu_memory_mb": peak_memory,
        "raw_max_fixed_drift": raw_drift,
        "output_to_reference": json.dumps(output_to_reference.tolist(), separators=(",", ":")),
        "output_sha256": _tensor_output_hash(output),
        "sdf_path": _relative(sdf_path),
    }


def _failure_record(
    *, method: str, material: dict[str, Any], seed: int, steps: int, stage: str, error: Exception
) -> dict[str, Any]:
    rollout_index = (
        seed - FORMAL_SEED_BASE - material["test_index"] * ROLLOUTS_PER_CASE
        if stage == "formal"
        else -1
    )
    return {
        "experiment_id": EXPERIMENT_ID,
        "method": method,
        "stage": stage,
        "test_index": material["test_index"],
        "system_id": material["system_id"],
        "rollout_index": rollout_index,
        "seed": seed,
        "steps": steps,
        "batch_size": 1,
        "run_status": "model_error",
        "error": f"{type(error).__name__}: {error}",
        "runtime_seconds": 0.0,
        "peak_gpu_memory_mb": 0.0,
        "raw_max_fixed_drift": float("nan"),
        "output_to_reference": "",
        "output_sha256": "",
        "sdf_path": "",
    }


def _run_method(
    *,
    method: str,
    cases: list[dict[str, Any]],
    seeds_by_case: dict[str, list[int]],
    steps: int,
    device: str,
    stage: str,
    output_root: Path,
    records_path: Path,
) -> list[dict[str, Any]]:
    systems = _system_map()
    records: list[dict[str, Any]] = _read_csv(records_path)
    completed = {
        (record["method"], record["system_id"], int(record["seed"]))
        for record in records
    }
    if method == "coordinate_only":
        model, vocabs = _load_coordinate_components(steps, device)
        components = (model, vocabs)
    else:
        components = _load_official_components(steps, device)

    total_expected = sum(len(values) for values in seeds_by_case.values())
    for case in cases:
        material = _load_material(case, systems[case["system_id"]])
        for seed in seeds_by_case[case["system_id"]]:
            key = (method, case["system_id"], seed)
            if key in completed:
                continue
            try:
                if method == "coordinate_only":
                    record = _run_coordinate_one(
                        model=components[0],
                        vocabs=components[1],
                        material=material,
                        seed=seed,
                        steps=steps,
                        device=device,
                        stage=stage,
                        output_root=output_root / method,
                    )
                else:
                    record = _run_official_one(
                        args=components[0],
                        model=components[1],
                        transform=components[2],
                        interpolant=components[3],
                        material=material,
                        seed=seed,
                        steps=steps,
                        device=device,
                        stage=stage,
                        output_root=output_root / method,
                    )
            except torch.OutOfMemoryError:
                if device.startswith("cuda"):
                    torch.cuda.empty_cache()
                raise
            except Exception as error:
                record = _failure_record(
                    method=method,
                    material=material,
                    seed=seed,
                    steps=steps,
                    stage=stage,
                    error=error,
                )
            records.append(record)
            records.sort(
                key=lambda item: (
                    item["method"],
                    int(item["test_index"]),
                    int(item["seed"]),
                )
            )
            _write_csv(records_path, records, SAMPLING_FIELDS)
            completed.add(key)
            method_done = sum(record["method"] == method for record in records)
            print(
                f"{method}: completed {method_done}/{total_expected}; "
                f"{case['system_id']} seed {seed}",
                flush=True,
            )

    del components
    if device.startswith("cuda"):
        torch.cuda.empty_cache()
    return records


def _smoke_cases(cases: list[dict[str, Any]], steps: int) -> list[dict[str, Any]]:
    ordered = sorted(cases, key=lambda case: (case["editable_atom_count"], case["test_index"]))
    if steps == 100:
        return [ordered[len(ordered) // 2]]
    positions = np.linspace(0, len(ordered) - 1, 5).round().astype(int)
    return [ordered[index] for index in positions]


def run_smoke(steps: int, device: str) -> dict[str, Any]:
    if steps not in {5, 100}:
        raise ValueError("Frozen smoke supports only 5 or 100 steps")
    run_preflight()
    cases = _load_cases()
    selected = _smoke_cases(cases, steps)
    seeds_by_case = {
        case["system_id"]: [424242 + int(case["test_index"])] for case in selected
    }
    smoke_root = OUTPUT_DIR / "smoke" / f"steps_{steps}"
    records_path = smoke_root / "sampling_runs.csv"
    for method in ("coordinate_only", "official_inpainting"):
        _run_method(
            method=method,
            cases=selected,
            seeds_by_case=seeds_by_case,
            steps=steps,
            device=device,
            stage="smoke",
            output_root=smoke_root / "repairs",
            records_path=records_path,
        )
    records = _read_csv(records_path)
    expected = len(selected) * 2
    technical_ok = bool(
        len(records) == expected
        and all(record["run_status"] == "completed" for record in records)
        and all(record["sdf_path"] for record in records)
    )
    summary = {
        "experiment_id": EXPERIMENT_ID,
        "status": "passed" if technical_ok else "failed",
        "steps": steps,
        "technical_ok": technical_ok,
        "case_ids": [case["system_id"] for case in selected],
        "records": len(records),
    }
    _write_json(smoke_root / "summary.json", summary)
    if not technical_ok:
        raise RuntimeError(f"{steps}-step paired smoke failed")
    print(json.dumps(summary, indent=2))
    return summary


def _require_smokes() -> None:
    for steps in (5, 100):
        path = OUTPUT_DIR / "smoke" / f"steps_{steps}" / "summary.json"
        if not path.is_file():
            raise FileNotFoundError(f"Required {steps}-step smoke has not run")
        summary = json.loads(path.read_text(encoding="utf-8"))
        if not summary.get("technical_ok"):
            raise RuntimeError(f"Required {steps}-step smoke did not pass")


def run_formal(method: str, device: str) -> dict[str, Any]:
    run_preflight()
    _require_smokes()
    cases = _load_cases()
    seeds_by_case = {
        case["system_id"]: [
            FORMAL_SEED_BASE + int(case["test_index"]) * ROLLOUTS_PER_CASE + rollout
            for rollout in range(ROLLOUTS_PER_CASE)
        ]
        for case in cases
    }
    methods = (
        ("coordinate_only", "official_inpainting")
        if method == "all"
        else (method,)
    )
    records_path = OUTPUT_DIR / "sampling_runs.csv"
    for selected_method in methods:
        _run_method(
            method=selected_method,
            cases=cases,
            seeds_by_case=seeds_by_case,
            steps=FORMAL_STEPS,
            device=device,
            stage="formal",
            output_root=OUTPUT_DIR / "repairs",
            records_path=records_path,
        )
    records = _read_csv(records_path)
    counts = Counter(record["method"] for record in records)
    expected_per_method = len(cases) * ROLLOUTS_PER_CASE
    inventory_complete = all(
        counts[selected] == expected_per_method
        for selected in ("coordinate_only", "official_inpainting")
    )
    result = {
        "experiment_id": EXPERIMENT_ID,
        "status": "sampling_complete" if inventory_complete else "method_complete",
        "eligible_cases": len(cases),
        "expected_per_method": expected_per_method,
        "method_counts": dict(counts),
    }
    _write_json(OUTPUT_DIR / "sampling_summary.json", result)
    print(json.dumps(result, indent=2))
    return result


def _as_bool(value: Any) -> bool:
    return value is True or str(value).lower() == "true"


def _as_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _numeric_summary(values: list[float]) -> dict[str, Any]:
    finite = np.asarray([value for value in values if math.isfinite(value)], dtype=float)
    if not finite.size:
        return {
            "evaluable": 0,
            "mean": None,
            "std": None,
            "median": None,
            "q1": None,
            "q3": None,
        }
    return {
        "evaluable": int(finite.size),
        "mean": float(np.mean(finite)),
        "std": float(np.std(finite)),
        "median": float(np.median(finite)),
        "q1": float(np.quantile(finite, 0.25)),
        "q3": float(np.quantile(finite, 0.75)),
    }


def _endpoint_funnel(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    stages = [
        ("attempted", lambda record: True),
        (
            "model_completed",
            lambda record: record.get("run_status", "completed") == "completed",
        ),
        ("valid", lambda record: _as_bool(record["valid"])),
        (
            "fully_connected_valid",
            lambda record: _as_bool(record["fully_connected_valid"]),
        ),
        ("condition_match", lambda record: _as_bool(record["condition_match"])),
        ("pb_valid", lambda record: _as_bool(record["pb_valid"])),
        ("same_molecule", lambda record: _as_bool(record["same_molecule"])),
        ("fixed_atoms_ok", lambda record: _as_bool(record["fixed_atoms_ok"])),
        ("fixed_bonds_ok", lambda record: _as_bool(record["fixed_bonds_ok"])),
        ("fixed_coords_ok", lambda record: _as_bool(record["fixed_coords_ok"])),
    ]
    active = list(records)
    previous = len(active)
    funnel = []
    for name, predicate in stages:
        if name != "attempted":
            active = [record for record in active if predicate(record)]
        remaining = len(active)
        funnel.append(
            {
                "stage": name,
                "remaining": remaining,
                "eliminated_at_stage": 0 if name == "attempted" else previous - remaining,
                "rate_of_attempts": remaining / len(records) if records else None,
            }
        )
        previous = remaining
    return funnel


def _summarize_method(records: list[dict[str, Any]]) -> dict[str, Any]:
    attempts = len(records)
    fields = [
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
        "strict_success",
        "strain_evaluable",
    ]
    result = {
        "attempts": attempts,
        "run_status": dict(
            Counter(record.get("run_status", "completed") for record in records)
        ),
        "funnel": _endpoint_funnel(records),
    }
    result.update(
        {
            field: {
                "count": sum(_as_bool(record[field]) for record in records),
                "rate": (
                    sum(_as_bool(record[field]) for record in records) / attempts
                    if attempts
                    else None
                ),
            }
            for field in fields
        }
    )
    result["strain_energy_kcal_mol"] = _numeric_summary(
        [_as_float(record["strain_energy_kcal_mol"]) for record in records]
    )
    result["editable_rmsd_to_good"] = _numeric_summary(
        [_as_float(record["editable_rmsd_to_good"]) for record in records]
    )
    result["all_atom_rmsd_to_good"] = _numeric_summary(
        [_as_float(record["all_atom_rmsd_to_good"]) for record in records]
    )
    result["pb_subtests"] = {
        field: {
            "count": sum(_as_bool(record[field]) for record in records),
            "rate": (
                sum(_as_bool(record[field]) for record in records) / attempts
                if attempts
                else None
            ),
        }
        for field in PB_FIELDS
    }
    return result


def _paired_bootstrap(
    records: list[dict[str, Any]], case_ids: list[str], metric: str
) -> dict[str, Any]:
    differences = []
    per_case = {}
    for system_id in case_ids:
        selected = [record for record in records if record["system_id"] == system_id]
        coordinate = np.mean(
            [
                _as_bool(record[metric])
                for record in selected
                if record["method"] == "coordinate_only"
            ]
        )
        official = np.mean(
            [
                _as_bool(record[metric])
                for record in selected
                if record["method"] == "official_inpainting"
            ]
        )
        difference = float(coordinate - official)
        differences.append(difference)
        per_case[system_id] = {
            "coordinate_only": float(coordinate),
            "official_inpainting": float(official),
            "difference": difference,
        }
    values = np.asarray(differences, dtype=float)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    indices = rng.integers(0, values.size, size=(BOOTSTRAP_REPLICATES, values.size))
    replicates = values[indices].mean(axis=1)
    lower, upper = np.quantile(replicates, [0.025, 0.975])
    if lower > 0:
        conclusion = "coordinate_only_better"
    elif upper < 0:
        conclusion = "official_inpainting_better"
    else:
        conclusion = "inconclusive"
    return {
        "metric": metric,
        "unit": "complex",
        "coordinate_minus_official": float(values.mean()),
        "ci95": [float(lower), float(upper)],
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "conclusion": conclusion,
        "per_case": per_case,
    }


def _paired_contingency(records: list[dict[str, Any]], metric: str) -> dict[str, int]:
    indexed = {
        (record["method"], record["system_id"], int(record["seed"])): _as_bool(record[metric])
        for record in records
    }
    pairs = {(record["system_id"], int(record["seed"])) for record in records}
    counts = Counter()
    for system_id, seed in pairs:
        coordinate = indexed[("coordinate_only", system_id, seed)]
        official = indexed[("official_inpainting", system_id, seed)]
        if coordinate and official:
            counts["both_pass"] += 1
        elif coordinate:
            counts["coordinate_only_pass"] += 1
        elif official:
            counts["official_inpainting_pass"] += 1
        else:
            counts["both_fail"] += 1
    return dict(counts)


def _editable_size_strata(
    records: list[dict[str, Any]], cases: list[dict[str, Any]], metric: str
) -> dict[str, Any]:
    case_size = {case["system_id"]: int(case["editable_atom_count"]) for case in cases}
    bins = [("1-3", 1, 3), ("4-6", 4, 6), ("7-12", 7, 12), ("13+", 13, math.inf)]
    output = {}
    for label, lower, upper in bins:
        ids = {sid for sid, size in case_size.items() if lower <= size <= upper}
        selected = [record for record in records if record["system_id"] in ids]
        output[label] = {
            "complexes": len(ids),
            "coordinate_only": {
                "successes": sum(
                    _as_bool(record[metric])
                    for record in selected
                    if record["method"] == "coordinate_only"
                ),
                "attempts": sum(
                    record["method"] == "coordinate_only" for record in selected
                ),
            },
            "official_inpainting": {
                "successes": sum(
                    _as_bool(record[metric])
                    for record in selected
                    if record["method"] == "official_inpainting"
                ),
                "attempts": sum(
                    record["method"] == "official_inpainting"
                    for record in selected
                ),
            },
        }
    return output


def _evaluate_baselines(
    cases: list[dict[str, Any]], systems: dict[str, Any], buster: PoseBusters
) -> list[dict[str, Any]]:
    records = []
    for case in cases:
        material = _load_material(case, systems[case["system_id"]])
        identity = np.arange(material["good_mol"].GetNumAtoms())
        for method, mol in (("g_good", material["good_mol"]), ("g_bad", material["bad_mol"])):
            metrics = evaluate_output(
                mol_pred=mol,
                material=material,
                output_to_reference=identity,
                raw_fixed_drift=0.0,
                run_completed=True,
                buster=buster,
            )
            records.append(
                {
                    "experiment_id": EXPERIMENT_ID,
                    "method": method,
                    "test_index": case["test_index"],
                    "system_id": case["system_id"],
                    **metrics,
                }
            )
    return records


def _write_case_rates(
    records: list[dict[str, Any]], cases: list[dict[str, Any]]
) -> None:
    output = []
    for case in cases:
        selected = [
            record for record in records if record["system_id"] == case["system_id"]
        ]
        coordinate = [
            record for record in selected if record["method"] == "coordinate_only"
        ]
        inpainting = [
            record
            for record in selected
            if record["method"] == "official_inpainting"
        ]
        coordinate_official = sum(
            _as_bool(record["official_quality_success"]) for record in coordinate
        )
        inpainting_official = sum(
            _as_bool(record["official_quality_success"]) for record in inpainting
        )
        coordinate_strict = sum(
            _as_bool(record["strict_success"]) for record in coordinate
        )
        inpainting_strict = sum(
            _as_bool(record["strict_success"]) for record in inpainting
        )
        output.append(
            {
                "experiment_id": EXPERIMENT_ID,
                "test_index": case["test_index"],
                "system_id": case["system_id"],
                "editable_atom_count": case["editable_atom_count"],
                "num_pairwise_clashes": case["num_pairwise_clashes"],
                "min_relative_distance": case["min_relative_distance"],
                "coordinate_official_successes": coordinate_official,
                "coordinate_official_rate": coordinate_official / len(coordinate),
                "inpainting_official_successes": inpainting_official,
                "inpainting_official_rate": inpainting_official / len(inpainting),
                "coordinate_strict_successes": coordinate_strict,
                "coordinate_strict_rate": coordinate_strict / len(coordinate),
                "inpainting_strict_successes": inpainting_strict,
                "inpainting_strict_rate": inpainting_strict / len(inpainting),
            }
        )
    _write_csv(OUTPUT_DIR / "case_rates.csv", output, CASE_RATE_FIELDS)


def _validate_sampling_inventory(
    sampling: list[dict[str, Any]], cases: list[dict[str, Any]]
) -> None:
    expected = {
        (
            method,
            case["system_id"],
            str(
                FORMAL_SEED_BASE
                + int(case["test_index"]) * ROLLOUTS_PER_CASE
                + rollout
            ),
        )
        for method in ("coordinate_only", "official_inpainting")
        for case in cases
        for rollout in range(ROLLOUTS_PER_CASE)
    }
    observed = [
        (record["method"], record["system_id"], record["seed"])
        for record in sampling
    ]
    if len(observed) != len(set(observed)):
        raise RuntimeError("Sampling inventory contains duplicate keys")
    observed_set = set(observed)
    if observed_set != expected:
        missing = sorted(expected - observed_set)[:5]
        unexpected = sorted(observed_set - expected)[:5]
        raise RuntimeError(
            "Sampling inventory does not match the frozen paired key set: "
            f"missing={missing}, unexpected={unexpected}"
        )


def _write_output_manifest(sampling: list[dict[str, Any]]) -> None:
    paths = {
        OUTPUT_DIR / "preflight.json",
        OUTPUT_DIR / "input_manifest.json",
        OUTPUT_DIR / "construction.csv",
        OUTPUT_DIR / "construction_summary.json",
        OUTPUT_DIR / "cases.json",
        OUTPUT_DIR / "sampling_runs.csv",
        OUTPUT_DIR / "sampling_summary.json",
        OUTPUT_DIR / "runs.csv",
        OUTPUT_DIR / "baselines.csv",
        OUTPUT_DIR / "case_rates.csv",
        OUTPUT_DIR / "summary.json",
    }
    paths.update(
        REPO_ROOT / record["sdf_path"]
        for record in sampling
        if record["run_status"] == "completed" and record["sdf_path"]
    )
    missing = sorted(_relative(path) for path in paths if not path.is_file())
    if missing:
        raise FileNotFoundError(f"Formal output manifest is missing files: {missing[:5]}")
    entries = [
        {
            "path": _relative(path),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in sorted(paths)
    ]
    _write_json(
        OUTPUT_DIR / "output_manifest.json",
        {
            "experiment_id": EXPERIMENT_ID,
            "entries": entries,
            "entry_count": len(entries),
        },
    )


def run_evaluate() -> dict[str, Any]:
    cases = _load_cases()
    sampling = _read_csv(OUTPUT_DIR / "sampling_runs.csv")
    expected_per_method = len(cases) * ROLLOUTS_PER_CASE
    method_counts = Counter(record["method"] for record in sampling)
    expected_counts = {
        "coordinate_only": expected_per_method,
        "official_inpainting": expected_per_method,
    }
    if dict(method_counts) != expected_counts:
        raise RuntimeError(
            f"Sampling inventory is incomplete: {dict(method_counts)} != {expected_counts}"
        )
    _validate_sampling_inventory(sampling, cases)

    systems = _system_map()
    materials = {
        case["system_id"]: _load_material(case, systems[case["system_id"]])
        for case in cases
    }
    buster = PoseBusters(config="dock")
    evaluated = []
    total = len(sampling)
    for index, record in enumerate(sampling, start=1):
        material = materials[record["system_id"]]
        mol_pred = None
        if record["run_status"] == "completed" and record["sdf_path"]:
            try:
                mol_pred = load_sdf(REPO_ROOT / record["sdf_path"], sanitize=False)
            except Exception:
                mol_pred = None
        output_to_reference = (
            np.asarray(json.loads(record["output_to_reference"]), dtype=int)
            if record["output_to_reference"]
            else np.asarray([], dtype=int)
        )
        metrics = evaluate_output(
            mol_pred=mol_pred,
            material=material,
            output_to_reference=output_to_reference,
            raw_fixed_drift=_as_float(record["raw_max_fixed_drift"]),
            run_completed=record["run_status"] == "completed",
            buster=buster,
        )
        evaluated.append({**record, **metrics})
        if index % 50 == 0 or index == total:
            _write_csv(OUTPUT_DIR / "runs.csv", evaluated, RUN_FIELDS)
            print(f"Evaluated {index}/{total}", flush=True)

    baselines = _evaluate_baselines(cases, systems, buster)
    _write_csv(OUTPUT_DIR / "baselines.csv", baselines, BASELINE_FIELDS)
    _write_case_rates(evaluated, cases)
    construction = _read_csv(OUTPUT_DIR / "construction.csv")
    methods = {
        method: _summarize_method(
            [record for record in evaluated if record["method"] == method]
        )
        for method in ("coordinate_only", "official_inpainting")
    }
    case_ids = [case["system_id"] for case in cases]
    strict_comparison = _paired_bootstrap(evaluated, case_ids, "strict_success")
    official_comparison = _paired_bootstrap(
        evaluated, case_ids, "official_quality_success"
    )
    summary = {
        "experiment_id": EXPERIMENT_ID,
        "status": "completed",
        "protocol": {
            "official_test_entries": 225,
            "eligible_cases": len(cases),
            "rollouts_per_case": ROLLOUTS_PER_CASE,
            "steps": FORMAL_STEPS,
            "checkpoint_sha256": EXPECTED_CHECKPOINT_SHA256,
            "posebusters_config_sha256": sha256(POSEBUSTERS_CONFIG),
            "run_py_sha256": sha256(SCRIPT_DIR / "run.py"),
            "protocol_py_sha256": sha256(SCRIPT_DIR / "protocol.py"),
            "repository_head_at_execution": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
            ).strip(),
        },
        "construction_funnel": dict(
            Counter(record["status"] for record in construction)
        ),
        "methods": methods,
        "primary_strict_comparison": strict_comparison,
        "secondary_official_quality_comparison": official_comparison,
        "paired_rollout_contingency": {
            "strict_success": _paired_contingency(evaluated, "strict_success"),
            "official_quality_success": _paired_contingency(
                evaluated, "official_quality_success"
            ),
        },
        "strict_by_editable_size": _editable_size_strata(
            evaluated, cases, "strict_success"
        ),
        "official_quality_by_editable_size": _editable_size_strata(
            evaluated, cases, "official_quality_success"
        ),
        "baselines": {
            method: _summarize_method(
                [record for record in baselines if record["method"] == method]
            )
            for method in ("g_good", "g_bad")
        },
    }
    _write_json(OUTPUT_DIR / "summary.json", summary)
    _write_output_manifest(sampling)
    print(json.dumps(_json_safe(summary), indent=2, ensure_ascii=False))
    return summary


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command", choices=("preflight", "construct", "smoke", "formal", "evaluate")
    )
    parser.add_argument("--steps", type=int, default=5)
    parser.add_argument("--device", default="cuda")
    parser.add_argument(
        "--method",
        choices=("coordinate_only", "official_inpainting", "all"),
        default="all",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    if args.command == "preflight":
        run_preflight()
    elif args.command == "construct":
        run_construct()
    elif args.command == "smoke":
        run_smoke(args.steps, args.device)
    elif args.command == "formal":
        run_formal(args.method, args.device)
    else:
        run_evaluate()
    return 0


if __name__ == "__main__":
    sys.exit(main())
