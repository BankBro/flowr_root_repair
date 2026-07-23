#!/usr/bin/env python
"""Reevaluate frozen repair artifacts with official FLOWR.ROOT metrics."""

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import posebusters
from posebusters import PoseBusters
from rdkit import Chem, rdBase

from flowr.gen.utils import check_substructure_match
from flowr.util import rdkit as smolRD
from flowr.util.metrics import evaluate_strain


EXPERIMENT_ID = "20260723-03-official-metrics-reevaluation"
COORDINATE_EXPERIMENT_ID = "20260723-01-oracle-mask-repair-pilot"
OFFICIAL_EXPERIMENT_ID = "20260723-02-official-fragment-inpainting-control"
ROLLOUTS_PER_CASE = 10
SEED_BASE = 2026072300

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
COORDINATE_OUTPUT_DIR = REPO_ROOT / "outputs" / COORDINATE_EXPERIMENT_ID
OFFICIAL_OUTPUT_DIR = REPO_ROOT / "outputs" / OFFICIAL_EXPERIMENT_ID
OUTPUT_DIR = REPO_ROOT / "outputs" / EXPERIMENT_ID
CASES_PATH = (
    REPO_ROOT / "experiments" / COORDINATE_EXPERIMENT_ID / "cases.json"
)
RAW_VAL_DIR = REPO_ROOT / ".data" / "datasets" / "spindr" / "data_prepared" / "val"
POSEBUSTERS_CONFIG = REPO_ROOT / "posebusters" / "config" / "dock.yml"

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

RUN_FIELDS = [
    "experiment_id",
    "method",
    "case_index",
    "system_id",
    "rollout_index",
    "seed",
    "source_sdf",
    "source_sha256",
    "load_ok",
    "valid",
    "fully_connected_valid",
    "condition_evaluable",
    "condition_match",
    "pb_evaluable",
    "pb_valid",
    "conditioned_pb_valid",
    "strain_evaluable",
    "strain_energy_kcal_mol",
    *PB_FIELDS,
    "error",
    "runtime_seconds",
]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def _write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RUN_FIELDS, lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow(
                {field: _csv_value(record.get(field, "")) for field in RUN_FIELDS}
            )
    temporary.replace(path)


def _load_cases() -> list[dict[str, Any]]:
    payload = json.loads(CASES_PATH.read_text(encoding="utf-8"))
    cases = payload.get("cases", [])
    if len(cases) != 5:
        raise RuntimeError(f"Expected 5 frozen cases, found {len(cases)}")
    return cases


def _load_sdf(path: Path, *, sanitize: bool = False) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(
        str(path),
        removeHs=False,
        sanitize=sanitize,
        strictParsing=True,
    )
    mol = next((candidate for candidate in supplier if candidate is not None), None)
    if mol is None:
        raise ValueError(f"Could not load molecule from {path}")
    return mol


def _expected_specs() -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []
    for case in _load_cases():
        case_index = int(case["case_index"])
        system_id = str(case["system_id"])
        for rollout_index in range(ROLLOUTS_PER_CASE):
            seed = SEED_BASE + 10 * case_index + rollout_index
            specs.extend(
                [
                    {
                        "method": "coordinate_only",
                        "case": case,
                        "case_index": case_index,
                        "system_id": system_id,
                        "rollout_index": rollout_index,
                        "seed": seed,
                        "sdf_path": COORDINATE_OUTPUT_DIR
                        / "repairs"
                        / system_id
                        / f"seed_{seed}.sdf",
                    },
                    {
                        "method": "official_inpainting",
                        "case": case,
                        "case_index": case_index,
                        "system_id": system_id,
                        "rollout_index": rollout_index,
                        "seed": seed,
                        "sdf_path": OFFICIAL_OUTPUT_DIR
                        / "official_repairs"
                        / system_id
                        / f"seed_{seed}.sdf",
                    },
                ]
            )
    return specs


def official_pb_valid(values: dict[str, Any]) -> bool:
    """Match FLOWR.ROOT's DataFrame all(axis=1) PB-valid aggregation."""
    return bool(all(bool(values.get(field, False)) for field in PB_FIELDS))


def official_condition_match(
    generated: Chem.Mol,
    reference_bad: Chem.Mol,
    editable_atom_indices: list[int],
) -> bool:
    """Apply the official substructure-inpainting condition matcher."""
    generated_copy = Chem.Mol(generated)
    reference_copy = Chem.Mol(reference_bad)
    Chem.SanitizeMol(generated_copy)
    Chem.SanitizeMol(reference_copy)
    return bool(
        check_substructure_match(
            generated_copy,
            reference_copy,
            inpainting_mode="substructure_inpainting",
            substructure_query=list(editable_atom_indices),
        )
    )


def strain_statistics(values: list[float]) -> dict[str, Any]:
    finite = np.asarray([value for value in values if math.isfinite(value)], dtype=float)
    if finite.size == 0:
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


def summarize_records(records: list[dict[str, Any]]) -> dict[str, Any]:
    attempts = len(records)
    valid = sum(bool(record["valid"]) for record in records)
    fully_connected = sum(
        bool(record["fully_connected_valid"]) for record in records
    )
    condition_evaluable = sum(
        bool(record["condition_evaluable"]) for record in records
    )
    condition_matches = sum(bool(record["condition_match"]) for record in records)
    pb_valid = sum(bool(record["pb_valid"]) for record in records)
    conditioned_pb_valid = sum(
        bool(record["conditioned_pb_valid"]) for record in records
    )

    def _rate(numerator: int, denominator: int) -> float | None:
        return numerator / denominator if denominator else None

    return {
        "attempts": attempts,
        "valid": {"count": valid, "rate": _rate(valid, attempts)},
        "fully_connected_valid": {
            "count": fully_connected,
            "rate": _rate(fully_connected, attempts),
        },
        "condition_match": {
            "count": condition_matches,
            "rate_of_attempts": _rate(condition_matches, attempts),
            "evaluable": condition_evaluable,
            "rate_of_evaluable": _rate(condition_matches, condition_evaluable),
        },
        "pb_valid": {"count": pb_valid, "rate": _rate(pb_valid, attempts)},
        "conditioned_pb_valid": {
            "count": conditioned_pb_valid,
            "rate": _rate(conditioned_pb_valid, attempts),
        },
        "pb_subtests": {
            field: {
                "count": sum(bool(record[field]) for record in records),
                "rate": _rate(
                    sum(bool(record[field]) for record in records), attempts
                ),
            }
            for field in PB_FIELDS
        },
        "strain_energy_kcal_mol": strain_statistics(
            [float(record["strain_energy_kcal_mol"]) for record in records]
        ),
    }


def paired_counts(
    records: list[dict[str, Any]], metric: str
) -> dict[str, int]:
    indexed = {
        (str(record["method"]), str(record["system_id"]), int(record["seed"])): bool(
            record[metric]
        )
        for record in records
    }
    pairs = {
        (str(record["system_id"]), int(record["seed"]))
        for record in records
    }
    counts = Counter()
    for system_id, seed in sorted(pairs):
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
    return {
        key: counts[key]
        for key in (
            "both_pass",
            "coordinate_only_pass",
            "official_inpainting_pass",
            "both_fail",
        )
    }


def _pb_results(buster: PoseBusters, mol: Chem.Mol, pdb_path: Path) -> dict[str, bool]:
    frame = buster.bust([mol], None, str(pdb_path))
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


def _evaluate_one(
    *,
    method: str,
    case: dict[str, Any],
    sdf_path: Path,
    buster: PoseBusters,
    rollout_index: int | None,
    seed: int | None,
) -> dict[str, Any]:
    started = time.perf_counter()
    errors: list[str] = []
    record: dict[str, Any] = {
        "experiment_id": EXPERIMENT_ID,
        "method": method,
        "case_index": int(case["case_index"]),
        "system_id": str(case["system_id"]),
        "rollout_index": "" if rollout_index is None else rollout_index,
        "seed": "" if seed is None else seed,
        "source_sdf": _relative(sdf_path),
        "source_sha256": _sha256(sdf_path),
        "load_ok": False,
        "valid": False,
        "fully_connected_valid": False,
        "condition_evaluable": False,
        "condition_match": False,
        "pb_evaluable": False,
        "pb_valid": False,
        "conditioned_pb_valid": False,
        "strain_evaluable": False,
        "strain_energy_kcal_mol": float("nan"),
        **{field: False for field in PB_FIELDS},
    }
    try:
        mol = _load_sdf(sdf_path, sanitize=False)
        record["load_ok"] = True
    except Exception as exc:
        errors.append(f"load: {type(exc).__name__}: {exc}")
        record["error"] = "; ".join(errors)
        record["runtime_seconds"] = time.perf_counter() - started
        return record

    record["valid"] = bool(smolRD.mol_is_valid(mol, connected=False))
    record["fully_connected_valid"] = bool(
        smolRD.mol_is_valid(mol, connected=True)
    )

    reference_bad_path = (
        COORDINATE_OUTPUT_DIR
        / "inputs"
        / f"{case['system_id']}__bad.sdf"
    )
    pdb_path = RAW_VAL_DIR / f"{case['system_id']}.pdb"

    if record["fully_connected_valid"]:
        record["condition_evaluable"] = True
        try:
            reference_bad = _load_sdf(reference_bad_path, sanitize=False)
            record["condition_match"] = official_condition_match(
                mol,
                reference_bad,
                list(case["editable_atom_indices"]),
            )
        except Exception as exc:
            record["condition_evaluable"] = False
            errors.append(f"condition: {type(exc).__name__}: {exc}")

        try:
            strain_mol = Chem.Mol(mol)
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
            record["strain_energy_kcal_mol"] = strain
            record["strain_evaluable"] = math.isfinite(strain)
        except Exception as exc:
            errors.append(f"strain: {type(exc).__name__}: {exc}")

    try:
        pb_values = _pb_results(buster, mol, pdb_path)
        record.update(pb_values)
        record["pb_evaluable"] = True
        record["pb_valid"] = official_pb_valid(pb_values)
    except Exception as exc:
        errors.append(f"posebusters: {type(exc).__name__}: {exc}")

    record["conditioned_pb_valid"] = bool(
        record["fully_connected_valid"]
        and record["condition_match"]
        and record["pb_valid"]
    )
    record["error"] = "; ".join(errors)
    record["runtime_seconds"] = time.perf_counter() - started
    return record


def _build_manifest() -> dict[str, Any]:
    entries: list[dict[str, Any]] = []
    for spec in _expected_specs():
        path = Path(spec["sdf_path"])
        entries.append(
            {
                "kind": "generated_result",
                "method": spec["method"],
                "case_index": spec["case_index"],
                "system_id": spec["system_id"],
                "rollout_index": spec["rollout_index"],
                "seed": spec["seed"],
                "path": _relative(path),
                "sha256": _sha256(path),
                "bytes": path.stat().st_size,
            }
        )
    for case in _load_cases():
        system_id = str(case["system_id"])
        for pose in ("good", "bad"):
            path = COORDINATE_OUTPUT_DIR / "inputs" / f"{system_id}__{pose}.sdf"
            entries.append(
                {
                    "kind": f"g_{pose}",
                    "case_index": int(case["case_index"]),
                    "system_id": system_id,
                    "path": _relative(path),
                    "sha256": _sha256(path),
                    "bytes": path.stat().st_size,
                }
            )
        pdb_path = RAW_VAL_DIR / f"{system_id}.pdb"
        entries.append(
            {
                "kind": "protein_pdb",
                "case_index": int(case["case_index"]),
                "system_id": system_id,
                "path": _relative(pdb_path),
                "sha256": _sha256(pdb_path),
                "bytes": pdb_path.stat().st_size,
            }
        )
    return {"experiment_id": EXPERIMENT_ID, "entries": entries}


def run_preflight() -> dict[str, Any]:
    specs = _expected_specs()
    if len(specs) != 100:
        raise RuntimeError(f"Expected 100 generated artifacts, found {len(specs)}")
    method_counts = Counter(str(spec["method"]) for spec in specs)
    if method_counts != {"coordinate_only": 50, "official_inpainting": 50}:
        raise RuntimeError(f"Unexpected method counts: {method_counts}")

    keys = [
        (str(spec["method"]), str(spec["system_id"]), int(spec["seed"]))
        for spec in specs
    ]
    if len(keys) != len(set(keys)):
        raise RuntimeError("Duplicate method/system/seed input keys")

    expected_pairs = Counter(
        (str(spec["system_id"]), int(spec["seed"])) for spec in specs
    )
    if set(expected_pairs.values()) != {2}:
        raise RuntimeError("The two methods do not cover identical system/seed pairs")

    required_paths = [Path(spec["sdf_path"]) for spec in specs]
    for case in _load_cases():
        system_id = str(case["system_id"])
        required_paths.extend(
            [
                COORDINATE_OUTPUT_DIR / "inputs" / f"{system_id}__good.sdf",
                COORDINATE_OUTPUT_DIR / "inputs" / f"{system_id}__bad.sdf",
                RAW_VAL_DIR / f"{system_id}.pdb",
            ]
        )
    missing = [_relative(path) for path in required_paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing frozen inputs: {missing}")

    sdf_paths = [path for path in required_paths if path.suffix == ".sdf"]
    for path in sdf_paths:
        _load_sdf(path, sanitize=False)

    manifest = _build_manifest()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    _write_json(OUTPUT_DIR / "input_manifest.json", manifest)
    preflight = {
        "experiment_id": EXPERIMENT_ID,
        "status": "passed",
        "generated_artifacts": len(specs),
        "coordinate_only_artifacts": method_counts["coordinate_only"],
        "official_inpainting_artifacts": method_counts["official_inpainting"],
        "baseline_sdfs": 10,
        "protein_pdbs": 5,
        "cases_sha256": _sha256(CASES_PATH),
        "posebusters_config": _relative(POSEBUSTERS_CONFIG),
        "posebusters_config_sha256": _sha256(POSEBUSTERS_CONFIG),
        "input_manifest_sha256": _sha256(OUTPUT_DIR / "input_manifest.json"),
    }
    _write_json(OUTPUT_DIR / "preflight.json", preflight)
    return preflight


def _evaluate_specs(specs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    buster = PoseBusters(config="dock")
    records: list[dict[str, Any]] = []
    total = len(specs)
    for index, spec in enumerate(specs, start=1):
        records.append(
            _evaluate_one(
                method=str(spec["method"]),
                case=spec["case"],
                sdf_path=Path(spec["sdf_path"]),
                buster=buster,
                rollout_index=int(spec["rollout_index"]),
                seed=int(spec["seed"]),
            )
        )
        if index == total or index % 10 == 0:
            print(f"Evaluated {index}/{total}", flush=True)
    return records


def run_smoke() -> dict[str, Any]:
    run_preflight()
    first_case_specs = [
        spec
        for spec in _expected_specs()
        if int(spec["case_index"]) == 0 and int(spec["rollout_index"]) == 0
    ]
    if len(first_case_specs) != 2:
        raise RuntimeError("Smoke requires one artifact from each method")
    records = _evaluate_specs(first_case_specs)
    smoke_dir = OUTPUT_DIR / "smoke"
    _write_csv(smoke_dir / "runs.csv", records)
    summary = {
        "experiment_id": EXPERIMENT_ID,
        "status": "passed",
        "records": len(records),
        "methods": {
            record["method"]: {
                "valid": record["valid"],
                "condition_match": record["condition_match"],
                "pb_valid": record["pb_valid"],
                "strain_energy_kcal_mol": record["strain_energy_kcal_mol"],
            }
            for record in records
        },
    }
    _write_json(smoke_dir / "summary.json", summary)
    return summary


def _evaluate_baselines() -> list[dict[str, Any]]:
    buster = PoseBusters(config="dock")
    records: list[dict[str, Any]] = []
    for case in _load_cases():
        system_id = str(case["system_id"])
        for pose in ("good", "bad"):
            records.append(
                _evaluate_one(
                    method=f"g_{pose}",
                    case=case,
                    sdf_path=COORDINATE_OUTPUT_DIR
                    / "inputs"
                    / f"{system_id}__{pose}.sdf",
                    buster=buster,
                    rollout_index=None,
                    seed=None,
                )
            )
    return records


def _version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        if name == "posebusters":
            return str(getattr(posebusters, "__version__", "unknown"))
        return "unknown"


def _git_head() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
    ).strip()


def _build_summary(
    records: list[dict[str, Any]], baseline_records: list[dict[str, Any]]
) -> dict[str, Any]:
    methods: dict[str, Any] = {}
    for method in ("coordinate_only", "official_inpainting"):
        selected = [record for record in records if record["method"] == method]
        if len(selected) != 50:
            raise RuntimeError(f"Expected 50 {method} records, found {len(selected)}")
        methods[method] = summarize_records(selected)
        methods[method]["per_case"] = {
            str(case["system_id"]): summarize_records(
                [
                    record
                    for record in selected
                    if record["system_id"] == case["system_id"]
                ]
            )
            for case in _load_cases()
        }

    baselines = {
        method: summarize_records(
            [record for record in baseline_records if record["method"] == method]
        )
        for method in ("g_good", "g_bad")
    }
    return {
        "experiment_id": EXPERIMENT_ID,
        "status": "completed",
        "protocol": {
            "attempts_per_method": 50,
            "posebusters_config": _relative(POSEBUSTERS_CONFIG),
            "posebusters_config_sha256": _sha256(POSEBUSTERS_CONFIG),
            "posebusters_version": _version("posebusters"),
            "rdkit_version": rdBase.rdkitVersion,
            "strain_force_field": "MMFF94s",
            "strain_max_iterations": 500,
            "condition_mode": "substructure_inpainting",
            "condition_reference": "G_bad",
            "repository_head_at_execution": _git_head(),
        },
        "methods": methods,
        "paired": {
            "pb_valid": paired_counts(records, "pb_valid"),
            "conditioned_pb_valid": paired_counts(
                records, "conditioned_pb_valid"
            ),
        },
        "baselines": baselines,
    }


def run_formal() -> dict[str, Any]:
    run_preflight()
    records = _evaluate_specs(_expected_specs())
    baseline_records = _evaluate_baselines()
    if len(records) != 100 or len(baseline_records) != 10:
        raise RuntimeError("Formal evaluation produced an incomplete result set")
    summary = _build_summary(records, baseline_records)
    _write_csv(OUTPUT_DIR / "runs.csv", records)
    _write_csv(OUTPUT_DIR / "baselines.csv", baseline_records)
    _write_json(OUTPUT_DIR / "summary.json", summary)
    return summary


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("preflight", "smoke", "formal"))
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    if args.command == "preflight":
        result = run_preflight()
    elif args.command == "smoke":
        result = run_smoke()
    else:
        result = run_formal()
    print(json.dumps(_json_safe(result), indent=2, ensure_ascii=False, allow_nan=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
