import importlib.util
from pathlib import Path

from rdkit import Chem


RUNNER_PATH = (
    Path(__file__).resolve().parents[2]
    / "experiments"
    / "20260723-03-official-metrics-reevaluation"
    / "run.py"
)
SPEC = importlib.util.spec_from_file_location(
    "official_metrics_reevaluation_runner", RUNNER_PATH
)
runner = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(runner)


def _record(**updates):
    record = {
        "valid": True,
        "fully_connected_valid": True,
        "condition_evaluable": True,
        "condition_match": True,
        "pb_valid": True,
        "conditioned_pb_valid": True,
        "strain_energy_kcal_mol": 5.0,
        **{field: True for field in runner.PB_FIELDS},
    }
    record.update(updates)
    return record


def test_official_condition_match_allows_editable_region_to_change():
    reference = Chem.MolFromSmiles("CCO")
    changed_editable = Chem.MolFromSmiles("CCN")
    missing_fixed_substructure = Chem.MolFromSmiles("CO")

    assert runner.official_condition_match(changed_editable, reference, [2])
    assert not runner.official_condition_match(
        missing_fixed_substructure, reference, [2]
    )


def test_official_pb_valid_requires_every_configured_field():
    values = {field: True for field in runner.PB_FIELDS}
    assert runner.official_pb_valid(values)

    values["bond_angles"] = False
    assert not runner.official_pb_valid(values)


def test_summary_preserves_attempt_denominator_and_missing_strain():
    records = [
        _record(),
        _record(
            valid=False,
            fully_connected_valid=False,
            condition_evaluable=False,
            condition_match=False,
            pb_valid=False,
            conditioned_pb_valid=False,
            strain_energy_kcal_mol=float("nan"),
            **{field: False for field in runner.PB_FIELDS},
        ),
    ]

    summary = runner.summarize_records(records)
    assert summary["attempts"] == 2
    assert summary["valid"] == {"count": 1, "rate": 0.5}
    assert summary["condition_match"]["evaluable"] == 1
    assert summary["condition_match"]["rate_of_attempts"] == 0.5
    assert summary["strain_energy_kcal_mol"]["evaluable"] == 1


def test_paired_counts_uses_matching_system_and_seed():
    records = []
    values = [
        ("a", 1, True, True),
        ("a", 2, True, False),
        ("b", 3, False, True),
        ("b", 4, False, False),
    ]
    for system_id, seed, coordinate, official in values:
        records.extend(
            [
                {
                    "method": "coordinate_only",
                    "system_id": system_id,
                    "seed": seed,
                    "pb_valid": coordinate,
                },
                {
                    "method": "official_inpainting",
                    "system_id": system_id,
                    "seed": seed,
                    "pb_valid": official,
                },
            ]
        )

    assert runner.paired_counts(records, "pb_valid") == {
        "both_pass": 1,
        "coordinate_only_pass": 1,
        "official_inpainting_pass": 1,
        "both_fail": 1,
    }


def test_frozen_inventory_has_two_methods_for_every_pair():
    specs = runner._expected_specs()
    assert len(specs) == 100
    pairs = {}
    for spec in specs:
        key = (spec["system_id"], spec["seed"])
        pairs.setdefault(key, set()).add(spec["method"])
        assert spec["sdf_path"].is_file()
    assert len(pairs) == 50
    assert all(
        methods == {"coordinate_only", "official_inpainting"}
        for methods in pairs.values()
    )
