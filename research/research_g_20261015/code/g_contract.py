"""G-only sidecar validation; never alters native hypothesis packages.

No physics is run here. Asset checks establish presence/key identity only.
"""
import json
import math
from pathlib import Path, PurePosixPath

VERSION = "GSearchEpisode/v1"
TARGET = "LEGACY_FROZEN_A30_24x24_scalar_TM"
PACKING = "A30_SOURCE_MAJOR_RE32_IM32/v1"
GAIN = "FIXED_ONE/v1"


class ContractError(ValueError):
    pass


def _require(condition, message):
    if not condition:
        raise ContractError(message)


def _number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _online_keys(value):
    if isinstance(value, dict):
        for key, item in value.items():
            normal = key.lower().replace("_", "").replace("-", "")
            _require(not any(word in normal for word in ("truth", "oracle", "holecount", "topologylabel")),
                     "FORBIDDEN_ONLINE_INPUT:" + key)
            _online_keys(item)
    elif isinstance(value, list):
        for item in value:
            _online_keys(item)


def _asset(reference, root, check):
    if isinstance(reference, str):
        name, key = reference, None
    else:
        _require(isinstance(reference, dict) and set(reference) == {"path", "array_key"}, "ASSET_REFERENCE_REQUIRED")
        name, key = reference["path"], reference["array_key"]
        _require(isinstance(key, str) and bool(key.strip()), "ARRAY_KEY_REQUIRED")
    _require(isinstance(name, str) and bool(name.strip()) and "\\" not in name,
             "RELATIVE_ASSET_PATH_REQUIRED")
    path = PurePosixPath(name)
    _require(not path.is_absolute() and ".." not in path.parts and ":" not in name,
             "RELATIVE_ASSET_PATH_REQUIRED")
    if check:
        _require(root is not None, "ASSET_ROOT_REQUIRED")
        base = Path(root).resolve()
        resolved = (base / name).resolve()
        _require(resolved.is_relative_to(base) and resolved.is_file(), "ASSET_MISSING_OR_OUTSIDE_ROOT:" + name)
        if key is not None:
            import numpy as np
            try:
                with np.load(resolved, allow_pickle=False) as arrays:
                    _require(key in arrays.files, "ARRAY_KEY_MISSING:" + key)
            except ContractError:
                raise
            except Exception as exc:
                raise ContractError("ARRAY_ASSET_UNREADABLE:" + name) from exc


def _diagnostic(value, name):
    _require(isinstance(value, dict) and set(value) == {"value", "reason"}, "DIAGNOSTIC_REQUIRED:" + name)
    if value["value"] is None:
        _require(isinstance(value["reason"], str) and bool(value["reason"].strip()), "MISSING_DIAGNOSTIC_REASON:" + name)
    else:
        _require(_number(value["value"]) and value["reason"] is None, "INVALID_DIAGNOSTIC:" + name)


def validate_episode(episode, asset_root=None, check_assets=False):
    """Validate v1 online sidecar, returning True or raising ContractError."""
    _require(isinstance(episode, dict), "EPISODE_OBJECT_REQUIRED")
    _online_keys(episode)
    required = {"version", "episode_id", "native_package_refs", "target", "observations", "candidate",
                "validation", "structure", "history", "probability", "cost"}
    _require(set(episode) == required, "EPISODE_FIELDS_MISMATCH")
    _require(episode["version"] == VERSION, "WRONG_EPISODE_VERSION")
    _require(isinstance(episode["episode_id"], str) and bool(episode["episode_id"].strip()), "EPISODE_ID_REQUIRED")
    refs = episode["native_package_refs"]
    _require(isinstance(refs, list) and bool(refs), "NATIVE_PACKAGE_REFS_REQUIRED")
    for ref in refs:
        _asset(ref, asset_root, check_assets)
    _require(episode["target"] == {"id": TARGET, "packing_version": PACKING, "gain_version": GAIN, "gain": 1},
             "TARGET_PACKING_GAIN_IDENTITY_MISMATCH")
    _require(type(episode["target"]["gain"]) in (int, float), "NUMERIC_GAIN_REQUIRED")
    observations = episode["observations"]
    _require(isinstance(observations, list) and bool(observations), "OBSERVATIONS_REQUIRED")
    ids = []
    for observation in observations:
        _require(isinstance(observation, dict) and set(observation) == {"observation_id", "data_path"}, "OBSERVATION_FIELDS_MISMATCH")
        identifier = observation["observation_id"]
        _require(isinstance(identifier, str) and bool(identifier.strip()), "OBSERVATION_ID_REQUIRED")
        ids.append(identifier)
        _asset(observation["data_path"], asset_root, check_assets)
    _require(len(set(ids)) == len(ids), "DUPLICATE_OBSERVATION_ID")
    candidate = episode["candidate"]
    paths = {"raw_field_path", "corrected_field_path", "verified_field_path", "prediction_path"}
    _require(isinstance(candidate, dict) and set(candidate) == paths | {"candidate_id"}, "CANDIDATE_FIELDS_MISMATCH")
    _require(isinstance(candidate["candidate_id"], str) and bool(candidate["candidate_id"].strip()), "CANDIDATE_ID_REQUIRED")
    for key in paths:
        _asset(candidate[key], asset_root, check_assets)
    validation = episode["validation"]
    _require(isinstance(validation, dict) and set(validation) == {"solver_diagnostics", "constitutive_diagnostics", "full_residual", "material_legal", "full_validation_receipt"}, "VALIDATION_FIELDS_MISMATCH")
    for key in ("solver_diagnostics", "constitutive_diagnostics", "full_residual"):
        _diagnostic(validation[key], key)
    _require(type(validation["material_legal"]) is bool, "MATERIAL_LEGAL_BOOLEAN_REQUIRED")
    _asset(validation["full_validation_receipt"], asset_root, check_assets)
    for key in ("structure", "history", "cost"):
        _require(isinstance(episode[key], dict), "OBJECT_REQUIRED:" + key)
    probability = episode["probability"]
    mass = {"posterior_weight", "represented_mass", "unmodeled_mass"}
    _require(isinstance(probability, dict) and mass <= set(probability) <= mass | {"proposal_density"}, "PROBABILITY_FIELDS_MISMATCH")
    _require(all(probability[key] is None for key in mass), "UNCALIBRATED_MASS_MUST_BE_NULL")
    if "proposal_density" in probability:
        density = probability["proposal_density"]
        _require(isinstance(density, dict) and set(density) == {"value", "interpretation"}, "PROPOSAL_DENSITY_FIELDS_MISMATCH")
        _require(density["interpretation"] == "PROPOSAL_ONLY" and _number(density["value"]) and density["value"] >= 0,
                 "PROPOSAL_DENSITY_IS_NOT_POSTERIOR")
    return True
