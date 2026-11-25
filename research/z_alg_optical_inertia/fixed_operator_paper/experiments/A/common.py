"""Immutable contracts, output formats and stored-operator numerical helpers."""
import hashlib
import json
import os
from pathlib import Path
from time import perf_counter

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                 "MKL_NUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
from scipy.linalg import solve

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
CONTRACT = json.loads((HERE / "PREREGISTRATION.json").read_text())
REGION = np.array(CONTRACT["cohort"]["local_region_indices_row_major_zero_based"])
BGS = CONTRACT["background_caps"]
RATIOS = CONTRACT["width_ratios"]
RESULTS = HERE / "results"
PUBLIC = RESULTS / "acquisition.npz"
PRIVATE = RESULTS / "sealed_truth.npz"
SEAL = RESULTS / "ACQUISITION_SEAL.json"


def digest(path):
    return hashlib.blake2b(path.read_bytes(), digest_size=32).hexdigest()


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(".part")
    temp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def save(path, **arrays):
    temp = path.with_suffix(".part")
    with temp.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    temp.replace(path)


def load(path):
    with np.load(path, allow_pickle=False) as stored:
        return dict(stored)


def check_preregistration():
    record = json.loads((HERE / "PREREGISTRATION_RECEIPT.json").read_text())
    assert record["digest"] == digest(HERE / record["file"])
    return record


def require_p0():
    gate = PAPER / "results/p0.json"
    record = json.loads(gate.read_text())
    assert record.get("pass") is True or record.get("passed") is True or record.get("status") in (
        "PASS", "PASSED", "pass", "passed"
    ), "P0 pass is a prerequisite for prospective decisions"
    return digest(gate)


def star(matrix):
    return matrix.conj().T


def herm(matrix):
    return (matrix + star(matrix)) * 0.5


def forward(green, sensor, incident, material):
    n = len(material)
    system = np.eye(n) - material[:, None] * green
    current = solve(system, material[:, None] * incident, check_finite=False)
    return sensor @ current, current
