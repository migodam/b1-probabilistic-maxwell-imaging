"""Independent four-material registry for the prospective mechanism ablation."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
if not __package__:
    sys.path.insert(0, str(HERE.parents[2]))
from research.z_alg_xudong_cie.holdout import data as data_contract
from research.z_alg_xudong_cie.holdout import generation as generator

FREEZE = HERE / "FREEZE.json"
DATA = HERE / "acquisition"


def read_freeze():
    return json.loads(FREEZE.read_text(encoding="utf-8"))


def configure_registry():
    """Inject only an independent anonymous cohort into the current process."""
    record = read_freeze()
    data_contract.OBJECT_IDS = tuple(record["objects"])
    data_contract.APERTURES = tuple(record["apertures"])
    data_contract.NOISE_TAGS = {.02: "02"}
    return record


def generate():
    record = configure_registry()
    keys = record["seeds"]
    generator.OBJECT_IDS = tuple(record["objects"])
    generator.NOISE_TAGS = {.02: "02"}
    generator.SRAP_SEED = keys["srap"]
    generator.NEW_SEED = keys["perforated"]
    generator.AMPLITUDE_SEED = keys["amplitude"]
    generator.ORDER_SEED = keys["permutation"]
    generator.NOISE_SEED = keys["noise"]
    recipe = [("srap",0,0), ("srap",4,1), ("new",0,0), ("new",1,0)]
    order = np.random.default_rng(keys["permutation"]).permutation(len(recipe))
    generator.target_description = lambda: {
        name: recipe[index] for name, index in zip(generator.OBJECT_IDS, order, strict=True)
    }
    generator.generation(DATA, FREEZE)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()
    with threadpool_limits(1):
        if args.list:
            configure_registry()
            print(json.dumps(list(data_contract.planned_cases()), indent=2))
        else:
            generate()
