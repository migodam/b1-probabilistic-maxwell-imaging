"""Run the three sealed mechanism arms without editing the original solver."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import sys

from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
if not __package__:
    sys.path.insert(0, str(HERE.parents[2]))
from research.z_alg_xudong_cie.mechanism.design import DATA, configure_registry, read_freeze
from research.z_alg_xudong_cie.holdout import data as data_contract
from research.z_alg_xudong_cie.benchmark import run as benchmark

import cross_engine
from profile_law import ProfileLaw


def configure_worker():
    configure_registry()

    def load(name):
        case, public = data_contract.get_case(name, DATA)
        _, native = data_contract.get_native_case(name, DATA)
        return case, public, native

    benchmark.load = load


def run_observation(task):
    name, index, options, output, methods = task
    configure_worker()
    original_law = cross_engine.PixelLaw
    generator = benchmark.generators
    order = methods[index%len(methods):] + methods[:index%len(methods)]
    results = []

    def extension(model,p,case,method,initial,upper,settings,measured):
        if method == "profile":
            return generator(model,p,case,"pixel",initial,upper,{"pixel":settings["profile"]},measured)
        return generator(model,p,case,method,initial,upper,settings,measured)

    try:
        benchmark.generators = extension
        with threadpool_limits(1):
            for arm in order:
                out = output/f"{name}__{arm}"
                if out.with_suffix(".json").is_file() and out.with_suffix(".npz").is_file():
                    continue
                cross_engine.PixelLaw = ProfileLaw if arm == "profile" else original_law
                receipt = benchmark.Benchmark(name,arm,options).run(output)
                results.append(receipt)
                print(*receipt,flush=True)
    finally:
        cross_engine.PixelLaw, benchmark.generators = original_law, generator
    return results


def main(output):
    freeze = configure_registry()
    timings = freeze["online"]
    methods = tuple(freeze["methods"])
    options = dict(
        synthetic_seconds=timings["wall_seconds"],
        synthetic_audit_seconds=timings["native_audit_reservation_seconds"],
        shortlist=timings["shortlist"], method_settings=freeze["methods"],
    )
    cases = [row["case"] for row in data_contract.planned_cases()]
    tasks = [(name,i,options,output,methods) for i,name in enumerate(cases)]
    with ProcessPoolExecutor(timings["workers"]) as pool:
        for f in as_completed([pool.submit(run_observation,job) for job in tasks]):
            f.result()


if __name__ == "__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output",type=Path,default=HERE/"predictions")
    args=ap.parse_args()
    main(args.output)
