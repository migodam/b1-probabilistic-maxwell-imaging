"""Stage public vector changes, execute unknown-support search, then score offline."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import vector_change as vector
import vector_change_subspace as subspace
import numpy as np


HERE = Path(__file__).resolve().parents[1]
PUBLIC_KEYS = {"case", "array", "background", "reference_y"} | {
    f"{mode}_{kind}" for mode in vector.DESIGN["conditions"] for kind in ("observation", "radius")
}
DICTIONARY = "dictionary"
SUBSPACE = "two_sided"
ONLINE_CONFIG = HERE / "configs/vector_change_blind_online.json"


def write_json(path, payload):
    temp = path.with_suffix(".tmp")
    temp.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temp.replace(path)


def case_keys():
    return [(f"{array}_{parent:02d}", array) for parent in vector.DESIGN["parents"]
            for array in vector.DESIGN["arrays"]]


def prepare(source, public, truth):
    public.mkdir(parents=True, exist_ok=True)
    truth.mkdir(parents=True, exist_ok=True)
    index = []
    for case, array in case_keys():
        with np.load(source / f"{case}.npz", allow_pickle=False) as original:
            payload = {"case": np.asarray(case), "array": np.asarray(array),
                       "background": original["background"], "reference_y": original["reference_y"]}
            clean = original["clean_differential"]
            radii = {"clean": 0.0, "diff10": 0.1 * np.linalg.norm(clean),
                     "full10": 0.1 * np.linalg.norm(original["target_y"])}
            for mode in vector.DESIGN["conditions"]:
                payload[f"{mode}_observation"] = original[f"{mode}_observation"]
                payload[f"{mode}_radius"] = np.asarray(radii[mode])
            if set(payload) != PUBLIC_KEYS:
                raise ValueError(f"Unexpected public fields: {case}")
            np.savez_compressed(public / f"{case}.npz", **payload)
            np.savez_compressed(truth / f"{case}.npz", true_sites=original["true_sites"],
                                true_change=original["true_change"], clean_differential=clean,
                                target_material=original["target"])
            index.append(dict(case=case, array=array, file=f"public/{case}.npz"))
    write_json(public.parent / "public_index.json", dict(files=index, fields=sorted(PUBLIC_KEYS)))


def _direct_model(operator):
    os.environ.setdefault("SRAP_DEVICE", "cpu")
    physics = vector._load_module("hc_cie_vector_blind_physics", operator)
    physics.CD = physics.torch.complex128
    return physics, {array: physics.Model(12, array) for array in vector.DESIGN["arrays"]}


def online(public, predictions, operator):
    if os.environ.get(vector.CONFIG_ENV) != str(ONLINE_CONFIG):
        raise RuntimeError("Use guarded-online with the independent generation-free configuration")
    if (predictions / "prediction_seal.json").exists():
        raise FileExistsError(f"Already sealed: {predictions}")
    predictions.mkdir(parents=True, exist_ok=True)
    started, c0 = time.perf_counter(), time.process_time()
    physics, models = _direct_model(operator)
    rows, setups = [], []
    for case, array in case_keys():
        model = models[array]
        path = public / f"{case}.npz"
        wall, cpu = time.perf_counter(), time.process_time()
        off = model.count.copy()
        with np.load(path, allow_pickle=False) as observed:
            if set(observed.files) != PUBLIC_KEYS or str(observed["case"]) != case or str(observed["array"]) != array:
                raise ValueError(f"Invalid public schema: {path}")
            background = observed["background"]
            reference = observed["reference_y"]
            again, current, state_error = vector.forward(physics, model, background)
            reference_error = float(np.linalg.norm(again - reference) / np.linalg.norm(reference))
            if reference_error > 5e-10:
                raise ValueError(f"Reference reforward disagrees: {case}, {reference_error:g}")
            o_map, u, dictionary, adjoint_error = vector.build_dictionary(physics, model, background, current)
            setups.append(dict(case=case, reference_reforward_relative=reference_error,
                               state_error=state_error, adjoint_error=adjoint_error,
                               counters=dict(model.count - off), wall_s=time.perf_counter()-wall,
                               cpu_s=time.process_time()-cpu))
            for mode in vector.DESIGN["conditions"]:
                # Load only the current condition's observed ports and public noise radius.
                measured = observed[f"{mode}_observation"]
                radius = float(observed[f"{mode}_radius"])
                for method in ((DICTIONARY, SUBSPACE) if mode in subspace.CONFIG["conditions"] else (DICTIONARY,)):
                    wall, cpu = time.perf_counter(), time.process_time()
                    off = model.count.copy()
                    if method == DICTIONARY:
                        winner, prediction, shortlist, audited, loss, state_error = vector.search_dictionary(
                            physics, model, background, reference, measured, o_map, u, dictionary)
                        sites, changes = winner["support"], winner["change"]
                        audited_losses = [c["independent_absolute_audit_loss"] for c in audited]
                    else:
                        winner, prediction, shortlist, _left, _right, _singular, _selected = subspace.search(
                            physics, model, background, reference, measured, o_map, u, dictionary)
                        sites, changes = winner["sites"], winner["change"]
                        loss, state_error = winner["independent_loss"], winner["independent_state_error"]
                        audited_losses = [c["independent_loss"] for c in shortlist if "independent_loss" in c]
                    m = background.copy()
                    m[sites] += changes
                    name = f"{case}__{mode}__{method}.npz"
                    np.savez_compressed(predictions / name, sites=np.asarray(sites), changes=np.asarray(changes),
                                        material=m, differential=prediction)
                    rows.append(dict(case=case, array=array, mode=mode, method=method, file=name,
                                     full_candidate_losses=audited_losses, chosen_loss=loss,
                                     chosen_state_error=state_error, noise_radius=radius,
                                     forward_selection_uses_public_ports=True, counters=dict(model.count - off),
                                     local_evaluations=sum(c["local_evaluations"] for c in shortlist),
                                     wall_s=time.perf_counter()-wall, cpu_s=time.process_time()-cpu))
            print("public replay sealed", case, flush=True)
    write_json(predictions / "prediction_seal.json", dict(rows=rows, setups=setups,
               wall_s=time.perf_counter()-started, cpu_s=time.process_time()-c0,
               sealed_epoch_s=time.time(), operator=str(operator), inputs_only="public NPZ and native dyadic operator"))


def guarded_online(public, predictions, operator, denied):
    sandbox = shutil.which("sandbox-exec")
    if sandbox is None:
        raise RuntimeError("This filesystem-enforced replay requires macOS sandbox-exec")
    blocked = [path.resolve() for path in denied]
    profile = "(version 1)\n(allow default)\n" + "\n".join(
        f"(deny file-read* ({'literal' if path.is_file() else 'subpath'} {json.dumps(str(path))}))"
        for path in blocked)
    for root in blocked:
        if public.resolve().is_relative_to(root) or predictions.resolve().is_relative_to(root):
            raise ValueError(f"Public input and predictions must be outside protected root: {root}")
    launcher = [sandbox, "-p", profile, sys.executable, "-B"]
    test = "import errno,os,sys\ntry: fd=os.open(sys.argv[1],os.O_RDONLY)\nexcept OSError as ex: sys.exit(0 if ex.errno in (errno.EACCES, errno.EPERM) else 1)\nelse: os.close(fd);sys.exit(1)"
    probes = []
    for root in blocked:
        file = root if root.is_file() else next((f for f in sorted(root.rglob("*.npz")) if f.is_file()), None)
        if file is None:
            continue
        subprocess.run([*launcher, "-c", test, str(file)], check=True)
        probes.append(str(file))
    subprocess.run([*launcher, str(Path(__file__).resolve()), "online", "--public", str(public),
                    "--predictions", str(predictions), "--operator", str(operator)], check=True,
                   env={**os.environ, vector.CONFIG_ENV: str(ONLINE_CONFIG)})
    write_json(predictions / "access_receipt.json", dict(profile=profile, confirmed_denied_open=probes,
                protected_descriptions=["offline truth", "original mixed/oracle archive", "original stored measurements",
                                        "original stored predictions", "historical sealed vector truth",
                                        "original generation parameters and seeds"],
                generation_free_online_config=str(ONLINE_CONFIG),
                filesystem_profile_applied_to="the full online Python subprocess"))


def score(public, truth, predictions, original):
    access = json.loads((predictions / "access_receipt.json").read_text())
    if not access["confirmed_denied_open"]:
        raise ValueError("This replay has no verified access restrictions")
    seal = json.loads((predictions / "prediction_seal.json").read_text())
    expected = len(case_keys()) * (len(vector.DESIGN["conditions"]) + len(subspace.CONFIG["conditions"]))
    if len(seal["rows"]) != expected or len({j["file"] for j in seal["rows"]}) != expected:
        raise ValueError("Online predictions are not complete")
    if not all((predictions / j["file"]).is_file() for j in seal["rows"]):
        raise FileNotFoundError("At least one sealed prediction file is missing")
    results = []
    with np.load(original / "two_sided.npz", allow_pickle=False) as historical_subspace:
        for row in seal["rows"]:
            case, mode, method = (row[k] for k in ("case", "mode", "method"))
            with (np.load(predictions / row["file"], allow_pickle=False) as new,
                  np.load(public / f"{case}.npz", allow_pickle=False) as measured,
                  np.load(truth / f"{case}.npz", allow_pickle=False) as gt,
                  np.load(original / f"{case}.npz", allow_pickle=False) as baseline):
                if method == DICTIONARY:
                    old = [baseline[f"{mode}_unknown_sites"], baseline[f"{mode}_unknown_change"],
                           baseline[f"{mode}_unknown_prediction"]]
                else:
                    prefix = f"{case}_{mode}"
                    old = [historical_subspace[f"{prefix}_proposed"], historical_subspace[f"{prefix}_proposed_change"],
                           historical_subspace[f"{prefix}_prediction"]]
                site_equal = bool(np.array_equal(new["sites"], old[0]))
                norms = vector.score_change(new["sites"], new["changes"], gt["true_sites"], gt["true_change"],
                                            gt["clean_differential"], new["differential"],
                                            measured[f"{mode}_observation"], row["noise_radius"])
                results.append(dict(case=case, array=row["array"], mode=mode, method=method,
                                    support_equal_to_initial_run=site_equal,
                                    largest_amplitude_difference_from_initial=float(np.max(abs(new["changes"]-old[1]))),
                                    port_difference_from_initial=float(np.linalg.norm(new["differential"]-old[2])), **norms))
    scored_epoch = time.time()
    write_json(predictions / "offline_scores.json", dict(rows=results, scored_epoch_s=scored_epoch,
               online_seal_epoch_s=seal["sealed_epoch_s"], all_predictions_precede_scoring=scored_epoch>=seal["sealed_epoch_s"],
               all_supports_agree=all(r["support_equal_to_initial_run"] for r in results),
               max_amplitude_difference=max(r["largest_amplitude_difference_from_initial"] for r in results),
               max_port_difference=max(r["port_difference_from_initial"] for r in results)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "online", "guarded-online", "score"))
    parser.add_argument("--workspace", type=Path, default=Path.cwd())
    parser.add_argument("--public", type=Path, default=HERE / "results/vector_change_blind/public")
    parser.add_argument("--truth", type=Path, default=HERE / "results/vector_change_blind/offline_truth")
    parser.add_argument("--predictions", type=Path, default=HERE / "results/vector_change_blind/sealed")
    parser.add_argument("--original", type=Path, default=HERE / "results/vector_change")
    parser.add_argument("--operator", type=Path)
    args = parser.parse_args()
    srap = args.workspace.resolve() / "research/srap_cie_array_study_20261121/vector"
    operator = (args.operator or srap / "physics.py").resolve()
    public, truth, predictions, original = (getattr(args, key).resolve() for key in
                                          ("public", "truth", "predictions", "original"))
    if args.stage == "prepare":
        prepare(original, public, truth)
    elif args.stage == "online":
        online(public, predictions, operator)
    elif args.stage == "guarded-online":
        guarded_online(public, predictions, operator,
                       [truth, original, srap / "data", srap / "predictions", srap / "sealed_truth",
                        HERE / "configs/vector_change.json"])
    else:
        score(public, truth, predictions, original)


if __name__ == "__main__":
    main()
