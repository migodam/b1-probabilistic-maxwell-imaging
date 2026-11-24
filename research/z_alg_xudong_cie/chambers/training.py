"""Training-only ports and mature-engine wiring for the frozen chamber splits."""

import argparse
from dataclasses import dataclass
import importlib.util
import json
from pathlib import Path
import sys
from time import perf_counter, process_time
from zipfile import ZipFile

import numpy as np
from scipy.fft import ifft2
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESEARCH = ROOT / "research"
OPTICAL = RESEARCH / "z_alg_optical_inertia"
LEGACY = RESEARCH / "limited_aperture_routed_dr_20261120/.runs"
MANITOBA = ROOT / ".artifacts/z_alg_validation/manitoba"
MANITOBA_TRAIN = ROOT / ".artifacts/z_alg_xudong_cie/chambers/manitoba_e_3p5GHz_even_tx_train.npz"
R8 = RESEARCH / "delegated/max_struct_r1/r8_worktree"
MANITOBA_FREQUENCY = 3.5


@dataclass(frozen=True)
class Case:
    id: str
    acquisition: str
    target: str
    phase: str
    frequency_ghz: float
    side: int
    width: float
    train_sources: tuple
    held_sources: tuple

    @property
    def receiver_count(self):
        return 360 if self.acquisition == "Fresnel" else 24

    @property
    def receivers_per_source(self):
        return 241 if self.acquisition == "Fresnel" else 23


CASES = {
    case.id: case for case in (
        Case("fresnel_twin_8GHz_half180_tx", "Fresnel", "FoamTwinDielTM", "development",
             8.0, 60, 0.20, tuple(range(1, 10)), tuple(range(10, 19))),
        Case("fresnel_int_6GHz_half180_tx", "Fresnel", "FoamDielIntTM", "confirmation",
             6.0, 60, 0.20, tuple(range(1, 5)), tuple(range(5, 9))),
        Case("fresnel_ext_6GHz_half180_tx", "Fresnel", "FoamDielExtTM", "confirmation",
             6.0, 60, 0.20, tuple(range(1, 5)), tuple(range(5, 9))),
        Case("manitoba_e_3p5GHz_even_tx", "Manitoba", "polyethylene_e", "confirmation",
             MANITOBA_FREQUENCY, 40, 0.14, tuple(range(2, 25, 2)), tuple(range(1, 24, 2))),
    )
}


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import the archived physical adapter at {path}.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def timers():
    return perf_counter(), process_time()


def elapsed(start):
    return dict(zip(("wall_seconds", "cpu_seconds"),
                    (new - old for new, old in zip(timers(), start, strict=True))))


def prepared_path(case):
    if case.acquisition == "Fresnel":
        return LEGACY / "public" / f"{case.id}.npz"
    return MANITOBA_TRAIN


def calibration_radius_manitoba(physical, fi, k, rx_ids, tx_zero, archived_y, target):
    """Use only calibration objects on the same trained channels."""
    with ZipFile(MANITOBA / "Tomographic_Dataset_1.zip") as raw:
        def channels(key):
            full, _ = physical.load_file(raw, physical.FILE_MAP[key])
            return full[fi, rx_ids.T, tx_zero[None, :]]

        empty = channels("empty_first")
        calibration = channels("cal_metal_3p5in") - empty
        if np.any(abs(calibration) == 0):
            raise ArithmeticError("A selected Manitoba channel has zero calibrator signal.")
        cal_model = physical.pec_scattered(k, .04445)[rx_ids.T, tx_zero[None, :]]
        held_model = physical.pec_scattered(k, .0254)[rx_ids.T, tx_zero[None, :]]
        gain = cal_model / calibration
        drift = gain * (channels("empty_last") - empty)
        steel_residual = gain * (channels("held_metal_2in") - empty) - held_model
        replay = gain * (channels(target) - empty)
        discrepancy = float(np.max(abs(replay - archived_y)))
        if discrepancy > 1e-13 * float(np.max(abs(archived_y))):
            raise ArithmeticError("Prepared Manitoba training data do not reproduce the fixed calibration.")
    drift_norm, hold_norm = float(np.linalg.norm(drift)), float(np.linalg.norm(steel_residual))
    return dict(calibrated_empty_drift=drift_norm, held_steel_residual=hold_norm,
                empirical_radius=drift_norm + hold_norm,
                gain_formula="PEC 3.5-inch prediction / (3.5-inch metal measurement - EMPTY1)",
                archived_train_calibration_max_discrepancy=discrepancy,
                scatter_calibration_fit_uses_target=False,
                radius_status="calibration-only sensitivity scale, no bound on target/model error")


def prepare_manitoba_training():
    """Materialize the standard training-only file from the existing local archive."""
    began = timers()
    case = CASES["manitoba_e_3p5GHz_even_tx"]
    physical = load_module(OPTICAL / "validation/manitoba/physical.py", "manitoba_calibration_reference")
    fi = int(np.flatnonzero(physical.FREQUENCIES == case.frequency_ghz)[0])
    tx = np.array(case.train_sources) - 1
    rx_ids = np.asarray([np.flatnonzero(physical.MASK[:, s]) for s in tx])
    # Only the public training file is exposed to the optimization engine.
    with np.load(MANITOBA / "prepared_private.npz", allow_pickle=False) as archive:
        incident = archive[f"f{fi}_B"][:, tx]
        sensor = archive[f"f{fi}_S"]
        symbol = archive[f"f{fi}_spectrum"]
        y = archive[f"f{fi}_{case.target}"][rx_ids.T, tx[None, :]]
    centers = physical.points(case.side)
    k = 2 * np.pi * case.frequency_ghz * 1e9 / 299792458.0
    calibration = calibration_radius_manitoba(physical, fi, k, rx_ids, tx, y, case.target)
    legacy_style_meta = dict(
        case=case.id, object=case.target, side=case.side, width=case.width, mode="even_tx_no_diagonal",
        suite="MEASURED_CHAMBER", frequency=case.frequency_ghz * 1e9, train_sources=case.train_sources,
        held_sources=case.held_sources, train_channels=y.size,
        held_channels=len(case.held_sources) * case.receivers_per_source,
        source_specific=True, ragged_rx=False, phase=case.phase, calibration=calibration,
    )
    MANITOBA_TRAIN.parent.mkdir(parents=True, exist_ok=True)
    pending = MANITOBA_TRAIN.with_suffix(".pending")
    with pending.open("wb") as stream:
        np.savez_compressed(
            stream, E=incident, S=sensor, Y=y, rx_ids=rx_ids, valid_mask=np.ones(y.shape, dtype=bool),
            spectrum=symbol, centers=centers, radius=(case.width / case.side) / np.sqrt(np.pi),
            k=k, delta=calibration["empirical_radius"], case_json=np.array(json.dumps(legacy_style_meta)),
        )
    pending.replace(MANITOBA_TRAIN)
    return dict(path=str(MANITOBA_TRAIN.relative_to(ROOT)), **calibration, **elapsed(began))


def load_training(case_id):
    """Return only online arrays, per-source receive masks, and calibration metadata."""
    case = CASES[case_id]
    path = prepared_path(case)
    names = ("E", "S", "Y", "rx_ids", "centers", "k", "radius", "delta")
    with np.load(path, allow_pickle=False) as source:
        bundle = {key: source[key] for key in names}
        meta = json.loads(source["case_json"].item())
        if "spectrum" in source.files:
            bundle["spectrum"] = source["spectrum"]
    y, e, s, rx = (bundle[k] for k in ("Y", "E", "S", "rx_ids"))
    source_count, state_count = len(case.train_sources), case.side**2
    if (e.shape != (state_count, source_count) or
            s.shape != (case.receiver_count, state_count) or
            y.shape != (case.receivers_per_source, source_count) or
            rx.shape != y.T.shape or bundle["centers"].shape != (state_count, 2) or
            any(not np.isfinite(bundle[k]).all() for k in ("Y", "E", "S", "centers", "k", "radius", "delta")) or
            float(bundle["delta"]) <= 0 or rx.min() < 0 or rx.max() >= s.shape[0] or
            any(len(np.unique(channel)) != len(channel) for channel in rx)):
        raise ValueError("The online chamber data do not match the frozen geometry or channel layout.")
    if (tuple(meta["train_sources"]) != case.train_sources or
            tuple(meta["held_sources"]) != case.held_sources):
        raise ValueError("The archive source ordering differs from the frozen split.")
    full_mask = np.zeros((s.shape[0], source_count), dtype=bool)
    full_mask[rx.T, np.arange(source_count)[None, :]] = True
    if case.acquisition == "Manitoba" and np.any(full_mask[np.array(case.train_sources)-1,
                                                            np.arange(source_count)]):
        raise ValueError("A missing Manitoba diagonal channel was included.")
    return dict(
        **bundle,
        B=e, E0=e, grid_shape=(case.side, case.side), valid_mask=np.ones(y.shape, dtype=bool),
        training_rx_mask=full_mask, train_port_ids=np.array(case.train_sources),
        fft_geometry=dict(coords=bundle["centers"], radius=float(bundle["radius"]), k=float(bundle["k"])),
        metadata=meta,
    )


def get_case(case_id):
    """Use the same (case metadata, public arrays) convention as the B3 runner."""
    public = load_training(case_id)
    metadata = public.pop("metadata")
    return metadata, public


def mature_imports():
    if str(R8) not in sys.path:
        sys.path.insert(0, str(R8))
    from research.research_max_struct_r8.src.fresnel_adapter import (
        BoundedMatrixFreeModel, SourceReceiverMatrixFreeModel,
    )
    return SourceReceiverMatrixFreeModel, BoundedMatrixFreeModel


def make_model(bundle, *, prior=True, packed_ports=True):
    """Wire the historical sourcewise matrix-free model to either finite operator."""
    packed_engine, full_engine = mature_imports()
    if packed_ports:
        source_indices = np.arange(bundle["E0"].shape[1])
        model = packed_engine(bundle, source_indices, bundle["rx_ids"])
        if not prior:
            model.reg = type(model.reg)(bundle["grid_shape"], smooth=0.0, tv=0.0)
    else:
        model = full_engine(bundle, prior=prior)
    if "spectrum" in bundle:
        symbol = bundle["spectrum"]
        model.green.spectrum = symbol
        model.green.pad = symbol.shape
        model.green.diagonal = ifft2(symbol)[0, 0]
    return model


def preparation_receipt(output):
    began = timers()
    calibration = prepare_manitoba_training()
    rows = {}
    for key, case in CASES.items():
        bundle = load_training(key)
        rows[key] = dict(
            phase=case.phase, acquisition=case.acquisition, source_file=str(prepared_path(case).relative_to(ROOT)),
            target=case.target, frequency_ghz=case.frequency_ghz, native_side=case.side,
            train_tx_one_based=case.train_sources, held_tx_one_based=case.held_sources,
            train_channels=bundle["Y"].size,
            held_channels_by_acquisition_contract=len(case.held_sources) * case.receivers_per_source,
            observed_train_norm=float(np.linalg.norm(bundle["Y"])), empirical_radius=float(bundle["delta"]),
            empirical_relative_radius=float(bundle["delta"] / np.linalg.norm(bundle["Y"])),
            separate_held_receivers_on_training_tx=False,
        )
    report = dict(
        scope="training assets and split layout only; held dielectric scores and materials untouched",
        cases=rows, manitoba_train_calibrator=calibration, material_forward_solves=0, **elapsed(began),
    )
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    for key, row in rows.items():
        print(key, row["train_channels"], row["empirical_relative_radius"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, default=HERE / "training_receipt.json")
    args = parser.parse_args()
    with threadpool_limits(1):
        preparation_receipt(args.receipt)
