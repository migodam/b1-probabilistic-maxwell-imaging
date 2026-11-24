"""Independent native material forward and held-transmitter scoring."""

import argparse
import json
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from training import (CASES, HERE, LEGACY, MANITOBA, OPTICAL, elapsed, load_module,
                      load_training, make_model, timers)


def validate_phase(case, freeze_record):
    if case.phase == "confirmation" and (freeze_record is None or not Path(freeze_record).is_file()):
        raise ValueError("Freeze the candidate algorithms and reference the freeze record before confirmation scoring.")


def held_channels(case, training, *, freeze_record):
    validate_phase(case, freeze_record)
    if case.acquisition == "Fresnel":
        path = LEGACY / "heldout" / f"{case.id}.npz"
        with np.load(path, allow_pickle=False) as held:
            if not np.array_equal(held["S"], training["S"]):
                raise ValueError("The Fresnel public and held receiver operators do not match.")
            incident, observed = held["E"], held["Y"]
    else:
        path = MANITOBA / "prepared_private.npz"
        tx = np.asarray(case.held_sources) - 1
        fi = 1  # 3.5 GHz in the original EIGOR frequency grid.
        with np.load(path, allow_pickle=False) as held:
            full_incident = held[f"f{fi}_B"]
            if not np.array_equal(full_incident[:, np.asarray(case.train_sources) - 1], training["E"]):
                raise ValueError("The Manitoba public and offline incident columns do not match.")
            incident = full_incident[:, tx]
            observed = held[f"f{fi}_{case.target}"][:, tx]
        if not np.array_equal(np.isfinite(observed), np.arange(case.receiver_count)[:, None] != tx):
            raise ValueError("The Manitoba held-out acquisition mask differs from the frozen split.")
    s, n = len(case.held_sources), case.side**2
    if (observed.shape != (case.receiver_count, s) or incident.shape != (n, s)
            or not np.isfinite(incident).all()):
        raise ValueError("Offline channel shapes differ from the original acquisition.")
    mask = np.isfinite(observed)
    if np.any(mask.sum(axis=0) != case.receivers_per_source):
        raise ValueError("The independent receiver sets disagree with the frozen acquisition.")
    rx_ids = np.asarray([np.flatnonzero(mask[:, s]) for s in range(s)])
    return incident, observed[rx_ids.T, np.arange(s)], rx_ids


def material_array(value, case):
    m = np.asarray(value)
    if m.shape == (case.side, case.side):
        m = m.ravel()
    if (m.shape != (case.side**2,) or np.iscomplexobj(m) or not np.isfinite(m).all()
            or np.any(m < 0)):
        raise ValueError("The candidate must be a full native-grid real nonnegative contrast array.")
    return m.astype(float, copy=False)


def residual_statistics(predicted, observed, source_ids):
    error = predicted - observed
    norm = np.linalg.norm
    amplitude = float(norm(observed))
    if amplitude <= 0:
        raise ValueError("Cannot normalize a prediction against a zero measured field.")
    err = float(norm(error))
    per_tx = []
    for i, source_id in enumerate(source_ids):
        e, q = float(norm(error[:, i])), float(norm(observed[:, i]))
        per_tx.append(dict(tx_one_based=int(source_id), channels=error.shape[0],
                           absolute_Frobenius=e, relative_to_measured=e / q if q else None))
    return dict(channels=error.size, measured_Frobenius=amplitude, absolute_Frobenius=err,
                relative_to_measured=err / amplitude, complex_RMSE=err / np.sqrt(error.size),
                per_tx=per_tx)


def nominal_metrics(case, centers, material):
    if case.acquisition == "Fresnel":
        physical = load_module(OPTICAL / "fresnel/physical.py", "fresnel_nominal_reference")
        nominal = physical.nominal_contrast(case.target, centers, case.side)
        level, provenance = .225, "historic R8/Geffrin Figure 5 nominal axes and disk cross-sections"
        nominal_label = "unregistered nominal material; nominal contrast 0.45 foam, 2 plastic"
    else:
        physical = load_module(OPTICAL / "validation/manitoba/physical.py", "manitoba_nominal_reference")
        nominal = physical.model_contrast(case.target, case.side)
        level, provenance = .65, "archived Manitoba physical.py, Gilmore et al. Figure 6b nominal E section"
        nominal_label = "unregistered nominal material; contrast 1.3 transferred from 5 GHz literature value"
    nominal = np.asarray(nominal, dtype=float)
    interior = nominal >= level
    air = nominal == 0
    recovered = material >= level
    diff = material - nominal
    support_union = np.count_nonzero(interior | recovered)
    return dict(
        scope="offline comparison to a fixed unregistered nominal reference; not registered truth",
        provenance=provenance, nominal_label=nominal_label, support_threshold=level,
        interior_cells=int(interior.sum()), air_cells=int(air.sum()),
        nominal_interior_relative_l2=float(np.linalg.norm(diff[interior]) / np.linalg.norm(nominal[interior])),
        nominal_interior_RMSE=float(np.sqrt(np.mean(diff[interior]**2))),
        nominal_air_RMSE=float(np.sqrt(np.mean(material[air]**2))),
        nominal_support_IoU=float(np.count_nonzero(interior & recovered) / support_union),
        nominal_full_grid_relative_l2=float(np.linalg.norm(diff) / np.linalg.norm(nominal)),
        material_min=float(np.min(material)), material_max=float(np.max(material)),
    ), nominal


def full_port_prediction(bundle, incident, material):
    model = make_model(dict(bundle, E0=incident, E=incident, B=incident), prior=False, packed_ports=False)
    state = model.state(material)
    residual = float(state["state_relative_residual"])
    if not (np.isfinite(state["Y"]).all() and residual <= 2e-8):
        raise ArithmeticError(f"Fresh native VIE scoring solve failed its state residual check ({residual}).")
    return state["Y"], dict(state_relative_residual=residual, operation_counts=dict(model.cost),
                            original_source_rhs=incident.shape[1])


def score_material(case_id, value, *, freeze_record=None, map_path=None):
    case = CASES[case_id]
    validate_phase(case, freeze_record)
    total = timers()
    m = material_array(value, case)
    bundle = load_training(case_id)
    held_e, held_y, held_rx = held_channels(case, bundle, freeze_record=freeze_record)
    train_tx = np.arange(len(case.train_sources))
    held_tx = np.arange(len(case.held_sources))
    started = timers()
    predicted, native = full_port_prediction(bundle, np.column_stack((bundle["E"], held_e)), m)
    native.update(elapsed(started))
    t_predicted = predicted[bundle["rx_ids"].T, train_tx]
    h_predicted = predicted[held_rx.T, held_tx + len(train_tx)]
    train = residual_statistics(t_predicted, bundle["Y"], case.train_sources)
    held = residual_statistics(h_predicted, held_y, case.held_sources)
    nominal, reference = nominal_metrics(case, bundle["centers"], m)
    train.update(empirical_radius=float(bundle["delta"]),
                 empirical_radius_margin=float(bundle["delta"] - train["absolute_Frobenius"]),
                 radius_status="sensitivity check only; no hard calibration/model discrepancy guarantee")
    if map_path is not None:
        write_sparse_maps(map_path, case, bundle, m, reference, t_predicted, h_predicted, held_y, held_rx)
    return dict(
        scope="fresh physical solve using the stored native geometry and a frozen material",
        case_id=case_id, phase=case.phase, native_side=case.side,
        freeze_record=str(Path(freeze_record).absolute()) if freeze_record is not None else None,
        train=train, held_transmit=held,
        held_receive_only=None,
        held_receive_only_reason="all measured receivers on training transmitters were used in the training split",
        nominal_reference=nominal, native_forward=native, **elapsed(total),
    )


def write_sparse_maps(path, case, bundle, material, reference, train, held, held_y, held_rx):
    """Local sparse absolute residual maps, with the missing channels as NaN."""
    train_map = np.full(bundle["training_rx_mask"].shape, np.nan)
    held_map = np.full((case.receiver_count, len(case.held_sources)), np.nan)
    train_map[bundle["rx_ids"].T, np.arange(train.shape[1])] = abs(train - bundle["Y"])
    held_map[held_rx.T, np.arange(held.shape[1])] = abs(held - held_y)
    with Path(path).open("wb") as stream:
        np.savez_compressed(stream, m=material, nominal=reference, centers=bundle["centers"],
                            train_absolute_port_error=train_map, held_absolute_port_error=held_map,
                            train_tx=np.array(case.train_sources), held_tx=np.array(case.held_sources))


def load_material(path, key, index):
    if path.suffix == ".npy":
        value = np.load(path, allow_pickle=False)
    else:
        with np.load(path, allow_pickle=False) as result:
            value = result[key]
    return value if index is None else value[index]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_id", choices=CASES)
    parser.add_argument("material", type=Path)
    parser.add_argument("--array", default="m")
    parser.add_argument("--index", type=int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--freeze-record", type=Path)
    parser.add_argument("--local-sparse-maps", type=Path)
    args = parser.parse_args()
    with threadpool_limits(1):
        validate_phase(CASES[args.case_id], args.freeze_record)
        value = load_material(args.material, args.array, args.index)
        report = score_material(args.case_id, value, freeze_record=args.freeze_record,
                                map_path=args.local_sparse_maps)
        report.update(material_path=str(args.material.absolute()),
                      material_array=args.array if args.material.suffix != ".npy" else None,
                      material_array_index=args.index)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
