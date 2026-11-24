#!/usr/bin/env python3
"""Replay fixed finite-port extinction on frozen synthetic and chamber masks."""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

HERE = Path(__file__).resolve().parent
BORN = HERE.parent
sys.path.insert(0, str(BORN))
import extinction as fixed

from flint import arb, ctx
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

STUDY = fixed.STUDY
SYNTHETIC = STUDY / "synthetic/public"
FRESNEL = STUDY / "fresnel"
MANITOBA = STUDY / "validation/manitoba"


def source_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def clock():
    return perf_counter(), process_time()


def elapsed(start):
    return dict(wall_seconds=perf_counter()-start[0], cpu_seconds=process_time()-start[1])


def norm_squared(x):
    return sum((abs(x[i,j])**2 for i in range(x.nrows()) for j in range(x.ncols())), arb(0))


def pair(zw, y):
    yb = fixed.mat(y)
    value = zw.conjugate().transpose() * yb
    return sum((value[i,i].real for i in range(value.nrows())), arb(0)), norm_squared(yb)


def readout(zoom, norm, y, groups):
    value, yy = arb(0), arb(0)
    for z, (receivers, sources) in zip(zoom, groups, strict=True):
        contribution, squares = pair(z, y[np.ix_(receivers, sources)])
        value += contribution
        yy += squares
    quotient = value / norm
    signed_interval = [
        float(np.nextafter(float(quotient.lower()), -np.inf)),
        float(np.nextafter(float(quotient.upper()), np.inf)),
    ]
    lower = max(0., signed_interval[0])
    field_norm = yy.sqrt()
    relative = max(0., float(np.nextafter(float((quotient/field_norm).lower()), -np.inf)))
    return dict(signed_over_norm_interval=signed_interval, distance_lower=lower,
                selected_field_norm=float(field_norm), lower_relative_to_selected_field=relative)


def extinction_setup(s, b, groups, cut=.01):
    """Build exactly the truncated-SVD inverse action on each real receive mask."""
    timings = {}
    t = clock()
    qs = np.zeros(len(b), complex)
    proposals, defects, ranks = [], [], []
    eigs_by_mask = {}
    for rx, tx in groups:
        a, fields = s[rx], b[:,tx]
        key = tuple(rx)
        if key not in eigs_by_mask:
            # SS* is small and shares the nonzero left singular vectors of S.
            vals, u = eigh(a @ a.conj().T, check_finite=False)
            keep = vals >= cut**2 * vals[-1]
            eigs_by_mask[key] = vals[keep], u[:,keep]
        lam, u = eigs_by_mask[key]
        rank = len(lam)
        c = u @ ((u.conj().T @ (a @ fields))/lam[:,None])
        fitted = a.conj().T @ c
        qs += np.sum(fitted * fields.conj(), axis=1)
        defects.append(float(np.linalg.norm(fitted - fields)**2))
        ranks.append(rank)
        proposals.append(c)
    count = sum(len(r) * len(t) for r,t in groups)
    lower_re = float(qs.real.min())
    row = dict(channels=count, group_count=len(groups), ranks=ranks, min_re=lower_re,
               max_re=float(qs.real.max()), sum_relative_field_defect=
                   float(np.sqrt(sum(defects))/np.linalg.norm(b)), eta=None)
    timings["geometry_weight"] = elapsed(t)
    if lower_re <= 0:
        row.update(status="no_phase_under_fixed_prescription", valid=False, timing=timings)
        return row, None, None, None
    eta = max(0., float(np.max(-qs.imag/qs.real))) + 2e-8
    weights = [(1j-eta) * c for c in proposals]
    t = clock()
    exact = [fixed.mat(z) for z in weights]
    norm = sum((norm_squared(z) for z in exact), arb(0)).sqrt()
    # The nonnegative image is shared by all source blocks. Sum each pixel's
    # projected weights across blocks before testing the common cone inequality.
    peak = None
    for lo in range(0, len(b), 384):
        sl = slice(lo,lo+384)
        cells = [arb(0) for _ in b[sl]]
        for (rx,tx), z in zip(groups, exact, strict=True):
            projected = fixed.mat(s[rx,sl].conj().T) * z
            incident = fixed.mat(b[sl,tx].conj())
            for i in range(projected.nrows()):
                cells[i] += sum(((projected[i,j]*incident[i,j]).real
                                 for j in range(len(tx))), arb(0))
        for value in cells:
            upper = value.upper()
            peak = upper if peak is None else peak.max(upper)
    timings["interval_geometry"] = elapsed(t)
    valid = bool(peak < 0)
    row.update(status="valid" if valid else "pixel_sign_not_proved", valid=valid, eta=eta,
               max_cell_upper=float(np.nextafter(float(peak),np.inf)),
               weight_norm_upper=float(np.nextafter(float(norm.upper()),np.inf)), timing=timings)
    return row, exact, norm, weights


def synthetic():
    t = clock()
    with np.load(SYNTHETIC/"fine_geometry.npz", allow_pickle=False) as arr:
        sensor, incident = arr["PW360_10x20__S"], arr["PW360_10x20__B"]
    with np.load(BORN.parent/"array/observations/arrays.npz", allow_pickle=False) as other:
        if not np.array_equal(sensor, other["S20"]) or not np.array_equal(incident, other["B10"]):
            raise ArithmeticError("The original and array extension disagree on s10p20")
    with np.load(BORN/"fixed_port_weights.npz", allow_pickle=False) as arr:
        weight = arr["w_10_20"]
    z, norm, peak = fixed.geometry_check(sensor, incident, weight)
    if not peak < 0:
        raise ArithmeticError("The saved original-array weighting fails its geometry test")
    groups = [(np.arange(len(sensor)), np.arange(incident.shape[1]))]
    saved_controls = json.loads((BORN/"controls/results.json").read_text())["observations"]
    weak = []
    with np.load(BORN/"controls/arrays.npz", allow_pickle=False) as arr:
        for row in saved_controls:
            ident, radius = row["observation"], row["radius"]
            now = clock()
            values = readout([z], norm, arr[ident+"__Y"], groups)
            weak.append(dict(id=ident, radius=radius, noise_pct=100*row["nominal_noise_fraction"],
                             normalized_by_radius=values["distance_lower"]/radius,
                             excluded=values["distance_lower"] > radius, **values, **elapsed(now)))
    with (BORN/"extinction.csv").open() as f:
        rows = list(csv.DictReader(f))
    strong = []
    for n in (10,16,24,32):
        for p in (20,32):
            selected = [r for r in rows if int(r["sources"])==n and int(r["receivers"])==p]
            ratios = [float(r["d_relative"]) for r in selected]
            strong.append(dict(sources=n,receivers=p,views=len(selected),all_excluded=
                               all(r["excludes"]=="True" for r in selected),
                               minimum_lower_relative=min(ratios),maximum_lower_relative=max(ratios)))
    return dict(weak=weak,strong_by_geometry=strong,strong_all_geometry_views=len(rows),
                original_geometry=dict(max_cell_upper=float(peak), weight_norm_upper=float(norm.upper())),
                **elapsed(t))


def laboratory_row(observation, freq, y, groups, audit, exact, norm, scales):
    row = dict(id=observation, frequency_GHz=freq, geometry_status=audit["status"],
               geometry_valid=audit["valid"], channels=audit["channels"], empirical_scales=scales)
    if audit["valid"]:
        t = clock()
        row.update(readout(exact,norm,y,groups))
        d = row["lower_relative_to_selected_field"]
        row["margin_vs_each_scale"] = {key: None if v is None else d-v for key,v in scales.items()}
        row["exceeds_both_empirical_scales"] = all(v is not None and d>v for v in scales.values())
        row["readout_timing"] = elapsed(t)
    return row


def manitoba():
    start = clock()
    physics = source_module("extinction_manitoba_physical", MANITOBA/"physical.py")
    rows, geometries, saved = [], [], {}
    data = json.loads((MANITOBA/"physical_audit.json").read_text())
    records = {(r["frequency_GHz"], r["target"]): r for r in data["cases"]}
    with np.load(physics.PREPARED, allow_pickle=False) as old:
        for fi, f in enumerate(physics.FREQUENCIES):
            s,b = old[f"f{fi}_S"],old[f"f{fi}_B"]
            groups = physics.GROUPS
            proof, z, norm, weights = extinction_setup(s,b,groups)
            geoid = f"manitoba_{f:g}GHz"
            geometries.append(dict(id=geoid,frequency_GHz=float(f),**proof))
            if weights:
                saved.update({f"{geoid}_group_{k}": w for k,w in enumerate(weights)})
            for target in physics.TARGETS:
                y = old[f"f{fi}_{target}"]
                rec = records[(f,target)]
                scales = {k: float(rec[v]) for k,v in (
                    ("two_metal_references","two_reference_columns_difference_relative"),
                    ("empty_drift_plus_held_metal","drift_plus_holdout_over_scatter"),
                )}
                rows.append(laboratory_row(rec["id"],float(f),y,groups,proof,z,norm,scales))
            print(geoid, proof["status"], proof["min_re"], proof["eta"], flush=True)
    return dict(geometries=geometries, observations=rows, **elapsed(start)),saved


def fresnel():
    start = clock()
    physics = source_module("extinction_fresnel_physical", FRESNEL/"physical.py")
    rows, geometries, saved = [], [], {}
    with (STUDY/"results/fresnel/ERROR_LEDGER.csv").open() as f:
        records = {(r["case"], int(r["GHz"])): r for r in csv.DictReader(f)}
    for freq in physics.FREQUENCIES:
        # Call the primary geometry, calibration and source-axis implementation.
        # Its volume kernel is unused in this linear forward-map test.
        geo = physics.operators(freq)
        geo.pop("G")
        s = geo["S"]
        for name in physics.NAMES:
            case = physics.load_case(name,freq,geo,model_audit=False)
            b,y = case["B"]*case["norm"],case["Y"]*case["norm"]
            mask = case["mask"]
            groups = [(np.flatnonzero(mask[:,i]),np.array([i])) for i in range(mask.shape[1])]
            proof, z, norm, weights = extinction_setup(s,b,groups)
            geoid = f"fresnel_{name}_{freq:g}GHz"
            geometries.append(dict(id=geoid,frequency_GHz=int(freq),**proof))
            if weights:
                saved.update({f"{geoid}_group_{k}": w for k,w in enumerate(weights)})
            row = records[(name,freq)]
            empty = float(row["calibration_empty_rms_transfer_all_scatter"])
            gain = float(row["diagnostic_gain_deviation_from_unit"])
            scales = dict(held_empty_rms_transfer=empty,
                          rms_and_diagnostic_complex_gain=(empty+gain)/(1-gain) if gain<1 else None)
            # Verify source column and receiver-mask identity with the archived numerical views.
            with np.load(STUDY/f"results/fresnel/{name}_{freq}GHz_ports.npz",allow_pickle=False) as old:
                assert np.array_equal(mask,old["mask"])
                assert np.allclose(y[mask],old["raw_scattered"][mask],rtol=5e-16,atol=0)
            rows.append(laboratory_row(case["record"]["id"],freq,y,groups,proof,z,norm,scales))
            print(geoid, proof["status"], proof["min_re"], proof["eta"], flush=True)
    return dict(geometries=geometries,observations=rows,**elapsed(start)),saved


def main():
    ctx.prec, ctx.threads=80,1
    run=clock()
    ss=synthetic()
    print("weak view distance/r",[(r["id"],r["normalized_by_radius"]) for r in ss["weak"]],flush=True)
    m,manitoba_weights=manitoba()
    f,fresnel_weights=fresnel()
    np.savez_compressed(HERE/"masked_geometry_weights.npz",**manitoba_weights,**fresnel_weights)
    full=dict(synthetic=ss,manitoba=m,fresnel=f,precision_bits=ctx.prec,**elapsed(run))
    (HERE/"results.json").write_text(json.dumps(full,allow_nan=False,indent=2)+"\n")


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
