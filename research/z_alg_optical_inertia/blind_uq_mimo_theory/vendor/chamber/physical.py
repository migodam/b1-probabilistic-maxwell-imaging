"""Calibrate the July 2010 University of Manitoba EIGOR tomography session."""
import argparse
import io
import json
from pathlib import Path
from time import perf_counter, process_time
from zipfile import ZipFile

import numpy as np
from scipy.fft import fft2, ifft2
from scipy.linalg import solve
from scipy.special import hankel1, jv
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
WORK = ROOT / ".artifacts/z_alg_validation/manitoba"
RAW = WORK / "Tomographic_Dataset_1.zip"
PREPARED = WORK / "prepared_private.npz"
FREQUENCIES = np.arange(3., 6.01, .5)
SIDE, WIDTH, RADIUS = 40, .14, .144
# Matrices use receiver rows and transmitter columns. This assignment uses
# the paper's port distances and one arbitrary common azimuth origin.
PORT_ANGLES = np.arange(24)*np.pi/12
GROUPS = ((np.arange(1,24,2),np.arange(0,24,2)),
          (np.arange(0,24,2),np.arange(1,24,2)))
MASK = ~np.eye(24, dtype=bool)
USED = np.zeros((24,24),bool)
for _rx,_tx in GROUPS:
    USED[np.ix_(_rx,_tx)] = True
FILE_MAP = {
    "empty_first": "tomographic_incident.txt",
    "empty_last": "tomographic_incident_2.txt",
    "cal_metal_3p5in": "tomographic_metal_3p5_inch_diam.txt",
    "held_metal_2in": "tomographic_metal_2p0_inch_diam.txt",
    "nylon_large": "tomographic_nylon_4p0_inch_diam.txt",
    "nylon_small": "tomographic_nylon_1p5_inch_diam.txt",
    "nylon_pair": "tomographic_nylon_cysl_0p5cm.txt",
    "polyethylene_e": "tomographic_e_phantom.txt",
}
TARGETS = ("nylon_large", "nylon_small", "nylon_pair", "polyethylene_e")


def timer():
    return perf_counter(), process_time()


def elapsed(start):
    return dict(zip(("wall_seconds", "cpu_seconds"),
                    (x-y for x,y in zip(timer(),start,strict=True))))


def load_file(archive, filename):
    lines = archive.read(filename).decode("latin-1").splitlines()
    numeric = "\n".join(line for line in lines if line and not line.startswith("//"))
    values = np.loadtxt(io.StringIO(numeric))
    freqs = np.searchsorted(FREQUENCIES,values[:,0])
    source = values[:,1].astype(int)-1
    receiver = values[:,2].astype(int)-1
    linear = np.ravel_multi_index((freqs, receiver, source), (len(FREQUENCIES),24,24))
    if (values.shape != (len(FREQUENCIES)*24*23,7) or
            not np.array_equal(FREQUENCIES[freqs],values[:,0]) or
            np.any(values[:,3:5]) or len(np.unique(linear)) != len(linear)):
        raise ValueError(f"Malformed port records in {filename}")
    fields = np.full((len(FREQUENCIES),24,24), np.nan+1j*np.nan)
    fields[freqs,receiver,source] = values[:,5]-1j*values[:,6]
    if not np.array_equal(np.isfinite(fields),np.broadcast_to(MASK,fields.shape)):
        raise ValueError(f"Unexpected port mask in {filename}")
    return fields, lines[0][2:].strip()


def points(side):
    axis=(np.arange(side)+.5)*WIDTH/side-WIDTH/2
    yy,xx=np.meshgrid(axis,axis,indexing="ij")
    return np.column_stack((xx.ravel(),yy.ravel()))


def operator(frequency, side):
    k = 2*np.pi*frequency*1.e9 / 299792458.
    centers = points(side)
    locations = RADIUS * np.column_stack((np.cos(PORT_ANGLES),np.sin(PORT_ANGLES)))
    spacing = WIDTH / side
    radius = spacing / np.sqrt(np.pi)
    amp = .5j*np.pi*k*radius*jv(1,k*radius)
    diag = .5j*np.pi*k*radius*hankel1(1,k*radius)-1
    offsets = spacing*np.minimum(np.arange(2*side),2*side-np.arange(2*side))
    d = np.hypot(offsets[:,None],offsets)
    d[0,0]=1.
    kernel = amp*hankel1(0,k*d)
    kernel[0,0]=diag
    rr = np.linalg.norm(locations[:,None] - centers,axis=2)
    sensor = amp*hankel1(0,k*rr)
    incident = .25j*hankel1(0,k*rr.T)
    return k, fft2(kernel), sensor, incident


def pec_scattered(k, radius, order=48):
    n = np.arange(order+1)
    exterior = hankel1(n, k*RADIUS)
    surface = hankel1(n,k*radius)
    mode = -0.25j * exterior * (exterior/surface) * jv(n,k*radius)
    mode[1:] *= 2
    delta=PORT_ANGLES[:,None]-PORT_ANGLES
    return np.cos(delta[:,:,None]*n) @ mode


def model_contrast(target, side):
    centers = points(side)
    step = WIDTH/side
    offsets = (np.arange(8)+.5)/8-.5
    fraction = np.zeros(len(centers))
    radius_small = .01905
    for dx in offsets:
        for dy in offsets:
            x,y = (centers + step*np.array([dx,dy])).T
            if target == "nylon_large":
                support = np.hypot(x,y) < .0508
            elif target == "nylon_small":
                support = np.hypot(x,y) < radius_small
            elif target == "nylon_pair":
                center_distance = 2*radius_small+.005
                support = ((np.hypot(x-center_distance/2,y) < radius_small) |
                           (np.hypot(x+center_distance/2,y) < radius_small))
            else:
                x,y = 100*x+5.4,100*y+3.85
                support = (0<x)&(x<10.8)&(0<y)&(y<7.7)
                for lo,hi,depth in ((2,4,3.2),(5.5,7,3.2),(8.5,9.7,2.4)):
                    support &= ~((lo<x)&(x<hi)&(y>7.7-depth))
                for lo,hi,depth in ((.8,1.6,1.6),(2.4,3.2,1.6),(4,5,2.4),
                                    (6,7.2,3.2),(8.4,9.6,3.2)):
                    support &= ~((lo<x)&(x<hi)&(y<depth))
            fraction += support / 64
    # Reference spectra held fixed across 3–6 GHz; no material setpoint
    # is used to calibrate or to propose certificate multipliers.
    return fraction * (1.3 if target == "polyethylene_e" else (2.03 + .03j))


def nominal_forward(contrast, side, symbol, sensor, incident):
    active = np.flatnonzero(contrast)
    iy,ix=np.divmod(active,side)
    kernel=ifft2(symbol)
    green=kernel[(iy[:,None]-iy)%(2*side),(ix[:,None]-ix)%(2*side)]
    m=contrast[active]
    a=-m[:,None]*green
    a.flat[::len(a)+1]+=1.
    rhs=m[:,None]*incident[active]
    x=solve(a,rhs,assume_a="gen",check_finite=False,overwrite_a=False)
    residual=np.linalg.norm(a@x-rhs)/np.linalg.norm(rhs)
    return sensor[:,active]@x, dict(active_cells=len(active), rhs_columns=x.shape[1],
                                    linear_relative_residual=float(residual))


def norm(x,mask=USED):
    return float(np.linalg.norm(x[mask]))


def estimate_eta(calibrated, model):
    # Held-out aggregate gain diagnostic only: a single gain from neither
    # this target nor this fit is used as part of calibration.
    good=USED & (abs(model) > 0)
    simple=np.vdot(model[good],calibrated[good])/np.vdot(model[good],model[good])
    eta=abs(calibrated[good]/model[good]-1)
    return dict(global_gain_deviation=float(abs(simple-1)),
                per_channel_gain_deviation_median=float(np.median(eta)),
                per_channel_gain_deviation_max=float(eta.max()))


def main(raw, physical_output, prepared):
    start=timer()
    archive_members={}
    measurements={}
    with ZipFile(raw) as archive:
        for key,filename in FILE_MAP.items():
            values, stamp=load_file(archive,filename)
            measurements[key]=values
            archive_members[key]=dict(member=filename,date=stamp,
                                      receiver_transmitter_frequency_rows=int(np.isfinite(values).sum()))
    rows, cases, arrays = [], [], {}
    nominal_materials = {n:{side:model_contrast(n,side) for side in (SIDE,2*SIDE)} for n in TARGETS}
    for fi,frequency in enumerate(FREQUENCIES):
        frequency_start=timer()
        k,symbol,sensor,incident=operator(frequency,SIDE)
        _,fine_symbol,fine_sensor,fine_incident=operator(frequency,2*SIDE)
        raw_empty = measurements["empty_first"][fi]
        raw_late = measurements["empty_last"][fi]
        cal_delta = measurements["cal_metal_3p5in"][fi] - raw_empty
        held_delta = measurements["held_metal_2in"][fi] - raw_empty
        cal_model = pec_scattered(k,.04445)
        held_model = pec_scattered(k,.0254)
        if not np.all(abs(cal_delta[MASK]) > 0):
            raise ArithmeticError("Zero calibration scatter on an acquired channel")
        gain=np.full_like(raw_empty,np.nan+1j*np.nan)
        gain[MASK]=cal_model[MASK]/cal_delta[MASK]
        gain_small=np.full_like(gain,np.nan+1j*np.nan)
        gain_late=np.full_like(gain,np.nan+1j*np.nan)
        small_delta = held_delta[MASK]
        late_delta = (measurements["cal_metal_3p5in"][fi] - raw_late)[MASK]
        if min(abs(small_delta).min(),abs(late_delta).min()) <= 0:
            raise ArithmeticError("A diagnostic calibration denominator is exactly zero.")
        gain_small[MASK]=held_model[MASK]/small_delta
        gain_late[MASK]=cal_model[MASK]/late_delta
        held = gain*held_delta
        drift = gain*(raw_late-raw_empty)
        held_err=norm(held-held_model)
        drift_abs=norm(drift)
        # An opposite convention is logged on the metal holdout only.
        # It never feeds the chosen gain or dielectric certificates.
        sign_control=np.full_like(held,np.nan+1j*np.nan)
        sign_control[MASK]=cal_model[MASK]*(held_delta[MASK]/cal_delta[MASK]).conj()
        cal_series=norm(cal_model-pec_scattered(k,.04445,36))
        hold_series=norm(held_model-pec_scattered(k,.0254,36))
        row=dict(frequency_GHz=float(frequency),mask_channels=int(USED.sum()),
                 all_observed_channels=int(MASK.sum()), raw_empty_drift=norm(raw_late-raw_empty),
                 raw_empty_drift_rel=norm(raw_late-raw_empty)/norm(raw_empty),
                 calibrated_empty_drift=drift_abs,
                 calibrated_empty_drift_all_observed=norm(drift,MASK),
                 calibrated_held_metal_residual=held_err,
                 calibrated_held_metal_residual_all_observed=norm(held-held_model,MASK),
                 held_metal_measured_norm=norm(held),
                 held_metal_reference_norm=norm(held_model),
                 held_metal_relative_to_measured=held_err/norm(held),
                 opposite_sign_held_residual=norm(sign_control-held_model),
                 analytic_mode_tail_delta=cal_series+hold_series,
                 raw_calibrator_minimum_scattered_amplitude=float(abs(cal_delta[USED]).min()),
                 raw_alternate_calibrator_minimum_scattered_amplitude=float(abs(held_delta[USED]).min()),
                 absolute_calibration_gain_quantiles=np.quantile(abs(gain[USED]),[0,.25,.5,.75,1]).tolist(),
                 absolute_alternate_gain_quantiles=np.quantile(abs(gain_small[USED]),[0,.25,.5,.75,1]).tolist(),
                 **estimate_eta(held,held_model))
        arrays[f"f{fi}_spectrum"],arrays[f"f{fi}_S"],arrays[f"f{fi}_B"]=symbol,sensor,incident
        for target in TARGETS:
            calibrated=gain*(measurements[target][fi]-raw_empty)
            alternate=gain_small*(measurements[target][fi]-raw_empty)
            after=gain_late*(measurements[target][fi]-raw_late)
            selected=norm(calibrated)
            if selected<=0:
                raise ArithmeticError("Blank dielectric scattering field")
            coarse,cmeta=nominal_forward(nominal_materials[target][SIDE],SIDE,symbol,sensor,incident)
            fine,fmeta=nominal_forward(nominal_materials[target][2*SIDE],2*SIDE,
                                      fine_symbol,fine_sensor,fine_incident)
            obskey=f"f{fi}_{target}"
            arrays[obskey]=calibrated
            cases.append(dict(id=f"manitoba_eigor2010/{target}/{frequency:g}GHz", target=target,
                              frequency_GHz=float(frequency), field_norm=selected,
                              all_port_field_norm=norm(calibrated,MASK),
                              measurement_energy_kept=(selected/norm(calibrated,MASK))**2,
                              drift_abs=drift_abs,drift_over_scatter=drift_abs/selected,
                              held_metal_residual_abs=held_err,held_metal_residual_over_scatter=held_err/selected,
                              drift_plus_holdout_over_scatter=(drift_abs+held_err)/selected,
                              two_reference_columns_difference_abs=norm(calibrated-alternate),
                              two_reference_columns_difference_relative=norm(calibrated-alternate)/selected,
                              refit_before_after_empty_difference_abs=norm(calibrated-after),
                              refit_before_after_empty_difference_relative=norm(calibrated-after)/selected,
                              two_reference_columns_difference_all_observed=norm(calibrated-alternate,MASK),
                              refit_before_after_empty_difference_all_observed=norm(calibrated-after,MASK),
                              nominal40_field_error_abs=norm(coarse-calibrated),
                              nominal80_field_error_abs=norm(fine-calibrated),
                              nominal40_rel=norm(coarse-calibrated)/selected,
                              nominal80_rel=norm(fine-calibrated)/selected,
                              nominal_grid_difference_over_scatter=norm(coarse-fine)/selected,
                              nominal80_predicted_norm=norm(fine), nominal_model_fit_gain_used=False,
                              nominal_state_solves=[dict(side=SIDE,**cmeta),dict(side=2*SIDE,**fmeta)]))
        row.update(elapsed(frequency_start))
        rows.append(row)
        print(frequency,"drift", drift_abs,"holdout",row["held_metal_relative_to_measured"],
              "opposite sign",row["opposite_sign_held_residual"]/norm(held), flush=True)
    np.savez_compressed(prepared,**arrays)
    report=dict(source_archive=raw.name,archive_size_bytes=raw.stat().st_size,members=archive_members,
                used_port_pairs=int(USED.sum()), raw_port_pairs_per_frequency=int(MASK.sum()),
                groups=[dict(receive_ids=(rx+1).tolist(),transmit_ids=(tx+1).tolist()) for rx,tx in GROUPS],
                frequency_audits=rows,cases=cases,nominal_full_wave_factorizations=2*len(cases),
                nominal_full_wave_columns=48*len(cases), **elapsed(start))
    physical_output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw",type=Path,default=RAW)
    parser.add_argument("--prepared",type=Path,default=PREPARED)
    parser.add_argument("--output",type=Path,default=HERE/"physical_audit.json")
    args=parser.parse_args()
    with threadpool_limits(1):
        main(args.raw,args.output,args.prepared)
