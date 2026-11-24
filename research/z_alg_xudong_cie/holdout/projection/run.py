"""Offline orthogonal projection of the sealed native materials into the 32² lift."""

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from flint import arb, ctx
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[3]
STUDY = HERE.parents[1]
sys.path.insert(0, str(WORKSPACE))
from research.z_alg_xudong_cie.holdout.data import (
    APERTURES, DEFAULT_ROOT, INVERSE_SIDE, NATIVE_SIDE, OBJECT_IDS, elapsed, get_native_model,
    observation_id, public_directory, public_geometry, require_freeze, timers, write_json,
)
from research.z_alg_xudong_cie.holdout.physics import NativeForward, NativeForwardError, nearest_lift
from research.z_alg_xudong_cie.holdout.scoring import material_scores, relative_error

DEFAULT_SEAL = STUDY / 'benchmark/ONLINE_SEAL.json'
NOISE = (.02, .10)


def layout_check(root):
    """Confirm row-major y,x tiles and the same coordinate convention as the frozen lift."""
    with np.load(public_geometry(root, INVERSE_SIDE), allow_pickle=False) as f:
        coarse = f['centers'].reshape(INVERSE_SIDE, INVERSE_SIDE, 2)
    with np.load(public_geometry(root, NATIVE_SIDE), allow_pickle=False) as f:
        native = f['centers'].reshape(NATIVE_SIDE, NATIVE_SIDE, 2)
    a = np.mean(native.reshape(INVERSE_SIDE, 3, INVERSE_SIDE, 3, 2), axis=(1, 3))
    error = float(np.max(np.abs(a - coarse)))
    for centers in (coarse, native):
        if not (np.all(np.diff(centers[0, :, 0]) > 0) and np.all(np.diff(centers[:, 0, 1]) > 0)
                and np.array_equal(centers[:, :, 0], np.broadcast_to(centers[0, :, 0], centers.shape[:2]))
                and np.array_equal(centers[:, :, 1], np.broadcast_to(centers[:, :1, 1], centers.shape[:2]))):
            raise ValueError('Unexpected native cell order.')
    if error > 1e-14:
        raise ValueError('The fine tile centers do not agree with the frozen coarse lift.')
    return dict(reshape_order='C: each row runs along x; rows run along y', native_block=[3, 3],
                maximum_coarse_vs_block_coordinate_difference=error)


def interval_endpoints(x):
    return float(np.nextafter(float(x.lower()), -np.inf)), float(np.nextafter(float(x.upper()), np.inf))


def floor(material, binary_mean):
    """Outward L2 floor for the exact mathematical mean of saved binary64 truth."""
    tiles = material.reshape(32, 3, 32, 3).transpose(0, 2, 1, 3).reshape(-1, 9)
    numerator, denominator, roundoff = arb(0), arb(0), arb(0)
    for row, midpoint in zip(tiles, binary_mean, strict=True):
        values = [arb(float(x)) for x in row]
        mean = sum(values, arb(0))/9
        numerator += sum(((x-mean)**2 for x in values), arb(0))
        denominator += sum((x*x for x in values), arb(0))
        roundoff = max(roundoff, abs(mean-arb(float(midpoint))).upper())
    return interval_endpoints((numerator/denominator).sqrt()), interval_endpoints(roundoff)[1]


def union_ports(root):
    named = {}
    for ap in APERTURES:
        p = get_native_model(ap, root)
        named[ap] = p['S'], p['E']
    with np.load(root/'offline/hidden_geometry96.npz', allow_pickle=False) as hidden:
        named['held'] = hidden['S'], hidden['B']
    sr = np.cumsum([s.shape[0] for s,b in named.values()])
    bs = np.cumsum([b.shape[1] for s,b in named.values()])
    unique_s, rx = np.unique(np.concatenate([s for s,b in named.values()], axis=0), axis=0, return_inverse=True)
    unique_b, tx = np.unique(np.concatenate([b.T for s,b in named.values()], axis=0), axis=0, return_inverse=True)
    receiver = dict(zip(named, np.split(rx, sr[:-1]), strict=True))
    source = dict(zip(named, np.split(tx, bs[:-1]), strict=True))
    return p['spectrum'], unique_s, unique_b.T, receiver, source


def get_observation(root, name):
    with np.load(public_directory(root)/'observations'/f'{name}.npz', allow_pickle=False) as f:
        return f['Y'], float(f['delta']), float(f['raw_radius'])


def write_csv(dst, rows):
    cols = list(dict.fromkeys(k for row in rows for k in row))
    with dst.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=cols)
        writer.writeheader()
        writer.writerows(rows)


def certify_projection(root, output, chosen, material, geometry, receiver, source):
    """Independently bound one concrete finite-state projection after scoring."""
    start = timers()
    sys.path.insert(0, str(STUDY/'minpeak_native'))
    from verify import ExactFourierGeometry, ExactGeometry, ball, norm, verified_ports
    from nonvacuity import hi, lo

    ctx.prec, ctx.threads = 128, 1
    spec, s, b = geometry
    exact = ExactFourierGeometry(public_geometry(root, NATIVE_SIDE))
    # Use precisely the same union receiver and incident bases as the floating replay.
    exact.S, exact.B = s, b
    fft = ExactGeometry(s, b, spec, [])
    prep = elapsed(start)
    ports, error, details = verified_ports(exact, fft, material)
    results = []
    oid = chosen['object_id']
    for ap in APERTURES:
        rr, cc = receiver[ap], source[ap]
        from flint import acb_mat
        value = acb_mat(len(rr), len(cc), [ports[int(i), int(j)] for i in rr for j in cc])
        for noise in NOISE:
            name = observation_id(oid, ap, noise)
            y, r, _ = get_observation(root, name)
            delta = norm(value-ball(y))
            lower = max(0., lo(delta.lower()-error))
            upper = hi(delta.upper()+error)
            status = 'inside' if upper <= r else ('outside' if lower > r else 'unresolved')
            results.append(dict(observation_id=name, residual_lower=lower, residual_upper=upper,
                                padded_radius=r, strict_margin_lower=lo(arb(r)-arb(upper)),
                                status=status))
    payload = dict(method='Independent active-state verified inverse and original Fourier FFT defect',
                   object_id=oid, candidate_archive='projections.npz', candidate_key=f'{oid}__m32',
                   selection_rule='After completing the entire offline batch, smallest numerical residual-to-radius ratio',
                   numerical_selection=chosen, geometry_preparation=prep,
                   validated_affine_nearest_lift=True, interval_bits=128,
                   observations=results, details=details, **elapsed(start))
    write_json(output/'strict_projection.json', payload)
    return payload


def run(root, output, seal, do_proof):
    start = timers()
    seal = require_freeze(seal)
    output.mkdir(parents=True, exist_ok=True)
    order = layout_check(root)
    spec, rs, bs, rx, tx = union_ports(root)
    native = NativeForward(spec, NATIVE_SIDE)
    geometry = dict(unique_receivers=rs.shape[0], unique_illuminations=bs.shape[1],
                    shared_setup=elapsed(start), **order)
    ctx.prec, ctx.threads = 128, 1
    rows, fit, saved = [], [], {}
    for oid in OBJECT_IDS:
        case_start = timers()
        with np.load(root/'offline/objects'/f'{oid}.npz', allow_pickle=False) as p:
            truth, reference_coarse = p['m_native'], p['m_coarse']
            clean = {ap: p[f'clean_{ap}'] for ap in APERTURES}
            clean['held'] = p['clean_held']
        m32 = truth.reshape(32, 3, 32, 3).mean(axis=(1, 3)).ravel()
        m96 = nearest_lift(m32)
        f, bar_roundoff = floor(truth, m32)
        metrics = material_scores(m96, truth)
        discrepancy = float(np.linalg.norm(m32-reference_coarse))
        row = dict(object_id=oid, exact_lift_class=True,
                   actual_projection_relative_L2=metrics['rel_L2'], relative_L2_floor_lower=f[0],
                   relative_L2_floor_upper=f[1], maximum_block_mean_roundoff_upper=bar_roundoff,
                   difference_from_saved_coarse_L2=discrepancy,
                   difference_from_saved_coarse_Linf=float(np.max(np.abs(m32-reference_coarse))),
                   inside_L2=metrics['inside_L2'], background_L2=metrics['background_L2'],
                   half_IoU=metrics['half_IoU'], peak_contrast=metrics['peak_contrast'],
                   true_peak_contrast=metrics['truth_peak_contrast'])
        saved[oid+'__m32'] = m32
        saved[oid+'__m96'] = m96
        try:
            response, work = native.forward(m96, bs, rs)
        except NativeForwardError as exc:
            rows.append({**row, 'forward_status':'unresolved','forward_error':str(exc),'native_cost':exc.cost,
                         **elapsed(case_start)})
            print(oid,'unresolved',str(exc),flush=True)
            continue
        projections = {ap:response[np.ix_(rx[ap],tx[ap])] for ap in rx}
        for ap, prediction in projections.items():
            saved[f'{oid}__{ap}_prediction'] = prediction
            clean_error = relative_error(prediction,clean[ap])
            tag = ap.lower().removesuffix('_10x20')
            row[f'{tag}_relative_clean_port_error'] = clean_error['relative_to_measured']
            row[f'{tag}_absolute_clean_port_error'] = clean_error['absolute_Frobenius']
            if ap=='held':
                continue
            for noise in NOISE:
                name = observation_id(oid,ap,noise)
                y, radius, raw_radius = get_observation(root,name)
                scores = relative_error(prediction,y)
                fit.append(dict(object_id=oid, observation_id=name, aperture=ap, noise_fraction=noise,
                                padded_radius=radius, raw_radius=raw_radius, **scores,
                                residual_over_radius=scores['absolute_Frobenius']/radius,
                                numerical_inside_ball=scores['absolute_Frobenius']<=radius))
        row.update(forward_status='scored', native_forward_rhs=work['forward_rhs'],
                   native_state_relative_residual=work['fine_state_residual'],
                   native_fft_green_columns=work['native_fft_green_columns'],
                   native_forward_wall_seconds=work['wall_seconds'], native_forward_cpu_seconds=work['cpu_seconds'],
                   formulation=work['method'], native_dense_LU=work['dense_LU'],
                   native_cost=work,**elapsed(case_start))
        rows.append(row)
        current_ratios = [r['residual_over_radius'] for r in fit if r['object_id']==oid]
        print(oid, 'floor', f[0], 'ratios', current_ratios, 'held error',row['held_relative_clean_port_error'],
              'state defect', row['native_state_relative_residual'],flush=True)
    np.savez_compressed(output/'projections.npz',**saved)
    write_csv(output/'objects.csv',[{k:v for k,v in row.items() if k!='native_cost'} for row in rows])
    write_csv(output/'fitting.csv',fit)
    numeric = dict(scope='Offline block L2 projection, performed only after the entire B3 online seal',
                   online_seal=str(seal), geometry=geometry, objects=rows, fitting=fit,
                   num_numerical_inside=sum(x['numerical_inside_ball'] for x in fit),
                   batch_including_material_floors=elapsed(start))
    write_json(output/'projection_receipt.json',numeric)
    if len(fit)!=4*len(OBJECT_IDS):
        raise ArithmeticError('Some projected materials lack full native port evaluation.')
    low=min(fit,key=lambda row:row['residual_over_radius'])
    print('NUMERICAL BATCH', len(fit),numeric['num_numerical_inside'],'best',low,flush=True)
    if do_proof:
        strict = certify_projection(root,output,low,saved[low['object_id']+'__m96'],(spec,rs,bs),rx,tx)
        print('STRICT WITNESS',strict['object_id'],[(r['observation_id'],r['status'],r['strict_margin_lower'])
                                                  for r in strict['observations']],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=DEFAULT_ROOT)
    p.add_argument('--output',type=Path,default=HERE)
    p.add_argument('--online-seal',type=Path,default=DEFAULT_SEAL)
    p.add_argument('--no-witness',action='store_true')
    args=p.parse_args()
    with threadpool_limits(1):
        run(args.root,args.output,args.online_seal,not args.no_witness)
