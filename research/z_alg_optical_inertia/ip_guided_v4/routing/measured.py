#!/usr/bin/env python3
"""Apply nonnegative Born fitting to the actual independently calibrated instruments."""
import csv
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from routing_run import HERE, STUDY, ROOT, BORN, clock, elapsed, module, read, write


def csv_rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def born_candidate(s,b,y,train,block):
    started = clock()
    side = int(np.sqrt(s.shape[1]))
    reduced = side//block
    selected = []
    measurement = []
    for i in range(b.shape[1]):
        indices = train[:,i]
        if not np.any(indices):
            continue
        a = s[indices]*b[:,i]
        selected.append(a.reshape(-1,reduced,block,reduced,block).sum(axis=(2,4)).reshape(-1,reduced**2))
        measurement.append(y[indices,i])
    a = np.concatenate(selected)
    data = np.concatenate(measurement)
    h = (a.conj().T@a).real
    lip = float(np.linalg.eigvalsh(h)[-1])
    common = elapsed(started)
    f = (a.conj().T@data).real
    r = .001*lip
    t = clock()
    z = x = np.zeros(len(h))
    alpha = 1.
    for _ in range(120):
        x_new = np.clip(z-(h@z-f+r*z)/(lip+r),0.,4.)
        a_new = .5*(1+np.sqrt(1+4*alpha**2))
        z = x_new+(alpha-1)/a_new*(x_new-x)
        alpha,x = a_new,x_new
    m = x.reshape(reduced,reduced).repeat(block,axis=0).repeat(block,axis=1).ravel()
    ports = s@(m[:,None]*b)
    return m,ports,dict(geometry_and_matrix=common,fit_and_prediction=elapsed(t),iterations=120,
                        native_grid=side,inversion_grid=reduced,material='real nonnegative [0,4]')


def fit(out,ident,target,freq,y,s,b,train,hold_tx,hold_rx,bound,levels,prep,block):
    if np.any(train&hold_tx) or np.any(train&hold_rx) or np.any(hold_tx&hold_rx):
        raise ArithmeticError(f'Non-disjoint fit and holdout port roles {ident}')
    m,p,meta = born_candidate(s,b,y,train,block)
    name = ident.replace('/','__')
    path = out/(name+'.npz')
    np.savez_compressed(path,material=m,complex_born_prediction=p,train_mask=train,
                        held_transmitter_mask=hold_tx,held_receiver_mask=hold_rx)
    with np.load(path,allow_pickle=False) as sealed:
        sealed_prediction = sealed['complex_born_prediction']
    groups = {}
    for label,idx in [('train',train),('held_source',hold_tx),('held_receiver',hold_rx)]:
        if np.any(idx):
            groups[label] = dict(channels=int(idx.sum()),relative_rms=float(
                np.linalg.norm(sealed_prediction[idx]-y[idx])/np.linalg.norm(y[idx])))
    result = dict(id=ident,object_id=target,frequency_GHz=freq,
                  previous_fixed_born_certificate=bound,previous_conditional_peak_margins=levels,
                  sigma=None,strict_loaded_target_radius=None,material_NRMSE=None,
                  fit_provenance=meta,calibration_and_operator=prep,port_fit_metrics=groups)
    write(out/(name+'.json'),result)
    print(ident,{k:round(v['relative_rms'],3) for k,v in groups.items()},flush=True)
    return result


def fresnel(out,bounds):
    source = module('routing_measured_fresnel',STUDY/'fresnel/physical.py')
    historically_verified = csv_rows(STUDY/'results/fresnel/SENSITIVITY.csv')
    rows = []
    for fr in source.FREQUENCIES:
        t = clock()
        op = source.operators(fr)
        op.pop('G')
        setup = elapsed(t)
        for target in source.NAMES:
            t = clock()
            d = source.load_case(target,fr,op,model_audit=False)
            with np.load(source.RESULTS/f'{target}_{fr}GHz_ports.npz') as f:
                htx,hrx,raw,fr_mask = (f[k] for k in ('nominal_held_tx','nominal_held_rx','raw_scattered','mask'))
            if not np.array_equal(fr_mask,d['mask']) or np.max(np.abs(raw[fr_mask]-d['Y'][fr_mask]*d['norm']))>1.e-13:
                raise ArithmeticError('Raw measured channel re-open does not match calibrated adapter.')
            train = d['mask'] & ~htx & ~hrx
            ident = d['record']['id']
            cs = []
            for r in historically_verified:
                if r['case']==target and float(r['GHz'])==fr and r['lower_outward']:
                    cs.append(dict(cap=float(r['cap_chi']),material_family=r['family'],
                                   unit='fraction of original actually observed Frobenius norm',
                                   distance_lower=float(r['lower_outward']),
                                   remaining_proxy_margin=float(r['bcrit_outward_empty_rms_and_diagnostic_gain'])
                                   if r['bcrit_outward_empty_rms_and_diagnostic_gain'] else None))
            rows.append(fit(out,ident,target,fr,d['Y'],op['S'],d['B'],train,htx,hrx,bounds[ident],cs,
                            dict(shared_per_frequency=setup,read_and_empty_calibration=elapsed(t)),2))
    return rows


def manitoba(out,bounds):
    archive = STUDY/'validation/manitoba'
    source = module('routing_measured_manitoba',archive/'physical.py')
    levels = csv_rows(STUDY/'validation/transfer/symmetric_margins.csv')
    result = []
    with np.load(source.PREPARED,allow_pickle=False) as a:
        for i,fr in enumerate(source.FREQUENCIES):
            t = clock()
            s,b = a[f'f{i}_S'],a[f'f{i}_B']
            prior_port_mask = source.USED
            htx = source.MASK & (np.arange(s.shape[0])%6==0)[None,:]
            hrx = source.MASK & ~htx & ~prior_port_mask
            train = prior_port_mask & ~htx
            for target in source.TARGETS:
                ident = f'manitoba_eigor2010/{target}/{fr:g}GHz'
                cs = [dict(cap=float(r['cap']),material_family=r['model'],
                           lower_distance=float(r['lower_distance']),
                           remaining_proxy_relative_margin=float(r['relative_margin_lower']))
                      for r in levels if r['target']==target and float(r['GHz'])==fr]
                result.append(fit(out,ident,target,float(fr),a[f'f{i}_{target}'],s,b,train,htx,hrx,bounds[ident],cs,
                                  dict(reopen_prepared_calibration=elapsed(t),calibration_paid_in_archived_receipt=
                                       'validation/manitoba/physical_audit.json'),2))
    return result


def main():
    t = clock()
    base = read(BORN/'ports/results.json')
    out = HERE/'laboratory'
    out.mkdir(exist_ok=True)
    for group,entry in [('manitoba',manitoba),('fresnel',fresnel)]:
        current = clock()
        bounds = {r['id']:r for r in base[group]['observations']}
        result = entry(out,bounds)
        write(out/(group+'__seal.json'),dict(model='fixed empty-field calibration and Born row operator; empirical loading envelopes',
                                           score_use='masked acquired receivers only; target held-out values read after predictions sealed',
                                           independent_count=len({r['object_id'] for r in result}),
                                           rows=[r['id'] for r in result],**elapsed(current)))
    write(out/'all__seal.json',dict(total_physical_objects=7,**elapsed(t)))


if __name__=='__main__':
    with threadpool_limits(1):
        main()
