#!/usr/bin/env python3
"""Score sealed predictions against material and independent synthetic ports."""
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from routing_run import (HERE, STUDY, ROOT, RESEARCH, BORN, LOSS, clock, elapsed, full_audit,
                         read, write, offset_geometry, native_suite, other_suite, module)
from cie import more_observations


def material_metrics(m, reference, floor):
    # The threshold represents the generating background, which may differ from the solver floor.
    support = reference>floor+1.e-10
    delta = m-reference
    value = dict(full_relative_l2=float(np.linalg.norm(delta)/np.linalg.norm(reference)),
                 inside_relative_l2=float(np.linalg.norm(delta[support])/np.linalg.norm(reference[support])),
                 background_rms=float(np.sqrt(np.mean(delta[~support]**2))) if np.any(~support) else None,
                 peak=float(m.max()),peak_error=float(m.max()-reference.max()))
    return value


def predict(geo,m):
    y = np.zeros((len(geo.s),geo.b.shape[1]),complex)
    return full_audit(geo,y,m)[1]


def score_case(geo,obs,candidates,truth,clean=None):
    start = clock()
    base = obs.object_id
    if clean is None:
        clean = predict(geo,truth)
    if geo.frequency is None or '_born__' in obs.id:
        held = None
    else:
        held = offset_geometry(geo,3.75 if geo.side==96 else 1.8)
    clean_held = predict(held,truth) if held else None
    arrays = dict()
    if clean_held is not None:
        arrays['reference_held_ports'] = clean_held
    result = dict(id=obs.id,object_id=base,ground_truth_peak=float(truth.max()),
                  clean_train_fit_to_noisy_over_radius=float(np.linalg.norm(clean-obs.y)/obs.radius),
                  holdout_kind='half-step offset simulated transmit and receive angles, used after online seal' if held else None,
                  arms={})
    background = float(truth.min()) if obs.cohort=='r14_ring7' else geo.lo
    result['reference_background_contrast'] = background
    for label,m,p in candidates:
        metric = material_metrics(m,truth,background)
        metric['clean_train_error']=float(np.linalg.norm(p-clean)/np.linalg.norm(clean))
        if held:
            try:
                p_held = predict(held,m)
            except ArithmeticError as exc:
                metric.update(held_error=None,held_failure=str(exc))
            else:
                metric['held_error']=float(np.linalg.norm(p_held-clean_held)/np.linalg.norm(clean_held))
                arrays[label+'_held_ports']=p_held
        else:
            metric['held_error']=None
        result['arms'][label]=metric
    return dict(**result,**elapsed(start)),arrays


def candidate_paths(obs, online, cie_path):
    with np.load(online/(obs.id+'.npz'),allow_pickle=False) as a:
        if not np.array_equal(a['Y'],obs.y):
            raise ArithmeticError(f'{obs.id} sealed observation changed')
        for label in ('born_only','full_zero','born_then_full'):
            if label+'_material' in a.files:
                yield label,a[label+'_material'],a[label+'_ports']
    if cie_path:
        for label in ('cie_from_zero','cie_from_born'):
            f = cie_path/(obs.id+'__'+label+'.npz')
            if not f.exists():
                continue
            with np.load(f,allow_pickle=False) as a:
                if not np.array_equal(a['Y'],obs.y):
                    raise ArithmeticError(f'{obs.id} sealed CIE observation changed')
                m = [k for k in a.files if k.endswith('_material')]
                for key in m:
                    yield label,a[key],a[key.replace('_material','_ports')]


def reference_loader(suite):
    if suite=='native':
        w = RESEARCH/'srap_cie_array_study_20261121'
        sr = module('routing_data2d',w/'data2d.py')
        axis = (np.arange(96)+.5)/96*2-1
        yy,xx = np.meshgrid(axis,axis,indexing='ij')
        points = np.c_[xx.ravel(),yy.ravel()]
        labels = np.random.default_rng(sr.FAMILY_SEED).permutation(np.repeat(np.arange(6),4))
        def load(obs):
            if obs.cohort=='optfresh':
                with np.load(STUDY/'synthetic/evaluation'/(obs.object_id+'.npz')) as data:
                    return data['m_fine'],None
            idx = int(obs.object_id.split('_')[-1])
            m = sr.random_material(points,idx,labels[idx-1])
            f = w/'data/truth'/(obs.id+'.npz')
            with np.load(f,allow_pickle=False) as data:
                saved_coarse = data['m']
            if np.max(np.abs(sr.area_average(m,96,32)-saved_coarse))>5e-14:
                raise ArithmeticError(f'{obs.id} regenerated native material disagrees with archive')
            return m,None
        return load
    if suite=='ring':
        path = RESEARCH/('delegated/max_struct_r1/r14_worktree/research/research_max_struct_r14/'
                         'sealed_truth/independent24')
        def load(obs):
            with np.load(path/(obs.object_id+'.npz'),allow_pickle=False) as data:
                return data['material_forward'],data['Y_clean']
        return load
    if suite=='loss':
        def load(obs):
            i = int(obs.object_id.split('_')[-1])
            with np.load(LOSS/f'finite_{int(i>2)}_inputs.npz',allow_pickle=False) as a:
                return a[obs.object_id],None
        return load
    if suite=='weak':
        def load(obs):
            with np.load(BORN/'controls/arrays.npz',allow_pickle=False) as a:
                name = obs.object_id
                gen = 'born' if '_born__' in obs.id else 'clean'
                return a[name+'__chi'].ravel(),a[name+'__'+gen]
        return load
    if suite=='mild_mid':
        with np.load(HERE/'mild_mid_offline_reference.npz',allow_pickle=False) as a:
            ids = dict(zip(a['ids'],a['materials'],strict=True))
        return lambda obs:(ids[obs.id],None)
    raise ValueError(suite)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--suite',choices=['native','ring','loss','weak','mild_mid','diagnostic'],required=True)
    parser.add_argument('--online',type=Path,default=HERE/'online')
    parser.add_argument('--out',type=Path,default=HERE/'offline')
    parser.add_argument('--with-cie',action='store_true')
    parser.add_argument('--diagnostic-arm')
    args = parser.parse_args()
    out = args.out
    out.mkdir(parents=True,exist_ok=True)
    start = clock()
    if args.suite=='diagnostic':
        source = HERE/'cie_diagnostic'
        cases = more_observations()
        load = reference_loader('native')
        for geo,obs in cases:
            receipt = read(source/(obs.id+'.json'))
            options = []
            for tag,branch in receipt['methods'].items():
                if args.diagnostic_arm and tag!=args.diagnostic_arm:
                    continue
                with np.load(source/(obs.id+'__'+tag+'.npz'),allow_pickle=False) as a:
                    for step in [k for k in a.files if k.endswith('_material')]:
                        label = tag+'_'+step.replace('_material','')
                        p = step.replace('_material','_ports')
                        if p in a.files:
                            options.append((label,a[step],a[p]))
            result,ports = score_case(geo,obs,options,*load(obs))
            if args.diagnostic_arm:
                prior = read(out/(obs.id+'.json'))
                prior['arms'].update(result['arms'])
                result = prior
                with np.load(out/(obs.id+'_held_ports.npz'),allow_pickle=False) as previous:
                    ports = {**dict(previous),**ports}
            write(out/(obs.id+'.json'),result)
            np.savez_compressed(out/(obs.id+'_held_ports.npz'), **ports)
    else:
        seal = read(args.online/(args.suite+'__seal.json'))
        cases = native_suite() if args.suite=='native' else other_suite(args.suite)
        cie = HERE/'cie_principal_common_box' if args.with_cie else None
        load = reference_loader(args.suite)
        completed = []
        for geo,obs in cases:
            if obs.id not in seal['ids']:
                raise ArithmeticError(f'{obs.id} missing from online run seal')
            result,ports = score_case(geo,obs,candidate_paths(obs,args.online,cie),*load(obs))
            write(out/(obs.id+'.json'),result)
            if ports:
                np.savez_compressed(out/(obs.id+'_held_ports.npz'), **ports)
            print(obs.id,[round(a['full_relative_l2'],4) for a in result['arms'].values()],flush=True)
            completed.append(obs.id)
        write(out/(args.suite+'__seal.json'),dict(observations=completed,**elapsed(start)))


if __name__=='__main__':
    with threadpool_limits(1):
        main()
