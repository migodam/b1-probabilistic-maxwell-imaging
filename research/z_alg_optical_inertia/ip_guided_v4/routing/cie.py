#!/usr/bin/env python3
"""Inherited native Dual-CIE with equal-width current-basis and contraction controls."""
from argparse import ArgumentParser
import csv
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
from scipy.linalg import qr, svd, subspace_angles
from threadpoolctl import threadpool_limits

from routing_run import (HERE, STUDY, clock, elapsed, full_audit, native_suite, Geometry, Observation,
                 read, write)
from cie_seed import Adapter, Regularizer
import engine

BASE = STUDY/'ip_unified/numerics/native'
ARRAY = STUDY/'physical_revision/array'
ORIGINAL_FOURIER = engine.FourierCurrents


def more_observations():
    path = ARRAY/'observations'
    with np.load(STUDY/'synthetic/public/fine_geometry.npz',allow_pickle=False) as data:
        spec = data['spectrum']
    with np.load(path/'arrays.npz',allow_pickle=False) as data:
        geo = Geometry(data['S32'],data['B32'],spec,96,2.,400e6)
    for r in read(path/'manifest.json')['observations']:
        if r['object_id'] not in ('optfresh_01','optfresh_05'):
            continue
        if (r['sources'],r['receivers'],r['noise_pct'])!=(32,32,2):
            continue
        f = path/r['file']
        with np.load(f,allow_pickle=False) as data:
            yield geo,Observation(r['id'],r['object_id'],'optfresh','native_diagnostic',str(f),
                                  data['Y'],r['padded_radius'],2.)


def basis_type(active, keep_arrays):
    class SwappedFourier(ORIGINAL_FOURIER):
        def __init__(self,shape,mf,groups,cost):
            t = clock()
            super().__init__(shape,mf,groups,cost)
            available = self.k
            replace = 0 if active is None else active.shape[1]
            f = np.stack([np.minimum(self.frequencies[:,k],shape[k]-self.frequencies[:,k])
                          for k in range(2)],axis=1)
            order = np.argsort(np.sum(f*f,axis=1),kind='stable')
            self.outside = order[available-replace:]
            self.inside = order[:available-replace]
            self.directions = []
            for ids,op in self.groups:
                if replace or keep_arrays:
                    projected = op.apply(np.eye(self.k,dtype=complex))
                    if keep_arrays:
                        keep_arrays['standard_U'] = projected.copy()
                if replace:
                    qr_keep, _, _ = qr(projected[:,self.inside],mode='economic',pivoting=True,check_finite=False)
                    extra = active-op.Vd@(op.Vd.conj().T@active)
                    extra -= qr_keep@(qr_keep.conj().T@extra)
                    angle = np.degrees(subspace_angles(projected,active))
                    u,s,_ = svd(extra,full_matrices=False,check_finite=False)
                    if len(s)!=replace or np.min(s)<1e-10:
                        raise ArithmeticError(f'Equal-width swap loses rank: {s}')
                    self.directions.append((ids,u))
                    cost['subspace_swap_qr_columns'] += len(self.inside)
                    cost['subspace_swap_replace_columns'] += replace
                    keep_arrays['swap_raw_complement_singular_values'] = s
                    keep_arrays['swap_full_fourier_angles_deg'] = angle
                    new_basis = projected.copy()
                    new_basis[:,self.outside] = u
                    keep_arrays['swapped_U'] = new_basis
                    keep_arrays['removed_frequencies'] = self.frequencies[self.outside]
            cost['subspace_swap_wall_seconds'] += elapsed(t)['wall_seconds'] if replace else 0.

        def apply(self, coefficients):
            if not self.directions:
                return super().apply(coefficients)
            base = coefficients.copy()
            base[self.outside] = 0
            j = super().apply(base)
            for ids,u in self.directions:
                j[:,ids] += u@coefficients[np.ix_(self.outside,ids)]
            return j

        def adjoint(self,j):
            v = super().adjoint(j)
            for ids,u in self.directions:
                v[np.ix_(self.outside,ids)] = u.conj().T@j[:,ids]
            return v

    return SwappedFourier


def current_metrics(current, v, lam, center):
    projection = v.conj().T@current
    centered = v.conj().T@(current-center)
    return dict(raw_negative_mode_projection_norm=float(np.linalg.norm(projection)),
                raw_negative_mode_projection_energy_fraction=float(np.linalg.norm(projection)**2/
                    max(np.linalg.norm(current)**2,1e-300)),
                weighted_centered_negative_coordinate_energy=float(np.sum(abs(lam)[:,None]*abs(centered)**2)),
                raw_negative_projection_energy_by_basis_column=np.sum(abs(projection)**2,axis=1).tolist())


def branch(geo, obs, tag, initial, rank, steps, beta, active=None, diag=None, mf=5,
           upper=20., full_schedule=False):
    start = clock()
    if geo.dense:
        green = SimpleNamespace(apply=lambda x: geo.g@x,adjoint=lambda x: geo.g.conj().T@x)
        model = SimpleNamespace(clock=engine.legacy.Clock(120.),E=geo.b,n=len(geo.b),
                                ns=geo.b.shape[1],shape=(geo.side,geo.side),ports=engine.legacy.Ports(geo.s),
                                green=green,reg=Regularizer((geo.side,geo.side)))
    else:
        model = Adapter(geo.g,geo.s,geo.b,obs.y,120.)
    state, history, basis, auxiliary = {}, [], {}, {}
    counts = model.clock.counts
    if diag:
        # First finite-dimensional current model used by the unchanged CIE loop.
        jd, groups = engine.dpcs(model,obs.y,rank,False)
        basis['Jd'] = jd
        basis['Vd'] = groups[0][1]

    def checkpoint(mark, step, stage, values):
        if mark == 'before':
            # A copy of what the unmodified solver currently represents, not a new fit.
            state['current'] = values['J'].copy()

    if diag:
        model.checkpoint_hook = checkpoint
    modes, betas = ([5,8,10],[8.,2.,.5]) if full_schedule else ([mf],[beta])
    schedule = dict(public_schedule=dict(modes=modes,betas=betas,rank=rank))
    schedule_end = 'completed first native Dual-CIE stage'
    engine.FourierCurrents = ORIGINAL_FOURIER if not diag else basis_type(active,basis)
    from research.research_max_struct_r2.src import cie_global
    original_readout = cie_global._baseline_readout

    def limited_readout(j,q,m,beta_material,source_weights,gain,box,*args):
        actual_box = box[0],np.minimum(box[1],upper)
        return original_readout(j,q,m,beta_material,source_weights,gain,actual_box,*args)

    if upper!=20.:
        cie_global._baseline_readout = limited_readout
    try:
        for m,count,stage,at_end in engine.run_cie_arm(model,obs.y,schedule,'dual',initial,obs.radius):
            if diag and (count in {4,8,16,steps} or at_end):
                auxiliary[str(count)] = m.copy()
            if count>=steps:
                schedule_end = 'fixed early-stage PR-current iteration limit'
                break
    except engine.legacy.TimeLimit:
        schedule_end = 'historical solver wall limit'
        if 'm' not in locals():
            return dict(failure=schedule_end,counts=dict(counts),**elapsed(start)),{}
    finally:
        engine.FourierCurrents = ORIGINAL_FOURIER
        cie_global._baseline_readout = original_readout
    work = elapsed(start)
    raw = np.asarray(m).copy()
    public_geo = Geometry(geo.s,geo.b,geo.g,geo.side,geo.width,geo.frequency,dense=geo.dense,hi=upper)
    check_points = auxiliary if diag else {str(count):raw}
    arrays = {k:v for k,v in basis.items() if (tag=='standard' or k.startswith(('swap','removed')))}
    v,lam,center = diag if diag else (None,None,None)
    for k,proposal in check_points.items():
        t = clock()
        arrays[k+'_material'] = proposal
        try:
            result,p,j = full_audit(public_geo,obs.y,proposal)
        except ArithmeticError as exc:
            recorded = dict(step=int(k),native=dict(failure=str(exc),**elapsed(t)))
        else:
            recorded = dict(step=int(k),native=result)
            arrays[k+'_ports'] = p
            if diag and k==str(count):
                arrays[k+'_physical_current'] = j
            if diag:
                recorded['native_mode_contents'] = current_metrics(j,v,lam,center)
        if diag:
            if k==str(count):
                clipped = np.minimum(proposal,1.)
                clip_result,p_clip,j_clip = full_audit(public_geo,obs.y,clipped)
                recorded['clipped_to_cap_one'] = dict(audit=clip_result,**current_metrics(j_clip,v,lam,center))
                arrays.update(cap1_material=clipped,cap1_ports=p_clip,cap1_current=j_clip)
                if 'current' in state:
                    recorded['auxiliary_mode_contents'] = current_metrics(state['current'],v,lam,center)
                    arrays['auxiliary_current'] = state['current']
        history.append(dict(**recorded,**elapsed(t)))
    return dict(arm=tag,method='inherited SRAP sourcewise Dual-CIE',step=count,stop=schedule_end,
                nominal_fourier_columns=(2*mf)**2,rank=rank,beta=beta,upper_cap=upper,
                full_schedule=full_schedule,counts=dict(counts),
                trace=model.trace_route, auxiliary_only=work,audits=history,**elapsed(start)),arrays


def diagnostic(out, objects, extra):
    with np.load(BASE/'u100_best.npz') as data:
        v,lam,centers = [data[k] for k in ('negative_vec','negative_eig','center')]
    weights = read(BASE/'u100_best.json')['profile']
    # Archive K eigensystem is in the K/scale representation.
    eig_scale = weights['scale']
    req = {r['id']:r['best'] for r in read(BASE/'required_currents.json')['rows']}
    with (ARRAY/'array_scan.csv').open() as stream:
        bounds = [r for r in csv.DictReader(stream) if r['strictly_excludes']=='True']
    for geo,obs in more_observations():
        if objects and obs.object_id not in objects:
            continue
        index = int(obs.object_id.split('_')[-1])-1
        info = (v,lam*eig_scale,centers[index])
        rows = [r for r in bounds if r['id']==obs.id]
        cap = max(float(r['tested_cap']) for r in rows)
        if cap!=.75:
            raise ValueError(f'Expected strict .75 exclusion on matching geometry: {obs.id} {cap}')
        package = dict(id=obs.id,object_id=obs.object_id,required_current=req[obs.id],
                       beta_peak_floor=cap,peak_floor_reuse_source='physical_revision/array/array_scan.csv',
                       n_rx=32,n_tx=32,noise_pct=2,forward_grid=96,radius=obs.radius,methods={})
        methods = [('standard',8.,None),('equal_width_negative_swap',8.,v),('beta_inverse_excluded_peak',1./cap,None)]
        if index%2==0:
            methods[0],methods[1] = methods[1],methods[0]
        if extra:
            package = read(out/(obs.id+'.json'))
            methods = [('beta2_no_certificate',2.,None)]
        for tag,beta,active in methods:
            result,raw = branch(geo,obs,tag,np.zeros(len(geo.b)),20,16,beta,active,info)
            package['methods'][tag]=result
            np.savez_compressed(out/(obs.id+'__'+tag+'.npz'),Y=obs.y,**raw)
            audit = result['audits'][-1]['native']
            status = audit['residual']/obs.radius if 'residual' in audit else audit['failure']
            print(obs.id,tag,result['step'],status,flush=True)
        write(out/(obs.id+'.json'),package)


def principal(out, routing, box4, objects, steps, full_schedule):
    for geo,obs in native_suite(only=objects):
        if not (routing/(obs.id+'.json')).exists():
            print(f'Skipping pending Born seat {obs.id}',flush=True)
            continue
        cases = []
        with np.load(routing/(obs.id+'.npz'),allow_pickle=False) as a:
            arms = [('cie_from_zero',np.zeros(len(geo.b))),('cie_from_born',a['born_only_material'])]
            if sum(map(ord,obs.id))%2==0:
                arms = arms[::-1]
            for tag,initial in arms:
                result,raw = branch(geo,obs,tag,initial,15,steps,8.,upper=4. if box4 else 20.,
                                    full_schedule=full_schedule)
                cases.append(result)
                np.savez_compressed(out/(obs.id+'__'+tag+'.npz'),Y=obs.y,**raw)
                print(obs.id,tag,result.get('step'),result.get('audits',[{}])[-1].get('native',{}).get('residual',-1)/obs.radius,
                      flush=True)
        write(out/(obs.id+'.json'),dict(id=obs.id,object_id=obs.object_id,radius=obs.radius,arms=cases))


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--branch',choices=['diagnostic','principal'],required=True)
    parser.add_argument('--routing',type=Path,default=HERE/'online')
    parser.add_argument('--out',type=Path)
    parser.add_argument('--objects',nargs='+')
    parser.add_argument('--box4',action='store_true')
    parser.add_argument('--steps',type=int,default=16)
    parser.add_argument('--full-schedule',action='store_true')
    parser.add_argument('--diagnostic-extra',action='store_true')
    args = parser.parse_args()
    out = args.out or HERE/('cie_'+args.branch)
    out.mkdir(parents=True,exist_ok=True)
    if args.branch=='diagnostic':
        diagnostic(out,args.objects,args.diagnostic_extra)
    else:
        principal(out,args.routing,args.box4,args.objects,args.steps,args.full_schedule)


if __name__ == '__main__':
    with threadpool_limits(1):
        main()
