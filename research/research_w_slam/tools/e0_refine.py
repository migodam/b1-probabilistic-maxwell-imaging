#!/usr/bin/env python3
"""Frozen numerical repair01; paid new reference solves, reused primary means."""
from pathlib import Path
from dataclasses import replace
import gc
import json
import time
import numpy as np
from wslam.accounting import BudgetBook,json_default
from wslam.world import World
from wslam.forward import Evaluator,volume_extraction
from wslam.fem_backend import PMLConfig
from wslam.physics import pack

ROOT=Path(__file__).resolve().parents[1]


def white(a):return float(np.linalg.norm(pack(a)))


def from_metadata(m):
    t=m['target_intrinsic']
    return World(m['support'],m['epsinf'],m['ell'],m['target_pose'],
                 radius=t['radius'],contrast=t['contrast'],target_loss=t['loss'])


def main():
    if not (ROOT/'docs/NUMERICAL_REPAIR_01.json').exists():
        raise RuntimeError('Numerical repair plan must exist before execution')
    prior=ROOT/'results/e0/qualification'
    old=json.loads((prior/'RESULTS.json').read_text())
    cases=[r for r in old['rows'] if 'point_grid_white' in r]
    if len(cases)!=23:raise RuntimeError('Expected exactly23 original fixed field cases')
    out=ROOT/'results/e0/refinement01_executable';out.mkdir(parents=True,exist_ok=False)
    book=BudgetBook();rows=[];started=time.perf_counter()
    def record(row):
        rows.append(row)
        (out/'RESULTS.json').write_text(json.dumps(dict(status='METRICS_ONLY_PARENT_REVIEW_REQUIRED',
            rows=rows,counts=book.totals('E0E1'),wall_seconds=time.perf_counter()-started,
            old_failed_costs_preserved=True,overall_gate='NOT_ASSIGNED'),
            indent=2,default=json_default,allow_nan=False))
        print(json.dumps(row,default=json_default,allow_nan=False),flush=True)
    def saved(label,ev):
        np.savez_compressed(out/(label+'.npz'),data=ev.data,volume_data=volume_extraction(ev),
                            tx=ev.backend.tx,rx=ev.backend.rx)
    for row in cases:
        label=row['case'];world=from_metadata(row['world']);freq=row['frequency']
        path=prior/(label+'_independent.npz')
        with np.load(path,allow_pickle=False) as z:primary=z['data'].copy()
        solver=Evaluator(book,'E0E1',h=.17,n=384,grid_refinement=3)
        reference=solver.evaluate(world,freq)
        saved(label+'_reference',reference)
        delta=white(primary-reference.data)
        record(dict(case=label,frequency=freq,world=row['world'],
            original_primary=row['primary'],original_independent=row['independent'],
            original_point_grid_white=row['point_grid_white'],
            new_primary=dict(source=str(path.relative_to(ROOT)),h=.17,n=256,grid_refinement=2,
                             reused_paid_means=True,new_rhs=0,new_LU=0,**{
                                 k:row['independent'][k] for k in ('ndofs','residual','point_volume_white')}),
            independent_reference=dict(h=.17,n=384,grid_refinement=3,ndofs=reference.backend.ndofs,
                residual=reference.state.residual,
                point_volume_white=white(volume_extraction(reference)-reference.data)),
            point_grid_white=delta,point_grid_exceeds_point1=delta>.1))
        solver.clear();del reference;gc.collect()
    # Fixed worst initial strong case, selected from initial grid error only.
    worst=max(cases,key=lambda r:r['point_grid_white'])
    world=from_metadata(worst['world']);freq=worst['frequency']
    with np.load(prior/(worst['case']+'_independent.npz'),allow_pickle=False) as z:primary=z['data'].copy()
    for label,pm in [('strength6',PMLConfig(strength=6)),('thickness1p25',PMLConfig(thickness=1.25))]:
        solver=Evaluator(book,'E0E1',h=.17,n=256,grid_refinement=2,pml=pm)
        ev=solver.evaluate(world,freq);saved('pml_'+label,ev)
        delta=white(ev.data-primary)
        record(dict(case='pml_'+label,original_case=worst['case'],frequency=freq,
            point_difference_white=delta,exceeds_point1=delta>.1,residual=ev.state.residual,
            pml=pm.as_dict(),baseline='reused original independent data'))
        solver.clear();del ev;gc.collect()
    base=World(np.ones((2,2),bool),np.array([[3.,4.],[3.,4.]]),.08,(.2,2.))
    solver=Evaluator(book,'E0E1',h=.17,n=256,grid_refinement=2)
    ev=solver.evaluate(base,1.,jacobian=True,kind='left_right')
    np.savez_compressed(out/'JACOBIAN.npz',data=ev.data,jacobian=ev.jacobian,names=np.asarray(ev.names))
    def shifted(name,step):
        if name.startswith('pose_'):
            p=base.pose.copy();p[0 if name=='pose_x' else 1]+=step
            return replace(base,pose=p)
        return replace(base,radius=base.radius+step)
    for name,steps in [('pose_x',[.002,.001,.0005]),('pose_y',[.002,.001,.0005]),('target_radius',[.0005])]:
        col=ev.names.index(name)
        for step in steps:
            fdsolver=Evaluator(book,'E0E1',h=.17,n=256,grid_refinement=2)
            plus=fdsolver.evaluate(shifted(name,step));minus=fdsolver.evaluate(shifted(name,-step))
            fd=(plus.data-minus.data).ravel()/(2*step);j=ev.jacobian[:,col]
            np.savez_compressed(out/f'FD_{name}_{step}.npz',plus=plus.data,minus=minus.data,
                                derivative=fd,jacobian=j)
            record(dict(case='central_FD_new_primary',parameter=name,step=step,
                relative_error=float(np.linalg.norm(fd-j)/max(np.linalg.norm(fd),np.finfo(float).tiny)),
                derivative_difference_white=white(fd-j),
                maximum_residual=max(plus.state.residual,minus.state.residual),
                grid=dict(h=.17,n=256,grid_refinement=2)))
            fdsolver.clear();del plus,minus;gc.collect()
    solver.clear();del ev;gc.collect()
    record(dict(case='driver_complete',overall_gate='NOT_ASSIGNED_PARENT_REVIEW',
        carried_evidence=dict(source_Mie_and_real_VJP='results/e0/qualification/RESULTS.json',
                             coefficient_FD='results/e0/qualification/RESULTS.json'),
        scope='numerical repair only; no E1/A0 performance execution or contract rewrite'))


if __name__=='__main__':main()
