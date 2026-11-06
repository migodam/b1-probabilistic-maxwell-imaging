"""Offline finite-step feasibility and fixed-shape low-amplitude controls.

GT is used only to measure representation ceilings; these are not online proposals.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import csv,json,sys,time
import numpy as np
from threadpoolctl import threadpool_limits
HERE=Path(__file__).resolve().parents[1];REPO=HERE.parents[1]
sys.path.insert(0,str(REPO))
from research.research_z_alg_r4_diversity.src.physics import Model
OLD=REPO/'research/research_z_alg_r4_diversity'
OUT=HERE/'results/finite'

def writecsv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)

def shrink(m,h):
    a=1.
    if np.any(h>0):a=min(a,float(np.min((20-m[h>0])/h[h>0])))
    if np.any(h<0):a=min(a,float(np.min((.02-m[h<0])/h[h<0])))
    a=max(0.,a)
    if a<1:a=float(np.nextafter(a,0.))
    return a*h,a

def main():
    started=time.perf_counter();OUT.mkdir(exist_ok=True)
    public=dict(np.load(OLD/'public_data/operators.npz',allow_pickle=False));model=Model(public)
    meta={r['task_id']:r for r in json.loads((OLD/'RESULTS_SUMMARY.json').read_text())['rows']}
    cases=list(csv.DictReader((HERE/'results/e1/cases.csv').open()))
    cases=[r for r in cases if r['method']=='native_trf' and float(r['amplitude'])==.4]
    ref=np.load(HERE/'results/e1/public_reference_geometry.npz');refbasis=ref['J_V'][:,-64:]
    finite=[];materials=[];indices=[]
    for case in cases:
        tid=case['task_id'];cache=np.load(HERE/'results/e1'/f'{tid}.npz',allow_pickle=False)
        m=cache['anchor'];err=cache['error_OFFLINE'];gt=m+err;Y=cache['Y'];yn=np.linalg.norm(Y)
        base=model.forward(m,label='finite_anchor');J=model.jacobian(base).copy();y0=base['y'].copy()
        before=float(np.linalg.norm(err)/24);eps=1e-4 if case['noise']=='clean' else 1e-3
        dirs={name:cache[name] for name in ('J_weakest_32','J_weakest_64','current_leakage_weakest_32','current_leakage_weakest_64','threshold_0.01')}
        dirs['reference_J_weakest_64']=refbasis@(refbasis.T@err)
        for name,p in dirs.items():
            prms=float(np.linalg.norm(p)/24)
            for radius in (.01,.03,.1):
                h=radius*p/prms if prms>0 else np.zeros(576)
                violation=int(np.count_nonzero((m+h<.02)|(m+h>20)))
                step,alpha=shrink(m,h);trial=model.forward(m+step,label='offline_projected_GT_ray')
                dy=trial['y']-y0;lin=J@step;packed=model.pack(dy)
                item=dict(task_id=tid,shape=case['shape'],noise=case['noise'],basis=name,requested_RMS=radius,
                    ray_scale=alpha,actual_RMS=float(np.linalg.norm(step)/24),unshrunk_box_violations=violation,
                    before_RMSE=before,after_RMSE=float(np.linalg.norm(trial['m']-gt)/24),
                    relative_data_drift=float(np.linalg.norm(dy)/yn),linear_relative_data_drift=float(np.linalg.norm(lin)/yn),
                    nonlinear_remainder_relative=float(np.linalg.norm(packed-lin)/yn),
                    anchor_data_residual=float(np.linalg.norm(y0-Y)/yn),trial_data_residual=float(np.linalg.norm(trial['y']-Y)/yn),
                    epsilon=eps,within_drift_budget=bool(np.linalg.norm(dy)/yn<=eps),classification='OFFLINE_GT_DIRECTED_REPRESENTATION_PROBE')
                finite.append(item);materials.append(trial['m'].copy());indices.append(tid+'|'+name+'|'+str(radius))
    writecsv(OUT/'finite_steps.csv',finite)
    np.savez_compressed(OUT/'finite_materials_OFFLINE.npz',key=np.array(indices),m=np.array(materials))
    counts_finite=dict(model.counts);wall_finite=time.perf_counter()-started
    controls=[];control_arrays={}
    fixed=model.forward(np.full(576,.1),label='public_background_control');model.jacobian(fixed)
    yref=fixed['y'].copy();jc=fixed['Jcomplex'].copy()
    for shape in ('digit7','ring'):
        original=next(r for r in meta.values() if r['backend']=='cpu' and r['method']=='native_trf' and r['shape']==shape and r['amplitude']==.4 and r['noise']=='clean')
        truth=np.load(OLD/'truth_after_freeze'/f"{original['observation_id']}.npz")
        for amplitude in (.1,.2,.4):
            m=original['background']+amplitude*truth['alpha'].ravel();state=model.forward(m,label='fixed_shape_low_control')
            y=state['y'];born=model.S@(m[:,None]*model.E0)
            linear=yref+(jc@(m-.1)).reshape(12,48).T
            controls.append(dict(shape=shape,amplitude=amplitude,background=original['background'],
                Born_full_relative_difference=float(np.linalg.norm(born-y)/np.linalg.norm(y)),
                background_linearization_relative_difference=float(np.linalg.norm(linear-y)/np.linalg.norm(y)),
                classification='DEVELOPMENT_FIXED_DISCLOSED_SHAPE_NO_NEW_INVERSION'))
            key=shape+'_'+str(amplitude);control_arrays[key+'_m']=m;control_arrays[key+'_Y']=y
    writecsv(OUT/'low_amplitude_controls.csv',controls);np.savez_compressed(OUT/'low_amplitude_controls.npz',**control_arrays)
    cost=dict(finite=counts_finite,controls={k:int(v-counts_finite.get(k,0)) for k,v in model.counts.items() if v-counts_finite.get(k,0)},total=dict(model.counts),finite_wall_seconds=wall_finite,total_wall_seconds=time.perf_counter()-started)
    (OUT/'RUN_SUMMARY.json').write_text(json.dumps(dict(status='COMPLETE',native_cases=len(cases),finite_trials=len(finite),controls=len(controls),cost=cost),indent=2))
    print(json.dumps(dict(native_cases=len(cases),finite_trials=len(finite),controls=len(controls),cost=cost)))

if __name__=='__main__':
    with threadpool_limits(limits=1):main()
