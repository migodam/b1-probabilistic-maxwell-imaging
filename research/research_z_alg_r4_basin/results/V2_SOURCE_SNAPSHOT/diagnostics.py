#!/usr/bin/env python3
"""OFFLINE_DIAGNOSTIC: frozen digit truth used only for basin interpretation.
No optimization, no new blind claim, no formula accuracy or duplicate audit test.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import sys,json,time,csv,importlib.util
from pathlib import Path
import numpy as np
from scipy.linalg import svd,eigvalsh
from scipy.linalg.lapack import get_lapack_funcs
from threadpoolctl import threadpool_limits,threadpool_info
ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT.parent/'research_z_alg_r3_digits'
spec=importlib.util.spec_from_file_location('frozen_digits_physics',OLD/'src/physics.py')
physics=importlib.util.module_from_spec(spec);spec.loader.exec_module(physics)
CASES=[('DEV7','dev_v2','dev','dev_2a5c221387b1'),('holdout7','holdout','holdout','obs_c3604cf3555d'),('holdout9','holdout','holdout','obs_42c48e795766'),('holdout8','holdout','holdout','obs_a63df7f20566')]
METHODS=['fullj_native','fullj_reciprocal','globaldual']
def clean(v):
    if isinstance(v,np.ndarray):return clean(v.tolist())
    if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    if isinstance(v,np.generic):return clean(v.item())
    if isinstance(v,float) and not np.isfinite(v):return None
    return v

def write(p,obj):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(clean(obj),indent=2,allow_nan=False)+'\n')
def cosine(a,b):
    n=np.linalg.norm(a)*np.linalg.norm(b)
    return float(np.dot(a,b)/n) if n>1e-24 else None

def geometry(model,m,target,gt,fixed_h,label,reference=None):
    begun,cpu=time.perf_counter(),time.process_time();counts=model.counts.copy()
    state=model.forward(m,label=label);J=model.jacobian(state)
    r=model.pack(state['y']-target);scale=max(np.linalg.norm(target)**2,1e-300)
    with model.operation('DIAGNOSTIC_FULL_REAL_J_SVD',diagnostic_SVD=1):
        _,s,Vt=svd(J,full_matrices=False,check_finite=False,lapack_driver='gesdd')
    V=Vt[:32].T;energy=s*s;energy/=max(energy.sum(),1e-300)
    with model.operation('DIAGNOSTIC_SOURCE_GRAM_SPECTRUM',diagnostic_Gram_spectra=1):
        gram_eig=eigvalsh(state['E'].conj().T@state['E'],check_finite=False)
    psi,gram_grad,gram_details=model.smallgram(state,target,gradient=True)
    A=model.identity-m[:,None]*model.G;anorm=float(np.linalg.norm(A,1))
    with model.operation('DIAGNOSTIC_ZGECON',resolvent_condition_estimators=1):
        rcond,info=get_lapack_funcs('gecon',(A,))(model.factor[0],anorm,norm='1')
    grad=J.T@r/scale;to_gt=gt-m
    row=dict(label=label,status='OFFLINE_DIAGNOSTIC',relative_residual=float(np.linalg.norm(r)/np.sqrt(scale)),material_RMSE_to_GT=float(np.sqrt(np.mean(to_gt**2))),
        J_singular_values=s,J_sigma_max=float(s[0]),J_sigma_32=float(s[31]),J_sigma_min=float(s[-1]),J_stable_rank=float(np.sum(s*s)/(s[0]*s[0])),
        J_energy_top32=float(energy[:32].sum()),J_count_above_relative_1e_minus6=int(np.count_nonzero(s>s[0]*1e-6)),
        rank_scope='Numerical spectrum only, relative threshold count is not exact rank or identifiability certificate',
        source_Gram_eigenvalues=gram_eig,source_Gram_condition=float(gram_eig[-1]/gram_eig[0]),
        resolvent_A='I-diag(m)G',resolvent_matrix_norm1=anorm,resolvent_rcond1_estimate=float(rcond),
        resolvent_norm1_estimate=float(1/(anorm*rcond)),resolvent_gecon_info=int(info),
        estimator_internal_RHS='OPAQUE_LAPACK_ESTIMATOR_WORK_NOT_EXPOSED_NOT_ZERO',
        physical_gradient_norm=float(np.linalg.norm(grad)),physical_descent_GT_cosine=cosine(-grad,to_gt),
        smallgram_psi=float(psi),smallgram_descent_GT_cosine=cosine(-gram_grad,to_gt),
        physical_and_smallgram_descent_cosine=cosine(grad,gram_grad),
        GT_direction_top32_energy_fraction=float(np.linalg.norm(V.T@to_gt)**2/np.linalg.norm(to_gt)**2) if np.linalg.norm(to_gt)>1e-12 else None,
        field_norm=float(np.linalg.norm(state['E'])),state_relative_residual=state['state_relative_residual'],smallgram_details=gram_details)
    if reference is not None:
        with model.operation('DIAGNOSTIC_TOP32_ROTATION',diagnostic_subspace_rotations=1):
            cc=svd(reference.T@V,compute_uv=False,check_finite=False)
        angles=np.degrees(np.arccos(np.clip(cc,0,1)))
        row.update(top32_principal_angles_degrees=angles,top32_max_angle_degrees=float(angles.max()),top32_mean_sin_angle=float(np.sqrt(np.mean(1-cc*cc))))
    else:
        row.update(top32_principal_angles_degrees=np.zeros(32),top32_max_angle_degrees=0.,top32_mean_sin_angle=0.,rotation_reference_identity=True)
    probes={}
    negative=-grad
    hgrad=negative*(.1*np.sqrt(model.n)/max(np.linalg.norm(negative),1e-300))
    for name,h in [('fixed_GT_direction_delta_t_0.05',fixed_h),('fixed_physical_gradient_RMS_0.1',hgrad)]:
        h=np.clip(m+h,.5,20)-m
        trial=model.forward(m+h,label=label+'_'+name)
        dr=model.pack(trial['y']-state['y']);jh=J@h
        rn=model.pack(trial['y']-target)
        pred=float(.5*(r@r-(r+jh)@(r+jh))/scale)
        actual=float(.5*(r@r-rn@rn)/scale)
        probes[name]=dict(step_RMS=float(np.sqrt(np.mean(h*h))),true_remainder_norm=float(np.linalg.norm(dr-jh)),
            remainder_relative_to_linear=float(np.linalg.norm(dr-jh)/max(np.linalg.norm(jh),1e-300)),
            predicted_decrease=pred,actual_decrease=actual,actual_over_predicted=float(actual/pred) if abs(pred)>1e-20 else None,
            predicted_descent=bool(pred>0),actual_descent=bool(actual>0),remainder_definition='F(m+h)-F(m)-J(m)h')
    row['fixed_step_probes']=probes
    row['cost_counts']={k:model.counts[k]-counts[k] for k in model.counts if model.counts[k]!=counts[k]}
    row['wall_seconds']=time.perf_counter()-begun;row['cpu_seconds']=time.process_time()-cpu
    return row,V

def main():
    out=ROOT/'diagnostics';out.mkdir(exist_ok=True,parents=True)
    begun,cpu=time.perf_counter(),time.process_time()
    op=dict(np.load(OLD/'public_data/operators.npz'));model=physics.Model(op)
    allrows=[];traces=[]
    with threadpool_limits(limits=1):
        for name,truthsplit,publicsplit,oid in CASES:
            gt=np.load(OLD/f'truth_after_freeze/{truthsplit}/{oid}.npz')['m']
            target=np.load(OLD/f'public_data/{publicsplit}/{oid}.npz')['Y']
            m0=np.ones(model.n);delta=gt-m0;reference=None;rows=[]
            for t in [0,.2,.4,.6,.8,1.]:
                row,V=geometry(model,m0+t*delta,target,gt,.05*delta,f'{name}_path_t{t:.1f}',reference)
                row.update(case=name,observation_id=oid,kind='GT_INTERPOLATION_PATH',t=t,truth_usage='OFFLINE_DIAGNOSTIC_ONLY_NOT_AN_OPTIMIZER_TRAJECTORY')
                if reference is None:reference=V
                rows.append(row);allrows.append(row)
                write(out/f'{name}_path_partial.json',rows)
            for method in METHODS:
                file=OLD/f'results/{truthsplit}/{method}/{oid}.json'
                frozen=json.loads(file.read_text());m=np.asarray(frozen['top1']['m'])
                row,V=geometry(model,m,target,gt,.05*delta,f'{name}_endpoint_{method}',reference)
                row.update(case=name,observation_id=oid,kind='FROZEN_TRUE_OPTIMIZER_ENDPOINT',method=method,source=str(file.relative_to(ROOT.parent)),
                    historical_status=frozen['status'],historical_wall_seconds=frozen['wall_seconds'],historical_accounting=frozen.get('accounting'),
                    original_Y_residual_scope='Nonzero residual is an approximate fit, never another solution of original Y')
                rows.append(row);allrows.append(row)
                trace=frozen.get('best_trace',[])
                traces.append(dict(case=name,method=method,kind='FROZEN_ACTUAL_BEST_TRACE',source=str(file.relative_to(ROOT.parent)),
                    trace_has_material_vectors=False,interpretation='Recorded best-objective trajectory; no interpolation or unrecorded state reconstruction',best_trace=trace,stage_events=[{k:v for k,v in x.items() if k not in ['initial_material','initial_coarse_material']} for x in frozen.get('stages',[])]))
            write(out/f'{name}_mechanism.json',rows)
            print(json.dumps(dict(case=name,path_residuals=[round(x['relative_residual'],6) for x in rows[:6]],gradient_GT_alignment=[x['physical_descent_GT_cosine'] for x in rows[:6]],path_gram_conditions=[x['source_Gram_condition'] for x in rows[:6]],wall_seconds=time.perf_counter()-begun)),flush=True)
        costs=dict(counts=dict(model.counts),operation_wall_seconds=dict(model.timings),operation_cpu_seconds=dict(model.cpu_timings),wall_seconds=time.perf_counter()-begun,cpu_seconds=time.process_time()-cpu,threadpools=threadpool_info(),resolvent_internal_RHS='OPAQUE_LAPACK_ESTIMATOR_WORK_NOT_EXPOSED_NOT_ZERO',scope='New offline diagnostic cost only; frozen historical solver cost is separately retained in endpoint rows')
        write(out/'mechanism_all.json',allrows);write(out/'frozen_actual_traces.json',traces);write(out/'COST_LEDGER.json',costs)
        columns=['case','kind','t','method','relative_residual','material_RMSE_to_GT','J_sigma_max','J_sigma_32','J_sigma_min','J_stable_rank','J_energy_top32','source_Gram_condition','resolvent_norm1_estimate','physical_descent_GT_cosine','smallgram_descent_GT_cosine','physical_and_smallgram_descent_cosine','top32_max_angle_degrees','top32_mean_sin_angle','wall_seconds','cpu_seconds']
        with (out/'mechanism_table.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows([{k:clean(r.get(k)) for k in columns} for r in allrows])
    print(json.dumps(clean(costs)),flush=True)
if __name__=='__main__':main()
