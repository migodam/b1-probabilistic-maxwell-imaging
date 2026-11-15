"""One exposed public case, three fixed updates: original GPG wiring only.

Also checks raw G/S rank headers on the 24 original public inputs, without
physical solves. Anchor transport uses the one already-paid initial LU of the
finite fixture; it never creates an additional anchor forward or new material.
"""
from __future__ import annotations
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import argparse,csv,importlib.util,json,sys,time
sys.dont_write_bytecode=True
from pathlib import Path


def write(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,required=True);ap.add_argument('--legacy-repo',type=Path,required=True);ap.add_argument('--out',type=Path);a=ap.parse_args();s=a.root.resolve();W=a.legacy_repo.resolve();R=s.parent;out=(a.out or s/'audit/wiring_validation_v2').resolve();out.mkdir(exist_ok=False)
    wall=time.perf_counter();start=time.process_time();write(out/'PROTOCOL_BEFORE_ACTIONS.json',dict(truth_access='NONE',header_inputs='Only original exposed24 public G/S/sigma; no forward for header check',fixture_case='r14_new_0000',fixed_native_updates_per_side=3,comparators=['Original R13 run_full_current/R11.prepare_full','Private internal_galerkin factory over same original code'],anchor_fixture='One prepare-only anchor Galerkin using new-side cached initialLU, before finite native trajectory; no additionalLU/forward',original_data_or_algorithm_modified=False,all_validation_work_charged=True))
    paths=[s/'src',s/'runtime',R/'runtime',W]
    sys.path[:0]=[str(p) for p in paths if p.is_dir()]
    import numpy as np
    from observable19_factory import svd_basis_headers,bind_runner
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    from research.research_max_struct_r8.src import benchmark_solver
    from research.research_max_struct_r10.src.twod_adapter import noise_ball_rank
    from observable19_factory import original
    run_full_current=original.run_full_current
    headers=[];headerCPU=time.process_time();spectra={}
    for i in range(24):
        cid=f'r14_new_{i:04d}';path=R/'data_public/independent24'/(cid+'.npz')
        with np.load(path,allow_pickle=False) as z:G=z['G'].copy();S=z['S'].copy();sigma=float(z['sigma'])
        _,info=svd_basis_headers(G,S)
        for key in ('singular_G','singular_S'):spectra[cid+'__'+key]=info.pop(key)
        headers.append(dict(case_id=cid,source=str(path),G_shape=list(G.shape),S_shape=list(S.shape),sigma=sigma,**info,
            internal_galerkin_header_qualified=info['original_G_threshold_0p1_rank']==19,
            internal_right_header_qualified=info['numerical_G_rank']>=19,
            observation_header_qualified=info['numerical_S_rank']>=19,
            anchor_header='Not inferred from raw S numerical rank; required initialLU/rank checked by actual anchor factory',new_physical_LU=0,new_complex_RHS=0))
    write(out/'PUBLIC_RANK_HEADERS.json',headers)
    with (out/'PUBLIC_RANK_HEADERS.csv').open('x',newline='') as f:w=csv.DictWriter(f,fieldnames=list(headers[0]));w.writeheader();w.writerows(headers)
    np.savez_compressed(out/'PUBLIC_RANK_SPECTRA.npz',**spectra)
    header_cost=dict(G_SVD=24,unscaled_S_SVD=24,CPU_seconds=time.process_time()-headerCPU,physical_LU=0,physical_complex_RHS=0)
    fixture_path=R/'data_public/independent24/r14_new_0000.npz'
    with np.load(fixture_path,allow_pickle=False) as z:public={k:z[k].copy() for k in z.files}
    m0=public.pop('m_init');Y=public['Y'];sigma=float(public['sigma']);decision=noise_ball_rank(public['S'],Y,sigma,out/'RECEIVER_DPCS.npz')
    trajectories={};models={};contexts={};results={}
    for side in ('original','fixed19'):
        directory=out/side;directory.mkdir();model=MaxwellModel(public);models[side]=model
        model.state(m0)  # Exactly the C2 paid initial state; subsequent warm/audits reuse cache.
        common=dict(host=.1,material_box=[.02,20.],dpcs_rank=decision['rank'])
        cfg=benchmark_solver.arm_config('a_plus',common,Y,model.E,[dict(mf=12,maxiter=1500)])
        cfg['stages']=[dict(cfg['stages'][0],maxiter=3)]
        cfg.update(fft_geometry=dict(coords=public['coords'],radius=float(public['cell_radius']),k=float(public['k'])),regularizer=lambda m,model=model:model.reg.evaluate(m)[:2],aplus_regularizer=model.reg,physical_sigma=float(np.linalg.norm(Y)),receipt_path=str(directory/'AUXILIARY_EXECUTION.json'),initialize_physical_current=True,coefficient_tol=0.,miniter=50)
        captured=[]
        def capture(row,captured=captured):captured.append({k:v.copy() if isinstance(v,np.ndarray) else v for k,v in row.items()})
        cfg['iteration_callback']=capture
        if side=='fixed19':
            # Test the declared cached-LU interface without an anchor inversion run.
            anchor=out/'anchor_prepare_only';anchor.mkdir();anchor_run,anchor_context=bind_runner('g19_anchor_observation_galerkin',anchor_m0=m0)
            initial_object=model._current;factor=model._physics.factor;before=dict(model.cost)
            anchor_run.prepare_for_validation(model,dict(cfg,receipt_path=str(anchor/'AUXILIARY_EXECUTION.json')))
            write(out/'ANCHOR_CACHE_INTERFACE.json',dict(status='CACHED_INITIAL_LU_USED',same_initial_state_object=model._current is initial_object,same_initial_factor_object=model._physics.factor is factor,material_unchanged=bool(np.array_equal(model._current['m'],m0)),physical_cost_delta={k:model.cost.get(k,0)-before.get(k,0) for k in set(model.cost)|set(before)},setup=anchor_context['setup'],qualification='Prepare-only current-defect-to-port operator; no material Jacobian, no anchor native/inverse run'))
            run,context=bind_runner('g19_internal_galerkin',anchor_m0=m0);contexts[side]=context
        else:run=run_full_current
        results[side]=run(model,Y,m0,cfg,mode='restricted',audit_every=0,g19=True);trajectories[side]=captured
        np.savez_compressed(directory/'THREE_FIXED_UPDATES.npz',m=np.asarray([r['m'] for r in captured]),J=np.asarray([r['J'] for r in captured]),coefficients=np.asarray([r['coefficients'] for r in captured]))
    comparisons={};all_equal=True
    for key in ('m','J','coefficients'):
        x=np.asarray([r[key] for r in trajectories['original']]);y=np.asarray([r[key] for r in trajectories['fixed19']]);same=bool(np.array_equal(x,y));all_equal&=same;comparisons[key]=dict(bitwise_equal=same,max_abs_difference=float(np.max(np.abs(x-y))))
    with np.load(out/'original/DENOMINATOR_BASIS.npz',allow_pickle=False) as z:U0=z['U'].copy();B0=z['B'].copy()
    with np.load(out/'fixed19/DENOMINATOR_BASIS.npz',allow_pickle=False) as z:U1=z['U'].copy();B1=z['B'].copy()
    for key,x,y in [('U',U0,U1),('B',B0,B1),('K_applied_to_original_E0',B0@(U0.conj().T@public['E0']),B1@(U1.conj().T@public['E0']))]:comparisons[key]=dict(bitwise_equal=bool(np.array_equal(x,y)),max_abs_difference=float(np.max(np.abs(x-y))))
    dependencies={k:importlib.util.find_spec(k).origin for k in ('field_current','material_events')}
    receipt=dict(status='FINITE_WIRING_MATCH' if all_equal else 'WIRING_MISMATCH_REQUIRES_REVIEW',truth_access='NONE',fixed_case='r14_new_0000',fixed_updates={k:len(v) for k,v in trajectories.items()},comparisons=comparisons,public_header_cost=header_cost,header_cases=24,all_original_G_strong_rank19=all(r['internal_galerkin_header_qualified'] for r in headers),all_raw_observation_rank_at_least19=all(r['observation_header_qualified'] for r in headers),public_G_supported_ranks=sorted({r['numerical_G_rank'] for r in headers}),public_S_supported_ranks=sorted({r['numerical_S_rank'] for r in headers}),finite_fixture_physical_cost={k:dict(m.cost) for k,m in models.items()},validation_lowrank_probe_actions=2,validation_lowrank_probe_complex_RHS=2*public['E0'].shape[1],fixture_lowrank_small_solve_cost={k:{name:value for name,value in result['cost'].items() if name.startswith(('denominator_','lowrank_','matrix_material_'))} for k,result in results.items()},resolved_dependencies=dependencies,CPU_seconds=time.process_time()-start,whole_process_CPU_seconds_at_receipt=time.process_time(),wall_seconds=time.perf_counter()-wall,remote_launched=False,full_recovery_experiment=False,rank_lowering=False)
    write(out/'RECEIPT.json',receipt);print(json.dumps(receipt,indent=2))
    return 0 if all_equal else 1
if __name__=='__main__':raise SystemExit(main())
