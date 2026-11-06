"""E1: offline geometry on frozen anchors, no optimization or GT-selected basis."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time,csv,faulthandler
faulthandler.enable()
import numpy as np
from scipy.linalg import svd,lu_solve,eigh,eigvalsh
from threadpoolctl import threadpool_limits
HERE=Path(__file__).resolve().parents[1];REPO=HERE.parents[1]
sys.path.insert(0,str(REPO))
from research.research_z_alg_r4_diversity.src.physics import Model
from research.research_gn_search_space_r2.src.space import prepare_geometry,current_complement_from_cache,euclidean_orth,real_j_svd
OUT=HERE/'results/e1';OUT.mkdir(parents=True,exist_ok=True)
SUBSETS=([0],[0,4,8],[0,2,4,6,8,10],list(range(12)))
RANKS=(8,16,32,64,128,256,576)

def writecsv(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)

def jspectra(jcomplex,Y,label,case,model):
    spectra={};rows=[]
    for ix in SUBSETS:
        jc=jcomplex.reshape(12,48,576)[ix].reshape(-1,576)
        jr=np.stack((jc.real,jc.imag),axis=1).reshape(-1,576)
        with model.operation('E1_SUBSET_J_SVD',material_svd=1):s,v=real_j_svd(jr)
        scale=float(np.linalg.norm(Y[:,ix]));ss=s/scale
        rank=int(np.count_nonzero(s>1e-12*s[0]))
        spectra[len(ix)]=(s,v)
        for rms in (.01,.03,.1):
            for eps in (1e-5,1e-4,1e-3):
                a=rms*np.sqrt(576)*ss/eps
                rows.append(dict(case=case,linearization=label,sources=len(ix),source_indices=json.dumps(ix),RMS=rms,epsilon=eps,
                                 observed_Y_norm=scale,d_weak=int(np.count_nonzero(a<=1)),d_soft=float(np.sum(1/(1+a*a))),
                                 numerical_rank=rank,numerical_null_dimension=576-rank,sigma_max=float(s[0]),sigma_min=float(s[-1])))
    return spectra,rows

def drift(a,b):
    cs=svd(a.T@b,compute_uv=False,check_finite=False,lapack_driver='gesvd');cs=np.clip(cs,0,1)
    return dict(projector_distance_fro=float(np.sqrt(max(0,2*a.shape[1]-2*np.sum(cs*cs)))),
                min_principal_angle_deg=float(np.degrees(np.arccos(cs[0]))),max_principal_angle_deg=float(np.degrees(np.arccos(cs[-1]))),
                mean_cosine_squared=float(np.mean(cs*cs)))

def main():
    start=time.perf_counter();old=REPO/'research/research_z_alg_r4_diversity'
    source=REPO/'research/research_gn_green_error_20261106/results/error_arrays.npz'
    meta={r['task_id']:r for r in json.loads((old/'RESULTS_SUMMARY.json').read_text())['rows']}
    arr=np.load(source,allow_pickle=False);allids=[str(s) for s in arr['task_id']]
    ids=[i for i,t in enumerate(allids) if meta[t]['amplitude']==.4 or (meta[t]['method']=='native_trf' and meta[t]['noise']=='clean' and meta[t]['shape'] in ('digit7','ring') and meta[t]['amplitude'] in (2.5,12.))]
    assert len(ids)==52,len(ids)
    public=dict(np.load(old/'public_data/operators.npz',allow_pickle=False));model=Model(public)
    U=current_complement_from_cache();print(json.dumps(dict(phase='fixed_cached',U_rank=U.shape[1],D_dim=576-U.shape[1])),flush=True)
    def previous(name):
        path=OUT/name
        return list(csv.DictReader(path.open())) if path.exists() else []
    dimensions=previous('dimensions.csv');coverage=previous('coverage.csv');drifts=previous('drift.csv');summary=previous('cases.csv');costs=previous('costs.csv')
    completed={r['task_id'] for r in summary}
    # Public fixed reference: entirely independent of all GT labels.
    state=model.forward(np.full(576,.1),label='public_mref')
    print('stage reference forward done',flush=True)
    if (OUT/'public_reference_geometry.npz').exists():
        z=np.load(OUT/'public_reference_geometry.npz')
        ref=dict(J_right_vectors=z['J_V'],J_singular_values=z['J_s'],leakage_vectors=z['leakage_vectors'],eigenvalues=z['leakage_values'],approx_basis=z['approx_basis'])
        model.jacobian(state)
    else:
        ref=prepare_geometry(model,state,current_complement=U)
    refJ=ref['J_right_vectors'][:,::-1];refG=ref['leakage_vectors']
    np.savez_compressed(OUT/'public_reference_geometry.npz',m=state['m'],J_s=ref['J_singular_values'],J_V=ref['J_right_vectors'],leakage_values=ref['eigenvalues'],leakage_vectors=refG,approx_basis=ref['approx_basis'],current_complement=U)
    # Exact reference and Born spectra differ only by linearization. Norm is Ycase.
    refjc=state['Jcomplex'].copy();bornjc=np.einsum('ri,is->sri',model.S,model.E0).reshape(576,576)
    gtcache={};refcache={};borncache={}
    for num,i in enumerate(ids):
        tic=time.perf_counter();tid=allids[i];r=meta[tid]
        if tid in completed:continue
        print(json.dumps(dict(phase='start_case',id=tid)),flush=True)
        before=dict(model.counts);Y=np.load(old/'public_data'/f"{r['observation_id']}.npz")['Y']
        gt=arr['truth_m'][i];m=arr['estimate_m'][i];e=gt-m
        obs=r['observation_id']
        if obs not in refcache:
            refcache[obs]=jspectra(refjc,Y,'public_mref',obs,model)
            borncache[obs]=jspectra(bornjc,Y,'Born',obs,model)
            dimensions.extend(refcache[obs][1]);dimensions.extend(borncache[obs][1])
        # Unique GT uses material values as offline labels; this state never enters solver.
        key=gt.tobytes()
        if key not in gtcache:
            print('stage GT forward '+tid,flush=True)
            gs=model.forward(gt,label='offline_GT');model.jacobian(gs)
            gtcache[key]=dict(jcomplex=gs['Jcomplex'].copy(),y=gs['y'].copy())
        print('stage GT spectra '+tid,flush=True)
        _,rows=jspectra(gtcache[key]['jcomplex'],Y,'offline_GT',tid,model);dimensions.extend(rows)
        print('stage anchor geometry '+tid,flush=True)
        state=model.forward(m,label='frozen_anchor');geo=prepare_geometry(model,state,current_complement=U)
        print('stage anchor spectra '+tid,flush=True)
        specs,rows=jspectra(state['Jcomplex'],Y,'anchor',tid,model);dimensions.extend(rows)
        # Offline projection is applied AFTER every basis is constructed from public state.
        bases={'J_weakest':geo['J_right_vectors'][:,::-1],'current_leakage_weakest':geo['leakage_vectors']}
        projections={};total=float(e@e)
        for name,raw in bases.items():
            for rank in RANKS:
                b=raw[:,:rank] if name=='J_weakest' else euclidean_orth(raw[:,:rank])
                p=b@(b.T@e);projmat=m+p
                lin=geo['W']@(p[:,None]*geo['Q']);leak=U.conj().T@lin
                entry=dict(task_id=tid,basis=name,rank=rank,error_coverage=float(p@p/total),
                           projected_RMS=float(np.linalg.norm(p)/24),remaining_RMS=float(np.linalg.norm(e-p)/24),
                           box_violations=int(np.count_nonzero((projmat<.02)|(projmat>20))),
                           box_min=float(projmat.min()),box_max=float(projmat.max()),
                           tangent_current_relative_leakage=float(np.linalg.norm(leak)**2/max(np.linalg.norm(lin)**2,1e-300)),
                           tangent_relative_data=float(np.linalg.norm(geo['J']@p)/np.linalg.norm(Y)))
                coverage.append(entry)
                if rank in (32,64):projections[f'{name}_{rank}']=p
        for threshold in (.01,.1):
            b=euclidean_orth(geo['leakage_vectors'][:,geo['eigenvalues']<=threshold]);p=b@(b.T@e);z=m+p
            lin=geo['W']@(p[:,None]*geo['Q'])
            coverage.append(dict(task_id=tid,basis=f'current_leakage_threshold_{threshold}',rank=b.shape[1],error_coverage=float(p@p/total),
                                 projected_RMS=float(np.linalg.norm(p)/24),remaining_RMS=float(np.linalg.norm(e-p)/24),box_violations=int(np.count_nonzero((z<.02)|(z>20))),box_min=float(z.min()),box_max=float(z.max()),tangent_current_relative_leakage=float(np.linalg.norm(U.conj().T@lin)**2/max(np.linalg.norm(lin)**2,1e-300)),tangent_relative_data=float(np.linalg.norm(geo['J']@p)/np.linalg.norm(Y))))
            projections[f'threshold_{threshold}']=p
        for name,raw,reference in [('J_weakest',bases['J_weakest'],refJ),('current_leakage_weakest',bases['current_leakage_weakest'],refG)]:
            for rank in (32,64):
                a=euclidean_orth(raw[:,:rank]);b=euclidean_orth(reference[:,:rank])
                drifts.append(dict(task_id=tid,basis=name,rank=rank,**drift(a,b)))
        details=dict(task_id=tid,observation_id=obs,method=r['method'],shape=r['shape'],amplitude=r['amplitude'],noise=r['noise'],
                     material_RMSE=float(np.linalg.norm(e)/24),**geo['diagnostics'],case_wall_seconds=time.perf_counter()-tic)
        summary.append(details)
        casecost={k:int(v-before.get(k,0)) for k,v in model.counts.items() if v-before.get(k,0)}
        costs.append(dict(task_id=tid,wall_seconds=time.perf_counter()-tic,**casecost))
        # Archive exact basis and offline projected directions for bounded nonlinear check.
        np.savez_compressed(OUT/f'{tid}.npz',anchor=m,error_OFFLINE=e,Y=Y,J_s=geo['J_singular_values'],J_V=geo['J_right_vectors'],
                            leakage_values=geo['eigenvalues'],leakage_vectors=geo['leakage_vectors'],approx_basis=geo['approx_basis'],
                            numerical_kernel=geo['numerical_kernel'],leakage_singular_values=geo['leakage_singular_values'],**projections)
        for name,rows0 in [('dimensions.csv',dimensions),('coverage.csv',coverage),('drift.csv',drifts),('cases.csv',summary),('costs.csv',costs)]:writecsv(OUT/name,rows0)
        (OUT/'progress.json').write_text(json.dumps(dict(completed=num+1,total=len(ids),last=tid,elapsed_seconds=time.perf_counter()-start,counts=dict(model.counts)),indent=2))
        print(json.dumps(dict(case=num+1,total=len(ids),id=tid,approx_rank=geo['approx_basis'].shape[1],kernel=geo['diagnostics']['numerical_material_kernel_dimension'],seconds=time.perf_counter()-tic)),flush=True)
        if os.environ.get('GN_E1_ONECASE')=='1':
            (OUT/f'process_cost_{num+1:02d}.json').write_text(json.dumps(dict(case=tid,counts=dict(model.counts),wall_seconds=time.perf_counter()-start),indent=2))
            return
    report=dict(status='COMPLETE',cases=len(ids),unique_GT=len(gtcache),reference_material=.1,source_subsets=SUBSETS,rank_rtol=1e-12,
                total_wall_seconds=time.perf_counter()-start,counts=dict(model.counts),operation_wall_seconds=dict(model.timings),operation_cpu_seconds=dict(model.cpu_timings),
                fixed_current_complement_rank=U.shape[1],fixed_current_dark_dimension=576-U.shape[1],
                gt_exposed=True,exposure_scope='E1 offline geometry and scoring only; src/space.py GT-free',
                finite_nonlinear_projected_GT_evaluations='NOT_RUN by E1; directions saved for parent')
    (OUT/'RUN_SUMMARY.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)

def finalize_existing():
    P=OUT
    from collections import Counter
    def read(name):return list(csv.DictReader(open(P/name)))
    def write(name,rows):
     keys=list(dict.fromkeys(k for r in rows for k in r))
     with open(P/name,'w',newline='') as f:
      w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    rows=read('dimensions.csv');dedup={tuple(r[k] for k in ('case','linearization','sources','RMS','epsilon')):r for r in rows};write('dimensions.csv',list(dedup.values()))
    cases=read('cases.csv');integrity=[];r=[]
    for row in cases:
     a=np.load(P/(row['task_id']+'.npz'))
     finite=all(np.all(np.isfinite(a[k])) for k in a.files)
     v=a['J_V'];b=a['approx_basis'];s=a['leakage_singular_values']
     r.append(s[-1]/s[0]);row.update(leakage_sigma_min_relative=s[-1]/s[0],leakage_sigma_max=s[0],leakage_sigma_min=s[-1])
     integrity.append(dict(task_id=row['task_id'],all_saved_arrays_finite=finite,J_orth_error=float(np.linalg.norm(v.T@v-np.eye(576))),approx_orth_error=float(np.linalg.norm(b.T@b-np.eye(b.shape[1]))),eigen_residual=float(row['generalized_eigen_residual']),kernel_dimension=a['numerical_kernel'].shape[1],strict_sigma_min_relative=float(s[-1]/s[0])))
    write('cases.csv',cases);write('integrity.csv',integrity)
    assert len(cases)==52 and all(x['all_saved_arrays_finite'] and x['J_orth_error']<1e-10 and x['approx_orth_error']<1e-10 and x['eigen_residual']<1e-10 for x in integrity)
    counts=Counter()
    for row in read('costs.csv'):
     for k,v in row.items():
      if k not in ('task_id','wall_seconds') and v:counts[k]+=int(v)
    counts.update(dict(FULL_FORWARD=1,LU=1,actual_complex_RHS=636,source_RHS=12,current_tangent_RHS=576,receiver_transpose_RHS=48,EXACT_CURRENT_TANGENT_INVERSE=1,REAL_MATERIAL_CURRENT_GEOMETRY=1,current_geometry_builds=1,FULL_ANALYTIC_J=1,FullJ=1,material_J_builds=1,MATERIAL_J_SVD=1,material_svd=1))
    raw=np.load(P.parents[2]/'research_gn_green_error_20261106/results/error_arrays.npz')
    ids={x['task_id'] for x in cases};unique=len({raw['truth_m'][i].tobytes() for i,t in enumerate(raw['task_id']) if str(t) in ids})
    report=json.load(open(P/'RUN_SUMMARY.json'));report.update(unique_GT=unique,status='COMPLETE_WITH_NATIVE_RETRY_COST_GAPS',counts=counts,costs_status='LOWER_BOUND; successful-case counters plus first public reference. Retried prep, failed uncheckpointed native work and repeated reference prep add unquantified costs.',total_wall_seconds=None,successful_case_wall_seconds=sum(float(x['wall_seconds']) for x in read('costs.csv')),max_J_orthogonality_error=max(x['J_orth_error'] for x in integrity),max_approx_orthogonality_error=max(x['approx_orth_error'] for x in integrity),max_generalized_eigen_residual=max(x['eigen_residual'] for x in integrity),strict_leakage_sigma_min_relative_range=[min(r),max(r)],dimension_rows_before_dedup=len(rows),dimension_rows_after_dedup=len(dedup),runtime=dict(python='3.13.13',numpy=np.__version__,scipy=__import__('scipy').__version__,architecture=__import__('platform').machine(),blas='macOS Accelerate; threadpoolctl reports []',threads='VECLIB_MAXIMUM_THREADS=1,OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1'),native_failures=dict(initial='exit139 after14 saved cases',resume01='exit133 during redundant reference prep',resume02='segfault after30 saved cases, faulthandler pure Python writecsv',resume03='debug allocator corruption on second case during orth',resume04='debug allocator corruption during orth',resume05='explicitly interrupted after47 saved cases to honor process isolation',final='five remaining cases completed in fresh processes; no further fault'),integrity_status='PASS all52 saved arrays finite and stored/recomputed orthogonality/eigen residual below1e-10')
    report.update(operation_wall_seconds=None,operation_cpu_seconds=None,timing_status='Successful-case wall times retained; aggregate operation and total retry timings unavailable after native failures.')
    (P/'RUN_SUMMARY.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))

def supplemental_subsets():
    """One paid anchor per fresh process; original12-source results reused."""
    old=REPO/'research/research_z_alg_r4_diversity'
    completed={r['task_id'] for r in csv.DictReader((OUT/'current_subset_dimensions.csv').open())} if (OUT/'current_subset_dimensions.csv').exists() else set()
    cases=list(csv.DictReader((OUT/'cases.csv').open()))
    row=next((r for r in cases if r['task_id'] not in completed),None)
    if row is None:return
    tid=row['task_id'];tic=time.perf_counter();z=np.load(OUT/f'{tid}.npz')
    model=Model(dict(np.load(old/'public_data/operators.npz')));state=model.forward(z['anchor'],label='supplemental_subset_anchor')
    U=current_complement_from_cache()
    with model.operation('SUPPLEMENTAL_CURRENT_TANGENT_INVERSE',rhs=576,current_tangent_RHS=576):w=lu_solve(model.factor,model.identity,check_finite=False)
    uw=U.conj().T@w;wg=w.conj().T@w;new=[]
    for ix in SUBSETS[:-1]:
        Q=state['Q'][:,ix]
        with model.operation('SUPPLEMENTAL_SUBSET_GEOMETRY',current_subset_builds=1):
            a=np.concatenate([uw*Q[:,j][None,:] for j in range(len(ix))],axis=0)
            ar=np.concatenate((a.real,a.imag),axis=0)
            sa=svd(ar,compute_uv=False,check_finite=False,lapack_driver='gesvd')
            ht=np.real(wg*(Q@Q.conj().T).T);ht=(ht+ht.T)*.5
            hc=ar.T@ar
            mu=eigh(hc,ht,eigvals_only=True,check_finite=False,driver='gv')
            hte=eigvalsh(ht,check_finite=False,driver='evr')
        rank=int(np.count_nonzero(sa>sa[0]*1e-12))
        new.append(dict(task_id=tid,sources=len(ix),source_indices=json.dumps(ix),real_leakage_rows=ar.shape[0],numerical_rank=rank,numerical_kernel_dimension=576-rank,
                        rank_rtol=1e-12,smallest_computed_relative_singular=float(sa[-1]/sa[0]),approx_dimension_0p01=int(np.count_nonzero(mu<=.01)),approx_dimension_0p1=int(np.count_nonzero(mu<=.1)),min_mu=float(mu[0]),max_mu=float(mu[-1]),ht_min_eigenvalue=float(hte[0]),ht_condition=float(hte[-1]/hte[0]),geometry_source='NEW_PAID_W'))
    new.append(dict(task_id=tid,sources=12,source_indices=json.dumps(SUBSETS[-1]),real_leakage_rows=1608,numerical_rank=row['numerical_leakage_rank'],numerical_kernel_dimension=row['numerical_material_kernel_dimension'],rank_rtol=1e-12,smallest_computed_relative_singular=row['leakage_sigma_min_relative'],approx_dimension_0p01=row['material_approx_rank'],approx_dimension_0p1=row['material_sensitivity_rank'],min_mu=row['min_relative_leakage'],max_mu=row['max_relative_leakage'],ht_min_eigenvalue=row['ht_min_eigenvalue'],ht_condition=row['ht_condition'],geometry_source='EXISTING_12_SOURCE_CACHE'))
    path=OUT/'current_subset_dimensions.csv';rows=list(csv.DictReader(path.open())) if path.exists() else []
    writecsv(path,rows+new)
    (OUT/f'subset_cost_{tid}.json').write_text(json.dumps(dict(task_id=tid,wall_seconds=time.perf_counter()-tic,counts=dict(model.counts),timings=dict(model.timings),status='COMPLETE'),indent=2))
    print(json.dumps(dict(task_id=tid,subset_rows=new,costs=dict(model.counts))),flush=True)

def fixed_reference_coverage():
    z=np.load(OUT/'public_reference_geometry.npz');raw=z['leakage_vectors'];mu=z['leakage_values']
    rows=list(csv.DictReader((OUT/'coverage.csv').open()));rows=[r for r in rows if not r['basis'].startswith('fixed_reference_Green')]
    bases=[(f'fixed_reference_Green_weakest',rank,euclidean_orth(raw[:,:rank])) for rank in RANKS]
    bases += [(f'fixed_reference_Green_threshold_{t}',int(np.count_nonzero(mu<=t)),euclidean_orth(raw[:,mu<=t])) for t in (.01,.1)]
    for case in csv.DictReader((OUT/'cases.csv').open()):
        a=np.load(OUT/(case['task_id']+'.npz'));e=a['error_OFFLINE'];m=a['anchor'];total=float(e@e)
        for name,rank,b in bases:
            p=b@(b.T@e);z=m+p
            rows.append(dict(task_id=case['task_id'],basis=name,rank=rank,error_coverage=float(p@p/total),projected_RMS=float(np.linalg.norm(p)/24),remaining_RMS=float(np.linalg.norm(e-p)/24),box_violations=int(np.count_nonzero((z<.02)|(z>20))),box_min=float(z.min()),box_max=float(z.max()),tangent_current_relative_leakage='',tangent_relative_data=''))
    writecsv(OUT/'coverage.csv',rows)

if __name__=='__main__':
    with threadpool_limits(limits=1):
        if '--finalize' in sys.argv:finalize_existing()
        elif '--subsets' in sys.argv:supplemental_subsets()
        elif '--fixed-reference-coverage' in sys.argv:fixed_reference_coverage()
        else:main()
