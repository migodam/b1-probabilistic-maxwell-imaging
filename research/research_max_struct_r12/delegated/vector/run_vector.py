"""Reduced true-vector material experiments; explicit generation/solve/score phases."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import sys,json,time,argparse
import numpy as np
from scipy.optimize import minimize
HERE=Path(__file__).absolute().parent
REPO=HERE.parents[3]
sys.path.insert(0,str(REPO))
from research.research_max_struct_r8.src.vector_pr3d import VectorPR3D
from research.research_max_struct_r8.src.vector_cie3d import box_fraction

CASES=['mixed_coated_cube','connected_inclusions','mixed_coated_cube_crossgrid8']
METHODS=['CIE_sourcewise','Dual_sourcewise','CIE_global','Aplus_global','AP','DR']
CONTRACT=dict(status='FROZEN_BEFORE_GENERATION',inverse_side=6,width=.6,real_chi_bounds=[0,4],initial_chi=.25,noise_relative_complex_RMS=.01,seed=202611111212,auxiliary_iterations=120,projection_iterations=3000,projection_audit_stride=50,projection_CPU_cap_seconds=320,common_corrector_maxiter=120,prior='inherited R8 lam=2e-4 mean neighbour squared differences + l2=1e-6 mean chi squared',selection='minimum original R8 source-normalized physical Phi over paid candidates, including initial and common corrector callbacks',methods=METHODS,cases=CASES,qualification='reduced true3D dyadic CM+RR DDA; not literal 2019 paper reproduction; exploratory public development geometry')
def dump(path,x):path.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def model(side=6):
    m=VectorPR3D(ns=side,width=.6,lossy=False)
    return m
def geometry(m,name):
    p,h=m.points,m.h
    if name.startswith('mixed_coated_cube'):
        outer=box_fraction(p,h,[.015,-.015,.01],[.19,.17,.18]);inner=box_fraction(p,h,[.045,-.035,.025],[.095,.075,.09])
        return .8*outer+2.4*inner
    a=box_fraction(p,h,[-.145,-.07,-.025],[.075,.095,.115])
    b=box_fraction(p,h,[.14,.07,.04],[.08,.075,.1])
    bridge=box_fraction(p,h,[0,0,0],[.11,.035,.035])
    return np.minimum(4,2.8*a+1.5*b+1.1*bridge)
def generate():
    if (HERE/'CONTRACT.json').exists():raise FileExistsError('Frozen experiment already exists')
    dump(HERE/'CONTRACT.json',CONTRACT);pub=HERE/'data_public';truth=HERE/'sealed_truth';pub.mkdir();truth.mkdir()
    start=time.perf_counter();cpu=time.process_time();rows=[]
    for i,cid in enumerate(CASES):
        m=model(8 if cid.endswith('crossgrid8') else 6);x=geometry(m,cid);Y,J,E=m.physical(x)
        rng=np.random.default_rng(CONTRACT['seed']+i);sigma=.01*np.linalg.norm(Y)/np.sqrt(Y.size)
        N=sigma*(rng.standard_normal(Y.shape)+1j*rng.standard_normal(Y.shape))/np.sqrt(2)
        np.savez_compressed(pub/(cid+'.npz'),Y=Y+N,sigma=sigma,source_points=m.receivers,grid_shape=np.array([6]*3))
        inv=model();xtarget=geometry(inv,cid)
        np.savez_compressed(truth/(cid+'.npz'),chi=xtarget,generation_chi=x,generation_side=m.ns)
        rows.append(dict(case=cid,generation_side=m.ns,inverse_side=6,generator_counts=m.counts,setup=m.setup,sigma=float(sigma),observation_shape=list(Y.shape),truth_access='GENERATOR_ONLY'))
    dump(HERE/'GENERATION_RECEIPT.json',dict(status='COMPLETE',rows=rows,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start))
def solve(cid,method):
    if cid not in CASES or method not in METHODS:raise ValueError('not frozen')
    out=HERE/'RESULTS'/cid/method;out.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((HERE/'CONTRACT.json').read_text());Y=np.load(HERE/'data_public'/(cid+'.npz'))['Y'];cpu=time.process_time();wall=time.perf_counter();m=model()
    if method in METHODS[:4]:
        arrays,receipt=m.optimize_pr(Y,method,outer=cfg['auxiliary_iterations'],inner=12,corrector=cfg['common_corrector_maxiter'])
    else:
        sigma=float(np.load(HERE/'data_public'/(cid+'.npz'))['sigma'])
        arrays,receipt=projection(m,Y,method,cfg,out,sigma)
    np.savez_compressed(out/'prediction.npz',**arrays)
    receipt.update(status='COMPLETE',case_id=cid,method=method,whole_child_CPU_seconds_near_exit=time.process_time()-cpu,process_CPU_seconds_near_exit=time.process_time(),whole_child_wall_seconds=time.perf_counter()-wall,setup=m.setup,truth_access='NONE',selection=cfg['selection'])
    dump(out/'RECEIPT.json',receipt)
    print(json.dumps(dict(case=cid,method=method,cpu=receipt['whole_child_CPU_seconds_near_exit'],phi=receipt['physical_phi'])),flush=True)
def projection(m,Y,method,cfg,out,sigma):
    from scipy.stats import chi2
    from research.research_max_struct_r12.src.field_current import DenseFieldCurrentProjectors,PhysicalAudit,run_projection
    weights=1/np.sqrt(m.P*np.sum(abs(Y)**2,axis=0))
    delta=float(sigma*np.sqrt(chi2.ppf(.99,2*Y.size)/2))
    p=DenseFieldCurrentProjectors(m.G,m.inc,m.S,Y,delta,scale=6,lower=0,upper=4,components=3,public_metadata=dict(case_scope='vector216 reduced public DEV',noise_rule='one raw iid 99% Frobenius noise event; raw unwhitened measured data ball'))
    def audit(x):
        before=m.counts.copy();phi,df,J,E=m.physical(x,Y,False)
        return PhysicalAudit(float(phi),J,E,metrics=dict(data_residual=float(np.sqrt(df))),costs={k:v-before.get(k,0) for k,v in m.counts.items()})
    run=run_projection('dr' if method=='DR' else 'ap',p,np.full(m.N,.25),audit,steps=cfg['projection_iterations'],audit_every=cfg['projection_audit_stride'],cpu_budget_seconds=cfg['projection_CPU_cap_seconds'])
    run.save(out/'projection')
    best=run.material.copy();bestphi=run.objective;raw=best.copy();curve=[];trials=[]
    def callback(x):
        nonlocal best,bestphi
        phi,df,J,E=m.physical(x,Y,False);trials.append(x.copy());curve.append(dict(phi=float(phi),physical_LU=m.counts['physical_LU']))
        if phi<bestphi:bestphi=phi;best=x.copy()
    res=minimize(lambda x:m.physical(x,Y),best,jac=True,method='L-BFGS-B',bounds=[(0,4)]*m.N,callback=callback,options=dict(maxiter=cfg['common_corrector_maxiter'],ftol=1e-13,gtol=1e-9,maxls=30))
    callback(res.x);f,df,J,E=m.physical(best,Y,False)
    arrays=dict(initial=np.full(m.N,.25),raw=raw,prediction=best,J=J,E=E,trials=np.asarray(trials))
    receipt=dict(method=method,physical_phi=float(f),physical_data_residual=float(np.sqrt(df)),pre_corrector_phi=float(run.objective),counts=m.counts.copy(),projection_counts=run.costs,projection_physical_counts=run.physical_costs,projection_setup_cpu=run.setup_cpu_seconds,projection_cpu=run.cpu_seconds,projection_iterations=run.iterations,projection_stop=run.stop_reason,projection_delta=delta,optimizer=[dict(stage='common_long_corrector',nit=int(res.nit),nfev=int(res.nfev),success=bool(res.success),message=str(res.message))],curve=curve)
    return arrays,receipt
def score():
    paths=[HERE/'RESULTS'/c/k/'prediction.npz' for c in CASES for k in METHODS]
    paths += [HERE/'RESULTS'/c/k/'long_corrector/prediction.npz' for c in CASES for k in METHODS]
    if not all(p.exists() for p in paths):raise RuntimeError('All predictions must be sealed before GT scoring')
    dump(HERE/'PREDICTION_SEAL.json',dict(status='SEALED_BEFORE_TRUTH_READ',predictions=[str(p.relative_to(HERE)) for p in paths],selection='original physical Phi only'))
    prefix_rows=seal_prefixes()
    rows=[]
    for cid in CASES:
        truth=np.load(HERE/'sealed_truth'/(cid+'.npz'))['chi']
        for method in METHODS:
            out=HERE/'RESULTS'/cid/method;z=np.load(out/'prediction.npz');r=json.loads((out/'RECEIPT.json').read_text())
            row=dict(case_id=cid,method=method,NRMSE=float(np.linalg.norm(z['prediction'].real-truth)/np.linalg.norm(truth)),raw_NRMSE=float(np.linalg.norm(z['raw'].real-truth)/np.linalg.norm(truth)),physical_phi=r['physical_phi'],physical_data_residual=r['physical_data_residual'],CPU_seconds=r['whole_child_CPU_seconds_near_exit'],counts=r['counts'])
            rows.append(row)
    dump(HERE/'SCORES.json',dict(status='OFFLINE_SCORED_AFTER_ALL_PREDICTIONS_SEALED',rows=rows,new_physical_solves=0,truth_access='SCORER_ONLY'))
    longrows=[]
    for cid in CASES:
        truth=np.load(HERE/'sealed_truth'/(cid+'.npz'))['chi']
        for method in METHODS:
            base=HERE/'RESULTS'/cid/method;out=base/'long_corrector';z=np.load(out/'prediction.npz');r=json.loads((out/'RECEIPT.json').read_text());original=json.loads((base/'RECEIPT.json').read_text())
            longrows.append(dict(case_id=cid,method=method,NRMSE=float(np.linalg.norm(z['prediction']-truth)/np.linalg.norm(truth)),physical_phi=r['physical_phi'],physical_data_residual=r['physical_data_residual'],additional_CPU_seconds=r['CPU_seconds'],total_CPU_seconds=original['whole_child_CPU_seconds_near_exit']+r['CPU_seconds'],counts=r['counts'],optimizer=r['optimizer']))
    dump(HERE/'LONG_SCORES.json',dict(status='OFFLINE_SCORED_AFTER_ALL36_PREDICTIONS_SEALED',rows=longrows,new_physical_solves=0,truth_access='SCORER_ONLY'))
    with np.load(HERE/'PREFIX_PREDICTIONS.npz') as z:
        for r in prefix_rows:
            truth=np.load(HERE/'sealed_truth'/(r['case_id']+'.npz'))['chi']
            r['NRMSE']=float(np.linalg.norm(z[r['array_key']]-truth)/np.linalg.norm(truth))
    dump(HERE/'PREFIX_SCORES.json',dict(status='OFFLINE_SCORED_AFTER_PREFIX_SELECTION_SEALED',rows=prefix_rows,new_physical_solves=0))
def seal_prefixes():
    prefixes=json.loads((HERE/'PREFIX_CONTRACT.json').read_text())['CPU_prefix_seconds'];arrays={};rows=[]
    for cid in CASES:
        for method in METHODS:
            out=HERE/'RESULTS'/cid/method;r=json.loads((out/'RECEIPT.json').read_text());z=np.load(out/'prediction.npz');candidates=[]
            if method in ['AP','DR']:
                h=json.loads((out/'projection/history.json').read_text());setup=r['setup']['cpu_s']+r['projection_setup_cpu']
                for i,record in enumerate(h['history']):candidates.append((setup+record['cpu_seconds'],record['incumbent_objective'],'audit',i))
            else:
                setup=r['setup']['cpu_s']
                for i,record in enumerate(r['curve']):candidates.append((setup+record['cpu_s'],record['phi'],'trial',i))
            candidates.append((r['whole_child_CPU_seconds_near_exit'],r['physical_phi'],'final',0))
            for limit in prefixes:
                legal=[x for x in candidates if x[0]<=limit]
                if not legal:continue
                picked=min(legal,key=lambda x:x[1]);source,index=picked[2:]
                if source=='final':m=z['prediction'].real
                elif source=='trial':m=z['trials'][index].real
                else:
                    with np.load(out/'projection/states.npz') as states:m=states[f'audit_{index:05d}_incumbent_material'].copy()
                key=f'p{len(rows):04d}';arrays[key]=m;rows.append(dict(case_id=cid,method=method,CPU_limit_seconds=limit,candidate_CPU_seconds=picked[0],physical_phi=picked[1],source=source,index=index,array_key=key))
    np.savez_compressed(HERE/'PREFIX_PREDICTIONS.npz',**arrays)
    dump(HERE/'PREFIX_SELECTION_SEAL.json',dict(status='SEALED_BEFORE_TRUTH_READ',rows=rows))
    return rows
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('phase',choices=['generate','solve','score']);p.add_argument('--case');p.add_argument('--method');a=p.parse_args()
    if a.phase=='generate':generate()
    elif a.phase=='solve':solve(a.case,a.method)
    else:score()
