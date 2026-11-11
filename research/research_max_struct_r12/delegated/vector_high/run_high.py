"""Public-only high-chi true3D solves; unchanged R8 physical/PR primitives."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
from datetime import datetime,timezone
import sys,json,time,argparse
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2

HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3]
sys.path.insert(0,str(REPO))
from model_adapter import HighContrastVectorPR3D
from research.research_max_struct_r12.src.field_current import DenseFieldCurrentProjectors,PhysicalAudit,run_projection


def dump(path,value):
    def convert(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value,indent=2,default=convert,allow_nan=False)+'\n')


class CorrectorBudgetEnd(Exception):pass


def solve(cid,method):
    cfg=json.loads((HERE/'FROZEN_PROTOCOL.json').read_text())
    assert cid in cfg['cases'] and method in cfg['methods']
    out=HERE/'RESULTS'/cid/method;out.mkdir(parents=True,exist_ok=False)
    cpu,wall=time.process_time(),time.perf_counter()
    with np.load(HERE/'data_public'/f'{cid}.npz',allow_pickle=False) as z:
        assert not {'chi','truth','generation_chi','Y_clean','noise'}&set(z.files)
        Y=z['Y'].copy();sigma=float(z['sigma'])
        assert np.array_equal(z['material_bounds'],[0,8])
    model=HighContrastVectorPR3D(ns=6);model.prstate=None
    initial=np.full(model.N,.25);history=[];paid=[];best=None;phase='initial';auxiliary=None;projection_cost=None
    delta=float(sigma*np.sqrt(chi2.ppf(.99,2*Y.size)/2))
    def record(x,gradient=None):
        nonlocal best
        phi,df,J,E=model.physical(x,Y,False)
        if best is None or phi<best['phi']:
            best=dict(x=np.asarray(x).copy(),phi=float(phi),df=float(df),J=J.copy(),E=E.copy(),
                       gradient=None if gradient is None else gradient.copy())
        history.append(dict(index=len(history),phase=phase,cpu_seconds=time.process_time()-cpu,
            wall_seconds=time.perf_counter()-wall,phi=float(phi),incumbent_phi=best['phi'],
            physical_data_residual=float(np.sqrt(df)),counts=dict(model.counts)))
        paid.append(np.asarray(x).copy())
        return float(phi),float(df),J,E
    record(initial)
    if method in ('CIE_sourcewise','Dual_sourcewise'):
        phase='auxiliary';startcpu,startwall=time.process_time(),time.perf_counter();before=dict(model.counts)
        subspace=model.prepare(Y);logs=[];x=initial.copy()
        for it in range(cfg['baseline']['cycles']):
            mf=1 if it<60 else 2;beta=6 if it<60 else 1
            bc=3 if method=='Dual_sourcewise' and it>=60 else beta
            J=model.current_pr(x,Y,'CIE_TSOM',bc,mf,False)
            res=minimize(lambda z:model.material(z,J,'CIE_TSOM',beta),x,jac=True,
                method='L-BFGS-B',bounds=[model.material_bounds]*model.N,
                options=dict(maxiter=cfg['baseline']['material_inner'],ftol=1e-12,gtol=1e-8,maxls=25))
            x=res.x;record(x)
            logs.append(dict(cycle=it,mf=mf,beta_material=beta,beta_current=bc,
                nit=int(res.nit),nfev=int(res.nfev),message=str(res.message),success=bool(res.success)))
        auxiliary=dict(subspace=subspace,logs=logs,cpu_seconds=time.process_time()-startcpu,
            wall_seconds=time.perf_counter()-startwall,counts={k:v-before.get(k,0) for k,v in model.counts.items()})
        dump(out/'AUXILIARY.json',auxiliary)
    else:
        phase='projection';pcpu,pwall=time.process_time(),time.perf_counter()
        p=DenseFieldCurrentProjectors(model.G,model.inc,model.S,Y,delta,scale=6.,lower=0.,upper=8.,
            components=3,grouped_ports=True,
            public_metadata=dict(scope='true3D216chi648vector highcontrast DEV',noise='raw iid CN99 event'))
        def audit(x):
            before=dict(model.counts);phi,df,J,E=record(x)
            return PhysicalAudit(phi,J,E,metrics=dict(data_residual=float(np.sqrt(df)),
                prior=float(model.prior(x)[0])),costs={k:v-before.get(k,0) for k,v in model.counts.items()})
        run=run_projection(method.lower(),p,initial,audit,steps=cfg['projection']['steps'],
            audit_every=cfg['projection']['audit_every'],omega=1.,
            cpu_budget_seconds=max(0.,cfg['projection']['CPU_cap_seconds']-(time.process_time()-pcpu)))
        projection_cost=dict(setup_cpu_seconds=run.setup_cpu_seconds,setup_wall_seconds=run.setup_wall_seconds,
            iteration_cpu_seconds=run.cpu_seconds,iteration_wall_seconds=run.wall_seconds,
            counts_including_setup=dict(p.costs),physical_audit_counts=run.physical_costs,
            iterations=run.iterations,stop=run.stop_reason,cpu_before_archive=time.process_time()-pcpu,
            wall_before_archive=time.perf_counter()-pwall,delta=delta)
        acpu,awall=time.process_time(),time.perf_counter();run.save(out/'projection')
        projection_cost.update(archive_cpu_seconds=time.process_time()-acpu,
                               archive_wall_seconds=time.perf_counter()-awall)
    proposal_cpu=time.process_time()-cpu;raw=best['x'].copy();raw_phi=best['phi']
    raw_df=best['df'];phase='corrector';cstart,wstart=time.process_time(),time.perf_counter();accepted=0
    def objective(x):
        if time.process_time()-cstart>=cfg['corrector']['CPU_cap_seconds']:raise CorrectorBudgetEnd()
        value,g=model.physical(x,Y,True)
        record(x,g)
        return value,g
    def callback(x):
        nonlocal accepted
        accepted+=1
    try:
        result=minimize(objective,raw,jac=True,method='L-BFGS-B',bounds=[model.material_bounds]*model.N,
            callback=callback,options=dict(maxiter=cfg['corrector']['maxiter'],ftol=cfg['corrector']['ftol'],
                gtol=cfg['corrector']['gtol'],maxls=cfg['corrector']['maxls']))
        stop=str(result.message);success=bool(result.success)
    except CorrectorBudgetEnd:
        stop='CPU_CAP_60_SECONDS';success=False
    corrector=dict(stop=stop,success=success,nit=accepted,
        evaluations=sum(r['phase']=='corrector' for r in history),
        cpu_seconds=time.process_time()-cstart,wall_seconds=time.perf_counter()-wstart,
        allowance=cfg['corrector'])
    # Already-paid original-Phi incumbent, with complete physical fields.
    np.savez_compressed(out/'prediction.npz',initial=initial,raw=raw,prediction=best['x'],J=best['J'],E=best['E'],
                         gradient=best['gradient'] if best['gradient'] is not None else np.empty(0))
    np.savez_compressed(out/'paid_materials.npz',m=np.asarray(paid))
    dump(out/'history.json',history)
    receipt=dict(status='SEALED_PUBLIC_PREDICTION',case_id=cid,method=method,
        physical_phi=best['phi'],physical_data_residual=float(np.sqrt(best['df'])),
        prior=float(model.prior(best['x'])[0]),raw_phi=raw_phi,raw_data_residual=float(np.sqrt(raw_df)),
        counts=dict(model.counts),projection=projection_cost,auxiliary=auxiliary,
        corrector=corrector,setup=model.setup,proposal_cpu_seconds=proposal_cpu,
        cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        truth_access='NONE',scoring='NOT_RUN',material_bounds=list(model.material_bounds),
        delta=delta,selection=cfg['selection'],qualification=cfg['qualification'],
        backend='unchanged inherited R8 VectorPR3D/VectorCIE3D CM+RR TensorDDA',
        current_APCS_coordinates='inherited whiten_apcs=True projected IFFT range SVD',
        physical_value_cache=True,finalization_extra_physical_LU=0,
        cost_scope='modelsetup/proposal/physicalaudits/archive/commoncorrector/output; moduleimports excluded; nestedtimers overlap')
    dump(out/'RECEIPT.json',receipt)
    print(json.dumps(dict(case=cid,method=method,Phi=best['phi'],CPU=receipt['cpu_seconds'],
                          corrector=stop,projection_steps=(projection_cost or {}).get('iterations'))),flush=True)


def seal():
    cfg=json.loads((HERE/'FROZEN_PROTOCOL.json').read_text());rows=[];prefixes=[];arrays={}
    for cid in cfg['cases']:
        for method in cfg['methods']:
            out=HERE/'RESULTS'/cid/method;rec=json.loads((out/'RECEIPT.json').read_text())
            assert rec['status']=='SEALED_PUBLIC_PREDICTION' and (out/'prediction.npz').is_file()
            history=json.loads((out/'history.json').read_text())
            with np.load(out/'paid_materials.npz',allow_pickle=False) as z:mats=z['m']
            assert len(history)==len(mats)
            for limit in cfg['CPU_prefix_seconds']:
                eligible=[(r['phi'],i) for i,r in enumerate(history) if r['cpu_seconds']<=limit]
                if not eligible:continue
                value,index=min(eligible);key=f'prefix_{len(prefixes):04d}';arrays[key]=mats[index]
                prefixes.append(dict(case_id=cid,method=method,CPU_limit_seconds=limit,
                    candidate_CPU_seconds=history[index]['cpu_seconds'],physical_phi=value,
                    index=index,phase=history[index]['phase'],array_key=key,counts=history[index]['counts']))
            rows.append(dict(case_id=cid,method=method,prediction=str(out/'prediction.npz'),
                             physical_phi=rec['physical_phi'],receipt=str(out/'RECEIPT.json')))
    np.savez_compressed(HERE/'PREFIX_PREDICTIONS.npz',**arrays)
    dump(HERE/'PREFIX_SELECTION_SEAL.json',dict(status='SEALED_BEFORE_OFFLINE_GT_READ',rows=prefixes,
        selection='minimum already-paid original Phi before CPU budget, neverGT'))
    dump(HERE/'ALL12_PREDICTIONS_SEALED.json',dict(status='ALL12_PREDICTIONS_AND_PREFIX_CHOICES_SEALED',
        created_utc=datetime.now(timezone.utc).isoformat(),rows=rows,truth_access='NONE',hash_checks='NONE'))
    print('ALL12_PREDICTIONS_SEALED',flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('phase',choices=['run','solve','seal'])
    parser.add_argument('--case');parser.add_argument('--method');args=parser.parse_args()
    if args.phase=='solve':solve(args.case,args.method)
    elif args.phase=='seal':seal()
    else:
        cfg=json.loads((HERE/'FROZEN_PROTOCOL.json').read_text())
        for cid in cfg['cases']:
            for method in cfg['methods']:
                path=HERE/'RESULTS'/cid/method/'RECEIPT.json'
                if path.exists():
                    assert json.loads(path.read_text())['status']=='SEALED_PUBLIC_PREDICTION'
                    continue
                solve(cid,method)
        seal()


if __name__=='__main__':main()
