"""Finite public-data inference batch. No private truth or morphology access."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import argparse,json,time,sys,gzip,traceback
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from threadpoolctl import threadpool_limits
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
from src.physics import load_model
from src.full_maxwell_corrector import corrector
from src.cie_tsom_baseline import run_cie
from src.controlled_release import run_release
from src.physical_budget import install_budget,PhysicalBudgetStop


def clean(x):
    if isinstance(x,np.ndarray):return clean(x.tolist())
    if isinstance(x,np.generic):return clean(x.item())
    if isinstance(x,complex):return dict(real=x.real,imag=x.imag)
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items() if not callable(v)}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,float) and not np.isfinite(x):return None
    return x


def infer(task):
    case_id,method,start_index,wall=task['case_id'],task['method'],task['start'],task['wall']
    tag=task['tag']
    outdir=HERE/'results'/tag
    public_dir=Path(task.get('public_dir',HERE/'data/e1_public'))
    model,public=load_model(public_dir/f'{case_id}.npz')
    Y=public['Y']
    m0=np.load(public_dir/f'starts_{model.shape[0]}.npz')['starts'][start_index]
    before=Counter(model.cost)
    physics_budget=None
    if task.get('max_lu') is not None or task.get('max_rhs') is not None:
        physics_budget=install_budget(model,Y,task.get('max_lu'),task.get('max_rhs'))
    started=time.perf_counter()
    cfg=dict(box=(.02,20.),max_nfev=400,wall_seconds=wall)
    history=[]
    candidates=[]
    status='OK'
    proposal_wall=0.
    try:
        with threadpool_limits(limits=1):
            if method in ('native','reciprocal','data','feedback','globaldual'):
                out=corrector(model,Y,m0,{**cfg,'method':method,'beta':.2})
                candidates=out['candidates'];history=out['history']
                solver_cost=out['cost'];stop=out['stop']
            else:
                # Reserve the same full-material corrector allocation for every
                # auxiliary method. All candidates are scored using original data.
                proposal_budget=.65*wall
                if method.startswith('cie') or method in ('tsom','csi'):
                    beta_j=1. if method=='cie_dual' else .2
                    widths=([1,2] if model.shape[0]==4 else [1,2,3])
                    betas=([1.,.2] if len(widths)==2 else [1.,.2,.2])
                    stages=[dict(mf=w,beta=b,beta_current=(1. if method=='cie_dual' and i else b),maxiter=700)
                            for i,(w,b) in enumerate(zip(widths,betas))]
                    if method=='cie_matched':
                        stages=[dict(mf=2,beta=.2,maxiter=700)]
                    paper=method=='cie_paperstop'
                    ccfg=dict(stages=stages,method=('cie' if method.startswith('cie') else method),
                        dpcs_rtol=1e-3,material_box=(.02,20.),rho=1.,
                        coefficient_tol=(1e-3 if paper else 1e-8),miniter=(2 if paper else 50),
                        regularizer=lambda m:model.regularizer(m)[:2],physical_sigma=np.linalg.norm(Y),
                        max_wall_seconds=proposal_budget,material_inner_maxiter=60)
                    out=run_cie(model,Y,m0,ccfg)
                    candidates=out['candidates'];history=out['history'];stop=out['status']
                elif method=='material_seeds':
                    jj=model.jacobian(m0)/np.linalg.norm(Y)
                    r=model.pack(model.forward(m0)-Y)/np.linalg.norm(Y)
                    d=-(jj.T@r)
                    d*=2./max(np.linalg.norm(d)/np.sqrt(model.n),1e-14)
                    candidates=[dict(m=np.clip(m0+sgn*d,.02,20),tag=f'material_seed_{sgn}') for sgn in (1,-1)]
                    out=dict(cost=dict(material_seed_count=2));stop='TWO_PUBLIC_GRADIENT_SEEDS'
                else:
                    options={
                        'n1_fourier':dict(selector='fourier',release='hard',finite_seed=False),
                        'n1_energy':dict(selector='energy',release='hard',finite_seed=False),
                        'n1_random':dict(selector='random',release='hard',finite_seed=False),
                        'n1_normal':dict(selector='normal',release='hard',finite_seed=False),
                        'n1_task':dict(selector='material_step',release='hard',finite_seed=False),
                        'n1_soft':dict(selector='material_step',release='soft',finite_seed=False),
                        'n1_replace':dict(selector='material_step',release='soft',operation='replace',finite_seed=False),
                        'n1_finite':dict(selector='material_step',release='soft',finite_seed=True,seed_rule='task'),
                        'n1_finite_fourier':dict(selector='fourier',release='soft',finite_seed=True,seed_rule='fourier'),
                        'n1_finite_random':dict(selector='random',release='soft',finite_seed=True,seed_rule='random'),
                        'n1_finite_krylov':dict(selector='krylov',release='soft',finite_seed=True,seed_rule='task'),
                    }
                    rcfg=dict(chart=task.get('chart','csi'),beta=.2,rho=1.,state_weight=1/np.linalg.norm(model.E),data_scale=np.linalg.norm(Y),
                        box=(.02,20.),eta_target=.5,lm=1e-4,
                        backbone_receiver_atoms=task.get('backbone',1),backbone_receiver_rtol=1e-3,backbone_fourier_atoms=1,
                        group_atoms=1,pool_groups=8,finalist_jvps=2,episodes=1,
                        seed_material_rms=2.,seed_scale_tolerance=.25,seed_current_scale=np.linalg.norm(model.E),
                        seed_defect_scale=np.linalg.norm(model.E),seed_current_cap=3.,seed_defect_cap=3.,
                        bounded_seed_scale=True,seed_scale_trials=8,seed_min_material_rms=.02,
                        hold_steps=3,joint_steps=task.get('joint_steps',8),rho_schedule=[1.],max_wall_seconds=proposal_budget,
                        qualification_residual=.001,**options[method])
                    rcfg['execute_local']=not rcfg['finite_seed']
                    if task.get('affine_mode')=='native_tsom':
                        rcfg.update(affine_mode='native_tsom',native_mf=2,native_core_apcs=15,
                            state_weight=1.,state_source_whitener=1/np.linalg.norm(model.E,axis=0),
                            data_whitener=np.repeat(1/np.linalg.norm(Y,axis=0),2*model.nr),
                            physical_data_scale=np.linalg.norm(Y),local_gain_min=-np.inf,
                            operation='replace')
                    out=run_release(model,Y,m0,rcfg)
                    candidates=out['candidates'];history=out['history'];stop=out['stop']
                proposal_wall=time.perf_counter()-started
                solver_cost=out.get('cost',{})
                # The initial physical point is always an incumbent; useful seeds
                # cannot be erased simply because an auxiliary score was large.
                seeds=[c for c in candidates if np.linalg.norm(np.asarray(c['m'])-m0)>1e-8]
                seeds.sort(key=lambda c:model.physical_score(np.asarray(c['m']),Y))
                seeds=seeds[:2] or [dict(m=m0.copy(),tag='no_proposal_anchor_fallback')]
                corrected=[]
                for i,seed in enumerate(seeds):
                    remaining=wall-(time.perf_counter()-started)
                    if remaining<=0:break
                    c=corrector(model,Y,np.asarray(seed['m']),dict(box=(.02,20.),method='trf',max_nfev=60,
                        wall_seconds=min(.35*wall/len(seeds),remaining)))
                    corrected.append(dict(m=c['m'],tag=f'common_corrector_{i}',history=c['history'],cost=c['cost']))
                candidates=list(candidates)+corrected
                initial_score=model.physical_score(m0,Y)
                m=m0.copy();best=initial_score
                for c in candidates:
                    score=model.physical_score(np.asarray(c['m']),Y)
                    c['physical_score']=score
                    c['relative_residual']=float(np.linalg.norm(model.state(np.asarray(c['m']))['Y']-Y)/np.linalg.norm(Y))
                    if score<best:m=np.asarray(c['m']).copy();best=score
                out['m']=m
            m=np.asarray(out['m'])
            pred=model.forward(m)
            score=model.physical_score(m,Y)
            relative=float(np.linalg.norm(pred-Y)/np.linalg.norm(Y))
    except PhysicalBudgetStop as exc:
        status='PHYSICAL_BUDGET_STOP';stop=str(exc)
        endpoint=physics_budget['incumbent']
        m,pred=endpoint['m'],endpoint['prediction']
        score,relative=endpoint['physical_score'],endpoint['relative_residual']
        candidates=list(candidates)+physics_budget['candidates']
        history.append(dict(event='physical_budget_stop',reason=stop))
        solver_cost=dict(physical_caps=dict(max_lu=task.get('max_lu'),max_rhs=task.get('max_rhs')))
    except Exception as exc:
        status='FAILED_EXECUTION';stop=repr(exc);solver_cost={}
        history.append(dict(event='exception',error=traceback.format_exc()))
        m=m0.copy();pred=model.forward(m);score=model.physical_score(m,Y)
        relative=float(np.linalg.norm(pred-Y)/np.linalg.norm(Y))
    elapsed=time.perf_counter()-started
    label_method=method+('_chart_cie' if task.get('chart')=='cie' and method.startswith('n1_') else '')
    if task.get('affine_mode')=='native_tsom' and method.startswith('n1_'):
        label_method+='_native_tsom'
    run_id=f'{case_id}__s{start_index}__{label_method}__b{wall:g}'
    candidate_m=np.stack([np.asarray(c['m']) for c in candidates]) if candidates else np.empty((0,model.n))
    np.savez_compressed(outdir/(run_id+'.npz'),m=m,prediction=pred,candidates=candidate_m)
    cost=dict(Counter(model.cost)-before)
    row=dict(run_id=run_id,case_id=case_id,side=model.shape[0],layout=task['layout'],start=start_index,
        method=label_method,base_method=method,chart=task.get('chart','csi'),affine_mode=task.get('affine_mode','physical_anchor'),
        wall_budget=wall,status=status,stop=stop,wall_seconds=elapsed,
        proposal_wall=proposal_wall,physical_score=score,relative_residual=relative,
        candidate_count=len(candidate_m),model_cost=cost,solver_cost=clean(solver_cost))
    row['physical_caps']=dict(max_lu=task.get('max_lu'),max_rhs=task.get('max_rhs'))
    with gzip.open(outdir/(run_id+'.json.gz'),'wt') as f:
        json.dump(clean(dict(row=row,history=history,candidates=[{k:v for k,v in c.items() if k!='m'} for c in candidates])),f)
    return clean(row)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--tag',default='e1_initial')
    ap.add_argument('--methods',nargs='+',default=['native','feedback','cie_single','cie_dual','tsom','csi',
        'n1_fourier','n1_energy','n1_random','n1_normal','n1_task','n1_soft','n1_replace','n1_finite',
        'n1_finite_fourier','n1_finite_random','material_seeds'])
    ap.add_argument('--starts',nargs='+',type=int,default=[0])
    ap.add_argument('--sides',nargs='+',type=int,default=[4,6])
    ap.add_argument('--layouts',nargs='+',default=['uniform'])
    ap.add_argument('--wall',type=float,default=4.)
    ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--cases',nargs='+')
    ap.add_argument('--backbone',choices=['1','observable'],default='1')
    ap.add_argument('--joint-steps',type=int,default=8)
    ap.add_argument('--chart',choices=['csi','cie'],default='csi')
    ap.add_argument('--max-lu',type=int)
    ap.add_argument('--max-rhs',type=int)
    ap.add_argument('--affine-mode',choices=['physical_anchor','native_tsom'],default='physical_anchor')
    args=ap.parse_args()
    cases=json.loads((HERE/'data/e1_case_ids.json').read_text())
    cases=[c for c in cases if c['side'] in args.sides and c['layout'] in args.layouts and
           (not args.cases or c['case_id'] in args.cases)]
    tasks=[dict(**c,start=s,method=m,wall=args.wall,tag=args.tag,
                backbone=(1 if args.backbone=='1' else args.backbone),joint_steps=args.joint_steps,chart=args.chart,
                max_lu=args.max_lu,max_rhs=args.max_rhs,affine_mode=args.affine_mode)
           for c in cases for s in args.starts for m in args.methods]
    out=HERE/'results'/args.tag;out.mkdir(exist_ok=True,parents=True)
    plan=dict(status='FROZEN_BEFORE_PREDICTIONS',tasks=tasks,truth_access=False,
        common_physics='frozen R4 disk-pulse reassembled per grid',real_material_box=[.02,20],
        public_initializations='starts_side.npz, independent of GT',main_budget='common actual wall cap; auxiliary share65%, common corrector share35%',
        common_corrector='bounded analytic TRF, at most2 seeds,60 eval each; actual seed counts recorded',
        auxiliary_scope='N1 all selectors same pool, current real DOF=2*p per spatial complex atom',
        seed_constraint='RMS2 maximum; same current/E and defect/E caps3; bounded feasible scale search8 trials; actual scales reported')
    (out/'PLAN.json').write_text(json.dumps(plan,indent=2))
    started=time.perf_counter();rows=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(infer,t) for t in tasks]
        for future in as_completed(futures):
            row=future.result();rows.append(row)
            with (out/'rows.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
            print(json.dumps({k:row[k] for k in ('run_id','status','relative_residual','wall_seconds')}),flush=True)
    summary=dict(status='PREDICTIONS_CLOSED_TRUTH_NOT_READ',rows=rows,count=len(rows),
        batch_wall=time.perf_counter()-started,workers=args.workers,blas_threads_per_worker=1)
    (out/'SUMMARY.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps({k:v for k,v in summary.items() if k!='rows'}),flush=True)


if __name__=='__main__':main()
