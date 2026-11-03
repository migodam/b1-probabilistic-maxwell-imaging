"""Independent32-cell-side PUBLIC adapter; frozen24 adapter is never changed.

Only original MaxwellModel, run_onehop_cie/global_cie and R4 refine are used.
No generator, original mapping, disclosed material or scoring is imported/read.
"""
from pathlib import Path
from collections import Counter
import argparse,json,sys,time,traceback
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parent.parents[1]
sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r2.src.cie_global import run_cie as global_cie
from research.research_max_struct_r4.src.cie_onehop import run_onehop_cie
from research.research_max_struct_r4.experiments.cooperative_paths import refine

ALLOW={'G','E0','S','Y','grid_shape','shape','sigma','noise_radius','host','amplitude',
    'material_bounds','smooth_weight','tv_weight','tv_epsilon','prior_width','coords','sources','receivers','k','radius'}
STAGES=[dict(mf=3,beta=1.,maxiter=500),dict(mf=6,beta_material=.2,beta_current=1.,maxiter=500),
        dict(mf=10,beta_material=.2,beta_current=1.,maxiter=500)]


def metric(model,Y,m):
    state=model.state(m);residual=state['Y']-Y
    return dict(m=np.asarray(m).copy(),objective=float(.5*np.vdot(residual,residual).real/np.vdot(Y,Y).real+model.reg.evaluate(m)[0]),
                relative_residual=float(np.linalg.norm(residual)/np.linalg.norm(Y)))


def histories(history):
    result={}
    for key in sorted({key for row in history for key in row}):
        values=[row.get(key) for row in history]
        if all(x is None or isinstance(x,(int,float,bool,np.number)) for x in values):
            result['history_'+key]=np.array([np.nan if x is None else x for x in values])
        elif all(x is None or isinstance(x,str) for x in values):
            result['history_'+key]=np.array(['' if x is None else x for x in values],dtype=str)
    return result


def solve(public,method,config,backend='cpu'):
    wall,cpu=time.perf_counter(),time.process_time();model=None
    result=dict(status='RUNNING',classification='DISCLOSED_DEVELOPMENT_GRID_DIAGNOSTIC',
        method=method,truth_access='NONE_TO_INVERSION',configuration=config,trajectory_arrays={},failed_attempts=[],
        scientific_wallcap=None,operator_assembly='Generator cost separately recorded, public G/E0/S loaded here',
        rhs_semantics='physical source/adjoint RHS separate from onehop small-block matrix RHS/scalar divisions')
    try:
        if backend!='cpu' or method not in ('onehop','global'):raise ValueError('Only frozen CPU onehop/global')
        if set(public)!=ALLOW:raise ValueError('Public payload keys differ; refuse payload')
        if (config['stages']!=STAGES or config['corrector_maxiter']!=1000 or config['host']!=.1
            or config['material_box']!=[.02,20.] or config['sourcewise_cg'] is not False):
            raise ValueError('Frozen schedule/host/box/source metric mismatch')
        for key in ('max_wall_seconds','budget_check','max_iterations','receipt_path'):
            if key in config:raise ValueError('No scientific wallcap or hidden budget override')
        with threadpool_limits(limits=1):
            tick,ct=time.perf_counter(),time.process_time()
            model=MaxwellModel(public)
            if tuple(model.shape)!=(32,32) or (model.n,model.ns,model.nr)!=(1024,12,48):
                raise ValueError('Declared32 geometry/acquisition mismatch')
            for key,value in dict(host=.1,amplitude=12.,smooth_weight=1e-4,tv_weight=1e-4,tv_epsilon=.01,prior_width=19.5,radius=.0135).items():
                if not np.isclose(float(public[key]),value,rtol=0,atol=1e-15):raise ValueError('Public frozen parameter mismatch: '+key)
            result['setup']=dict(wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-ct,physical_cost=dict(model.cost))
            Y=public['Y'];m0=np.full(model.n,.1);initial=metric(model,Y,m0)
            if method=='onehop':
                from research.research_max_struct_r4.src.cie_onehop import OneHopAbsorption
                bound_cost=Counter();blocks=OneHopAbsorption(model.G,model.shape,bound_cost)
                result['whole_box_bounds']=[blocks.whole_box_bound((.02,20.),beta) for beta in (.2,1.)]
                result['bound_setup_counts']=dict(bound_cost)
            numerical=dict(config,method='cie',regularizer=lambda x:model.reg.evaluate(x)[:2],physical_sigma=float(np.linalg.norm(Y)))
            tick,ct,before=time.perf_counter(),time.process_time(),Counter(model.cost)
            if method=='onehop':aux=run_onehop_cie(model,Y,m0,numerical)
            else:aux=global_cie(model,Y,m0,numerical)
            result['auxiliary']=dict(status=aux['status'],wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-ct,
                operation_counts=aux['cost'],execution_receipt=aux.get('execution_receipt'),physical_cost=dict(Counter(model.cost)-before),
                stage_history=[row for row in aux['history'] if row.get('event')=='stage_end'])
            result['trajectory_arrays'].update(histories(aux['history']))
            scored=[metric(model,Y,c['m']) for c in aux['candidates']]
            seed=min([initial]+scored,key=lambda c:c['objective'])
            result['candidate_scores']=[{key:value for key,value in row.items() if key!='m'} for row in scored]
            result['trajectory_arrays'].update(candidate_m=np.array([c['m'] for c in scored]),
                candidate_original_Phi=np.array([c['objective'] for c in scored]),selected_before_corrector_m=seed['m'])
            tick,ct,before=time.perf_counter(),time.process_time(),Counter(model.cost)
            refined,corrector=refine(model,Y,seed['m'],maximum=1000)
            result['trajectory_arrays']['corrector_history']=np.array(corrector.pop('history'))
            result['corrector']=dict(corrector,cpu_seconds=time.process_time()-ct,physical_cost=dict(Counter(model.cost)-before))
            best=min([initial,seed,metric(model,Y,refined)],key=lambda c:c['objective']);state=model.state(best['m'])
            result['trajectory_arrays'].update(m=best['m'],Y_prediction=state['Y'],j=state['j'],Q=state['Q'])
            result.update(status='SEALED_DISCLOSED_DEVELOPMENT_PREDICTION',physical_objective=best['objective'],
                physical_relative_residual=best['relative_residual'],physical_whitened_residual=float(np.linalg.norm(state['Y']-Y)/float(public['sigma'])),
                noise_radius99=float(public['noise_radius']),linear_relative_residual=float(state['state_relative_residual']),
                selection='Original normalized Phi plus historical prior, host and true3stage endpoints then one common refine1000; retain incumbent')
    except Exception as exc:
        result.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc())
        result['failed_attempts'].append(dict(error=repr(exc),charged=True))
        if hasattr(exc,'cie_execution_receipt'):result['failed_auxiliary_execution']=exc.cie_execution_receipt
    finally:
        result.update(wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
            physical_model_cost=dict(model.cost) if model else {},physical_operation_wall=dict(model.timings) if model else {},
            physical_operation_cpu=dict(model._physics.cpu_timings) if model else {},
            cost_scope='adapter setup/initial/auxiliary/ranking/corrector/final state; caller must charge public read/output; imports excluded and not claimed free')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--method',required=True,choices=['onehop','global'])
    args=parser.parse_args();cfg=json.loads((ROOT/'config.json').read_text())
    output=ROOT/'results'/args.method
    if output.exists():raise FileExistsError('Preserve existing diagnostic result/attempt')
    output.mkdir(parents=True)
    wall,cpu=time.perf_counter(),time.process_time()
    with np.load(ROOT/'data_public/grid32_curved_bar.npz',allow_pickle=False) as archive:
        public={key:archive[key] for key in archive.files}
    result=solve(public,args.method,cfg['methods'][args.method])
    arrays=result.pop('trajectory_arrays')
    if arrays:np.savez_compressed(output/'prediction_and_trajectory.npz',**arrays)
    result.update(caller_wall_seconds=time.perf_counter()-wall,caller_cpu_seconds=time.process_time()-cpu,
        caller_scope='public load/solve/output serialization; imports and receipt own write excluded')
    (output/'RECEIPT.json').write_text(json.dumps(result,indent=2,default=lambda x:x.tolist(),allow_nan=False)+'\n')
    print(json.dumps({key:result.get(key) for key in ('status','method','physical_relative_residual','wall_seconds','physical_model_cost')}))
