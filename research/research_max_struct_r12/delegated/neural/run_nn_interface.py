"""Frozen R8 NN interface comparison, public-only prepare/run; separate score.

One process uses one BLAS/PyTorch thread. No training or physics replacement.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import argparse
from datetime import datetime, timezone
import json
import sys
import time
import traceback
import platform

HERE=Path(__file__).resolve().parent
STUDY=HERE.parents[1]
REPO=STUDY.parents[1]
R8REPO=REPO.parent/'r8_worktree'
R8=R8REPO/'research/research_max_struct_r8'
sys.path.insert(0,str(REPO))


def dump(path,value):
    import numpy as np
    def convert(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value,indent=2,default=convert,allow_nan=False)+'\n')


def load_arrays(path):
    import numpy as np
    with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}


def prepare():
    """Freeze IDs/contracts BEFORE outcomes, then actual features/inference."""
    import numpy as np
    from research.research_max_struct_r8.experiments.run_case import load_public
    manifest=json.loads((R8/'configs/warm_primary_manifest.json').read_text())
    ids=['r8_test_0001','r8_test_0004','r8_test_0022','r8_test_0027']
    tasks=[]
    for ident in ids:
        matched=[t for t in manifest['tasks'] if t['task_id']==ident+'__cie_single__nn']
        assert len(matched)==1,ident
        task=matched[0]
        files=[R8REPO/p for p in task['public_files']]
        assert all(p.is_file() for p in files)
        tasks.append(dict(case_id=ident,cohort=task['cohort'],public_files=files,
                          original_task=task,historical_receipt=R8REPO/task['output_dir']/'RECEIPT.json'))
    protocol=dict(status='FROZEN_BEFORE_FEATURE_INFERENCE_AND_ANY_NEW_OUTCOME',
        created_utc=datetime.now(timezone.utc).isoformat(),ids=ids,tasks=tasks,
        arms=['uniform_direct','nn_direct','nn_dr','nn_cie_single'],
        weights=REPO/'research/research_max_struct_r8/RESULTS/training/physics/weights.pt',
        features='Unchanged R8 FeatureOperator native ten-channel public transform',
        retraining=False,initial_Phi_rejection=False,box=[.02,20.],
        objective='Original MaxwellModel Phi=.5||F(m)-Y||²/||Y||² + original smooth/TV prior',
        prior=dict(smooth=1e-4,tv=1e-4,epsilon=.01),
        dr=dict(steps=3000,audit_every=50,scale=6.,omega=1.,grouped_ports=True,
                initialization='NN material lifted through complete unchanged physical solve',
                data_ball='raw iid circular-complex sigma, full 99% global event; no whitening'),
        cie=dict(backend='unchanged R8 benchmark_solver/cie_fft cie_single',
                 initialize_physical_current=False,dpcs_rtol=.001,
                 stages='mf3/6/10 each500 then mf12/1500; inherited PR/beta/readout'),
        corrector=dict(maxiter=3000,maxfun=12000,scale=1e6,ftol=1e-12,gtol=.01,maxls=30,
                        cpu_cap_seconds=120.,cap_scope='corrector only; proposals/setup/inference separate'),
        selection='Minimum already-paid original Phi including prior; never GT',
        scoring='All16 predictions and four initializations sealed before first private read',
        scope='Four fixed exposed R8 cases, including historicalNN0027; not independent confirmation',
        truth_access='NONE')
    dump(HERE/'FROZEN_PROTOCOL.json',protocol)
    wall,cpu=time.perf_counter(),time.process_time()
    import torch
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    from research.research_max_struct_r8.src.warm_features import FeatureOperator
    from research.research_max_struct_r8.src.warm_network import load_predictor
    w,c=time.perf_counter(),time.process_time()
    predict,checkpoint=load_predictor(protocol['weights'],device='cpu')
    model_load=dict(cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w,
                    feature_mode=checkpoint['feature_mode'],architecture=checkpoint['architecture'],
                    epoch=checkpoint.get('selected_epoch',checkpoint.get('epoch')))
    operators={};rows=[]
    (HERE/'public').mkdir(exist_ok=True)
    for task in tasks:
        ident=task['case_id'];public=load_public([str(p) for p in task['public_files']])
        key=str(task['public_files'][0]);setup_charged=key not in operators
        if setup_charged:operators[key]=FeatureOperator(public)
        op=operators[key]
        features,basics,feature_cost=op(public['Y'])
        w,c=time.perf_counter(),time.process_time()
        material=np.asarray(predict(features)[0],float).ravel()
        infer=dict(cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w)
        assert material.shape==(576,) and np.min(material)>=.02 and np.max(material)<=20.
        np.savez_compressed(HERE/'public'/f'{ident}.npz',**public)
        np.savez_compressed(HERE/'public'/f'{ident}_nn.npz',m_init=material,features=features,
                             born=basics['born'],dpcs=basics['dpcs'])
        # Stored historical init is public: compare as provenance, never select.
        historical=R8REPO/task['original_task']['initialization_file']
        old=load_arrays(historical)['m_init'].ravel()
        row=dict(case_id=ident,cohort=task['cohort'],source_count=public['E0'].shape[1],
                  receiver_count=public['S'].shape[0],feature_setup=op.setup,
                  feature_setup_charged_here=setup_charged,features=feature_cost,inference=infer,
                  historical_initializer=str(historical),
                  relative_initializer_difference=float(np.linalg.norm(material-old)/np.linalg.norm(old)),
                  material_min=float(material.min()),material_max=float(material.max()),truth_access='NONE')
        rows.append(row)
    dump(HERE/'PREPROCESS_RECEIPT.json',dict(status='FOUR_FROZEN_PHYSICS_NN_INITIALIZATIONS_SEALED',
         records=rows,model_load=model_load,cpu_seconds=time.process_time()-cpu,
         wall_seconds=time.perf_counter()-wall,truth_access='NONE',LU=0,complex_RHS=0,
         imports='PyTorch/feature/network import charged in total; model_load nested; shared once'))
    print('PREPARED four actual R8 public cases',flush=True)


class CPUBudgetEnd(Exception):pass


def execute(case_id,arm):
    import numpy as np
    from scipy.optimize import minimize
    from scipy.stats import chi2
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    from research.research_max_struct_r12.src.field_current import DenseFieldCurrentProjectors,PhysicalAudit,run_projection
    out=HERE/'runs'/f'{case_id}__{arm}';out.mkdir(parents=True,exist_ok=False)
    wall,cpu=time.perf_counter(),time.process_time()
    public=load_arrays(HERE/'public'/f'{case_id}.npz')
    Y=public['Y'];model=MaxwellModel(public)
    m0=(np.full(model.n,.1) if arm=='uniform_direct'
        else load_arrays(HERE/'public'/f'{case_id}_nn.npz')['m_init'])
    delta=float(public['sigma'])*np.sqrt(chi2.ppf(.99,2*Y.size)/2)
    history=[];materials=[];best=None;phase='initial';last_cost={};projection_cost=None;auxiliary_cost=None
    def record(m,state,value,gradient=None):
        nonlocal best
        value=float(value)
        if best is None or value<best['objective']:
            best=dict(m=m.copy(),J=state['j'].copy(),Q=state['Q'].copy(),
                      Y_prediction=state['Y'].copy(),objective=value,
                      gradient=None if gradient is None else gradient.copy())
        history.append(dict(evaluation=len(history),phase=phase,objective=value,incumbent=best['objective'],
            cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
            relative_residual=float(np.linalg.norm(state['Y']-Y)/np.linalg.norm(Y)),
            prior=float(model.reg.evaluate(m)[0]),physics_cost=dict(model.cost)))
        materials.append(m.copy())
    def audit(m):
        nonlocal last_cost
        state=model.state(m);value=model.physical_score(m,Y);record(m,state,value)
        now=dict(model.cost);increments={k:float(v-last_cost.get(k,0)) for k,v in now.items()};last_cost=now
        metrics=dict(relative_residual=history[-1]['relative_residual'],prior=history[-1]['prior'],
                      noise_feasible=bool(np.linalg.norm(state['Y']-Y)<=delta))
        return PhysicalAudit(value,state['j'],state['Q'],metrics,increments)
    audit(m0)
    np.savez_compressed(out/'initial.npz',m=m0,J=best['J'],Q=best['Q'],Y_prediction=best['Y_prediction'])
    initial_objective=best['objective']
    if arm=='nn_dr':
        phase='dr_proposal'
        pr=DenseFieldCurrentProjectors(model.G,model.E,model.S,Y,delta,scale=6.,lower=.02,upper=20.,grouped_ports=True)
        projection=run_projection('dr',pr,m0,audit,steps=3000,audit_every=50)
        projection_cost=dict(setup_cpu_seconds=projection.setup_cpu_seconds,
            setup_wall_seconds=projection.setup_wall_seconds,cpu_seconds=projection.cpu_seconds,
            wall_seconds=projection.wall_seconds,total_counts=dict(pr.costs),
            iterations=projection.iterations,stop_reason=projection.stop_reason)
        projection.save(out/'projection')
    elif arm=='nn_cie_single':
        phase='cie_proposal'
        from research.research_max_struct_r8.src.benchmark_solver import arm_config,run_auxiliary
        config=json.loads((REPO/'research/research_max_struct_r8/configs/solver.json').read_text())
        numerical=arm_config('cie_single',config['common'],Y,model.E,config['continuation_stages'])
        numerical.update(initialize_physical_current=False,
            fft_geometry=dict(coords=public.get('coords',public.get('centers')),radius=float(public['radius']),k=float(public['k'])),
            aplus_regularizer=model.reg,regularizer=lambda m:model.reg.evaluate(m)[:2],physical_sigma=float(np.linalg.norm(Y)))
        class ObservedModel:
            def __getattr__(self,name):return getattr(model,name)
            def state(self,m):
                state=model.state(m);record(m,state,model.physical_score(m,Y));return state
        aux=run_auxiliary(ObservedModel(),Y,m0,'cie_single',numerical)
        auxiliary_cost=dict(status=aux['status'],execution_receipt=aux['execution_receipt'],cost=aux['cost'],
                             history=aux['history'],candidate_scores=[{k:v for k,v in row.items() if k!='m'} for row in aux['candidates']])
        np.savez_compressed(out/'auxiliary_candidates.npz',m=np.array([r['m'] for r in aux['candidates']]))
        dump(out/'AUXILIARY.json',auxiliary_cost)
    proposal_cpu=time.process_time()-cpu
    proposal_wall=time.perf_counter()-wall
    precorrector={k:(v.copy() if isinstance(v,np.ndarray) else v) for k,v in best.items()}
    np.savez_compressed(out/'selected_before_corrector.npz',**{k:v for k,v in precorrector.items() if isinstance(v,np.ndarray)})
    phase='corrector';cstart,wstart=time.process_time(),time.perf_counter();accepted=0
    def objective(m):
        if time.process_time()-cstart>=120.:raise CPUBudgetEnd()
        value,g=model.data_value_gradient(m,Y)
        record(m,model._current,value,g)
        return 1e6*value,1e6*g
    def accepted_callback(m):
        nonlocal accepted
        accepted+=1
    try:
        result=minimize(objective,precorrector['m'],jac=True,method='L-BFGS-B',
            bounds=[(.02,20.)]*model.n,callback=accepted_callback,
            options=dict(maxiter=3000,maxfun=12000,ftol=1e-12,gtol=.01,maxls=30))
        termination=str(result.message);success=bool(result.success)
    except CPUBudgetEnd:
        termination='CPU_CAP_120_SECONDS';success=False
    corrector=dict(stop=termination,success=success,accepted_iterations=accepted,
        evaluations=sum(row['phase']=='corrector' for row in history),
        cpu_seconds=time.process_time()-cstart,wall_seconds=time.perf_counter()-wstart,
        allowance=dict(maxiter=3000,maxfun=12000,cpu_seconds=120.))
    # Export already-paid original-Phi incumbent; no finalization Maxwell solve.
    arrays={k:v for k,v in best.items() if isinstance(v,np.ndarray)}
    np.savez_compressed(out/'prediction.npz',**arrays)
    np.savez_compressed(out/'paid_materials.npz',m=np.asarray(materials))
    dump(out/'history.json',history)
    residual=float(np.linalg.norm(best['Y_prediction']-Y)/np.linalg.norm(Y))
    receipt=dict(status='SEALED_PUBLIC_PREDICTION',case_id=case_id,arm=arm,
        initial_objective=initial_objective,precorrector_objective=precorrector['objective'],
        objective=best['objective'],relative_residual=residual,delta=delta,
        noise_feasible=bool(np.linalg.norm(best['Y_prediction']-Y)<=delta),
        prior=float(model.reg.evaluate(best['m'])[0]),corrector=corrector,
        proposal_cpu_seconds=proposal_cpu,proposal_wall_seconds=proposal_wall,
        projection_cost=projection_cost,auxiliary_cost=auxiliary_cost,
        physical_cost=dict(model.cost),physical_operation_wall=dict(model.timings),
        cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        upstream_preprocess='PREPROCESS_RECEIPT.json; features/inference shared but attributable perNNarm',
        truth_access='NONE',scoring='NOT_RUN',device='CPU',threads=1,
        execution_platform=platform.platform(),python=sys.executable,
        selection='Minimum original Phi among initial/paidproposal/correctorevaluations including prior',
        finalization_physical_cost=0)
    dump(out/'RECEIPT.json',receipt)
    print(json.dumps(dict(case_id=case_id,arm=arm,Phi=best['objective'],cpu_seconds=receipt['cpu_seconds'],corrector=termination)),flush=True)
    return receipt


def run():
    from threadpoolctl import threadpool_limits
    protocol=json.loads((HERE/'FROZEN_PROTOCOL.json').read_text())
    with threadpool_limits(limits=1):
        for ident in protocol['ids']:
            for arm in protocol['arms']:
                out=HERE/'runs'/f'{ident}__{arm}'
                if (out/'RECEIPT.json').exists():
                    assert json.loads((out/'RECEIPT.json').read_text())['status']=='SEALED_PUBLIC_PREDICTION'
                    continue
                execute(ident,arm)
        entries=[]
        for ident in protocol['ids']:
            for arm in protocol['arms']:
                out=HERE/'runs'/f'{ident}__{arm}'
                receipt=json.loads((out/'RECEIPT.json').read_text())
                assert receipt['status']=='SEALED_PUBLIC_PREDICTION' and (out/'prediction.npz').is_file()
                entries.append(dict(case_id=ident,arm=arm,prediction=str(out/'prediction.npz'),receipt=str(out/'RECEIPT.json'),objective=receipt['objective']))
        dump(HERE/'ALL_PREDICTIONS_SEALED.json',dict(status='ALL16_PREDICTIONS_AND4_INITIALIZATIONS_SEALED',
             created_utc=datetime.now(timezone.utc).isoformat(),entries=entries,truth_access='NONE',
             initializations=[str(HERE/'public'/f'{i}_nn.npz') for i in protocol['ids']],hash_checks='NONE'))
    print('ALL16_PREDICTIONS_SEALED',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('phase',choices=['prepare','run']);args=parser.parse_args()
    prepare() if args.phase=='prepare' else run()
