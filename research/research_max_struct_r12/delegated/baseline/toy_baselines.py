"""R12 exact Pro toy bridge to mature R8 CIE/TSOM and R10 Green readout.

No GT access in prepare/run. Score is an explicit separate offline phase.
Existing outputs are never replaced. Run all arms before invoking score.
"""
from pathlib import Path
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[key]='1'
import argparse, importlib.util, json, sys, time
from collections import Counter
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import svd
from threadpoolctl import threadpool_limits, threadpool_info

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[3]
sys.path.insert(0,str(REPO))
SOURCE=REPO/'research/research_max_struct_r12/sources/pro/MAX_STRUCT_RESEARCH/feasibility'
CASES=('ring_high','split_shell_high','bridge_high','disk_low')
ARMS=('cie_single','cie_dual','global_dual','a_plus','green_strong_full')
SUITE='scaled'

def suite_path(stem,suffix=''):
    return HERE/(stem+('' if SUITE=='scaled' else '_clipped')+suffix)

def dump(path,obj):
    def convert(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False,default=convert)+'\n')

def toy_module():
    spec=importlib.util.spec_from_file_location('pro_exact_toy',SOURCE/'consensus_reflection_toy.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def prepare():
    out=HERE/'public';out.mkdir(exist_ok=False)
    cpu,wall=time.process_time(),time.perf_counter()
    xy,G,E,S=toy_module().setup()
    for case in CASES:
        # Only the named observation is accessed from the supplied source file.
        with np.load(SOURCE/(case+'.npz'),allow_pickle=False) as z:Y=z['Y'].copy()
        np.savez_compressed(out/(case+'.npz'),G=G,E0=E,S=S,Y=Y,grid_shape=[10,10],coords=xy)
    results=json.loads((SOURCE/'toy_results.json').read_text())
    noise={r['case']:dict(sigma=r['sigma_complex_rms'],delta=r['delta_CN99']) for r in results['cases']}
    dump(out/'NOISE.json',noise)
    dump(out/'PREPARATION.json',dict(truth_access='NONE',source='Pro setup() and existing Y only',
         cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall))

class ExactDenseGreen:
    """Injection point only: supplied exact G, with unchanged core PR/readouts."""
    G=None
    def __init__(self,coords,shape,radius,k,cost):
        self.cost=cost;self.spectrum=np.empty(0,complex)
    def apply(self,x,adjoint=False):
        self.cost['dense_G_adjoint_calls' if adjoint else 'dense_G_forward_calls']+=1
        self.cost['dense_G_complex_RHS']+=x.shape[1]
        return (self.G.conj().T if adjoint else self.G)@x
    def adjoint(self,x):return self.apply(x,True)

def numerical_config(arm,Y,E,rank):
    # Same mature beta schedule, iteration allowances and inner settings.
    # Spatial frequencies are scaled by 10/24 and rounded (half up), with
    # Nyquist mf5 full mask containing exactly 100 distinct frequencies.
    single=arm=='cie_single'
    masks=(1,3,4,5) if SUITE=='scaled' else (3,5,5,5)
    return dict(arm='a_plus' if arm=='green_strong_full' else arm,method='cie',
        suite=SUITE,
        material_box=[0.,12.],dpcs_rank=rank,sourcewise_cg=arm.startswith('cie_'),
        coefficient_tol=0.,miniter=50,material_inner_maxiter=60,
        material_inner_ftol=1e-12,material_inner_gtol=1e-8,
        data_weights=1/np.maximum(np.sum(abs(Y)**2,axis=0),1e-300),
        state_weights=1/np.maximum(np.sum(abs(E)**2,axis=0),1e-300),
        stages=[dict(mf=mf,beta_material=bm,beta_current=bm if single else 1.,maxiter=its)
                for mf,bm,its in zip(masks,(1.,.2,.2,.2),(500,500,500,1500))],
        aplus_mode='mobius',aplus_arc_candidates=True,aplus_inherited_gtol=True,
        aplus_sweeps=2,aplus_step_tol=1e-8)

def run(case,arm):
    out=suite_path('runs')/(case+'__'+arm);out.mkdir(parents=True,exist_ok=False)
    wall,cpu=time.perf_counter(),time.process_time()
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel,Regularizer
    from research.research_max_struct_r8.src import cie_fft
    from research.research_max_struct_r8.src.benchmark_solver import run_auxiliary
    from research.research_max_struct_r9.src.baseline_adapter import noise_ball_rank
    from research.research_max_struct_r10.src.lowrank_denominator import prepare_green_strong,matrix_material_readout
    from research.research_max_struct_r6.src.fft_tsom_ops import NativeFourierAPCS
    with np.load(HERE/'public'/(case+'.npz'),allow_pickle=False) as z:p={k:z[k] for k in z.files}
    noise=json.loads((HERE/'public/NOISE.json').read_text())[case]
    Y=p['Y'];model=MaxwellModel(p);model.reg=Regularizer(model.shape,smooth=0.,tv=0.)
    rank_decision=noise_ball_rank(model.S,Y,noise['sigma']);dump(out/'DPCS.json',rank_decision)
    if rank_decision['rank'] is None:raise RuntimeError(rank_decision['status'])
    cfg=numerical_config(arm,Y,model.E,rank_decision['rank'])
    dump(out/'FROZEN_CONFIG.json',cfg)
    cfg.update(fft_geometry=dict(coords=p['coords'],radius=1.,k=1.),
        physical_sigma=float(np.linalg.norm(Y)),
        regularizer=lambda m:model.reg.evaluate(m)[:2],aplus_regularizer=model.reg,
        receipt_path=str(out/'AUXILIARY_EXECUTION.json'))
    # Audit actual independent affine current rank, distinct from coefficients.
    _,_,vh=svd(model.S,full_matrices=False);vd=vh[:rank_decision['rank']].conj().T
    capacity=[]
    for stage in cfg['stages']:
        u=NativeFourierAPCS(model.shape,stage['mf'],vd)
        q=u.apply(np.eye(u.k,dtype=complex));s=svd(q,compute_uv=False)
        independent=int(np.count_nonzero(s>1e-10))
        capacity.append(dict(mf=stage['mf'],distinct_fourier_columns=u.k,
            independent_complex_APCS_per_source=independent,
            real_free_current_dof=2*independent*model.ns,
            fixed_affine_DPCS_rank=rank_decision['rank'],material_dof=model.n))
    dump(out/'CAPACITY.json',capacity)
    states=[]
    def paid(m,phase,value=None):
        st=model.state(m)
        value=model.physical_score(m,Y) if value is None else float(value)
        row=dict(cpu_seconds=time.process_time()-cpu,phase=phase,objective=value,
                 physical_cost=dict(model.cost),m=np.asarray(m).copy())
        states.append(row)
        return value
    stage_cpu=[]
    cfg['stage_callback']=lambda row:stage_cpu.append(dict(stage=row['stage'],cpu_seconds=time.process_time()-cpu))
    prep_cpu=time.process_time()
    if arm=='green_strong_full':
        cfg=prepare_green_strong(model,cfg)
        cfg['material_readout']=matrix_material_readout
    setup_cpu=time.process_time()-prep_cpu
    ExactDenseGreen.G=model.G
    original=cie_fft.DiskPulseToeplitzGreen;cie_fft.DiskPulseToeplitzGreen=ExactDenseGreen
    aux_wall,aux_cpu=time.perf_counter(),time.process_time()
    try:
        if arm=='green_strong_full':aux=cie_fft.run_cie_fft(model,Y,np.full(model.n,.1),cfg)
        else:aux=run_auxiliary(model,Y,np.full(model.n,.1),arm,cfg)
    finally:cie_fft.DiskPulseToeplitzGreen=original
    aux_seconds=time.process_time()-aux_cpu
    # Core validation endpoints are already paid. Record original-Phi candidates
    # with their physical validation time inferred from stage callbacks; do not
    # invoke a second forward simply to create a budget snapshot.
    for row in aux['candidates']:
        t=0. if row['stage']==-1 else stage_cpu[row['stage']]['cpu_seconds']
        states.append(dict(cpu_seconds=t,phase='auxiliary',objective=row['objective'],m=row['m'].copy()))
    chosen=min(aux['candidates'],key=lambda r:r['objective'])
    selected=chosen['m'].copy();best=[selected.copy(),chosen['objective']]
    correct_cpu=time.process_time();correct_wall=time.perf_counter()
    # R6 refine_observed numerical body, only toy bounds changed to [0,12].
    # Same original Phi and best paid trajectory policy; no regularizer here.
    history=[]
    def fun(m):
        value,g=model.data_value_gradient(m,Y);paid(m,'corrector',value)
        history.append([value,float(np.max(abs(g)))])
        if value<best[1]:best[:]=[m.copy(),value]
        return 1e6*value,1e6*g
    opt=minimize(fun,selected,method='L-BFGS-B',jac=True,bounds=[(0.,12.)]*model.n,
        options=dict(maxiter=3000,maxfun=12000,ftol=1e-12,gtol=.01,maxls=30))
    corrector=dict(status=str(opt.message),success=bool(opt.success),iterations=int(opt.nit),
        evaluations=int(opt.nfev),cpu_seconds=time.process_time()-correct_cpu,
        wall_seconds=time.perf_counter()-correct_wall,history=history)
    # Best physical current may differ from the final cached optimizer point;
    # this one explicitly charged export forward guarantees matching arrays.
    final=model.state(best[0]);paid(best[0],'export',best[1])
    np.savez_compressed(out/'prediction.npz',m=best[0],j=final['j'],Q=final['Q'],Y_prediction=final['Y'])
    np.savez_compressed(out/'PAID_STATES.npz',m=np.array([s['m'] for s in states]),
        cpu_seconds=np.array([s['cpu_seconds'] for s in states]),
        objective=np.array([s['objective'] for s in states]),phase=np.array([s['phase'] for s in states]))
    dump(out/'AUXILIARY_HISTORY.json',aux['history']);dump(out/'CORRECTOR.json',corrector)
    dump(out/'AUXILIARY_CANDIDATES.json',[{k:v for k,v in r.items() if k!='m'} for r in aux['candidates']])
    receipt=dict(status='SEALED_PUBLIC_PREDICTION',case=case,arm=arm,suite=SUITE,truth_access='NONE',
        material_dof=model.n,material_box=[0.,12.],prior='none; inherited Regularizer with smooth=tv=0',
        auxiliary_iterations=len(aux['history']),auxiliary_status=aux['status'],
        auxiliary_selected_stage=chosen['stage'],auxiliary_CPU_seconds=aux_seconds,
        custom_denominator_setup_CPU_seconds=setup_cpu,corrector= {k:v for k,v in corrector.items() if k!='history'},
        physical_relative_residual=float(np.linalg.norm(final['Y']-Y)/np.linalg.norm(Y)),
        original_Phi=best[1],physical_model_cost=dict(model.cost),
        auxiliary_execution=aux['execution_receipt'],
        limitations=['Domain-adapted mature baseline: 10x10 and box[0,12], not bitwise historical contract',
        'Fourier schedule spatially scaled from 24x24; current coefficient metric and PR rules retained',
        'Exact supplied dense G action replaces FFT geometry reconstruction; no FFT complexity claim',
        'Green strong .1 threshold inherited; its toy rank is not forced to historical19',
        'Legacy state_term_material field for custom Green readout is onehop diagnostic; use material_final_block_loss',
        'Fixed iteration allowances; actual work differs; no equal-cost superiority claim'],
        threads=threadpool_info(),wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu)
    dump(out/'RECEIPT.json',receipt)
    print(json.dumps({k:receipt[k] for k in ('case','arm','status','physical_relative_residual','cpu_seconds')}),flush=True)

def seal():
    rows=[]
    for case in CASES:
        for arm in ARMS:
            path=suite_path('runs')/(case+'__'+arm)
            r=json.loads((path/'RECEIPT.json').read_text())
            if r['status']!='SEALED_PUBLIC_PREDICTION':raise RuntimeError(path)
            rows.append(dict(case=case,arm=arm,prediction=str(path/'prediction.npz'),receipt=str(path/'RECEIPT.json'),
                         objective=r['original_Phi'],cpu_seconds=r['cpu_seconds']))
    out=suite_path('SEALED_COHORT','.json')
    if out.exists():raise FileExistsError(out)
    dump(out,dict(status='ALL20_PUBLIC_PREDICTIONS_SEALED_BEFORE_TRUTH',rows=rows))

def score():
    sealed=json.loads(suite_path('SEALED_COHORT','.json').read_text());rows=[];prefix=[]
    for case in CASES:
        subset=[r for r in sealed['rows'] if r['case']==case]
        common=min(r['cpu_seconds'] for r in subset)
        with np.load(SOURCE/(case+'.npz'),allow_pickle=False) as z:truth=z['truth'].copy()
        for row in subset:
            with np.load(row['prediction'],allow_pickle=False) as z:m=z['m'].copy()
            receipt=json.loads(Path(row['receipt']).read_text())
            rows.append(dict(case=case,arm=row['arm'],relative_material_l2=float(np.linalg.norm(m-truth)/np.linalg.norm(truth)),
                relative_data_l2=receipt['physical_relative_residual'],cpu_seconds=receipt['cpu_seconds'],
                auxiliary_CPU_seconds=receipt['auxiliary_CPU_seconds'],LU=receipt['physical_model_cost'].get('LU'),
                complex_RHS=receipt['physical_model_cost'].get('actual_complex_RHS')))
            with np.load(Path(row['prediction']).parent/'PAID_STATES.npz',allow_pickle=False) as z:
                ids=np.flatnonzero(z['cpu_seconds']<=common)
                ix=ids[np.argmin(z['objective'][ids])];mm=z['m'][ix]
                prefix.append(dict(case=case,arm=row['arm'],common_CPU_seconds=common,paid_state_CPU_seconds=float(z['cpu_seconds'][ix]),
                    phase=str(z['phase'][ix]),relative_material_l2=float(np.linalg.norm(mm-truth)/np.linalg.norm(truth)),
                    relative_data_l2=float(np.sqrt(2*z['objective'][ix]))))
    out=suite_path('OFFLINE_SCORES','.json')
    if out.exists():raise FileExistsError(out)
    dump(out,dict(truth_access='OFFLINE_ONLY_AFTER_ALL20_SEALED',endpoints=rows,
        common_paid_CPU_prefix=prefix,prefix_limitations='Snapshots only include stage endpoints and corrector evaluations; no GT selection. Prefix is observed paid-state comparison, not equal expenditure.'))
    print(json.dumps(rows),flush=True)

def main():
    global SUITE
    ap=argparse.ArgumentParser();ap.add_argument('phase',choices=('prepare','run','seal','score'))
    ap.add_argument('--suite',choices=('scaled','clipped'),default='scaled')
    ap.add_argument('--case',choices=CASES);ap.add_argument('--arm',choices=ARMS);args=ap.parse_args()
    SUITE=args.suite
    with threadpool_limits(limits=1):
        if args.phase=='prepare':prepare()
        elif args.phase=='run':run(args.case,args.arm)
        elif args.phase=='seal':seal()
        else:score()

if __name__=='__main__':main()
