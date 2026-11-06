"""Run complete physical NIE Austria contracts with declared grid adaptations.

Build and solve are separate invocations; solve never loads truth. Offline score
is a third invocation after every material/prediction has been written/sealed.
"""
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):os.environ[_k]='1'
import argparse,json,time,sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from research.research_max_struct_r8.src.paper_physics import acquisition,austria,fft_forward,make_model
from research.research_max_struct_r8.src.paper_cie_fft import run_cie_fft
from research.research_max_struct_r8.src.paper_nie import run_paper_tsom
from research.research_max_struct_r7.src.aplus_readout import mobius_material_readout
from research.research_max_struct_r8.src.paper_scores import material_scores

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/paper_benchmarks';ARMS=('cie_single','cie_dual','global_dual','a_material','a_plus','a_plus_sourcewise','paper_nie','paper_lsie','physical_lbfgs')
CASES=json.loads((ROOT/'delegated/paper_baselines/PAPER_CONTRACTS.json').read_text())['NIE2016']['cases']
def write(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    def default(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        raise TypeError(type(v).__name__)
    path.write_text(json.dumps(obj,indent=2,default=default,allow_nan=False))

def build(args):
    OUT.mkdir(parents=True,exist_ok=True);private=ROOT/'data/paper_benchmarks_private';private.mkdir(parents=True,exist_ok=True)
    for i,case in enumerate(CASES):
        if args.case and case['id']!=args.case:continue
        folder=OUT/case['id'];folder.mkdir(parents=True,exist_ok=True)
        t,c=time.perf_counter(),time.process_time();fwd=acquisition(args.forward,case['frequency_Hz'],dense=False)
        truth=austria(fwd['centers'],case['epsilon_r']);Y,receipt=fft_forward(fwd,truth)
        rng=np.random.default_rng(8600+i);noise=rng.normal(size=Y.shape)+1j*rng.normal(size=Y.shape)
        noise*=.1*np.linalg.norm(Y)/np.linalg.norm(noise);Y+=noise
        np.savez_compressed(folder/'PUBLIC_DATA.npz',Y=Y,frequency=case['frequency_Hz'],inverse_side=args.inverse,forward_side=args.forward)
        invcoords=acquisition(args.inverse,case['frequency_Hz'],dense=False)['centers']
        np.savez_compressed(private/(case['id']+'.npz'),m=austria(invcoords,case['epsilon_r']))
        write(folder/'DATA_RECEIPT.json',dict(status='COMPLETED',case=case,forward=receipt,forward_side=args.forward,inverse_side=args.inverse,noise_relative=.1,pulse='equal-area disk h/sqrt(pi), analytic disk self-term',wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,qualification='paper physical geometry and acquisition, adapted grids and disk-pulse quadrature; not exact published reproduction'))
        print('BUILD',case['id'],receipt['wall_seconds'],flush=True)

def solve(args):
    started=time.perf_counter();cpu=time.process_time()
    for case in CASES:
        if args.case and case['id']!=args.case:continue
        folder=OUT/case['id'];data=dict(np.load(folder/'PUBLIC_DATA.npz'));side=int(data['inverse_side']);Y=data['Y']
        assembly_t,assembly_c=time.perf_counter(),time.process_time();public=acquisition(side,case['frequency_Hz']);setup=dict(wall_seconds=time.perf_counter()-assembly_t,cpu_seconds=time.process_time()-assembly_c)
        for arm in ARMS:
            if args.arm and arm!=args.arm:continue
            dest=folder/arm
            if (dest/'PREDICTION_SEAL.json').exists():continue
            dest.mkdir(parents=True,exist_ok=True);wall,c=time.perf_counter(),time.process_time()
            model=make_model(public,prior=not arm.startswith('paper_'));model.paper_geometry=dict(coords=public['centers'],radius=public['radius'],k=public['k'])
            stages=[dict(mf=mf,beta_material=beta,beta_current=beta if arm=='cie_single' else case['beta'][0],maxiter=args.iterations,paper_delta=.01 if j==0 else .001,coefficient_tol=(.01 if j==0 else .001)**2) for j,(mf,beta) in enumerate(zip(case['mf'],case['beta']))]
            m0=np.zeros(model.n);history=[];aux={};candidate=[]
            try:
                if arm.startswith('paper_'):
                    result=run_paper_tsom(model,Y,m0,stages,lsie=arm=='paper_lsie',maxiter=args.iterations)
                    candidate=[s['m'] for s in result['stage_states']];history=result['history'];aux=result['receipt']
                elif arm!='physical_lbfgs':
                    actual='a_plus' if arm=='a_plus_sourcewise' else arm
                    cfg=dict(arm=actual,stages=stages,dpcs_rank=15,material_box=(0,20),coefficient_tol=0.,miniter=2,
                        sourcewise_cg=arm in ('cie_single','cie_dual','a_plus_sourcewise'),regularizer=lambda z:model.reg.evaluate(z)[:2],
                        fft_geometry=model.paper_geometry,material_inner_maxiter=60,cache_material_value_gradient=True,
                        physical_sigma=float(np.linalg.norm(Y)),aplus_regularizer=model.reg,aplus_sweeps=2,
                        aplus_arc_candidates=True)
                    if actual=='a_plus':cfg['material_readout']=mobius_material_readout
                    result=run_cie_fft(model,Y,m0,cfg);candidate=[s['m'] for s in result['stage_states']];history=result['history'];aux=result['execution_receipt']
                else:candidate=[m0]
                # Candidate selection by original data+prior only; no truth loaded.
                values=[model.physical_score(m,Y) for m in candidate];chosen=int(np.argmin(values));m=candidate[chosen]
                auxiliary_m=m.copy();auxiliary_cost=dict(model.cost);corrector=[]
                def fg(z):
                    v,g=model.data_value_gradient(z,Y);corrector.append(dict(objective=v,LU=model.cost.get('LU',0),cpu_seconds=time.process_time()-c));return v,g
                opt=minimize(fg,m,jac=True,method='L-BFGS-B',bounds=[(0,20)]*model.n,options=dict(maxiter=args.corrector,maxfun=4*args.corrector,ftol=1e-12,gtol=1e-8,maxls=30))
                final=opt.x if model.physical_score(opt.x,Y)<=model.physical_score(m,Y) else m
                state=model.state(final);np.savez_compressed(dest/'PREDICTION.npz',m=final,Y_prediction=state['Y'],m_auxiliary=auxiliary_m)
                write(dest/'HISTORY.json',dict(auxiliary=history,corrector=corrector))
                receipt=dict(status='COMPLETED',arm=arm,auxiliary=aux,setup_shared=setup,auxiliary_physical_cost=auxiliary_cost,physical_cost=dict(model.cost),wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-c,original_Phi=model.physical_score(final,Y),relative_data_residual=float(np.linalg.norm(state['Y']-Y)/np.linalg.norm(Y)),corrector=dict(iterations=int(opt.nit),evaluations=int(opt.nfev),stop=str(opt.message)),candidate_selection='original physical data plus declared prior, all stage endpoints',selected_stage=chosen,all_material_pixels_free=True,GT_loaded=False)
                write(dest/'RECEIPT.json',receipt);write(dest/'PREDICTION_SEAL.json',dict(status='SEALED_BEFORE_OFFLINE_SCORING',prediction_file='PREDICTION.npz',time_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),GT_loaded=False))
                print('SOLVED',case['id'],arm,round(receipt['wall_seconds'],2),receipt['relative_data_residual'],flush=True)
            except Exception as e:
                write(dest/'FAILURE.json',dict(status='FAILED',error=repr(e),physical_cost=dict(model.cost),wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-c));raise
    write(OUT/'RUN_TOTAL.json',dict(wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu))

def score(args):
    rows=[]
    for case in CASES:
        folder=OUT/case['id']
        for arm in ARMS:
            dest=folder/arm
            if not (dest/'PREDICTION_SEAL.json').exists():continue
            p=np.load(dest/'PREDICTION.npz');m=np.load(ROOT/'data/paper_benchmarks_private'/(case['id']+'.npz'))['m'];r=json.loads((dest/'RECEIPT.json').read_text());est=p['m'];support=m>0;pred=est>.5*(case['epsilon_r']-1)
            rows.append(dict(case=case['id'],arm=arm,NRMSE=float(np.linalg.norm(est-m)/np.linalg.norm(m)),auxiliary_NRMSE=float(np.linalg.norm(p['m_auxiliary']-m)/np.linalg.norm(m)),IoU=float(np.sum(support&pred)/np.sum(support|pred)),relative_data_residual=r['relative_data_residual'],CPU_seconds=r['cpu_seconds'],wall_seconds=r['wall_seconds'],LU=r['physical_cost'].get('LU',0),complex_RHS=r['physical_cost'].get('actual_complex_RHS',0),**material_scores(est,m,case['epsilon_r']-1)))
    write(OUT/'OFFLINE_SCORES.json',rows);print(json.dumps(rows,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['build','solve','score']);ap.add_argument('--case');ap.add_argument('--arm');ap.add_argument('--forward',type=int,default=64);ap.add_argument('--inverse',type=int,default=32);ap.add_argument('--iterations',type=int,default=300);ap.add_argument('--corrector',type=int,default=80);args=ap.parse_args();globals()[args.mode](args)
