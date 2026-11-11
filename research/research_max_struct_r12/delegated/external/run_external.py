"""Fresh external anchors with inherited R8 physics and original observations."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time,argparse,subprocess,resource
import numpy as np
from scipy.optimize import minimize
HERE=Path(__file__).absolute().parent;REPO=HERE.parents[3];sys.path.insert(0,str(REPO))
R8=REPO.parent/'r8_worktree/research/research_max_struct_r8'
from research.research_max_struct_r8.src.fresnel_adapter import prepare,load_file
from research.research_max_struct_r8.src.fresnel_solver import run as fresnel_run
from research.research_max_struct_r8.src.paper_physics import acquisition,make_model
from research.research_max_struct_r8.src.paper_cie_fft import run_cie_fft
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r12.src.field_current import DenseFieldCurrentProjectors,PhysicalAudit,run_projection
CASES=['dielTM_dec8f','twodielTM_8f','austria_eps4_400MHz'];METHODS=['direct','Dual','AP','DR']
CONTRACT=dict(status='FROZEN_BEFORE_EXTERNAL_INVERSION',cases=CASES,methods=METHODS,bounds=[0,20],Fresnel=dict(side=24,width_m=.16,frequency_GHz=2,initial=.1,illumination='existing outgoing multipole empty-target fit',masks='27training sources x39actualreceiver rows; heldout711channels never used in inversion/selection',delta_relative=.15,delta_status='MODEL_DISCREPANCY_DIAGNOSTIC; no known noise covariance and no CN certificate',aux_stages=[dict(mf=m,beta_material=b,beta_current=1,maxiter=i) for m,b,i in [(3,1,50),(6,.2,50),(10,.2,50),(12,.2,150)]]),Austria=dict(case='austria_eps4_400MHz',inverse_side=32,forward_side=150,initial=0,frequency_Hz=400e6,epsilon_r=4,material='chi=epsilon_r-1, real; unknown all1024cells',delta_relative=.12,noise='Inherited actual R8 observations with10% fixed-Frobenius noise; .12||noisyY|| data ball is a predetermined discrepancy diagnostic, not iid CN99',aux_stages=[dict(mf=m,beta_material=b,beta_current=8,maxiter=100,coefficient_tol=t*t) for m,b,t in [(5,8,.01),(8,2,.001),(10,.5,.001)]]),prior='R8 inherited smooth1e-4/TV1e-4/epsilon.01',projection=dict(scale=6,steps=5000,audit_every=50,omega=1,grouped_ports=True,CPU_cap_seconds=180),auxiliary_CPU_cap_seconds=180,common_corrector=dict(maxiter=1000,maxfun=4000,CPU_cap_seconds=60,ftol=1e-12,gtol=.01,scale=1e6,maxls=30),selection='minimum original training physical Phi including same prior; all paid line-search/accepted/stage/audit candidates',heldout_selection=False,no_raw_original_measurement_redistribution=True,truth_access='NONE during inversion/selection; Austria scorer only after all12 predictions sealed')
def dump(p,x):
    def cv(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        raise TypeError(type(v).__name__)
    p.write_text(json.dumps(x,indent=2,default=cv,allow_nan=False)+'\n')
class BudgetReached(BaseException):pass
def inputs(cid):
    if cid.startswith('austria'):
        p=R8/'RESULTS/paper_benchmarks'/cid/'PUBLIC_DATA.npz';z=np.load(p);public=acquisition(int(z['inverse_side']),float(z['frequency']));m=make_model(public);return m,z['Y'],dict(public=public,source_path=str(p),initial=0,delta_relative=.12)
    p=R8/'data_local/fresnel'/(cid+'.exp');d=prepare(p,frequency=2.,side=24,width=.16,illumination='multipole',matrix_free=False);return d['model'],d['Y_train'],dict(**d,source_path=str(p),initial=.1,delta_relative=.15)
def solve(cid,method):
    if cid not in CASES or method not in METHODS:raise ValueError('not frozen')
    out=HERE/'RESULTS'/cid/method;out.mkdir(parents=True,exist_ok=False);cpu=time.process_time();wall=time.perf_counter();m,Y,d=inputs(cid)
    cfg=json.loads((HERE/'CONTRACT.json').read_text());fresnel=not cid.startswith('austria');best=np.full(m.n,d['initial']);bestphi=m.physical_score(best,Y);history=[];candidates=[]
    original=m.physical_score
    def observe(x,label):
        nonlocal best,bestphi
        value=float(original(x,Y));history.append(dict(label=label,Phi=value,CPU_seconds=time.process_time()-cpu));candidates.append(np.array(x,copy=True))
        if value<bestphi:bestphi=value;best=x.copy()
        return value
    def wrapped(x,yy):return observe(x,'inherited_stage')
    m.physical_score=wrapped
    auxstart=time.process_time();auxreceipt={};projection=None
    def budget(limit,start):
        if time.process_time()-start>=limit:raise BudgetReached()
    def physical_opt(start,seconds,maxiter,label):
        nonlocal best,bestphi
        began=time.process_time()
        def fg(x):
            nonlocal best,bestphi
            budget(seconds,began);v,g=m.data_value_gradient(x,Y);history.append(dict(label=label,Phi=float(v),CPU_seconds=time.process_time()-cpu));candidates.append(x.copy())
            if v<bestphi:best=x.copy();bestphi=float(v)
            return v*1e6,g*1e6
        try:
            rr=minimize(fg,start,jac=True,method='L-BFGS-B',bounds=[(0,20)]*m.n,options=dict(maxiter=maxiter,maxfun=4*maxiter,ftol=1e-12,gtol=.01,maxls=30))
            return dict(success=bool(rr.success),nit=int(rr.nit),nfev=int(rr.nfev),message=str(rr.message),CPU_seconds=time.process_time()-began)
        except BudgetReached:return dict(success=False,message='PREDECLARED_CPU_CAP',CPU_seconds=time.process_time()-began)
    if method=='direct':auxreceipt=physical_opt(best,180,2000,'direct_auxiliary')
    elif method=='Dual':
        def reg(x):budget(180,auxstart);return m.reg.evaluate(x)[:2]
        try:
            if fresnel:
                stages=cfg['Fresnel']['aux_stages'];config=dict(host=.1,material_box=[0,20],dpcs_rtol=.001,coefficient_tol=0,miniter=50,material_inner_maxiter=60,material_inner_ftol=1e-12,material_inner_gtol=1e-8,regularizer=reg,stages=stages)
                aux=fresnel_run(m,Y,best,'cie_dual',config)
            else:
                public=d['public'];config=dict(arm='cie_dual',stages=cfg['Austria']['aux_stages'],dpcs_rank=15,material_box=(0,20),coefficient_tol=0,miniter=2,sourcewise_cg=True,regularizer=reg,fft_geometry=dict(coords=public['centers'],radius=public['radius'],k=public['k']),material_inner_maxiter=60,cache_material_value_gradient=True,physical_sigma=float(np.linalg.norm(Y)))
                aux=run_cie_fft(m,Y,best,config)
            auxreceipt=aux['execution_receipt'];observe(aux['m'],'auxiliary_selected')
        except BudgetReached:auxreceipt=dict(status='PREDECLARED_AUX_CPU_CAP',CPU_seconds=time.process_time()-auxstart,selection='Retain earlier paid physical stage candidates; partial inner readout not promoted')
    else:
        ports=list(np.moveaxis(m.Ss,2,0)) if fresnel else m.S
        projector=DenseFieldCurrentProjectors(m.G,m.E,ports,Y,d['delta_relative']*np.linalg.norm(Y),scale=6,lower=0,upper=20,grouped_ports=True,public_metadata=dict(noise_certificate=False,scope='external discrepancy diagnostic'))
        def audit(x):
            before=dict(m.cost);v=observe(x,'projection_audit');st=m.state(x)
            return PhysicalAudit(v,st['j'],st['Q'],costs={k:v-before.get(k,0) for k,v in m.cost.items()})
        projection=run_projection(method.lower(),projector,best,audit,steps=5000,audit_every=50,cpu_budget_seconds=max(0,180-(time.process_time()-auxstart)))
        projection.save(out/'projection');best=projection.material.copy();bestphi=projection.objective
        auxreceipt=dict(status=projection.stop_reason,iterations=projection.iterations,CPU_seconds=projection.cpu_seconds,setup_CPU_seconds=projection.setup_cpu_seconds,counts=projection.costs,setup_counts=projection.history[0].costs,delta=projector.delta)
    raw=best.copy();rawphi=bestphi;corrector=physical_opt(best,60,1000,'common_corrector');st=m.state(best)
    if fresnel:
        full=MaxwellModel(d['public']);prediction=full.forward(best);extra=dict(full.cost)
        np.savez_compressed(out/'derived.npz',m=best,raw=raw,Y_prediction=prediction,train_mask=d['train_mask'],withheld_mask=d['withheld_mask'],held_sources=d['held_sources'],centers=d['centers'],candidates=np.asarray(candidates))
    else:
        extra={};np.savez_compressed(out/'derived.npz',m=best,raw=raw,Y_prediction=st['Y'],centers=d['public']['centers'],candidates=np.asarray(candidates))
    dump(out/'RECEIPT.json',dict(status='COMPLETE',case_id=cid,method=method,source_path=d['source_path'],original_Phi=bestphi,raw_Phi=rawphi,train_relative_complex_residual=float(np.linalg.norm(st['Y']-Y)/np.linalg.norm(Y)),CPU_seconds=time.process_time()-cpu,process_CPU_seconds_near_exit=time.process_time(),wall_seconds=time.perf_counter()-wall,physical_counts=dict(m.cost),full_prediction_counts=extra,auxiliary=auxreceipt,corrector=corrector,history=history,calibration=d.get('calibration'),training_channels=int(Y.size),heldout_selection=False,truth_access='NONE',raw_original_saved=False))
    print(json.dumps(dict(case=cid,method=method,Phi=bestphi,CPU=time.process_time()-cpu)),flush=True)
def queue():
    if not (HERE/'CONTRACT.json').exists():dump(HERE/'CONTRACT.json',CONTRACT)
    rows=[]
    for cid in CASES:
        for method in METHODS:
            out=HERE/'RESULTS'/cid/method
            if (out/'RECEIPT.json').exists():continue
            before=resource.getrusage(resource.RUSAGE_CHILDREN);t=time.time();p=subprocess.run([sys.executable,__file__,'solve','--case',cid,'--method',method],capture_output=True,text=True);after=resource.getrusage(resource.RUSAGE_CHILDREN)
            rows.append(dict(case_id=cid,method=method,returncode=p.returncode,external_child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,wall_seconds=time.time()-t,stdout=p.stdout,stderr=p.stderr));dump(HERE/'QUEUE.json',dict(status='RUNNING' if not p.returncode else 'FAILED',rows=rows));print(p.stdout or p.stderr,flush=True)
            if p.returncode:sys.exit(p.returncode)
    dump(HERE/'QUEUE.json',dict(status='COMPLETE',rows=rows))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['queue','solve']);p.add_argument('--case');p.add_argument('--method');a=p.parse_args()
    if a.mode=='queue':queue()
    else:solve(a.case,a.method)
