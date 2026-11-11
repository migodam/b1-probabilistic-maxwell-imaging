"""Authorized matched32 Austria diagnostic; reuse one offlineGT clean prediction."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,copy,json,time,subprocess,resource,argparse
import numpy as np
HERE=Path(__file__).absolute().parent;sys.path.insert(0,str(HERE.parent/'external'))
import run_external as base
base.HERE=HERE;base.CASES=['austria_matched32'];base.METHODS=['direct','Dual','DR'];cid=base.CASES[0]
def matched_inputs(case):
    z=np.load(HERE/'data_public/OBSERVATIONS.npz');public=base.acquisition(32,400e6);m=base.make_model(public)
    return m,z['Y'],dict(public=public,source_path=str(HERE/'data_public/OBSERVATIONS.npz'),initial=0,delta_relative=.12)
base.inputs=matched_inputs
def generate():
    if (HERE/'CONTRACT.json').exists():raise FileExistsError('Do not repeat matched diagnostic')
    cfg=copy.deepcopy(base.CONTRACT);cfg.update(status='OFFLINE_PHYSICS_DIAGNOSTIC_FROZEN_BEFORE_MATCHED_RUNS',cases=base.CASES,methods=base.METHODS,qualification='Matched32 diagnostic to isolate original150-to32 model discrepancy, not originalpaper reproduction/freshconfirmation',generation_seed=2026111232,generation_source='Reuse single authorized offlineGT32 forward saved external_physics/OFFLINE_GT_COARSE_PREDICTION.npz; no redundant physical generation solve',noise_relative_fixed_norm=.1)
    cfg['Austria'].update(case=cid,inverse_side=32,forward_side=32,epsilon_r=4,noise='Fresh complex Gaussian direction, normalized to exact.1 normclean; fixedseed2026111232')
    for st in cfg['Austria']['aux_stages']:st['maxiter']=300
    HERE.mkdir(exist_ok=True);base.dump(HERE/'CONTRACT.json',cfg);(HERE/'data_public').mkdir();(HERE/'sealed_truth').mkdir()
    cpu=time.process_time();cleanpath=HERE.parent/'external_physics/OFFLINE_GT_COARSE_PREDICTION.npz';z=np.load(cleanpath);clean=z['Y_prediction'];truth=z['m_GT'];rng=np.random.default_rng(2026111232);noise=rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape);noise*=.1*np.linalg.norm(clean)/np.linalg.norm(noise);Y=clean+noise
    np.savez_compressed(HERE/'data_public/OBSERVATIONS.npz',Y=Y,frequency=400e6,inverse_side=32,forward_side=32)
    np.savez_compressed(HERE/'sealed_truth/GT.npz',m=truth,cleanY=clean)
    base.dump(HERE/'GENERATION_RECEIPT.json',dict(status='COMPLETE_MATCHED_OFFLINE_PHYSICS_DIAGNOSTIC',seed=2026111232,reused_clean_forward=str(cleanpath),new_physical_forward_count=0,new_LU=0,new_source_RHS=0,prior_authorized_clean_cost='external_physics/RESULTS.json:1LU20sourceRHS, accounted once',relative_noise_to_clean=float(np.linalg.norm(noise)/np.linalg.norm(clean)),coarse_GT_observed_relative_residual=float(np.linalg.norm(noise)/np.linalg.norm(Y)),CPU_seconds=time.process_time()-cpu,truth_access='GENERATOR_ONLY; solver reads publicY only'))
def queue():
    generate();rows=[]
    for method in base.METHODS:
        before=resource.getrusage(resource.RUSAGE_CHILDREN);t=time.time();p=subprocess.run([sys.executable,__file__,'--method',method],capture_output=True,text=True);after=resource.getrusage(resource.RUSAGE_CHILDREN)
        rows.append(dict(case_id=cid,method=method,returncode=p.returncode,external_child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,wall_seconds=time.time()-t,stdout=p.stdout,stderr=p.stderr));base.dump(HERE/'QUEUE.json',dict(status='RUNNING'if not p.returncode else'FAILED',rows=rows));print(p.stdout or p.stderr,flush=True)
        if p.returncode:sys.exit(p.returncode)
    base.dump(HERE/'QUEUE.json',dict(status='COMPLETE',rows=rows))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--method');a=p.parse_args()
    if a.method:base.solve(cid,a.method)
    else:queue()
