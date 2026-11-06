"""Root-authorized five-task amendment. No original24 overwrite or scoring.
Run `all` only after the original24 queue/children exit. Same public-only worker.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[key]='1'
import time
WALL,CPU=time.perf_counter(),time.process_time()
from pathlib import Path
import sys,json,subprocess,gc,argparse
ROOT=Path(__file__).resolve().parents[2];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits
from research.research_max_struct_r9.experiments import vector_extension as base
from research.research_max_struct_r8.src.vector_cie3d import box_fraction
CONFIG=ROOT/'configs/vector_extension_resolved_hollow.json';CASE='hollow_cube_resolved_control'
LOGS=ROOT/'delegated/3d_extension/resolved_hollow_queue'

def prepare(cfg):
    public=base.PUBLIC/CASE;private=base.TRUTH/CASE;public.mkdir(exist_ok=False);private.mkdir(exist_ok=False)
    contract=base.load(base.TRUTH/'RESOLVED_HOLLOW_CONTRACT.json');t,c=time.perf_counter(),time.process_time();states=[]
    for nside in (10,8):
        model=base.StrictVectorInverse(nside,.6,cfg,auxiliary=False)
        x=(contract['epsilon_r']-1)*(box_fraction(model.points,model.h,contract['center'],[contract['outer_side']/2]*3)-box_fraction(model.points,model.h,contract['center'],[contract['inner_side']/2]*3))
        st=model.state(x);states.append(dict(chi=x.copy(),Y=st['Y'].copy(),counts=dict(model.counts),setup=model.setup,geometry=base.fieldmeta(model)))
        del model,st;gc.collect()
    fine,coarse=states;clean=fine['Y'];sigma=.03*la.norm(clean)/np.sqrt(clean.size)
    rng=np.random.default_rng(contract['noise_seed']);noise=sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape));Y=clean+noise
    np.savez_compressed(public/'observations.npz',Y=Y)
    base.dump(public/'PUBLIC.json',dict(case_id=CASE,noise_model='iid_circular_complex',sigma=float(sigma),sigma_convention='E|n|²=sigma² commonrawcomplexobservation',expected_relative_Frobenius_noise=.03,acquisition=coarse['geometry'],noise_rule_scope='noise only; discretization mismatch not included'))
    np.savez_compressed(private/'truth.npz',chi=coarse['chi'],Y_clean=clean)
    row=dict(status='GENERATED',case_id=CASE,crossgrid_clean_relative_difference=float(la.norm(coarse['Y']-clean)/la.norm(clean)),
        actual_new_iid_noise_ratio=float(la.norm(noise)/la.norm(clean)),simulation_counts=dict(fine=fine['counts'],coarse=coarse['counts']),
        operator_setup=dict(fine=fine['setup'],coarse=coarse['setup']),CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,
        scope='root-authorized geometry resolution control; not a pure grid-refinement experiment; no original24 modification')
    base.dump(base.OUT/'generation'/f'{CASE}.json',row);print(json.dumps(row),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('command',choices=['all','worker']);p.add_argument('--arm');args=p.parse_args();cfg=base.load(CONFIG)
    base.CONFIG=CONFIG;base.ENTRY_CPU=CPU;base.ENTRY_WALL=WALL
    if args.command=='worker':return base.worker(cfg,CASE,args.arm)
    LOGS.mkdir(exist_ok=False);prepare(cfg);pending=list(cfg['arms']);active={};jobs=[];start=time.perf_counter()
    while pending or active:
        while pending and len(active)<2:
            arm=pending.pop(0);path=LOGS/f'{arm}.log';fh=path.open('x');command=[sys.executable,str(Path(__file__).resolve()),'worker','--arm',arm]
            child=subprocess.Popen(command,cwd=REPO,stdout=fh,stderr=subprocess.STDOUT);row=dict(arm=arm,case_id=CASE,pid=child.pid,command=command,started_utc=base.utc(),log=str(path.relative_to(REPO)))
            active[child]=(fh,time.perf_counter(),row);print(json.dumps(dict(event='START',**row)),flush=True)
        for child,(fh,t,row) in list(active.items()):
            code=child.poll()
            if code is None:continue
            fh.close();row.update(exit_code=code,ended_utc=base.utc(),wall_seconds=time.perf_counter()-t)
            receipt=base.load(base.OUT/CASE/row['arm']/'RECEIPT.json');row.update(status=receipt['status'],error=receipt.get('error'))
            jobs.append(row);del active[child];print(json.dumps(dict(event='END',**row)),flush=True)
        base.dump(LOGS/'RECEIPT.json',dict(status='RUNNING' if active or pending else 'COMPLETE',jobs=jobs,running=len(active),pending=len(pending),max_processes=2,wall_seconds=time.perf_counter()-start,truth_access='NONE',scoring='NOT_RUN'))
        if active:time.sleep(.5)
    base.dump(LOGS/'PARENT_COST.json',dict(parent_CPU_seconds=time.process_time()-CPU,parent_wall_seconds=time.perf_counter()-WALL,
        scope='includes imports/new generation/queue; children reported separately; no GT scoring'))
    return 0 if all(j['exit_code']==0 for j in jobs) else 1

if __name__=='__main__':
    with threadpool_limits(limits=1):raise SystemExit(main())
