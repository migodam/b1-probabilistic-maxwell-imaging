"""Postscore full native-Dual budget control; externalv1 endpoints immutable."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,copy,time,json,subprocess,resource,argparse
HERE=Path(__file__).absolute().parent;sys.path.insert(0,str(HERE.parent/'external'))
import run_external as base
base.HERE=HERE;base.METHODS=['Dual']
cfg=copy.deepcopy(base.CONTRACT);cfg['status']='FROZEN_POST_V1_SCORE_BEFORE_NEW_RUNS';cfg['methods']=['Dual'];cfg['reason']='Root requested full native baseline budget control after v1 scoring; exploratory postscore budget diagnostic, not new confirmation'
for st,it in zip(cfg['Fresnel']['aux_stages'],[500,500,500,1500]):st['maxiter']=it
for st in cfg['Austria']['aux_stages']:st['maxiter']=300
cfg['schedule_provenance']=dict(Fresnel='research/research_max_struct_r8/experiments/fresnel_invert.py stages full:500/500/500/1500 mf3/6/10/12; only lower bound0 air correction retained',Austria='research/research_max_struct_r8/experiments/paper_benchmark.py solve(): case mf5/8/10 beta8/2/.5 currentbeta8; parser default iterations300 per stage; DPCS15 literal inherited',caps='same180proposal+60correctorCPU as externalv1; no rank/radius/box/prior/start changes',selection='same paid original physical Phi, no GT/heldout tuning',no_raw_measurement_redistribution=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case');a=p.parse_args()
    if a.case:base.solve(a.case,'Dual');sys.exit(0)
    path=HERE/'CONTRACT.json'
    if path.exists():raise FileExistsError('Do not overwrite a completed budget control')
    base.dump(path,cfg);records=[]
    for cid in base.CASES:
        before=resource.getrusage(resource.RUSAGE_CHILDREN);t=time.time();p=subprocess.run([sys.executable,__file__,'--case',cid],capture_output=True,text=True);after=resource.getrusage(resource.RUSAGE_CHILDREN)
        records.append(dict(case_id=cid,method='Dual',returncode=p.returncode,external_child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,wall_seconds=time.time()-t,stdout=p.stdout,stderr=p.stderr));base.dump(HERE/'QUEUE.json',dict(status='RUNNING' if not p.returncode else 'FAILED',rows=records));print(p.stdout or p.stderr,flush=True)
        if p.returncode:sys.exit(p.returncode)
    base.dump(HERE/'QUEUE.json',dict(status='COMPLETE',rows=records))
