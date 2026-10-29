"""Explicit offline known-wall diagnostic. Requires legitimate method endpoints sealed."""
from pathlib import Path
import argparse,json,time,sys,resource
import os
for _thread_key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS']:os.environ[_thread_key]='1'
import numpy as np
from ordinary_controls import ROOT,target_control
from r3.physics import JointModel,default_config
from r3.accounting import BudgetBook,plain

def main():
    p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--repeat',type=int,default=0);p.add_argument('--legitimate-output',required=True,type=Path,action='append');p.add_argument('--out',required=True,type=Path);a=p.parse_args()
    for out in a.legitimate_output:
        if not (out/'ONLINE_SEALED.json').exists():raise RuntimeError('Legitimate outputs must seal before oracle truth read')
    book=BudgetBook(root=ROOT);start=time.process_time();startwall=time.time();before=book.snapshot();cfg=default_config();cfg.update(phase='R3_ORACLE');model=JointModel(cfg,book)
    # Only wall coefficients are used. True target values are never an initializer.
    wallfile=np.load(ROOT/'inputs/private'/a.case/'oracle_wall.npz');wall=np.r_[wallfile['eps'],wallfile['loss']];base=np.asarray(model.spec['initial'],float);base[:12]=wall
    online=ROOT/'inputs/online'/f'{a.case}_r{a.repeat}';y=np.load(online/'observation.npz')['y'];model.observation=y;a.out.mkdir(parents=True,exist_ok=True)
    if (a.out/'endpoint.npz').exists() or (a.out/'TRACE.jsonl').exists():raise FileExistsError('Control output already contains scientific artifacts')
    summary=target_control(model,y,a.out,base,'KNOWN_WALL_ORACLE');summary['diagnostic_only']=True;summary['truth_access']='wall coefficients only, after legitimate endpoints sealed'
    sealpath=a.out/'ONLINE_SEALED.json';seal=json.loads(sealpath.read_text());seal.update(truth_access='private wall coefficients only; no target truth',diagnostic_only=True,summary=summary);sealpath.write_text(json.dumps(seal,default=plain,indent=2)+'\n')
    (a.out/'ORACLE_SCOPE.json').write_text(json.dumps(summary,default=plain,indent=2)+'\n')
    after=book.snapshot();receipt=dict(method='KNOWN_WALL_ORACLE',pid=os.getpid(),source_root=str(ROOT/'frozen/src' if (ROOT/'frozen/src').exists() else ROOT/'src'),peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),process_cpu_seconds=time.process_time()-start,wall_seconds=time.time()-startwall,physical_cost={k:after[k]-before[k] for k in before});(a.out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');book.event('R3_ORACLE','control_job_completion',**{k:v for k,v in receipt.items() if k!='pid'})
if __name__=='__main__':main()
