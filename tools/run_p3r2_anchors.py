"""Truth-free native/probit physics contexts, separate from packed labels."""
from pathlib import Path
import sys,argparse,gc,json,subprocess,os
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT,Book,Target,save,write
from p3r2.inversion import reconstruct
from p3r2.anchors import probit_optimize

def main():
    ap=argparse.ArgumentParser();ap.add_argument('group',choices=['train','validation','locked']);ap.add_argument('--count',type=int,required=True);ap.add_argument('--start',type=int,default=0);ap.add_argument('--native',type=int,default=48);ap.add_argument('--probit',type=int,default=48);a=ap.parse_args()
    if a.group=='locked' and not (OUT/'METHODS_FROZEN_FOR_LOCKED.json').exists():raise RuntimeError('locked method freeze required')
    stage='L' if a.group=='locked' else 'A';job=f'anchors_{a.group}_{a.start}_{a.count}_{a.native}_{a.probit}'
    b=Book(stage,job,max_wall=14400);status='COMPLETE';error=None
    try:
        t=Target(b,keep_fields=True)
        for k in range(a.start,a.start+a.count):
            parent=f'R2_{a.group}_{k:05d}';dest=OUT/f'{stage}/contexts/{parent}.npz'
            if dest.exists():raise FileExistsError('context overwrite prohibited')
            y=t.observation(parent);folder=OUT/f'{stage}/contexts/details/{parent}';folder.mkdir(parents=True,exist_ok=True)
            native,report=reconstruct(t,y,parent,method='full_TRF',max_evaluations=a.native)
            save(folder/'NATIVE.npz',q=native.q,chi=native.state.chi,prediction=native.state.data)
            p,x,H,stationarity=probit_optimize(t,y,parent,folder/'probit',native,cap=a.probit)
            J=t.jacobian(p)
            save(dest,observation=y,q0=p.q,chi0=p.state.chi,J=J,prediction=p.state.data,residual=y-p.state.data,phi=t.phi(p,y),x=x,H=H)
            t.persist_cache(p,OUT/f'{stage}/contexts/cache/{parent}.npz')
            write(dest.with_suffix('.json'),dict(parent=parent,native_cap=a.native,probit_cap=a.probit,Phi=t.phi(p,y),native_stationarity=report['kkt'],probit_stationarity=stationarity,GT_read=False,online_target=t.cfg['target_id']))
            del p,native,J;t.backend()._transfers.clear();gc.collect()
            rss=int(subprocess.check_output(['ps','-o','rss=','-p',str(os.getpid())],text=True).strip())/1024
            if rss>2000:raise MemoryError('bounded online RSS')
            print(json.dumps(dict(parent=parent,completed=True,RSS_MiB=rss)),flush=True)
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:b.close(status,error)
if __name__=='__main__':main()
