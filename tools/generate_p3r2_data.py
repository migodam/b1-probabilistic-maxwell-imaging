"""Dedicated joint simulator, independent of online reconstruction/scoring."""
from pathlib import Path
import argparse,sys,json,gc,os,subprocess,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT,Book,Target,save,write

def main():
    ap=argparse.ArgumentParser();ap.add_argument('group',choices=['train','validation','locked']);ap.add_argument('--count',type=int,required=True);ap.add_argument('--start',type=int,default=0);a=ap.parse_args()
    if a.group=='locked' and not (OUT/'METHODS_FROZEN_FOR_LOCKED.json').exists():raise RuntimeError('freeze before new locked observations')
    stage='L' if a.group=='locked' else 'A';job=f'data_{a.group}_{a.start}_{a.count}'
    b=Book(stage,job,max_wall=3600);status='COMPLETE';error=None;start=time.perf_counter()
    try:
        t=Target(b);ys=[];chis=[];qs=[];ids=[];base={'train':0,'validation':1,'locked':2}[a.group]
        for k in range(a.start,a.start+a.count):
            parent=f'R2_{a.group}_{k:05d}';dest=OUT/f'inputs/online/{parent}.npz'
            if dest.exists():raise RuntimeError('no overwrite or implicit campaign replay')
            rng=np.random.default_rng(np.random.SeedSequence([2026102202,base,k,0]));q=t.prior.draw(rng)[0]
            nr=np.random.default_rng(np.random.SeedSequence([2026102202,base,k,1]));noise=t.sigma/np.sqrt(2)*(nr.normal(size=(32,8))+1j*nr.normal(size=(32,8)))
            p=t.evaluate(q,'TRUTH_GENERATION_'+parent);y=p.state.data+noise
            save(OUT/f'inputs/offline_only/{parent}.npz',q=q,chi=p.state.chi,prediction=p.state.data,noise=noise)
            save(dest,observation=y,target_id=np.array(t.cfg['target_id']),sigma=np.array(t.sigma))
            ys.append(y);chis.append(p.state.chi);qs.append(q);ids.append(parent);del p
            if (k-a.start+1)%256==0:
                t.backend()._transfers.clear();gc.collect();rss=int(subprocess.check_output(['ps','-o','rss=','-p',str(os.getpid())],text=True).strip())/1024
                print(json.dumps(dict(group=a.group,completed=k-a.start+1,RHS=8*(k-a.start+1),wall=time.perf_counter()-start,RSS_MiB=rss)),flush=True)
                if rss>2000:raise MemoryError('bounded simulator RSS')
        # Joint labels accessible only to supervised training/offline score.
        dest=OUT/f'{stage}/data/{a.group}_{a.start}_{a.count}.npz'
        save(dest,observation=np.stack(ys),chi=np.stack(chis),q=np.stack(qs),parent=np.array(ids))
        if a.start==0:save(OUT/f'{stage}/data/{a.group}.npz',observation=np.stack(ys),chi=np.stack(chis),q=np.stack(qs),parent=np.array(ids))
        write(b.jobdir/'DATA_MANIFEST.json',dict(group=a.group,parents=ids,seed=[2026102202,base],GT_scope='generation and training labels only',packed_path=str(dest),dedicated_process=True,elapsed=time.perf_counter()-start))
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:b.close(status,error)

if __name__=='__main__':main()
