"""Online-only paired anchor campaign. No offline truth path or generator import."""
import sys,json,time,gc,argparse
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from shared.physics import OUT,Book,Target,write,SIGMA
from track_c.anchors import run,CONFIGS

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--group',default='CDEV');ap.add_argument('--limit',type=int)
    ap.add_argument('--configs',nargs='+',default=list(CONFIGS));a=ap.parse_args()
    rows=json.loads((OUT/f'configs/parents_{a.group}.json').read_text())['parents']
    if a.limit:rows=rows[:a.limit]
    stage='C' if a.group=='CDEV' else 'D'
    for i,row in enumerate(rows):
        parent=row['parent']
        with np.load(OUT/'raw/online'/f'{parent}.npz') as z:y=z['observation'].copy();sigma=float(z['sigma'])
        if y.shape!=(32,8) or not np.isfinite(y).all():raise ValueError('observation contract')
        for config in a.configs:
            dest=OUT/f'track_c/contexts/{a.group}'/parent/config
            if (dest/'REPORT.json').exists():continue
            b=Book(stage,f'anchor_{parent}_{config}',parent,'ONLINE_INFERENCE',max_wall=120)
            status='COMPLETE';error=None;tick=time.perf_counter()
            try:
                target=Target(b,sigma=sigma)
                context,report=run(target,y,parent,dest,config)
                write(dest/'REPORT.json',dict(report,GT_read=False,sigma=sigma,end_to_end_anchor_wall_seconds=time.perf_counter()-tick))
                target.backend()._transfers.clear();del context,report,target;gc.collect()
            except BaseException as exc:status='FAILED';error=repr(exc);raise
            finally:b.close(status,error)
        print(json.dumps(dict(event='ANCHOR_PARENT_COMPLETE',group=a.group,completed=i+1,parent=parent,RHS=Book.totals()['RHS'])),flush=True)
    write(OUT/f'track_c/ANCHORS_COMPLETE_{a.group}.json',dict(status='COMPLETE',n=len(rows),configs=a.configs))

if __name__=='__main__':main()
