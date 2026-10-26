"""Truth-free original Native96 + Probit64. No offline evaluator imports."""
import argparse,json,time,gc,resource
import numpy as np
from physics_confirm import OUT,Book,Target,save,write,plain,read_online,reconstruct,probit_optimize

def main():
    ap=argparse.ArgumentParser();ap.add_argument('group',choices=['probe','A','B']);a=ap.parse_args()
    manifest=json.loads((OUT/f'configs/parents_{a.group}.json').read_text())
    for index,row in enumerate(manifest['parents']):
        parent=row['parent'];dest=OUT/'raw/contexts'/f'{parent}.npz'
        if dest.exists() and dest.with_suffix('.json').exists():continue
        if 3300000-Book.totals()['RHS']<6500:
            write(OUT/f'raw/RESOURCE_LIMITED_{a.group}.json',dict(parent=parent,completed=index,RHS=Book.totals()['RHS'],reason='reserve full-fidelity next-parent worst case'))
            raise RuntimeError('RESOURCE_LIMITED: no low-fidelity replacement')
        b=Book('anchors_'+parent,parent,'ONLINE_INFERENCE');tick=time.perf_counter();before=b.totals();status='COMPLETE';error=None
        try:
            y=read_online(OUT/'raw/online'/f'{parent}.npz');t=Target(b)
            folder=b.jobdir/'details';folder.mkdir(parents=True,exist_ok=True)
            native,nstat=reconstruct(t,y,parent,method='full_TRF',max_evaluations=96)
            save(folder/'NATIVE.npz',q=native.q,chi=native.state.chi,prediction=native.state.data)
            point,x,H,pstat=probit_optimize(t,y,parent,folder/'probit',native,cap=64)
            J=t.jacobian(point)
            save(dest,observation=y,q0=point.q,chi0=point.state.chi,J=J,prediction=point.state.data,
                 residual=y-point.state.data,phi=np.array(t.phi(point,y)),x=x,H=H,parent=np.array(parent),native_chi=native.state.chi,native_q=native.q)
            after=b.totals()
            write(dest.with_suffix('.json'),dict(parent=parent,group=a.group,GT_read=False,native_cap=96,probit_cap=64,
                  Phi=t.phi(point,y),native_stationarity=nstat['kkt'],native_report=nstat,probit_stationarity=pstat,
                  new_RHS=after['RHS']-before['RHS'],new_states=after['state_attempts']-before['state_attempts'],
                  new_LU_attempts=after['LU_attempts']-before['LU_attempts'],wall_seconds=time.perf_counter()-tick,
                  target_id=t.cfg['target_id'],optimizer_math='unchanged original R4 imported functions',cache_output='lossless fields/J; no persisted dense LU',
                  fallback='original optimizer best finite endpoint retained; no GT restart'))
            t.backend()._transfers.clear();del native,point,t,J,H;gc.collect()
        except BaseException as exc:
            status='FAILED';error=repr(exc);write(OUT/'raw/failures'/f'{parent}.json',dict(error=error,parent=parent,group=a.group));raise
        finally:b.close(status,error)
        print(json.dumps(dict(group=a.group,completed=index+1,parent=parent,RHS=Book.totals()['RHS'],parent_RHS=after['RHS']-before['RHS'],wall_seconds=time.perf_counter()-tick,peak_RSS_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))),flush=True)
    write(OUT/f'raw/ANCHORS_COMPLETE_{a.group}.json',dict(status='COMPLETE',count=len(manifest['parents']),group=a.group,GT_read=False))
if __name__=='__main__':main()
