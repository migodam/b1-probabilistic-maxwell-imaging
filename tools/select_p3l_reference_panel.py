from pathlib import Path
import sys,json,shutil,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target,write,pack

def main():
    if (OUT/'A/DEEP_REFERENCE_PARENTS.json').exists():raise RuntimeError('reference selection already frozen')
    b=Book('A','reference_panel_y_only_selection',max_wall=300);status='COMPLETE';error=None
    try:
        T=Target(b);P=T.prior;method=json.loads((OUT/'A/BASELINE_SELECTION.json').read_text())['method'];rows=[]
        for i in range(32):
            parent=f'DEV_{i:03d}';stem=OUT/f'A/anchor_runs/{parent}_{method}';r=json.loads(stem.with_suffix('.json').read_text())
            with np.load(stem.with_suffix('.npz')) as a:
                q=a['q'];D=P.derivative(q);transfer=a['receiver_T'];field=a['field']
                if not transfer.size:raise RuntimeError('expected paid receiver factor at anchor')
                C=np.einsum('qp,ps->qsp',transfer,field).reshape(256,576);Z=C@D
                J=np.column_stack([pack(Z[:,k].reshape(32,8))*np.sqrt(2)/T.sigma for k in range(66)])
            aux=P.decode(q)[1];f=aux['fraction'][P.mask]
            ev=np.linalg.eigvalsh(J[:,:64].T@J[:,:64]);condition=float(np.sqrt(max(ev[-1],1e-30)/max(ev[0],1e-20)))
            feat=[r['phi'],float(np.mean((f<1e-6)|(f>1-1e-6))),float(np.min(np.r_[q[64:],1-q[64:]])),float(np.log1p(condition))]
            rows.append(dict(parent=parent,features=feat,noise_ball_qualified=r['noise_ball_qualified'],KKT=r['kkt']['scaled_KKT'],objective=r['objective'],method=method))
            dest=OUT/f'A/anchors/{parent}.npz';dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(stem.with_suffix('.npz'),dest);write(dest.with_suffix('.json'),r)
        X=np.asarray([r['features'] for r in rows]);X=(X-X.mean(axis=0))/np.maximum(X.std(axis=0),1e-8)
        chosen=[0]
        while len(chosen)<16:
            dist=np.min(np.sum((X[:,None]-X[chosen][None])**2,axis=2),axis=1);dist[chosen]=-1;chosen.append(int(np.argmax(dist)))
        parents=[rows[k]['parent'] for k in chosen]
        write(OUT/'A/ANCHOR_PANEL.json',dict(parents=rows,truth_used=False,all32_in_risk=True))
        write(OUT/'A/DEEP_REFERENCE_PARENTS.json',dict(parents=parents,indices=chosen,selection='predeclared deterministic y-only maxmin standardized anchor features; first DEV_000',feature_names=['Phi','saturation_fraction','property_boundary_distance','log_condition_Jz'],truth_used=False,selected_before_posterior_results=True))
        print(json.dumps(dict(deep_parents=parents,qualified_anchors=sum(r['noise_ball_qualified'] for r in rows),stationary_anchors=sum(r['KKT'] is not None and r['KKT']<.001 for r in rows))),flush=True)
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:b.close(status,error)

if __name__=='__main__':main()
