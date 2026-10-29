"""Independent offline-only audit. Never imports the physics or scoring engines."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[k]='1'
from pathlib import Path
import json, csv, time, hashlib
import numpy as np

P=Path(__file__).resolve().parents[2]
A=Path(__file__).resolve().parent
F=P/'results/candidates_fresh10'
SIGMA=.00119211798514519
GATE=282.8737921622988
NEW=('FULL66_RML','SHAPE_TARGET_RML_BLOCK','BILATERAL_TARGET_RML_BLOCK')
def read(p): return json.loads(Path(p).read_text())
def arrays(p):
    with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')

def ledger(path,kind):
    rows=[json.loads(x) for x in path.read_text().splitlines()]
    errors=[]
    if [r['sequence'] for r in rows]!=list(range(1,len(rows)+1)):errors.append('sequence gap/reset')
    if kind=='transfer':
        rhs=sum(r.get('rhs',0) for r in rows);lu=sum(r.get('lu',0) for r in rows)
        if rhs!=rows[-1]['RHS256'] or lu!=rows[-1]['LU256']:errors.append('reservation sum mismatch')
        return dict(RHS256=rhs,LU256=lu,errors=errors,rows=len(rows),jobs=[dict(event=r['event'],utc=r['utc'],RHS256=r['RHS256'],LU256=r['LU256'],wall_seconds=r.get('wall_seconds')) for r in rows if r['event'] in ('BEGIN','END','END_REVISION2')])
    keys=('reduced_RHS256','reduced_LU256','full_RHS576','full_LU576','actual_forward_RHS256','actual_J_RHS256','actual_LU256','completed_forward_RHS256','completed_J_RHS256')
    totals={k:sum(r.get(k,0) for r in rows) for k in keys}
    for k,v in totals.items():
        if v!=rows[-1]['after'][k]:errors.append(k+' cumulative mismatch')
    if totals['reduced_RHS256']!=totals['actual_forward_RHS256']+totals['actual_J_RHS256']:errors.append('reserved vs actual RHS mismatch')
    if totals['reduced_LU256']!=totals['actual_LU256']:errors.append('reserved vs actual LU mismatch')
    phases={}
    for r in rows:
        phase=r.get('phase','DEV');b=phases.setdefault(phase,{k:0 for k in keys})
        for k in keys:b[k]+=r.get(k,0)
    return dict(**totals,errors=errors,rows=len(rows),phases=phases,final_wall_seconds=rows[-1].get('wall_seconds'),jobs=[dict(event=r['event'],utc=r['utc'],phase=r.get('phase'),wall_seconds=r.get('wall_seconds')) for r in rows if r['event'] in ('JOB_BEGIN','JOB_END')])

def transfer_nodes():
    result=[]
    for p in sorted((P/'results/transfer').glob('*_integration_nodes_revision2.csv')):
        rows=list(csv.DictReader(p.open()));r=read(p.with_name(p.name.replace('_integration_nodes_revision2.csv','_conditional_revision2.json')))
        for lev in ('0','1','2'):
            q=[x for x in rows if x['level']==lev];v=np.array([[float(x['delta_real']),float(x['delta_imag'])] for x in q]);w=np.array([float(x['integration_weight']) for x in q]);phi=np.array([float(x['phi']) for x in q]);a=np.log(w)-phi-np.log(.15);top=a.max();z=top+np.log(np.exp(a-top).sum());pn=np.exp(a-z);mean=pn@v;event=pn@np.column_stack((v[:,0]>1.3,v[:,1]>.105));stored=r['axis_references'][int(lev)] if lev!='2' else r['independent_local_reference']
            result.append(dict(parent=r['parent'],level=lev,nodes=len(q),positive_weights=bool(np.all(w>0)),all_nodes_in_original_prior=bool(np.all(v>=np.array([.8,.03])) and np.all(v<=np.array([1.8,.18]))),log_evidence=float(z),mean=mean.tolist(),events=event.tolist(),max_stored_error=max(abs(z-stored['log_evidence']),float(np.max(abs(mean-stored['amplitude_mean']))),float(np.max(abs(event-stored['event_probability']))))))
    return result

def main():
    ledgers=dict(transfer=ledger(P/'results/transfer/COST_LEDGER.jsonl','transfer'),DEV=ledger(P/'results/candidates/RESOURCE_LEDGER.jsonl','candidate'),fresh=ledger(F/'RESOURCE_LEDGER.jsonl','candidate'))
    write(A/'LEDGER_RECONCILIATION.json',dict(new_RHS=0,ledgers=ledgers,global_RHS256=sum(x.get('RHS256',x.get('reduced_RHS256',0)) for x in ledgers.values()),global_LU256=sum(x.get('LU256',x.get('reduced_LU256',0)) for x in ledgers.values())))
    write(A/'TRANSFER_NODE_RECOMPUTATION.json',dict(new_RHS=0,rows=transfer_nodes()))
    # Guard before even locating private files or opening score artifacts.
    if not (F/'PREDICTIONS_SEALED.json').exists() or not (F/'GT_SCORING_OPENED.json').exists():
        print('HOLD: complete prediction seal and GT scoring open required; truth untouched');return
    c=read(P/'config/candidate_fresh10_protocol.json');seal=read(F/'PREDICTIONS_SEALED.json');opened=read(F/'GT_SCORING_OPENED.json')
    assert len(seal['rows'])==480 and seal['GT_used'] is False
    assert opened['utc']>seal['utc']
    stats=read(F/'STATISTICS.json')
    errors=[];rows=[];endpoint_checks=[]
    for entry in c['parents']:
        parent=entry['parent'];y=arrays(F/'inputs/online'/f'{parent}.npz')['observation'];anchor=arrays(F/'parents'/parent/'anchor/ANCHOR.npz');truth=arrays(P/'private/candidate_fresh10/truth'/f'{parent}.npz')['material'];aphi=float(np.sum(abs(anchor['prediction']-y)**2)/SIGMA**2)
        for arm in c['arms']:
            rr=sorted([r for r in seal['rows'] if r['parent']==parent and r['arm']==arm],key=lambda r:r['k']);assert [r['k'] for r in rr]==list(range(8))
            loaded={}
            for r in rr:
                if not r.get('path'):continue
                z=arrays(r['path']);loaded[r['k']]=z;phi=float(np.sum(abs(z['prediction']-y)**2)/SIGMA**2);legal=bool(z['q'].shape==(66,) and np.isfinite(z['q']).all() and np.all(z['q'][64:]>=0) and np.all(z['q'][64:]<=1));qual=legal and phi<=GATE
                if arm not in NEW:qual=qual and r['radius_progress_qualified']
                err=abs(phi-r['Phi']);endpoint_checks.append(dict(parent=parent,arm=arm,k=r['k'],Phi=phi,stored_Phi_error=err,qualified=qual,stored_qualified=bool(r['qualified']),legal=legal))
                if err>1e-7 or qual!=bool(r['qualified']):errors.append([parent,arm,r['k'],'gate/Phi mismatch'])
            for K in (1,2,4,8):
                candidates=[(anchor['material'],aphi)] if aphi<=GATE else []
                for r in rr[:K]:
                    if r.get('path') and r['qualified']:candidates.append((loaded[r['k']]['material'],float(r['Phi'])))
                unique=[]
                for m,phi in candidates:
                    if all(np.sqrt(np.mean(abs(m-u[0])**2))>1e-9 for u in unique):unique.append((m,phi))
                if unique:
                    X=np.stack([x[0] for x in unique]);ph=np.array([x[1] for x in unique]);mse=np.mean(abs(X-truth)**2,axis=1);center=np.mean(X,axis=0);weights=np.exp(-ph+ph.min());weights/=weights.sum();dist=np.sqrt(np.mean(abs(X[:,None,:]-X[None,:,:])**2,axis=2));s=np.linalg.svd(np.concatenate(((X-center).real,(X-center).imag),axis=1)/24,compute_uv=False)
                    score=dict(oracle_MSE=float(mse.min()),oracle_RMS=float(np.sqrt(mse.min())),coverage_eps075=bool(mse.min()<=.075**2),uniform_action_MSE=float(np.mean(abs(center-truth)**2)),likelihood_heuristic_action_MSE=float(np.mean(abs(weights@X-truth)**2)),pair_diameter=float(dist.max()),mean_pair_distance=float(dist[np.triu_indices(len(X),1)].mean()) if len(X)>1 else 0.,RMS_dispersion=float(np.sqrt(np.mean(abs(X-center)**2))),cloud_rank=int(np.sum(s>max(s[0]*1e-8,1e-10))),unique_qualified=len(unique),action_fallback=False)
                else:score=dict(oracle_MSE=None,oracle_RMS=None,coverage_eps075=False,uniform_action_MSE=float(np.mean(abs(anchor['material']-truth)**2)),likelihood_heuristic_action_MSE=float(np.mean(abs(anchor['material']-truth)**2)),pair_diameter=None,mean_pair_distance=None,RMS_dispersion=None,cloud_rank=0,unique_qualified=0,action_fallback=True)
                raw=[loaded[r['k']]['material'] for r in rr[:K] if r.get('path')]
                score.update(anchor_MSE=float(np.mean(abs(anchor['material']-truth)**2)),qualified_proposals=sum(bool(r['qualified']) for r in rr[:K]),attempted_proposals=sum(bool(r.get('path')) for r in rr[:K]),qualified_rate=sum(bool(r['qualified']) for r in rr[:K])/K,anchor_selected_proposals=sum(bool(r.get('anchor_selected')) for r in rr[:K]),failure_count=sum(r['status']=='FAILED' for r in rr[:K]),right_censored_count=sum('CENSORED' in r['status'] for r in rr[:K]),all_selected_oracle_MSE_diagnostic=None if not raw else float(np.min(np.mean(abs(np.stack(raw)-truth)**2,axis=1))),RHS256=sum(r.get('cost',{}).get('reduced_RHS256',0) for r in rr[:K]),LU256=sum(r.get('cost',{}).get('reduced_LU256',0) for r in rr[:K]),online_wall_seconds=sum(r.get('wall_seconds',0) for r in rr[:K]),CPU_seconds=sum(r.get('CPU_seconds',0) for r in rr[:K]))
                row=dict(parent=parent,arm=arm,K=K,**score);rows.append(row);stored=next(x for x in stats['rows'] if x['parent']==parent and x['arm']==arm and x['K']==K)
                for key,value in score.items():
                    actual=stored[key]
                    if value is None or isinstance(value,(bool,int)):
                        ok=value==actual
                    else:ok=abs(value-actual)<1e-10
                    if not ok:errors.append([parent,arm,K,key,value,actual])
    write(A/'FRESH10_SCORE_RECOMPUTATION.json',dict(status='PASS' if not errors else 'FAIL',new_RHS=0,recomputed_utc=time.time(),prediction_seal_utc=seal['utc'],scoring_opened_utc=opened['utc'],rows=rows,endpoint_checks=endpoint_checks,errors=errors,all_10_parents_retained=True))
    print('Independent offline score audit:',len(rows),'cloud-prefix rows,',len(endpoint_checks),'endpoint checks,',len(errors),'errors')

if __name__=='__main__':main()
