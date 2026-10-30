"""Independent saved-response/J, trace and prepaid-budget audit. Never opens GT.

This is an offline array reader, not an engine. Chronological paid-state secants
include rejected trials and branch changes, explicitly marked as such.
"""
import argparse
import json
import fcntl
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adaptive_depth import finite_secant, projector_change, orthobasis
from offline_score import load_predictions, resolve, read, write, csv_write


def white(F,sigma=.00119211798514519):
    F=np.asarray(F)
    return np.concatenate((F.real.T,F.imag.T),axis=1).ravel()*np.sqrt(2)/sigma


def state(path):
    with np.load(path,allow_pickle=False) as z:
        if any('truth' in k.lower() or k.lower() in ('xgt','target_digit','target_occupancy') for k in z.files):
            raise ValueError('GT-bearing state rejected')
        keys=('x','xs','material','F','prediction','jac','J','dmaterial','residual','phi','Phi')
        return {k:np.asarray(z[k]).copy() for k in keys if k in z.files}


def physical_weak_basis(J,D,rank=8):
    """Data-only weak tangent in original ROI material metric, not engine frame."""
    M=np.r_[D.real,D.imag]/24
    _,s,vh=np.linalg.svd(M,full_matrices=False)
    keep=s>max(float(s[0])*1e-10,np.finfo(float).eps*max(M.shape)*float(s[0])) if len(s) else np.zeros(0,bool)
    if not keep.any():return np.zeros((M.shape[0],0))
    T=vh[keep].T/s[keep][None]
    _,_,kvh=np.linalg.svd(J@T,full_matrices=True)
    W=T@kvh[-min(rank,len(T.T)):].T
    return orthobasis(M@W)


def traces(path):
    path=Path(path)
    if not path.exists():return []
    if path.suffix=='.jsonl':return [json.loads(s) for s in path.read_text().splitlines() if s.strip()]
    obj=read(path)
    return obj if isinstance(obj,list) else obj.get('rows',obj.get('trace',[]))


def audit_budget(root):
    root=Path(root)
    with (root/'GLOBAL_BUDGET.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_SH)
        snapshot=read(root/'GLOBAL_BUDGET_STATE.json')
        rows=[json.loads(s) for s in (root/'COST_LEDGER.jsonl').read_text().splitlines() if s.strip()]
        fcntl.flock(lock,fcntl.LOCK_UN)
    reserved={k:0 for k in ('LU','RHS','RHS576')};ticket={};complete=set();by_order={};by_group={};issues=[]
    lastseq=0
    for row in rows:
        seq=row['sequence']
        if seq<=lastseq:issues.append('nonincreasing ledger sequence')
        lastseq=seq
        if row['event']=='RESERVE':
            r=row['request'];t=row['ticket']
            if t in ticket:issues.append('duplicate reservation ticket')
            ticket[t]=row
            if r['RHS576']>r['RHS']:issues.append('RHS576 subset exceeds total RHS')
            for k in reserved:reserved[k]+=int(r[k])
            group=row['group'];order=str(row.get('meta',{}).get('order','UNSPECIFIED'))
            for mapping,key in ((by_group,group),(by_order,order)):
                v=mapping.setdefault(key,dict(LU=0,RHS=0,RHS576=0,reservations=0))
                for k in reserved:v[k]+=int(r[k])
                v['reservations']+=1
        if row['event']=='COMPLETE':
            t=row['ticket']
            if t not in ticket:issues.append('completion has no preceding reservation')
            if t in complete:issues.append('duplicate completion ticket')
            complete.add(t)
    undercharged={k:reserved[k]>snapshot['used'][k] for k in reserved}
    overcharge={k:snapshot['used'][k]-reserved[k] for k in reserved}
    if any(undercharged.values()):issues.append('snapshot undercharged versus durable ledger')
    if any(snapshot['used'][k]>snapshot['limits'][k] for k in reserved):issues.append('global hard cap exceeded')
    return dict(status='PASS' if not issues else 'FAIL',issues=issues,limits=snapshot['limits'],snapshot_used=snapshot['used'],
                ledger_reserved=reserved,conservative_unmatched_snapshot_charge=overcharge,
                uncompleted_reservations=[t for t in ticket if t not in complete],by_order=by_order,by_group=by_group,
                snapshot_sequence=snapshot['sequence'],scope='reservation charges, including failures; completion timing is not charged again',new_physics_calls=0)


def diagnose(manifest_path,out,max_J_states=4):
    manifest,groups,cfg=load_predictions(manifest_path);base=Path(manifest_path).resolve().parent
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    secants=[];projectors=[];layers=[];events=[]
    for g in groups:
        ident={k:g[k] for k in ('dataset','parent','method')}
        folder=g['source_path'].parent if g['source_path'] else None
        if folder is None:continue
        tp=resolve(g['trace_path'],base) if g.get('trace_path') else folder/'TRACE.json'
        if not tp.exists():tp=folder/'TRACE.jsonl'
        tr=traces(tp)
        for row in tr:
            if row.get('phase')=='layer':events.append(dict(**ident,**row))
        for branch in sorted(set(r.get('branch',0) for r in tr)):
            br=[r for r in tr if r.get('branch',0)==branch];prev=0
            for lay in sorted(set(r['layer'] for r in br if 'layer' in r)):
                rr=[r for r in br if r.get('layer')==lay];end=max([r.get('cost_RHS',prev) for r in rr]+[prev])
                trials=[t for r in rr for h in r.get('correction',{}).get('history',[]) for t in h.get('trials',[])]
                defects=[t['secant_relative_defect'] for t in trials if 'secant_relative_defect' in t]
                layers.append(dict(**ident,branch=branch,layer=lay,RHS_layer=end-prev,RHS_prefix=end,
                     finite_secant_max=max(defects) if defects else None,finite_secant_count=len(defects),
                     accepted_trials=sum(t.get('status')=='ACCEPTED' for t in trials),
                     attempted_trials=len(trials),source='branch PaidEvaluator cumulative requested cost; final validation separately charged'))
                prev=end
        sp=resolve(g['states_path'],base) if g.get('states_path') else folder/'states'
        previous_J=None;previous_path=None;previous_basis=None
        paths=sorted(sp.glob('*.npz')) if sp.exists() else []
        Jpaths=[]
        for path in paths:
            with np.load(path,allow_pickle=False) as z:
                if 'jac' in z.files or 'J' in z.files:Jpaths.append(path)
        indices=np.unique(np.linspace(0,len(Jpaths)-1,min(max_J_states,len(Jpaths))).astype(int)) if Jpaths and max_J_states>0 else []
        sampled=[Jpaths[i] for i in indices]
        if paths and max_J_states>0 and paths[-1] not in sampled:sampled.append(paths[-1])
        for path in sorted(sampled):
            z=state(path);x=z.get('x');F=z.get('F',z.get('prediction'));J=z.get('jac',z.get('J'))
            if x is None or F is None:continue
            if previous_J is not None:
                row=finite_secant(previous_J['x'],x,white(previous_J['F']),white(F),previous_J['J'],J)
                row.update(**ident,from_state=previous_path.name,to_state=path.name,
                           path_semantics='uniform-index sample of at most4 paid J states plus final paid state; includes rejected trials and branch changes',
                           new_J_at_to=J is not None)
                secants.append(row)
            if J is not None:
                if 'dmaterial' in z:
                    B=physical_weak_basis(J,z['dmaterial'])
                    if previous_basis is not None:projectors.append(dict(**ident,from_state=previous_path.name,to_state=path.name,
                          metric='originalROI realified complex material /24; data-only weak8; distinct from prior-aware proposal frame',
                          **projector_change(previous_basis,B)))
                    previous_basis=B
                previous_J=dict(x=x,F=F,J=J);previous_path=path
    write(out/'FINITE_SECANTS.json',secants);write(out/'J_PROJECTOR_EVOLUTION.json',projectors)
    write(out/'FRAME_REFRESH_EVENTS.json',events);write(out/'LAYER_COSTS.json',layers)
    csv_write(out/'LAYER_COSTS.csv',layers)
    if manifest.get('budget_root'):write(out/'BUDGET_AUDIT.json',audit_budget(resolve(manifest['budget_root'],base)))
    if secants:
        fig,axs=plt.subplots(1,2,figsize=(11,4.5))
        for method in sorted(set(r['method'] for r in secants)):
            rr=[r for r in secants if r['method']==method]
            axs[0].plot(np.sort([r['secant_defect_relative'] for r in rr]),np.linspace(0,1,len(rr)),label=method)
            pp=[r['new_outside_old_fraction'] for r in projectors if r['method']==method]
            if pp:axs[1].plot(np.sort(pp),np.linspace(0,1,len(pp)),label=method)
        axs[0].set_xlabel('Finite whitened-response secant defect');axs[1].set_xlabel('Material weak8 projector novelty')
        for ax in axs:ax.set_ylabel('Empirical fraction of saved transitions');ax.grid(alpha=.2);ax.legend(fontsize=7)
        fig.suptitle('Saved paid states, including rejected trials; no statistical information-gain claim',fontsize=11)
        fig.tight_layout();fig.savefig(out/'SECANT_PROJECTOR_EVOLUTION.png',dpi=150);plt.close(fig)
    write(out/'ONLINE_DIAGNOSTIC_RECEIPT.json',dict(groups=len(groups),finite_secants=len(secants),projector_pairs=len(projectors),
          GT_files_opened=0,new_physics_calls=0,independent_J_state_sample_per_group=max_J_states,
          metrics_scope='finite actual states; sample for independent projector check, full traces for layer costs/defects; no between-node bounds'))
    return dict(secants=len(secants),projector_pairs=len(projectors),layer_cost_rows=len(layers))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest');p.add_argument('--out',required=True);p.add_argument('--budget-root');a=p.parse_args()
    if a.budget_root:
        write(Path(a.out)/'BUDGET_AUDIT.json',audit_budget(a.budget_root))
    else:print(json.dumps(diagnose(a.manifest,a.out)))
