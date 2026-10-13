"""Use17 where13 failed, preserving frozen models and all previous outputs."""
import os
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from tools.p2_cdf_partition import integrate
from src.b1_joint.reference import Reference
from src.b1_joint.response import Response
from src.b1_joint.accounting import dump
def main():
    out=ROOT/'research/joint_p2/results';audit=json.loads((out/'ALL_FRESH_REFERENCE_REFINEMENT.json').read_text());ids=audit['failed_original13_observations'];accepted=out/'FRESH_CONTINUOUS_REFERENCE_ACCEPTED.npz'
    if accepted.exists():
        loaded=np.load(accepted);arrays={k:loaded[k] for k in loaded.files};orders=arrays['reference_orders']
    else:
        base=np.load(out/'FRESH_CONTINUOUS_REFERENCE_PARTITION.npz');data=np.load(out/'FRESH_FINAL_OBSERVATIONS.npz');ref=Reference(Response(out/'COMMON_RESPONSE_BANK.npz'));cdf,cond,z,ep=integrate(ref,data['y'][ids,:32],data['sigma'][:32],base['grid'],17);arrays={k:base[k] for k in base.files};orders=np.full(len(arrays['mean']),13)
        for j,i in enumerate(ids):
            row=audit['rows'][i];arrays['CDF'][i]=cdf[j];arrays['conditional_CDF'][i]=cond[j];arrays['mean'][i]=row['order17_mean'];arrays['SD'][i]=row['order17_SD'];arrays['event_mass'][i]=row['order17_event_mass'];orders[i]=17
        arrays['reference_orders']=orders;np.savez_compressed(accepted,**arrays)
    detail=json.loads((out/'FRESH_CONTINUOUS_PARTITION_AUDIT.json').read_text());detail.update(final_reference='FRESH_CONTINUOUS_REFERENCE_ACCEPTED.npz',refined_observations=ids,mean_SD_event_refinement_all_observations=True,unresolved_after21=audit['unresolved_after21'],monotone_violation=float(max(0.,-np.diff(arrays['CDF'],axis=-1).min())),upper_error=float(abs(arrays['CDF'][:,:,-1]-1).max()),claim_scope='positiveCDFgrid all128; complete event/mean/SD13/17check all128, three upgraded17 and crosschecked21; CDF13/17 onlyfixed8 plusadaptive refinement3, no uniform/continuum certificate')
    dump(out/'ACCEPTED_REFERENCE_AUDIT.json',detail)
    C1=Path(os.environ.get('C1_P2',ROOT.parent.parent/'C1/c1-joint-p2'));summary=json.loads((C1/'p2/results/evaluation/EVALUATION_SUMMARY.json').read_text());rows=summary['rows'];p=arrays['event_mass'];metrics=[];rng=np.random.default_rng(26101329)
    for row in rows:
        i=row['observation_index'];prob=np.array(row['model_event_mass']);sig=p[i]>=.05;metrics.append(dict(method=row['method'],seed=row['seed'],observation=i,parent=int(row['parent_id']),event_MAE=float(abs(prob-p[i]).mean()),mass_recall=float(p[i][sig&(prob>=.01)].sum()/p[i][sig].sum()),count_recall=float(np.mean(prob[sig]>=.01)),false_exclusion=bool(np.any(sig&(prob<.005))),reference_order=int(orders[i])))
    aggregates=[]
    for method,seed in sorted(set((a['method'],a['seed']) for a in metrics)):
        a=[b for b in metrics if b['method']==method and b['seed']==seed];aggregates.append(dict(method=method,seed=seed,event_MAE=float(np.mean([b['event_MAE'] for b in a])),count_recall=float(np.mean([b['count_recall'] for b in a])),mass_recall=float(np.mean([b['mass_recall'] for b in a])),false_exclusion_count=sum(b['false_exclusion'] for b in a)))
    delta=[]
    for parent in range(64):
        f=[a['event_MAE'] for a in metrics if a['parent']==parent and a['method']=='full_distill'];d=[a['event_MAE'] for a in metrics if a['parent']==parent and a['method']=='diag_distill'];delta.append(np.mean(f)-np.mean(d))
    b=rng.choice(delta,(2000,64),replace=True).mean(1);dump(out/'REFINED_FRESH_EVENT_COMPARISON.json',dict(aggregates=aggregates,rows=metrics,full_minus_diagonal=dict(mean=float(np.mean(delta)),CI95=np.quantile(b,[.025,.975]).tolist()),models_unchanged=True,new_physics=0,old_metrics_retained=True,reference='all128 moment/event audit; upgraded3 rows17 crosschecked21',latent_KL='old13quadrature expectation retained empirical, not recomputed/material KL'))
if __name__=='__main__':main()
