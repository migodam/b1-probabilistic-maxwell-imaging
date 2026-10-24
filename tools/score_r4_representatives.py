#!/usr/bin/env python3
"""Separate-process offline scoring after actual inference/request/Book gates."""
from pathlib import Path
import argparse,json,os,sys,csv,time,resource
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.special import ndtri
from p3r4.representative_eval import check_representative_gates
from p3r4.evaluation import material_losses,paired_bootstrap,write_new_json
ACTIONS=('raw_mu','clipped_mu','ordinary_extra_Phi','representative_raw','representative_clipped','original_anchor')


def finite(v):
    if isinstance(v,dict):return {k:finite(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [finite(x) for x in v]
    if isinstance(v,np.ndarray):return finite(v.tolist())
    if isinstance(v,(float,np.floating)) and not np.isfinite(v):return None
    if isinstance(v,np.generic):return v.item()
    return v


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('physical','prepared-request','online-gate','receipt','truth','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();wall,cpu=time.monotonic(),time.process_time()
    if a.out.exists():raise FileExistsError('fresh score output required')
    prepared,run,gate=check_representative_gates(a.physical,a.prepared_request,a.online_gate,a.receipt)
    if 'offline_only' not in a.truth.resolve().parts:raise ValueError('isolated offline_only GT directory required')
    stage=ROOT/'research/p3l_r4_conditional_transport';cfg=json.loads((stage/'configs/TARGET.json').read_text());points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
    from p1_material.prior import MaterialPrior
    prior=MaterialPrior(cfg,points);roi=np.flatnonzero(prior.mask);req=prepared['request'];parents=req['parents'];n=len(parents)
    source=Path(prepared['source_predictions']['path'])
    with np.load(source,allow_pickle=False) as z:
        parent_source=z['parent'].astype(str).tolist();lookup={name:i for i,name in enumerate(parent_source)};mu={key:z[key].copy() for key in req['outputs']}
    if req['mu_layout']=='realimag_N2_24_24':mu={k:(v[:,0]+1j*v[:,1]).reshape(len(v),576) for k,v in mu.items()}
    qualified={r['parent']:r for r in prepared['qualified']};skips={r['parent']:r for r in prepared['skips']}
    metrics={key:np.full((n,len(ACTIONS)),np.nan) for key in ('full_mse','full_mae','roi_mse','roi_mae','Phi','probit_prior_energy','native_z_prior_energy','RMS_from_anchor')}
    rows=[];failures=[];perparent=[];a.out.mkdir(parents=True)
    for index,parent in enumerate(parents):
        predictions=np.full((len(ACTIONS),576),np.nan+1j*np.nan);metadata={};anchor_q=None;anchor_chi=None
        # First GT file access: all actual online/request/closedBook gates above
        # have already been checked in this independent process.
        truth_error=None
        try:
            with np.load(a.truth/f'{parent}.npz',allow_pickle=False) as z:truth=z['chi'];truth_q=z['q']
            if truth.shape!=(576,) or truth_q.shape!=(66,) or not np.iscomplexobj(truth) or not np.isfinite(truth).all() or not np.isfinite(truth_q).all():raise ValueError('invalid isolated truth')
            true_s=(cfg['field_mean']+prior.A[roi]@truth_q[:64])/prior.temperature
        except Exception as e:truth=None;true_s=None;truth_error=repr(e)
        item=qualified.get(parent);cache=Path(item['anchor_cache']) if item else stage/'results/contexts/cache'/f'{parent}.npz'
        try:
            with np.load(cache,allow_pickle=False) as z:
                anchor_q=z['q'];anchor_chi=z['chi'];prediction=z['prediction']
            with np.load(stage/'inputs/online'/f'{parent}.npz',allow_pickle=False) as z:y=z['observation']
            predictions[5]=anchor_chi;anchor_phi=float(np.linalg.norm(prediction-y)**2/cfg['noise']['sigma']**2)
            metadata['original_anchor']=dict(status='AVAILABLE',source_kind='original_paid_physics_cache',source=str(cache),Phi=anchor_phi,fallback=False,q=anchor_q)
        except Exception as e:metadata['original_anchor']=dict(status='NOT_EVALUATED',reason=repr(e),source_kind='original_paid_physics_cache')
        for key,action in [('prediction_raw','raw_mu'),('prediction_clipped','clipped_mu')]:
            if key in mu and parent in lookup:
                estimate=mu[key][lookup[parent]]
                if estimate.shape==(576,) and np.isfinite(estimate).all():predictions[ACTIONS.index(action)]=estimate;metadata[action]=dict(status='AVAILABLE',source_kind='frozen_online_risk_action',source=str(source),Phi=None,fallback=False)
                else:metadata[action]=dict(status='NOT_EVALUATED',reason='NONFINITE_OR_INVALID_FROZEN_MU')
            else:metadata[action]=dict(status='NOT_EVALUATED',reason='OUTPUT_NOT_DECLARED_OR_PARENT_MISSING')
        for arm,action in [('ordinary_likelihood_shared','ordinary_extra_Phi'),('prediction_raw','representative_raw'),('prediction_clipped','representative_clipped')]:
            case=a.physical/parent/arm
            try:
                if item is None:raise ValueError('UNQUALIFIED_SKIP:'+str(skips.get(parent)))
                if arm!='ordinary_likelihood_shared' and arm not in req['outputs']:raise ValueError('REPRESENTATIVE_NOT_REQUESTED')
                result=json.loads((case/'RESULTS.json').read_text())
                if result.get('GT_read') is not False or result.get('posterior_mean') is not False or result.get('posterior_sample') is not False:raise ValueError('online output role/GT contract invalid')
                with np.load(case/'FINAL.npz',allow_pickle=False) as z:
                    estimate=z['chi'];q=z['q'];x=z['x'];phi=float(z['Phi']);fallback=bool(z['fallback'])
                if estimate.shape!=(576,) or not np.isfinite(estimate).all() or not np.array_equal(prior.decode(q)[0],estimate):raise ValueError('final originaldecoder material mismatch')
                predictions[ACTIONS.index(action)]=estimate
                metadata[action]=dict(status='AVAILABLE',source_kind='new_paid_ordinary_control' if arm=='ordinary_likelihood_shared' else 'new_paid_Output3_representative',source=str(case/'FINAL.npz'),
                    Phi=phi,fallback=fallback,x=x,q=q,optimizer=result.get('optimizer'),stop_reason=result.get('stop_reason'),kkt=result.get('kkt'),actual_cost=result.get('actual_cost'))
            except Exception as e:metadata[action]=dict(status='NOT_EVALUATED',reason=repr(e),source_kind='new_paid_output')
        anchor_s=None if anchor_q is None else (cfg['field_mean']+prior.A[roi]@anchor_q[:64])/prior.temperature
        parent_metrics={};parent_records=[]
        for col,action in enumerate(ACTIONS):
            info=metadata[action];estimate=predictions[col]
            info['qualified_in_frozen_request']=parent in qualified
            if info['status']=='AVAILABLE' and truth is not None:
                loss=material_losses(estimate[None],truth[None],roi,None if anchor_s is None else anchor_s[None],true_s[None]);parent_metrics[action]={k:v[0] for k,v in loss.items()}
                for k in ('full_mse','full_mae','roi_mse','roi_mae'):metrics[k][index,col]=loss[k][0]
                status='EVALUATED_FALLBACK' if info.get('fallback') else 'EVALUATED'
            else:status='NOT_EVALUATED';failures.append(dict(parent=parent,action=action,reason=truth_error if truth is None else info.get('reason')))
            phi=info.get('Phi');metrics['Phi'][index,col]=np.nan if phi is None else phi
            q=info.get('q');x=info.get('x')
            if q is not None:
                metrics['native_z_prior_energy'][index,col]=.5*float(q[:64]@q[:64])
                if x is None and np.all((q[64:]>0)&(q[64:]<1)):x=np.r_[q[:64],ndtri(q[64:])]
                if x is not None:metrics['probit_prior_energy'][index,col]=.5*float(x@x)
            if anchor_chi is not None and np.isfinite(estimate).all():metrics['RMS_from_anchor'][index,col]=float(np.sqrt(np.mean(abs(estimate-anchor_chi)**2)))
            record=dict(parent=parent,action=action,status=status,fallback=info.get('fallback',False),qualified=parent in qualified,
                        **{k:metrics[k][index,col] for k in metrics},source_kind=info.get('source_kind'),source=info.get('source'),reason=info.get('reason',truth_error))
            rows.append(record);parent_records.append({**record,'phase_strata':parent_metrics.get(action),'online_optimizer':{k:info.get(k) for k in ('optimizer','stop_reason','kkt','actual_cost')}})
        perparent.append(dict(parent=parent,qualified_in_request=parent in qualified,predeclared_skip=skips.get(parent),actions=parent_records))
        write_new_json(a.out/f'{parent}.json',finite(perparent[-1]))
        np.savez_compressed(a.out/f'{parent}.npz',actions=np.array(ACTIONS),predictions=predictions,**{k:v[index] for k,v in metrics.items()},
                            **{action+'__'+key:np.asarray(value) for action,loss in parent_metrics.items() for key,value in loss.items()})
    comparisons=[]
    for representative,target in [('representative_raw','raw_mu'),('representative_clipped','clipped_mu'),('representative_raw','original_anchor'),('representative_clipped','original_anchor'),('ordinary_extra_Phi','original_anchor')]:
        x=metrics['full_mse'][:,ACTIONS.index(representative)];y=metrics['full_mse'][:,ACTIONS.index(target)];valid=np.isfinite(x)&np.isfinite(y)
        item=dict(representative=representative,reference=target,planned_parents=n,evaluable_pairs=int(valid.sum()),MSE_change_rep_minus_reference=(x-y).tolist())
        if valid.all() and n>1:item.update(status='EVALUATED_ALL_DECLARED_PARENTS',**paired_bootstrap(x,y))
        else:item.update(status='NOT_EVALUATED_INCOMPLETE_NO_PARENT_DROPPING',finite_pair_mean_change_descriptive=float(np.mean(x[valid]-y[valid])) if valid.any() else None)
        comparisons.append(item)
    summaries={action:{k:dict(planned_parents=n,evaluated_parents=int(np.isfinite(v[:,col]).sum()),all_parent_mean=float(v[:,col].mean()) if np.isfinite(v[:,col]).all() else None,
              finite_only_mean_descriptive=float(v[np.isfinite(v[:,col]),col].mean()) if np.isfinite(v[:,col]).any() else None) for k,v in metrics.items()} for col,action in enumerate(ACTIONS)}
    with (a.out/'PER_PARENT.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(finite(rows))
    np.savez_compressed(a.out/'PER_PARENT.npz',parent=np.array(parents),actions=np.array(ACTIONS),**metrics)
    write_new_json(a.out/'RESULTS.json',finite(dict(status='OFFLINE_REPRESENTATIVE_SCORE_COMPLETE',gate=gate,scoring_pid=os.getpid(),parents=perparent,summary=summaries,comparisons=comparisons,
        failures=failures,predeclared_skips=prepared['skips'],physical_cost=run.get('actual_cost'),physics_calls=0,training_updates=0,
        output_role='materialrisk actions plus constrained Output3 representatives/control; no posterior label',
        field_provenance='original physics cache, frozen NNmu, and new physical FINAL fields remain separate source kinds',
        wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024))))
    print(json.dumps(dict(status='OFFLINE_REPRESENTATIVE_SCORE_COMPLETE',parents=n,unevaluated_entries=len(failures))))

if __name__=='__main__':main()
