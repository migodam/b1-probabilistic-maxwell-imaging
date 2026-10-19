"""Isolated post-freeze evaluator. This is the only scoring entry reading fresh GT."""
from pathlib import Path
from functools import lru_cache
from collections import Counter,defaultdict
import sys,json,time,csv
import numpy as np
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from material_fiber.core import OUT,CAPS,write,save,mpack
from p1_material.prior import MaterialPrior,topology,structural_change

def lines(path):return [json.loads(x) for x in path.read_text().splitlines()] if path.exists() else []

@lru_cache(maxsize=3000)
def endpoint(key):
    path=OUT/key[0]/'endpoints'/f'{key}.npz'
    with np.load(path,allow_pickle=False) as a:return {k:a[k].copy() for k in a.files}

def main():
    if not (OUT/'ONLINE_CAMPAIGN_FROZEN.json').exists():raise RuntimeError('online methods not sealed')
    if (OUT/'OFFLINE_SCORING_CLOSED.json').exists():raise RuntimeError('preserve frozen scoring; use a new analysis version')
    start=time.perf_counter();cfg=json.loads((OUT/'configs/TARGET.json').read_text());P=MaterialPrior(cfg,json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']);parents=json.loads((OUT/'configs/PARENTS.json').read_text());anchors=json.loads((OUT/'A/ANCHORS.json').read_text());ar=json.loads((OUT/'A/SUMMARY.json').read_text())['parents'];base_rows={r['parent']:r for r in ar}
    candidates=[]
    for stage in ['B','C','D']:
        for r in lines(OUT/stage/'candidates.jsonl'):
            if r['parent'] in anchors and not r.get('post_scoring_baseline_audit',False):r['phase']=stage;candidates.append(r)
    metrics=[];packages={};selection={};gtmap={};risk=[];weak=[]
    def error(x,gt):
        d=x-gt;return dict(MAE=float(np.mean(abs(d))),RMS=float(np.linalg.norm(d)/24),Re_MAE=float(np.mean(abs(d.real))),Im_MAE=float(np.mean(abs(d.imag))),ROI_MAE=float(np.mean(abs(d[P.mask]))))
    def ordered(anchor,rr):
        # Online rule: anchor first; then material-space farthest-first; never GT.
        chosen=[anchor];available={r['key']:r for r in rr if r['data_feasible'] and r.get('DeltaPhi',1e9)<=2}
        keys=[]
        while available and len(chosen)<12:
            distance={k:min(P.material_rms(endpoint(k)['chi'],endpoint(c)['chi']) for c in chosen) for k in available}
            k=max(available,key=lambda k:(distance[k],-available[k]['Phi'],k))
            if distance[k]<.02-1e-10:break
            chosen.append(k);keys.append(k);del available[k]
        return chosen
    for pr in parents:
        parent=pr['parent'];base=endpoint(anchors[parent]);truthid=json.loads((OUT/f'offline_only/{parent}_truth.json').read_text())['key'];truth=endpoint(truthid);gt=truth['chi'];gtmap[parent]=truthid;bt=topology(P.decode(base['q'])[1]['fraction']);gtt=topology(P.decode(truth['q'])[1]['fraction']);be=error(base['chi'],gt)
        allrows=[r for r in candidates if r['parent']==parent];groups=sorted({r['phase']+':'+r['method'] for r in allrows})
        for group in groups+['ALL']:
            rr=allrows if group=='ALL' else [r for r in allrows if r['phase']+':'+r['method']==group]
            # Stress and pre-freeze duplicate versions remain in raw evidence, outside primary metrics.
            rr=[r for r in rr if r.get('version')!='stress']
            chosen=ordered(anchors[parent],rr);errs=[error(endpoint(k)['chi'],gt) for k in chosen];alternative=chosen[1:]
            structures=[structural_change(P,base['q'],endpoint(k)['q']) for k in alternative]
            record=dict(parent=parent,split=pr['split'],method=group,base_qualified=base_rows[parent]['data_feasible'],base=be,K_available=len(chosen),distinct_alternatives=len(alternative),structural_alternatives=sum(x['task_structure_pass'] for x in structures),topology_alternatives=sum(x['topology_before']['event']!=x['topology_after']['event'] for x in structures),candidate_keys=chosen,physics_rows=len(rr),data_feasible_rows=sum(r['data_feasible'] for r in rr),competitive_rows=sum(r['data_feasible'] and r.get('DeltaPhi',1e9)<=2 for r in rr),same_topology_as_GT=bt['event']==gtt['event'])
            for K in [1,2,4,8]:
                selected=errs[:K];record[f'BestOf{K}_MAE']=min(e['MAE'] for e in selected);record[f'BestOf{K}_RMS']=min(e['RMS'] for e in selected);record[f'Coverage{K}_RMS002']=any(e['RMS']<=.02 for e in selected)
            metrics.append(record)
            if group=='ALL':selection[parent]=chosen;risk.append(record)
        # Packages describe candidates/continuous fields, not probability events with inferred mass.
        selected=selection[parent];hy=[]
        with np.load(OUT/f'inputs/online/{parent}.npz') as a:y=a['observation'].copy()
        for i,k in enumerate(selected):
            ep=endpoint(k);struc=structural_change(P,base['q'],ep['q']);match=next((r for r in allrows if r['key']==k),base_rows[parent]);pred=ep['prediction'];phi=float(np.vdot(pred-y,pred-y).real/cfg['noise']['sigma']**2)
            hy.append(dict(candidate_id=f'candidate_{i}',endpoint_key=k,material_array=f'arrays/{k}.npz',representative='full passive chi field, epsilon_r=1+chi; 2D scalar TM only',event_definition='derived topology/structure descriptor, not a posterior partition',topology=topology(P.decode(ep['q'])[1]['fraction']),material_RMS_from_base=struc['material_RMS'],structure=struc,white_residual=float(np.sqrt(2*phi)),Phi=phi,DeltaPhi_from_fixed_base=phi-base_rows[parent]['Phi'],latent_logprior=P.logprior(ep['q']),posterior_mass=None,probability_status='UNESTIMATED_CANDIDATE_ONLY',data_feasible=bool(np.sqrt(2*phi)<=23.785449004057032),source_method=match.get('method'),forward_prediction=f'arrays/{k}.npz:prediction',full_Maxwell_numerical_status='saved exact endpoint; backward residual <=1e-10',offline_GT_not_used_for_selection=True))
            save(OUT/f'MaterialHypothesisPackage/arrays/{k}.npz',q=ep['q'],chi=ep['chi'],epsilon_r=1+ep['chi'],prediction=pred,fraction=P.decode(ep['q'])[1]['fraction'])
        pair=[]
        for i in range(len(selected)):
            for j in range(i):
                a,b=endpoint(selected[i]),endpoint(selected[j]);d=np.linalg.norm(np.r_[(a['prediction']-b['prediction']).real.ravel(),(a['prediction']-b['prediction']).imag.ravel()])*np.sqrt(2)/cfg['noise']['sigma']
                pair.append(dict(i=i,j=j,material_RMS=P.material_rms(a['chi'],b['chi']),whitened_mean_prediction_distance=float(d),fixed_pair_future_classification_error_interpretation='only meaningful for fixed means with same known Gaussian noise; not posterior odds'))
        package=dict(schema_version='P2-candidate-material-package-v1',observation_id=parent,target_identity=cfg['target_id'],support='original P1 correlated Gaussian64 x uniform2; passive effective material',hypotheses=hy,probability_semantics='candidate_only',posterior_samples=False,represented_mass=None,unshown_modeled_mass=None,unexplored_mass_status='not quantified; search is not exhaustive',outside_model_support_mass=None,model_adequacy_status='same-model synthetic benchmark only',pairwise_material_and_measurement_distances=pair,GT_not_in_online_payload=True)
        write(OUT/f'MaterialHypothesisPackage/{parent}.json',package);packages[parent]=package
        # Offline oracle displacement diagnosis using only the paid A receiver factor.
        with np.load(OUT/f'cache/{anchors[parent]}.npz') as a:T=a['receiver_T'].copy()
        if T.size:
            D=P.derivative(base['q']);M=np.r_[D.real,D.imag]/24;U,d,Vh=np.linalg.svd(M,full_matrices=False);keep=d>d[0]*1e-8;Z=Vh[keep].T/d[keep];Q=U[:,keep];C=np.einsum('qp,ps->qsp',T,base['field']).reshape(256,576);J=C@D
            from material_fiber.core import pack
            A=np.column_stack([pack(J[:,i].reshape(32,8))*np.sqrt(2)/cfg['noise']['sigma'] for i in range(66)]);K=A@Z;_,sk,Vk=np.linalg.svd(K,full_matrices=False);h=mpack(gt-base['chi']);coord=Q.T@h;weakmask=sk*.02<=1;energy=float(np.sum((Vk@coord)[weakmask]**2)/max(1e-30,h@h))
            weak.append(dict(parent=parent,split=pr['split'],weak_GT_energy_at_task002=energy,outside_tangent_energy=float(max(0,1-coord@coord/max(1e-30,h@h))),base_qualified=base_rows[parent]['data_feasible'],GT_prediction_distance=float(np.sqrt(2)*np.linalg.norm(truth['prediction']-base['prediction'])/cfg['noise']['sigma']),BestOf4_MAE_improvement=be['MAE']-risk[-1]['BestOf4_MAE'],ORACLE_GEOMETRY_DIAGNOSTIC=True))
    write(OUT/'OFFLINE_PARENT_METRICS.json',metrics);write(OUT/'OFFLINE_GT_WEAK_DISPLACEMENT.json',weak);write(OUT/'ONLINE_SELECTED_CANDIDATES.json',selection)
    with (OUT/'PARENT_METHOD_RESULTS.csv').open('w') as f:
        fields=['parent','split','method','base_qualified','K_available','distinct_alternatives','structural_alternatives','topology_alternatives','BestOf1_MAE','BestOf2_MAE','BestOf4_MAE','BestOf8_MAE','BestOf4_RMS','Coverage4_RMS002'];w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(metrics)
    held=[r for r in risk if r['split']=='heldout'];qualified=[r for r in held if r['base_qualified']]
    rng=np.random.default_rng(2026101999)
    def bootstrap(values):
        a=np.array(values,float)
        return list(np.quantile(a[rng.integers(len(a),size=(2000,len(a)))].mean(1),[.025,.975])) if len(a) else None
    method_summary=[]
    for group in sorted({r['method'] for r in metrics}):
        rr=[r for r in metrics if r['method']==group and r['split']=='heldout'];rrq=[r for r in rr if r['base_qualified']]
        if not rr:continue
        method_summary.append(dict(method=group,heldout_parents_observed=len(rr),qualified_anchor_parents=len(rrq),mean_distinct_alternatives=float(np.mean([r['distinct_alternatives'] for r in rr])),structural_alternatives=sum(r['structural_alternatives'] for r in rr),topology_alternatives=sum(r['topology_alternatives'] for r in rr),mean_BestOf4_MAE=float(np.mean([r['BestOf4_MAE'] for r in rr])),mean_base_MAE=float(np.mean([r['base']['MAE'] for r in rr])),MAE_improvement_CI=bootstrap([r['base']['MAE']-r['BestOf4_MAE'] for r in rr])))
    # Paired method comparisons are on the identical eight fresh held-out parents.
    paired=[]
    table={(r['parent'],r['method']):r for r in metrics}
    for competitor in ['B:full_J','B:hard_prior','B:soft_frozen','D:nonlearned']:
        for method in ['B:soft_dynamic','B:soft_adaptive','D:gaussian','D:flow']:
            pp=[r['parent'] for r in parents if r['split']=='heldout' and (r['parent'],competitor) in table and (r['parent'],method) in table]
            if pp:paired.append(dict(method=method,baseline=competitor,n=len(pp),BestOf4_MAE_delta=float(np.mean([table[p,method]['BestOf4_MAE']-table[p,competitor]['BestOf4_MAE'] for p in pp])),CI=bootstrap([table[p,method]['BestOf4_MAE']-table[p,competitor]['BestOf4_MAE'] for p in pp]),comparison='available anchor panel; not an all-eight intention-to-treat algorithm advantage claim'))
    write(OUT/'PAIRED_COMPARISONS.json',paired)
    costs=lines(OUT/'RESOURCE_LEDGER.jsonl');resource={}
    for stage,caps in CAPS.items():
        rr=[r for r in costs if r['stage']==stage];resource[stage]=dict(states=sum(r['amount'] for r in rr if r['event']=='STATE'),RHS=sum(r['amount'] for r in rr if r['event']=='RHS'),controlled_wall_seconds=sum(r['controlled_wall_seconds'] for r in rr if r['event']=='JOB_END'),physics_operation_wall_seconds=sum(r['wall_seconds'] for r in rr if r['event']=='PHYSICS_WALL' and r.get('additive',True)),LU=sum(r.get('counts',{}).get('LU_factorizations',0) for r in rr if r['event']=='ACTION'),full_native_forward_calls=sum(r.get('counts',{}).get('full_forward_calls',0) for r in rr if r['event']=='ACTION'),caps=caps)
    resource['total']={k:sum(resource[s][k] for s in CAPS) for k in ['states','RHS','controlled_wall_seconds','physics_operation_wall_seconds','LU','full_native_forward_calls']};resource['GPU_training_seconds']=0;resource['historical_26_states_976_RHS_readonly']=True
    for stage,caps in CAPS.items():
        if resource[stage]['states']>caps[0] or resource[stage]['RHS']>caps[1] or resource[stage]['controlled_wall_seconds']>caps[2]:raise RuntimeError('resource reconciliation mismatch')
    write(OUT/'RESOURCE_SUMMARY.json',resource)
    er=json.loads((OUT/'E/SUMMARY.json').read_text());acceleration=[]
    for method in sorted({r['method'] for r in er['rows']}):
        rr=[r for r in er['rows'] if r['method']==method];setup=sum(p['setup_seconds'][method] for p in er['setups'])
        acceleration.append(dict(method=method,queries=len(rr),fallbacks=sum(r['full_fallback'] for r in rr),total_RHS=sum(r['RHS'] for r in rr),median_cold_seconds=float(np.median([r['total_seconds'] for r in rr if r['repetition']==0])),median_warm_seconds=float(np.median([r['total_seconds'] for r in rr if r['repetition']>0])),setup_seconds=setup,amortized_seconds=(setup+sum(r['total_seconds'] for r in rr))/len(rr),max_white_error=max(r['white_prediction_error'] or 0 for r in rr),same_space_max=max([max(r.get('same_space_relative_errors',[0]) or [0]) for r in rr])))
    write(OUT/'ACCELERATION_COMPARISON.json',acceleration)
    chosen_demo=next((r['parent'] for r in parents if r['split']=='heldout' and len(selection[r['parent']])>1),next(r['parent'] for r in parents if len(selection[r['parent']])>1))
    wi=[r for r in weak if r['split']=='heldout' and r['base_qualified']]
    corr=float(spearmanr([r['weak_GT_energy_at_task002'] for r in wi],[r['BestOf4_MAE_improvement'] for r in wi]).statistic) if len(wi)>2 else None
    summary=dict(fresh_parents=len(parents),anchor_qualified=sum(r['data_feasible'] for r in ar),heldout_total=len(held),heldout_anchor_qualified=len(qualified),heldout_mean_base_MAE=float(np.mean([r['base']['MAE'] for r in held])),heldout_mean_BestOf4_MAE=float(np.mean([r['BestOf4_MAE'] for r in held])),heldout_BestOf4_improvement_CI=bootstrap([r['base']['MAE']-r['BestOf4_MAE'] for r in held]),heldout_distinct_alternatives=sum(r['distinct_alternatives'] for r in held),heldout_structural_alternatives=sum(r['structural_alternatives'] for r in held),heldout_topology_alternatives=sum(r['topology_alternatives'] for r in held),all_distinct_alternatives=sum(r['distinct_alternatives'] for r in risk),all_structural_alternatives=sum(r['structural_alternatives'] for r in risk),all_topology_alternatives=sum(r['topology_alternatives'] for r in risk),methods=method_summary,demo_parent=chosen_demo,weak_GT_energy_vs_MAE_gain_spearman=corr,posterior_mass=None,continuous_path_certificate=False,GT_used_online=False)
    write(OUT/'METRICS.json',summary);write(OUT/'OFFLINE_SCORING_CLOSED.json',dict(new_Maxwell_states=0,new_RHS=0,wall_seconds=time.perf_counter()-start,truth_reads='fresh saved dataset endpoints, isolated scoring after online freeze',online_algorithms_changed=False))
    print(json.dumps(summary,ensure_ascii=False))

if __name__=='__main__':main()
