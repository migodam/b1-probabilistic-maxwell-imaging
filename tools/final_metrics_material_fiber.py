from pathlib import Path
from collections import Counter
import sys,json,shutil
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tools'))
from material_fiber.core import OUT,CAPS,write,save
from p1_material.prior import MaterialPrior,structural_change
from analyze_material_fiber import endpoint,lines

def main():
    cfg=json.loads((OUT/'configs/TARGET.json').read_text());P=MaterialPrior(cfg,json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']);parents=json.loads((OUT/'configs/PARENTS.json').read_text());anchors=json.loads((OUT/'A/ANCHORS.json').read_text());primary=json.loads((OUT/'OFFLINE_PARENT_METRICS.json').read_text());b=lines(OUT/'B/candidates.jsonl');supp=[r for r in b if r.get('post_scoring_baseline_audit')];results=[];rng=np.random.default_rng(2026101998)
    for method in ['fullJ_sphere_only','hard_prior_sphere_only']:
        for r in parents:
            parent=r['parent'];a=endpoint(anchors[parent]);gtkey=json.loads((OUT/f'offline_only/{parent}_truth.json').read_text())['key'];gt=endpoint(gtkey)['chi'];rr={v['key']:v for v in supp if v['parent']==parent and v['method']==method and v['data_feasible'] and v['DeltaPhi']<=2};chosen=[anchors[parent]]
            while rr and len(chosen)<12:
                dist={k:min(P.material_rms(endpoint(k)['chi'],endpoint(j)['chi']) for j in chosen) for k in rr};k=max(dist,key=dist.get)
                if dist[k]<.02-1e-10:break
                chosen.append(k);del rr[k]
            err=[float(np.mean(abs(endpoint(k)['chi']-gt))) for k in chosen];st=[structural_change(P,a['q'],endpoint(k)['q']) for k in chosen[1:]]
            results.append(dict(parent=parent,split=r['split'],method=method,distinct_alternatives=len(chosen)-1,structural_alternatives=sum(x['task_structure_pass'] for x in st),base_MAE=err[0],BestOf4_MAE=min(err[:4]),candidate_keys=chosen,status='exposed supplementary control'))
    summary=[]
    for method in ['fullJ_sphere_only','hard_prior_sphere_only']:
        rr=[r for r in results if r['method']==method and r['split']=='heldout'];delta=[]
        for r in rr:
            dyn=next((v for v in primary if v['parent']==r['parent'] and v['method']=='B:soft_dynamic'),None)
            if dyn:delta.append(dyn['BestOf4_MAE']-r['BestOf4_MAE'])
        ci=np.quantile(np.array(delta)[rng.integers(len(delta),size=(2000,len(delta)))].mean(1),[.025,.975]).tolist() if delta else None
        summary.append(dict(method=method,heldout_total=8,qualified_anchors=4,heldout_distinct_alternatives=sum(r['distinct_alternatives'] for r in rr),heldout_structural_alternatives=sum(r['structural_alternatives'] for r in rr),all8_mean_base_MAE=float(np.mean([r['base_MAE'] for r in rr])),all8_mean_BestOf4_MAE=float(np.mean([r['BestOf4_MAE'] for r in rr])),dynamic_minus_baseline_MAE_CI=ci,comparison_n=len(delta),status='EXPOSED_AFTER_PRIMARY_SCORING'))
    write(OUT/'SUPPLEMENT_RESULTS.json',dict(parent_metrics=results,summary=summary,post_scoring=True,primary_METRICS_unchanged=True))
    costs=lines(OUT/'RESOURCE_LEDGER.jsonl');resource={}
    for stage,caps in CAPS.items():
        rr=[r for r in costs if r['stage']==stage];resource[stage]=dict(states=sum(r['amount'] for r in rr if r['event']=='STATE'),RHS=sum(r['amount'] for r in rr if r['event']=='RHS'),controlled_wall_seconds=sum(r['controlled_wall_seconds'] for r in rr if r['event']=='JOB_END'),physics_operation_wall_seconds=sum(r['wall_seconds'] for r in rr if r['event']=='PHYSICS_WALL' and r.get('additive',True)),LU=sum(r.get('counts',{}).get('LU_factorizations',0) for r in rr if r['event']=='ACTION'),full_native_forward_calls=sum(r.get('counts',{}).get('full_forward_calls',0) for r in rr if r['event']=='ACTION'),caps=caps)
        assert resource[stage]['states']<=caps[0] and resource[stage]['RHS']<=caps[1] and resource[stage]['controlled_wall_seconds']<=caps[2]
    resource['total']={k:sum(resource[s][k] for s in CAPS) for k in ['states','RHS','controlled_wall_seconds','physics_operation_wall_seconds','LU','full_native_forward_calls']};resource['GPU_training_seconds']=0;resource['old_26_states_976_RHS_unchanged']=True;resource['paid_physics_stage_status']='FINAL_CLOSED_NO_AUTOMATIC_REUSE_OF_REMAINDERS'
    write(OUT/'RESOURCE_SUMMARY_FINAL.json',resource)
    # Preserve original payloads; final payload never calls an incompatible base a feasible world.
    for source in list((OUT/'MaterialHypothesisPackage').glob('*.json')):
        package=json.loads(source.read_text());dest=OUT/'MaterialHypothesisPackage/FINAL'/source.name
        qualified=[h for h in package['hypotheses'] if h['data_feasible']]
        failed=[h for h in package['hypotheses'] if not h['data_feasible']]
        package.update(schema_version='P2-candidate-material-package-v1.1',hypotheses=qualified,failed_baseline=failed,package_status='CANDIDATE_ONLY' if qualified else 'NO_DATA_COMPATIBLE_HYPOTHESES',array_reference_base='MaterialHypothesisPackage root; not FINAL subfolder',search_is_exhaustive=False)
        write(dest,package)
    # Add all setup RHS to the evaluation table; query-only RHS would hide exact-port costs.
    acc=json.loads((OUT/'ACCELERATION_COMPARISON.json').read_text())
    setup={'full':0,'exact_T':4*(8+256),'goal_ROM':4*(8+32),'Krylov':4*8,'A17_exchange':4*(8+32)}
    for r in acc:r['setup_RHS']=setup[r['method']];r['total_including_setup_RHS']=r['total_RHS']+r['setup_RHS']
    assert sum(r['total_including_setup_RHS'] for r in acc)==resource['E']['RHS']
    write(OUT/'ACCELERATION_FULL_COST.json',acc)
    metrics=json.loads((OUT/'METRICS.json').read_text());topo=json.loads((OUT/'TOPOLOGY_STABILITY.json').read_text());decision=dict(overall='REDIRECT',local_finite_noise_exploration='GO_AS_EMPIRICAL_DIAGNOSTIC',soft_dynamic_superiority='NOT_ESTABLISHED',strong_FullJ_control='COMPLETED_EXPOSED',learned_explorer='NO_GO_CURRENT_PILOT',A17_X_integration='NO_GO_ON_TM24_PLATFORM',posterior_calibration='NOT_RUN_NO_REFERENCE',exact_local_fiber='NOT_ESTABLISHED',cross_topology_important_modes='NOT_ESTABLISHED',strong_TAP_algorithm_paper='NO_GO_CURRENT_EVIDENCE',reason='Task-scale finite-noise alternatives exist, but soft/dynamic costs exceed hard split; removing Full-J angular constraint yields strong competing solutions; no heldout topology mode evidence, learned gain or A17 speed benefit',support=dict(primary=metrics,supplement=summary,robust_topology_changes_all_three_thresholds=sum(x['structure_pass'] and all(x['topology_change_by_threshold'].values()) for x in topo)),resources=resource['total'])
    write(OUT/'GATE_DECISION.json',decision)
    print(json.dumps(dict(supplement=summary,resources=resource['total']),ensure_ascii=False,indent=2))
if __name__=='__main__':main()
