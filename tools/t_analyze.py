"""Analyze already paid Research T arrays; no solver or truth access."""
import json,csv,time
from pathlib import Path
import numpy as np
from research_t.target import OUT,pack,TAU
from research_t.budget import dump,Book
from research_t.structure import material_distances

def main():
    tic=time.perf_counter();allrows=[];versions=[];ranking=[];pairs=[]
    for v in [1,2,3]:
        rows=json.loads((OUT/f'mechanism_v{v}/SUMMARY.json').read_text())
        for r in rows:
            p=OUT/f'mechanism_v{v}'/r['case'];f=p/f'edit_{r["edit_id"]:02}'
            with np.load(p/'base.npz') as z:base={k:z[k].copy() for k in z.files}
            e=f/'best_valid.npz'
            if not e.exists():e=f/'final.npz'
            with np.load(e) as z:end={k:z[k].copy() for k in z.files}
            with np.load(f/'raw.npz') as z:raw=z['prediction'].copy()
            with np.load(f/'edit_geometry.npz') as z:core=z['core'].ravel()
            with np.load(OUT/'qualification'/r['case']/'observation.npz') as z:sigma=float(z['sigma'])
            pd=float(np.linalg.norm(pack(end['prediction']-base['prediction'])*np.sqrt(2)/sigma))
            rd=float(np.linalg.norm(pack(raw-base['prediction'])*np.sqrt(2)/sigma))
            x=dict(version=v,case=r['case'],edit_id=r['edit_id'],raw_residual=r['raw_validation']['residual_norm'],
                   corrected_residual=r['corrected_validation']['residual_norm'],valid=r['genuine_new_candidate'],
                   invalid_to_valid=r['invalid_to_valid'],linear_U=r['compensation_score']['linear_primal_upper'],
                   linear_L=r['compensation_score']['linear_dual_lower'],
                   normalized_material_L2=r['material_distance']['normalized_material_l2'],
                   phase_changed_cells=r['material_distance']['phase_changed_cells'],
                   core_new_air_cells=int(np.sum(core&(abs(base['chi'])>.1))),
                   prediction_distance_raw_to_base=rd,prediction_distance_corrected_to_base=pd,
                   pair_distance_reduction_fraction=1-pd/rd,
                   core_guard_preserved=r['structural_protection']['preserved'],
                   state_attempts=r['cost']['new_full_state_attempts'],RHS=r['cost']['complex_full_solver_RHS'])
            allrows.append(x)
        s=[r for r in allrows if r['version']==v];valid=[r for r in s if r['valid']]
        versions.append(dict(version=v,proposals=len(s),raw_valid=sum(r['raw_residual']<=TAU for r in s),
                             invalid_to_valid=sum(r['invalid_to_valid'] for r in s),
                             independent_parents_with_success=len({r['case'] for r in valid}),
                             one_step_U_above_tau_but_nonlinear_valid=sum(r['linear_U']>TAU and r['valid'] for r in s),
                             corrected_residual_range=[min(r['corrected_residual'] for r in s),max(r['corrected_residual'] for r in s)],
                             raw_residual_range=[min(r['raw_residual'] for r in s),max(r['raw_residual'] for r in s)],
                             valid_material_L2_range=[min(r['normalized_material_L2'] for r in valid),max(r['normalized_material_L2'] for r in valid)],
                             prediction_distance_to_base_range=[min(r['prediction_distance_corrected_to_base'] for r in valid),max(r['prediction_distance_corrected_to_base'] for r in valid)]))
        for case in sorted({r['case'] for r in s}):
            rr=[r for r in s if r['case']==case];plain=sorted(rr,key=lambda x:x['edit_id']);ranked=sorted(rr,key=lambda x:x['linear_U'])
            for name,order in [('frozen_order',plain),('U_order_replay',ranked)]:
                # U ranking requires ALL raw endpoints and initial32-adjoint
                # score states, not a free ordering oracle. Shared base omitted
                # identically here and declared in source stage summaries.
                states=len(rr) if name=='U_order_replay' else 0
                rhs=40*len(rr) if name=='U_order_replay' else 0
                found=0
                for t,r in enumerate(order):
                    states+=r['state_attempts']-(1 if name=='U_order_replay' else 0)
                    rhs+=r['RHS']-(40 if name=='U_order_replay' else 0)
                    found+=r['valid']
                    ranking.append(dict(version=v,case=case,method=name,rank=t+1,edit_id=r['edit_id'],valid_found=int(found),states=states,RHS=rhs,
                                        status='ACCOUNTED_OFFLINE_ORDER_REPLAY' if name=='U_order_replay' else 'OBSERVED_FROZEN_SEQUENCE'))
    for case in sorted({r['case'] for r in allrows}):
        candidates=[]
        for r in allrows:
            if r['case']!=case or not r['valid']:continue
            path=OUT/f'mechanism_v{r["version"]}'/case/f'edit_{r["edit_id"]:02}'/'best_valid.npz'
            with np.load(path) as z:candidates.append((r['version'],r['edit_id'],z['chi'].copy(),z['prediction'].copy()))
        with np.load(OUT/'qualification'/case/'observation.npz') as z:sigma=float(z['sigma'])
        for i,a in enumerate(candidates):
            for b in candidates[i+1:]:
                d=material_distances(a[2],b[2],mass=1/576)
                pairs.append(dict(case=case,a=f'v{a[0]}edit{a[1]}',b=f'v{b[0]}edit{b[1]}',
                                  material_L2=d['normalized_material_l2'],phase_cells=d['phase_changed_cells'],
                                  mean_prediction_whitened_distance=float(np.linalg.norm(pack(a[3]-b[3])*np.sqrt(2)/sigma))))
    for name,rows in [('PER_EDIT_METRICS',allrows),('RANKING_COST_REPLAY',ranking),('PAIRWISE_VALID_MATERIALS',pairs)]:
        with (OUT/f'{name}.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    q=json.loads((OUT/'qualification/SUMMARY.json').read_text())
    offline=json.loads((OUT/'offline/METRICS.json').read_text())['records'];bestk=[]
    for v in [1,2,3]:
        for case in sorted({r['case'] for r in allrows if r['version']==v}):
            base_role='qualified_minchange_base' if v==3 else 'qualified_base'
            base=next(r for r in offline if r['case']==case and r['phase']=='qualification' and r['role']==base_role)
            pool=[dict(id='base',MAE=base['metrics']['full_complex_MAE'])]
            for r in sorted([r for r in allrows if r['version']==v and r['case']==case and r['valid']],key=lambda r:r['edit_id']):
                scored=next(a for a in offline if a['case']==case and a['phase']==f'mechanism_v{v}' and a['role']=='actual_best_valid' and int(a['edit_id'])==r['edit_id'])
                pool.append(dict(id=f'edit{r["edit_id"]}',MAE=scored['metrics']['full_complex_MAE']))
            for k in [1,2,4,8]:
                bestk.append(dict(version=v,case=case,K_requested=k,K_available=min(k,len(pool)),base_MAE=pool[0]['MAE'],
                                  best_MAE=min(a['MAE'] for a in pool[:k]),offline_only=True,order='base then successful frozen edit_id order'))
    dump(OUT/'OFFLINE_BEST_OF_K.json',bestk)
    metrics=dict(target='A30-original-24x24-scalarTM-v1',candidate_weight_semantics='candidate_only; probability null',
                 qualified_conditions=sum(r['qualified'] for r in q),total_conditions=len(q),iterations=versions,
                 physical_qualification=json.loads((OUT/'actual_physics_qualification/RESULT.json').read_text()),
                 resources=Book(OUT).receipt(),truth_arrays_read_by_this_script=False,offline_truth_derived_metrics_read=True,
                 offline_best_of_K=bestk,archive=json.loads((OUT/'T_PHYSICS_HANDOFF/MATERIAL_ARCHIVE.json').read_text()),
                 comparison_limit='No multistart/QD comparison, no A17 incremental experiment. C score selection is explicitly charged replay, not independent online trial.',
                 analysis_wall_seconds=time.perf_counter()-tic)
    dump(OUT/'METRICS.json',metrics)
    dump(OUT/'RUN_CLOSED.json',dict(reason='192 full-state attempts exhausted',new_physics_allowed=False,
                                 states=metrics['resources']['new_full_state_attempts'],rhs=metrics['resources']['complex_full_solver_RHS'],
                                 unspendable_RHS_balance=6144-metrics['resources']['complex_full_solver_RHS'],
                                 requires='New explicit human authorization; do not reopen or reset this stage'))
    print(json.dumps(dict(qualified=metrics['qualified_conditions'],iterations=versions,resources={k:v for k,v in metrics['resources'].items() if k!='action_counts'})))
if __name__=='__main__':main()
