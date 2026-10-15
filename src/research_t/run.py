"""Frozen bounded Research T experiments; no truth imports or reads."""
from pathlib import Path
from dataclasses import asdict
import argparse,json,time,subprocess
import numpy as np
from .budget import Book,dump,Limit
from .target import OUT,CASES,TAU,Targets,validate,save_state,pack
from .corrector import CorrectorConfig,correct,linearize,bounded_compensation_score
from .structure import cavity_birth,protected_bounds,check_protected,material_distances,topology_summary

def costs(book,start):
    rows=book.rows()[start:];a=[r for r in rows if r['event']=='ACTION'];roles={}
    for r in rows:
        if r['event']=='RHS':roles[r['role']]=roles.get(r['role'],0)+r['amount']
    return dict(new_full_state_attempts=sum(r.get('amount',0) for r in rows if r['event']=='STATE'),
                complex_full_solver_RHS=sum(roles.values()),RHS_by_role=roles,
                action_counts={k:sum(r.get('counts',{}).get(k,0) for r in a) for k in sorted({k for r in a for k in r.get('counts',{})})},
                failed_actions=sum(r.get('status')=='FAIL' for r in a),ledger_rows=[start,len(book.rows())])

def persist_result(folder,result,y,sigma):
    folder.mkdir(parents=True,exist_ok=True)
    save_state(folder/'final.npz',result.state)
    for i,state in enumerate(result.attempted_states):save_state(folder/f'trial_{i:02}.npz',state)
    if result.best_valid_state is not None:save_state(folder/'best_valid.npz',result.best_valid_state)
    row=dict(diagnostics=result.diagnostics,path=result.path,validation=validate(result.state,y,sigma))
    dump(folder/'result.json',row);return row

def qualify(book):
    root=OUT/'qualification';root.mkdir(parents=True,exist_ok=True);targets=Targets(book);summary=[]
    for case in CASES:
        book.context=dict(phase='E0-Q',case=case,truth_access=False)
        start=len(book.rows());model,anchor,y,sigma=targets.case(case)
        folder=root/case;folder.mkdir(exist_ok=True);save_state(folder/'original_n3.npz',anchor)
        np.savez_compressed(folder/'observation.npz',observation=y,sigma=sigma)
        initial=validate(anchor,y,sigma)
        cache=linearize(model,anchor,sigma)
        shared=costs(book,start);arms=[];states={}
        for mode in ['full','multiscale_haar']:
            book.context.update(arm=mode);start=len(book.rows());tic=time.perf_counter()
            cfg=CorrectorConfig(mode=mode,max_steps=4,max_full_attempts=6,trust_radius=.6,damping=.001,
                                tau=TAU,compatibility_target_factor=.95)
            result=correct(model,anchor,y,sigma,book,cfg,cache=cache)
            row=persist_result(folder/mode,result,y,sigma)
            row.update(arm=mode,cost=costs(book,start),wall_seconds=time.perf_counter()-tic)
            arms.append(row);states[mode]=result.best_valid_state
            print(json.dumps(dict(phase='qualification',case=case,arm=mode,initial=initial['residual_norm'],
                                  final=row['validation']['residual_norm'],valid=row['validation']['status'],cost=row['cost'])),flush=True)
        legal=[a for a in arms if states[a['arm']] is not None and validate(states[a['arm']],y,sigma)['status']=='FULL_VALIDATED_NUMERICAL']
        winner=min(legal,key=lambda a:validate(states[a['arm']],y,sigma)['residual_norm']) if legal else None
        if winner:save_state(folder/'qualified_base.npz',states[winner['arm']])
        entry=dict(case=case,initial=initial,shared_setup=shared,arms=arms,selected_arm=None if winner is None else winner['arm'],
                   qualified=winner is not None,selection='lowest observed residual among qualified actual full endpoints; no truth',
                   selected_validation=None if winner is None else validate(states[winner['arm']],y,sigma))
        summary.append(entry);dump(root/'SUMMARY.json',summary)
    return summary

def mechanism(book):
    # These fixed locations and the ring-first case choice are frozen before Q.
    protocol=json.loads((OUT/'configs/MECHANISM_v1.json').read_text())
    q=json.loads((OUT/'qualification/SUMMARY.json').read_text());lookup={r['case']:r for r in q}
    targets=Targets(book);root=OUT/'mechanism_v1';root.mkdir(parents=True,exist_ok=True);summary=[]
    for parent in protocol['parents']:
        choices=[f'p{parent}_{bg}_a1_snr20' for bg in protocol['background_priority']]
        choices=[c for c in choices if lookup[c]['qualified']]
        if not choices:continue
        case=choices[0];book.context=dict(phase='finite_mechanism_v1',case=case,truth_access=False)
        start=len(book.rows());model,_,y,sigma=targets.case(case)
        # A new matching base factor/endpoint is paid, not silently reconstructed.
        with np.load(OUT/'qualification'/case/'qualified_base.npz') as z:chi=z['chi'].copy()
        base=model.full_state(chi,purpose='ResearchT_mechanism_base_refactor');bv=validate(base,y,sigma)
        if bv['status']!='FULL_VALIDATED_NUMERICAL':raise RuntimeError('FROZEN_BASE_REPLAY_FAILED')
        base_cost=costs(book,start);case_folder=root/case;case_folder.mkdir(exist_ok=True)
        save_state(case_folder/'base.npz',base)
        for idx,center in enumerate(protocol['centers_row_col']):
            book.context.update(edit=idx);start=len(book.rows());tic=time.perf_counter()
            folder=case_folder/f'edit_{idx:02}';folder.mkdir(exist_ok=True)
            edit=cavity_birth(base.chi,center_row_col=tuple(center),radius_cells=2,min_feature_radius_cells=2,grid_spacing=1/24)
            low,high=protected_bounds(edit.core,edit.guard);protected=dict(lower=low,upper=high)
            np.savez_compressed(folder/'edit_geometry.npz',core=edit.core,guard=edit.guard,lower=low,upper=high)
            raw=model.full_state(edit.chi.ravel(),purpose='ResearchT_fixed_structural_edit')
            save_state(folder/'raw.npz',raw);rv=validate(raw,y,sigma)
            cfg=CorrectorConfig(mode='full',max_steps=3,max_full_attempts=3,trust_radius=.6,damping=.001,tau=TAU)
            score=bounded_compensation_score(model,raw,y,sigma,cfg,protected)
            cache=score.pop('jacobian_cache');score.pop('linear_delta')
            # Scoring does not change the ordinary corrector or its starting point.
            result=correct(model,raw,y,sigma,book,cfg,protected,cache=cache)
            saved=persist_result(folder,result,y,sigma);end=result.best_valid_state or result.state
            fv=validate(end,y,sigma);dist=material_distances(base.chi,end.chi,mass=1/576)
            guard=check_protected(end.chi,edit.core,edit.guard)
            genuine=bool(fv['status']=='FULL_VALIDATED_NUMERICAL' and guard['preserved'] and dist['significant_phase_edit'] and dist['normalized_material_l2']>=.03)
            row=dict(case=case,edit_id=idx,edit_metadata=edit.metadata,base_validation=bv,shared_base_cost=base_cost,
                     raw_validation=rv,compensation_score=score,corrected_validation=fv,
                     structural_protection=guard,material_distance=dist,topology=topology_summary(end.chi),
                     genuine_new_candidate=genuine,invalid_to_valid=genuine and rv['status']=='REJECTED',
                     cost=costs(book,start),wall_seconds=time.perf_counter()-tic,corrector=saved)
            dump(folder/'SUMMARY.json',row);summary.append(row);dump(root/'SUMMARY.json',summary)
            if genuine:save_state(folder/'new_valid_candidate.npz',end)
            print(json.dumps(dict(phase='mechanism',case=case,edit=idx,raw=rv['residual_norm'],corrected=fv['residual_norm'],
                                  genuine=genuine,phase_cells=dist['phase_changed_cells'],cost=row['cost'])),flush=True)
    return summary

def main():
    p=argparse.ArgumentParser();p.add_argument('phase',choices=['qualification','mechanism']);args=p.parse_args()
    book=Book(OUT)
    try:
        with book.job(args.phase+'_v1'):
            qualify(book) if args.phase=='qualification' else mechanism(book)
    finally:dump(OUT/'RESOURCE_RECEIPT.json',book.receipt())

if __name__=='__main__':main()
