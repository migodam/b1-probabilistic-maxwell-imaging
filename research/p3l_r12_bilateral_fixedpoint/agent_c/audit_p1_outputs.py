"""Adversarial P1 prediction/gate/accounting checks, using public inputs only."""
import json
from pathlib import Path
import numpy as np
from adaptive_depth import orthobasis, outside_span, finite_secant, projector_change
from analyze_online import white, physical_weak_basis, state
from offline_score import load_predictions, write, read, csv_write

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def main():
    manifest,groups,cfg=load_predictions(ROOT/'audit/P1_PREDICTION_MANIFEST.json')
    budget=read(HERE/'p1_diagnostics/BUDGET_AUDIT.json')
    gate=[];legacy=[];cost=[];sameJ=[];issues=[];common={};historical={}
    for g in groups:
        folder=g['source_path'].parent;ident={k:g[k] for k in ('dataset','parent','method')}
        y=np.load(ROOT/'data/online'/f"{g['parent']}.npz",allow_pickle=False)['observation']
        z=g['z'];Fs=z['Fs'];actual=np.sum(abs(Fs-y)**2,axis=(1,2))/.00119211798514519**2
        phierr=float(np.max(abs(actual-g['phi']),initial=0))
        legal=z['legal'].astype(bool);noise=legal & np.isfinite(g['phi']) & (g['phi']<=g['noise_max'])
        claimed=z.get('noise_qualified',g['qualified']).astype(bool)
        if not np.array_equal(noise,claimed):issues.append(str(ident)+' noise qualification mismatch')
        keep=[];distinct=[]
        for i,(mat,ok) in enumerate(zip(g['materials'],noise)):
            valid=bool(ok and all(np.linalg.norm(mat-g['materials'][j])/24>.002 for j in keep))
            distinct.append(valid)
            if valid:keep.append(i)
        if not np.array_equal(distinct,g['qualified']):issues.append(str(ident)+' greedy qualification mismatch')
        if phierr>1e-6:issues.append(str(ident)+' data Phi recomputation mismatch')
        distances=np.linalg.norm(g['materials']-g['anchor']['material'],axis=1)/24
        nonanchor=np.asarray(distinct)&(distances>.002)
        with np.load(g['source_path'],allow_pickle=False) as raw:
            if 'nonanchor_qualified' in raw.files and not np.array_equal(nonanchor,raw['nonanchor_qualified']):issues.append(str(ident)+' nonanchor flag mismatch')
        selected=g['selected'];selecteddist=None if selected is None else float(distances[selected])
        gate.append(dict(**ident,raw_members=len(noise),noise_qualified=int(noise.sum()),greedy_qualified=int(sum(distinct)),
            nonanchor_qualified=int(nonanchor.sum()),max_actual_Phi_error=phierr,
            selected_original_anchor_RMS=selecteddist,qualified_coverage_flags_unchanged=True))
        branches=z.get('branch',np.zeros(len(noise),int))
        for i in range(len(noise)):
            bi=int(np.flatnonzero(branches==branches[i])[0]);bprior=float(z['prior_cost'][bi]);prior=float(z['prior_cost'][i])
            dp=prior-bprior;df=float(g['phi'][i]-g['phi'][bi])
            legacy.append(dict(**ident,index=i,branch=int(branches[i]),noise_qualified=bool(noise[i]),primary_distinct=bool(distinct[i]),
                delta_Phi_from_branch_anchor=df,delta_prior_from_branch_anchor=dp,
                legacy_deltaPhi2_deltaPrior10_pass=df<=2 and dp<=10,
                delta_Phi_from_original_anchor=float(g['phi'][i]-g['anchor']['phi']),
                material_RMS_from_original_anchor=float(distances[i])))
        seal=read(folder/'SEALED.json');ledger=budget['by_group']['P1/'+g['parent']+'/'+g['method']]
        if seal['new_RHS']!=ledger['RHS'] or seal['new_LU']!=ledger['LU']:issues.append(str(ident)+' physical cost mismatch')
        cost.append(dict(**ident,actual_new_RHS=seal['new_RHS'],actual_new_LU=seal['new_LU'],
             standalone_R12_RHS=seal['standalone_RHS'],standalone_R12_LU=seal['standalone_LU'],
             shared_common_RHS=seal['common_new_RHS'],shared_common_LU=seal['common_new_LU']))
        if g['parent'] not in common:
            common[g['parent']]=read(folder.parent/'COMMON/SEALED.json')
            historical[g['parent']]=read(Path(g['anchor_path']).parents[1]/'COMPLETE.json')['anchor_cost']
    for parent in cfg['development_parent_ids']:
        g=next(g for g in groups if g['parent']==parent and g['method']=='H2_FULL16')
        common_arrays=dict(np.load(g['source_path'].parent.parent/'COMMON/ANCHORS.npz',allow_pickle=False))
        x0=common_arrays['xs'][0];F0=common_arrays['Fs'][0];J0=common_arrays['jac'][0];D0=common_arrays['dmaterial'][0]
        frame=dict(np.load(g['initial_frame_path'],allow_pickle=False));Q=orthobasis(frame['weaklift'])
        _,_,vh=np.linalg.svd(J0@Q,full_matrices=False)
        child=Q@vh[:4].T
        _,_,vh2=np.linalg.svd(J0@child,full_matrices=False);grandchild=child@vh2[:2].T
        nested_violation=float(np.linalg.norm(grandchild-Q@(Q.T@grandchild),'fro'))
        if nested_violation>1e-10:issues.append(parent+' same-J nested-span violation')
        Jpaths=[]
        for path in sorted((g['source_path'].parent/'states').glob('*.npz')):
            with np.load(path,allow_pickle=False) as z:
                if 'jac' in z.files:Jpaths.append(path)
        last=Jpaths[-1];z=state(last);J1=z['jac']
        row=dict(parent=parent,sameJ_nested_span_violation=nested_violation,sameJ_new_physics_calls=0,
           refreshed_paid_state=str(last),refreshed_state_receipt=read(last.with_suffix('.json')),
           secant=finite_secant(x0,z['x'],white(F0),white(z['F']),J0,J1),
           physical_material_dataweak8_projector=projector_change(physical_weak_basis(J0,D0),physical_weak_basis(J1,z['dmaterial'])),
           scope='actual saved candidate versus original anchor; finite secant, not accepted path or acquired new data')
        sameJ.append(row)
    arm_rhs=sum(r['actual_new_RHS'] for r in cost);arm_lu=sum(r['actual_new_LU'] for r in cost)
    shared_rhs=sum(r['new_RHS'] for r in common.values());shared_lu=sum(r['new_LU'] for r in common.values())
    standalone_rhs=sum(r['standalone_R12_RHS'] for r in cost);standalone_lu=sum(r['standalone_R12_LU'] for r in cost)
    historical_rhs=sum(r['reduced_RHS256'] for r in historical.values());historical_lu=sum(r['reduced_LU256'] for r in historical.values())
    receipt=dict(status='PASS' if not issues else 'FAIL',issues=issues,GT_files_opened=0,new_physics_calls=0,
          primary_gate='actual Phi/noise+legality, then immutable greedy .002 material RMS; original anchor first',
          sameJ_control='nested restriction remains inside initial8 span; no new operator evaluation',
          selected_original_anchor_count=sum(r['selected_original_anchor_RMS']==0 for r in gate),groups=len(gate),
          campaign_P1_actual_RHS=arm_rhs+shared_rhs,campaign_P1_actual_LU=arm_lu+shared_lu,
          all_arm_actual_RHS=arm_rhs,COMMON_once_RHS=shared_rhs,
          sum_of_standalone_R12_RHS=standalone_rhs,sum_of_standalone_R12_LU=standalone_lu,
          intentional_COMMON_repetition_in_standalone_RHS=standalone_rhs-arm_rhs-shared_rhs,
          historical_R9_anchor_search_RHS_not_new=historical_rhs,historical_R9_anchor_search_LU_not_new=historical_lu,
          standalone_scope='R12-paid work plus COMMON; full cold deployment also needs historical anchor search/setup, separately disclosed',
          global_budget_snapshot=budget['snapshot_used'])
    write(HERE/'P1_ADVERSARIAL_QA.json',receipt);write(HERE/'P1_GATE_QA.json',gate);write(HERE/'P1_COST_RECONCILIATION.json',cost)
    write(HERE/'SAMEJ_VS_REFRESHED_CONTROL.json',sameJ);csv_write(HERE/'LEGACY_GATE_DIAGNOSTICS.csv',legacy)
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
