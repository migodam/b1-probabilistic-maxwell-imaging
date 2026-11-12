"""Post hoc mechanism audit of frozen MAX-STRUCT R12 paired evidence.

No optimization, truth-based selection, or old source modification. Uses saved
operators, physical states and exact saved DPCS bases. New work is a decomposition
of the initialization projection, rather than an algebraic error-floor test.
"""
from pathlib import Path
import os
import json, time, sys
import numpy as np
import pandas as pd
from scipy.linalg import svd
from scipy.stats import spearmanr, pearsonr, chi2
from threadpoolctl import threadpool_limits

ROOT=Path(os.environ.get('MAX_STRUCT_RUN_ROOT', str(Path(__file__).resolve().parent.parent)))
REPO=Path(os.environ.get('MAX_STRUCT_REPO', str(ROOT/'research_repo')))
BASE=REPO/'research/research_max_struct_r12/supplement_paired'
OUT=ROOT/'work/evidence_analysis'
OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(REPO))
from research.research_max_struct_r6.src.fft_tsom_ops import NativeFourierAPCS

def norm2(a): return float(np.vdot(a,a).real)

def main():
    t=time.process_time()
    scores=pd.read_csv(BASE/'RESULTS/offline_v2/SCORES.csv')
    defects=pd.read_csv(BASE/'RESULTS/offline_v2/G19_INITIAL_PROJECTION_DEFECT.csv')
    costs=pd.read_csv(BASE/'RESULTS/offline_v2/CHILD_COST_LEDGER.csv')
    primary=scores.loc[scores['view'].eq('primary')]
    initial=scores.loc[scores['view'].eq('initial') & scores['method'].eq('g19')]
    pivot=primary.pivot(index=['cohort','case_id'],columns='method',values='NRMSE').reset_index()
    joined=defects.merge(pivot,on=['cohort','case_id']).merge(initial[['cohort','case_id','NRMSE','residual_over_delta','lower_active','upper_active']],on=['cohort','case_id'])
    joined=joined.rename(columns={'NRMSE':'initial_NRMSE','residual_over_delta':'initial_residual_over_delta'})
    joined['g19_improvement']=joined['initial_NRMSE']-joined['g19']
    joined['dr_advantage']=joined['g19']-joined['dr']
    joined.to_csv(OUT/'projection_outcomes.csv',index=False)
    correlations=[]
    for group in ['ALL','AB','C']:
        z=joined if group=='ALL' else joined[joined.cohort.isin(list(group))]
        for predictor in ['relative_current_loss','constitutive_defect','relative_port_change']:
            for outcome in ['g19','g19_improvement','dr_advantage']:
                v=spearmanr(z[predictor],z[outcome])
                correlations.append(dict(group=group,n=len(z),predictor=predictor,outcome=outcome,spearman=float(v.statistic),p_descriptive_only=float(v.pvalue)))
    pd.DataFrame(correlations).to_csv(OUT/'projection_correlations_exploratory.csv',index=False)
    rows=[]
    for rec in defects.to_dict('records'):
        task=rec['task_id']; case=rec['case_id']; cohort=rec['cohort']
        d=BASE/'RESULTS'/task
        with np.load(BASE/'inputs'/f'{cohort}__{case}.npz',allow_pickle=False) as z:
            S=z['S']; G=z['G']; E=z['E0']; Y=z['Y']; shape=tuple(z['grid_shape']); m=z['m_init']
            delta=float(z['sigma'])*np.sqrt(chi2.ppf(.99,2*Y.size)/2)
        with np.load(d/'initial.npz',allow_pickle=False) as z: J=z['j']; Q=z['Q']
        with np.load(d/'RECEIVER_DPCS.npz',allow_pickle=False) as z: Vd=z['Vd']; JD=z['JD']; rank=int(z['receiver_rank'])
        U=NativeFourierAPCS(shape,3,Vd)
        UU=U.apply(np.eye(U.k,dtype=complex))
        Uo,s,_=svd(UU,full_matrices=False,check_finite=False)
        Uo=Uo[:,s>s[0]*1e-12]
        pdj=Vd@(Vd.conj().T@J)
        puj=Uo@(Uo.conj().T@J)
        fixed_jump=JD-pdj
        truncated=J-pdj-puj
        deltaJ=fixed_jump-truncated
        Jp=J+deltaJ; Qp=E+G@Jp
        denom=norm2(J); total=norm2(deltaJ)
        # Shared-real material defect component: orthogonal rowwise projection
        # of the constitutive residual onto directions h_i Q_i at the START.
        # This is a descriptive local component, not a proposed success predictor.
        residual=Jp-m[:,None]*Qp
        qenergy=np.sum(abs(Q)**2,axis=1)
        hm=np.real(np.sum(np.conj(Q)*residual,axis=1))/np.maximum(qenergy,1e-300)
        aligned=hm[:,None]*Q
        exec_record=json.loads((d/'G19_EXECUTION.json').read_text())
        external_cost=float(costs.loc[costs.task_id.eq(task),'external_whole_child_CPU_seconds'].iloc[0])
        count=exec_record.get('operation_counts',{})
        row=dict(task_id=task,cohort=cohort,case_id=case,receiver_rank=rank,initial_mf=3,initial_fourier_columns=U.k,
            final_mf_reached=12 if count.get('PR_CG_steps',0)>1500 else (10 if count.get('PR_CG_steps',0)>1000 else 6 if count.get('PR_CG_steps',0)>500 else 3),
            PR_steps=count.get('PR_CG_steps'),initial_current_norm=np.sqrt(denom),projection_loss=np.sqrt(total/denom),
            fixed_DPCS_jump_relative=np.sqrt(norm2(fixed_jump)/denom),fourier_truncation_relative=np.sqrt(norm2(truncated)/denom),
            fixed_DPCS_fraction=norm2(fixed_jump)/total,fourier_truncation_fraction=norm2(truncated)/total,
            fixed_DPCS_port_change_relative=np.sqrt(norm2(S@fixed_jump)/norm2(Y)),fourier_port_change_relative=np.sqrt(norm2(S@truncated)/norm2(Y)),
            total_port_change_relative=np.sqrt(norm2(S@deltaJ)/norm2(Y)),
            start_port_residual_relative=np.sqrt(norm2(S@J-Y)/norm2(Y)),
            projected_port_residual_relative=np.sqrt(norm2(S@Jp-Y)/norm2(Y)),
            shared_material_aligned_defect_fraction=norm2(aligned)/norm2(residual),
            whole_child_CPU_seconds=external_cost)
        # At eventual complete Fourier stage, only the fixed DPCS jump remains.
        for method in ['g19','dr','ap','direct']:
            for state in ['raw','primary']:
                with np.load(BASE/'RESULTS'/f'{cohort}__{case}__{method}'/f'{state}.npz',allow_pickle=False) as z:
                    j=z['j']; yp=z['Y_prediction']
                jump=JD-Vd@(Vd.conj().T@j)
                key=f'{method}_{state}'
                row[f'{key}_fixed_DPCS_current_relative']=np.sqrt(norm2(jump)/norm2(j))
                row[f'{key}_fixed_DPCS_port_over_delta']=np.sqrt(norm2(S@jump))/delta
                row[f'{key}_residual_over_delta']=np.sqrt(norm2(yp-Y))/delta
        rows.append(row)
    split=pd.DataFrame(rows)
    split.to_csv(OUT/'projection_split.csv',index=False)
    merged=joined.merge(split,on=['task_id','cohort','case_id'])
    merged.to_csv(OUT/'projection_split_outcomes.csv',index=False)
    # Corrector analysis excludes 20 bypass cases repeated under method names.
    corr=pd.read_csv(BASE/'RESULTS/offline_v2/CORRECTOR_DECOMPOSITION.csv')
    actual=costs[['cohort','case_id','method']]
    actual_corr=corr.merge(actual,on=['cohort','case_id','method'])
    actual_corr.to_csv(OUT/'actual_corrector_cases.csv',index=False)
    harms=actual_corr.loc[actual_corr.delta_NRMSE>1e-8]
    harms.to_csv(OUT/'corrector_harms.csv',index=False)
    summary={'frozen_commit':'ec3a37da6f3ee90beb61e753067d4f639f4118f5','projection_cases':len(split),
      'all_reached_full_fourier':bool((split.final_mf_reached==12).all()),
      'fourier_loss_fraction_min':float(split.fourier_truncation_fraction.min()),
      'fourier_loss_fraction_max':float(split.fourier_truncation_fraction.max()),
      'fourier_loss_fraction_median':float(split.fourier_truncation_fraction.median()),
      'actual_corrector_cases':len(actual_corr),'corrector_harms_cases':len(harms),
      'DR_corrector_harms_cases':int((harms.method=='dr').sum()),
      'CPU_seconds':time.process_time()-t,'new_forward_solves':0,
      'scope':'Retrospective exposed cases. Correlations are exploratory, not validation of a predictor.'}
    (OUT/'audit_summary.json').write_text(json.dumps(summary,indent=2))
    print(split[['case_id','projection_loss','fixed_DPCS_fraction','fourier_truncation_fraction','fixed_DPCS_port_change_relative','fourier_port_change_relative','PR_steps','dr_raw_fixed_DPCS_current_relative','dr_primary_fixed_DPCS_current_relative']].to_string(index=False))
    print(json.dumps(summary,indent=2))
    print('Corrector harms:')
    print(harms[['cohort','case_id','method','raw_NRMSE','corrected_NRMSE','delta_D','delta_R','raw_residual_over_delta','corrected_residual_over_delta']].to_string(index=False))

if __name__=='__main__':
    with threadpool_limits(limits=1): main()
