"""Aggregate frozen receipts; offline GT labels only, zero physical solves."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ[key]='1'
import json,csv
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parents[1]
def read(name):return json.loads((P/name).read_text())
def med(x):return float(np.median(x))
def savecsv(name,rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with (P/name).open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def main():
    cases=list(csv.DictReader((P/'results/e1/cases.csv').open()));low=[r for r in cases if float(r['amplitude'])==.4];ids={r['task_id']for r in low}
    dimensions=list(csv.DictReader((P/'results/e1/dimensions.csv').open()));scale_table=[]
    for radius in (.01,.03,.1):
        for eps in (1e-5,1e-4,1e-3):
            rs=[r for r in dimensions if r['case']in ids and r['linearization']=='anchor'and int(r['sources'])==12 and float(r['RMS'])==radius and float(r['epsilon'])==eps]
            scale_table.append(dict(RMS=radius,epsilon=eps,n=len(rs),dimension_median=med([float(r['d_weak'])for r in rs]),soft_dimension_median=med([float(r['d_soft'])for r in rs])))
    savecsv('results/SCALE_TOLERANCE_TABLE.csv',scale_table)
    pca=[]
    for label,rs in [('Native_clean_12',[r for r in low if r['method']=='native_trf'and r['noise']=='clean']),('Native_all24',[r for r in low if r['method']=='native_trf']),('all48',low)]:
        E=np.stack([np.load(P/'results/e1'/(r['task_id']+'.npz'))['error_OFFLINE']for r in rs]);E-=E.mean(0);s=np.linalg.svd(E,compute_uv=False);c=np.cumsum(s*s)/np.sum(s*s)
        pca.append(dict(cohort=label,n_rows=len(rs),independent_materials=12,centered=True,sample_rank=int((s>s[0]*1e-12).sum()),r90=int(np.searchsorted(c,.9)+1),r95=int(np.searchsorted(c,.95)+1),r99=int(np.searchsorted(c,.99)+1),singular_values=s.tolist(),scope='OFFLINE_SAMPLE_STATISTIC_NOT_PHYSICAL_FREEDOM'))
    (P/'results/SAMPLE_ERROR_PCA.json').write_text(json.dumps(pca,indent=2)+'\n')
    z=np.load(P.parent/'research_gn_green_error_20261106/results/fixed_spectra.npz');fixed=[]
    for k in ('in_s','ob_s','real_in_s','real_ob_s'):
        s=z[k];rank=int((s>s.max()*1e-12).sum());strong=int((s>=s.max()*.1).sum());fixed.append(dict(operator=k,domain_dimension=576,numerical_rank=rank,strong=strong,weak_nonzero=rank-strong,complement=576-rank,scope='REUSED_FIXED_SPECTRUM_NOT_MATERIAL_OBSERVABILITY'))
    savecsv('results/FIXED_OPERATOR_DIMENSIONS.csv',fixed)
    e2=read('results/e2/SCORES.json')['rows'];e3=read('results/e3/SCORES.json')['rows']
    comp=read('results/e3/COMPLETION.json');geo=read('results/GEOMETRY_SUMMARY.json')
    comparisons=[];flat=[];transverse=[];errors=[]
    for arm in ('continue_trf','haar_unfiltered','green_pullback','soft_jacobian'):
        rs=[r for r in e3 if r['arm']==arm]
        comparisons.append(dict(arm=arm,n=len(rs),top1_to_anchor_NRMSE_median=med([r['NRMSE']/r['anchor_metrics']['NRMSE']for r in rs]),bestK_to_anchor_NRMSE_median=med([r['oracle_bestK_NRMSE_OFFLINE']/r['anchor_metrics']['NRMSE']for r in rs]),top1_NRMSE_median=med([r['NRMSE']for r in rs]),top1_improves_count=sum(r['top1_delta_NRMSE'] < -1e-8 for r in rs),data_qualified=sum(r['strict_data_qualification']for r in rs),joint_qualified=sum(r['historical_r4_material_gate']for r in rs),wall_median=med([r['wall_seconds']for r in rs]),LU_median=med([r['costcounts']['LU']for r in rs]),RHS_median=med([r['costcounts']['actual_complex_RHS']for r in rs]),attempted_seeds=sum(r['attempted_seed_count']for r in rs),nonzero_seeds=sum(r['nonzero_seed_count_gt_1e_12']for r in rs)))
        qr=[];nq=0;nt=0
        for row in rs:
            record=read('results/e3/'+row['task_id']+'.json');rec=record.get('result',{});limit=.001 if row['noise']=='clean'else .003
            errs=[row['anchor_metrics']['NRMSE'],row['NRMSE']]
            for phase in ('before','after'):
                for meta,score in zip(rec.get('proposals',[]),row['candidate_'+phase+'_metrics']):
                    nt+=1
                    if meta[phase]['relative_residual']<=limit:errs.append(score['NRMSE']);nq+=1
            qr.append(min(errs)/row['anchor_metrics']['NRMSE'])
        comparisons[-1].update(qualified_bestK_to_anchor_median=med(qr),qualified_before_after_slots=nq,before_after_slots=nt)
    for r in e2+e3:
        row={k:r.get(k)for k in ('task_id','experiment','shape','noise','arm','reg_scale','stopping_profile','RMSE','NRMSE','IoU_Otsu','wall_seconds','cpu_seconds','oracle_bestK_NRMSE_OFFLINE','oracle_candidate_count','strict_data_qualification','historical_r4_material_gate','reused_from')}
        row.update(anchor_NRMSE=r['anchor_metrics']['NRMSE'],relative_residual=r['top1_online']['relative_residual'],LU=r['costcounts']['LU'],complex_RHS=r['costcounts']['actual_complex_RHS']);flat.append(row)
        if r['experiment']=='e2':
            with np.load(P/'results/e1'/(r['anchor_task_id']+'.npz'))as a,np.load(P/'results/e2'/(r['task_id']+'.npz'))as z:
                V=a['J_V'];s=a['J_s']/np.linalg.norm(a['Y']);strong=s*.03*24>1e-4;e=a['error_OFFLINE']-(z['top1']-a['anchor']);e0=a['error_OFFLINE']
                errors.append(dict(task_id=r['task_id'],shape=r['shape'],reg_scale=r['reg_scale'],stopping_profile=r.get('stopping_profile','frozen'),strong_dim=int(strong.sum()),before_strong_energy=float(np.linalg.norm(V[:,strong].T@e0)**2),after_strong_energy=float(np.linalg.norm(V[:,strong].T@e)**2),before_weak_energy=float(np.linalg.norm(V[:,~strong].T@e0)**2),after_weak_energy=float(np.linalg.norm(V[:,~strong].T@e)**2)))
        elif r['arm']!='continue_trf':
            with np.load(P/'results/e3'/(r['task_id']+'.npz'))as z:
                V=z['anchor_J_right_vectors'];s=z['anchor_J_singular_values'];eps=1e-4 if r['noise']=='clean'else 1e-3;strong=s*.03*24>eps;anchor=z['anchor']
                for k,(before,after)in enumerate(zip(z['proposal_before'],z['proposal_after'])):
                    h=before-anchor;g=after-anchor
                    bn=np.linalg.norm(V[:,~strong].T@h);an=np.linalg.norm(V[:,~strong].T@g)
                    bs=np.linalg.norm(V[:,strong].T@h);aas=np.linalg.norm(V[:,strong].T@g)
                    transverse.append(dict(task_id=r['task_id'],shape=r['shape'],noise=r['noise'],arm=r['arm'],candidate=k,seed_RMS=float(np.linalg.norm(h)/24),after_RMS=float(np.linalg.norm(g)/24),weak_before=float(bn),weak_after=float(an),weak_retention=float(an/bn)if bn>1e-10 else '',strong_before=float(bs),strong_after=float(aas),strong_retention=float(aas/bs)if bs>1e-10 else '',direction_cosine=float(h@g/(np.linalg.norm(h)*np.linalg.norm(g)))if np.linalg.norm(h)*np.linalg.norm(g)>1e-20 else ''))
    ledger=[]
    for name,d in comp['groups'].items():
        ledger.append(dict(stage=name,LU=d['LU'],complex_RHS=d['complex_RHS'],wall_seconds=d.get('solver_wall_seconds',d.get('wall_seconds')),CPU_seconds=d.get('cpu_seconds'),scope='sum_unique_receipts',status='MEASURED',notes='12 E3 clean continuation receipts reuse E2; counted only once'if name=='e3_main'else''))
    a=read('results/e1/RUN_SUMMARY.json');b=read('results/e1/SUPPLEMENTAL_SUMMARY.json');f=read('results/finite/RUN_SUMMARY.json')
    ledger += [dict(stage='e1_geometry',LU=a['counts']['LU'],complex_RHS=a['counts']['actual_complex_RHS'],wall_seconds=a['successful_case_wall_seconds'],scope='successful_case_wall_and_known_counts',status='LOWER_BOUND',notes='failed/redundant preparation not recoverable'),dict(stage='e1_current_subsets',LU=b['counts']['LU'],complex_RHS=b['counts']['actual_complex_RHS'],wall_seconds=b['sum_process_wall_seconds'],scope='sum_child_numerical_wall',status='MEASURED'),dict(stage='finite_and_amplitude_controls',LU=f['cost']['total']['LU'],complex_RHS=f['cost']['total']['actual_complex_RHS'],wall_seconds=f['cost']['total_wall_seconds'],scope='one_process_numerical_wall',status='MEASURED'),dict(stage='native_failed_uncheckpointed_attempts',LU='',complex_RHS='',wall_seconds='',scope='lost native process state',status='UNKNOWN_NOT_ZERO'),dict(stage='historical_fixed_G_S_spectra',LU=0,complex_RHS=0,wall_seconds='',scope='reused immutable offline artifact',status='REUSED_NOT_FREE_COLD_SETUP',notes='No repeated Green SVD; original setup timing not isolated; no end-to-end speedup claim')]
    summary=dict(status='COMPLETE_BOUNDED_DEVELOPMENT_STUDY',scope='Historical disclosed objects; no fresh blind confirmation; no NN training',independent_low_contrast_materials=12,e1_endpoints=52,source_subset_rows=208,e2_main_predictions=36,e3_predictions=96,tight_diagnostic_predictions=12,finite_GT_directed_trials=432,amplitude_controls=6,geometry=geo,search_comparison_24= comparisons,e2_grouped=read('results/e2/SCORES.json')['summary'],cost_known_lower_bound=dict(LU=sum(int(r['LU'])for r in ledger if r['LU']!=''),complex_RHS=sum(int(r['complex_RHS'])for r in ledger if r['complex_RHS']!='')),failure_costs='UNKNOWN native lost attempts; all 8 retained corrector failures included',decision=dict(small_weak_latent='NOT_SUPPORTED: first64 covers only about3percent',online_search_superiority='NO_GAIN_SHOWN',NN_training='NOT_RUN',full_TSOM_inverse_ablation='DEFERRED',transverse='EMPIRICAL_TRF_ONLY; local damped-GN theorem distinct'))
    (P/'results/RESULTS_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    for name,rows in [('SEARCH_COMPARISON.csv',comparisons),('RESULTS_BY_CASE.csv',flat),('TRANSVERSE_RETENTION.csv',transverse),('E2_ERROR_SUBSPACE_CHANGE.csv',errors),('COST_LEDGER.csv',ledger)]:savecsv('results/'+name,rows)
    print(json.dumps(dict(comparisons=comparisons,cost=summary['cost_known_lower_bound']),indent=2))
if __name__=='__main__':main()
