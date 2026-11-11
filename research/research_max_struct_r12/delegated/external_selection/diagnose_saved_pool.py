"""Posthoc exposed saved-pool selection diagnostic. No physical operator imports."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,time,json,numpy as np
HERE=Path(__file__).absolute().parent;REPO=HERE.parents[3];sys.path.insert(0,str(REPO))
from research.research_z_alg_r4_diversity.src.physics import Regularizer
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
CONFIG=dict(status='POSTHOC_EXPOSED_SELECTION_DIAGNOSTIC',prior_GT_exposure='externalv1 and fullDual already scored; this is not blinded or confirmatory',data_ball_relative_delta=.12,source='unchanged external/Austria predetermined raw-data discrepancy ball',selection_rule='Among actual saved paid candidates with reconstructed raw relative data residual <=.12, minimize exact inheritedR; deterministic first chronological record resolves ties',original_selection='Registered minimum original physical Phi final endpoint remains unchanged',oracle='OFFLINE reference only: minimum material NRMSE in actual saved pool; also report eligible-pool oracle. Never becomes an algorithm prediction',forbidden='No new forwards/inversion/ball/rank/prior/box tuning; no filling missing iterates',objective='Phi=.5||Ypred-Y||_F²/||Y||_F²+R(m); R uses same smoothedTV and quadratic neighbor prior smooth=tv=1e-4,epsilon=.01,32x32',tolerance='Eligibility uses exact recorded floating value residual<=.12, no extra tolerance; tiny negative reconstructed squaredresidual only clipped after gross-error check')
start=time.process_time();cid='austria_eps4_400MHz';reg=Regularizer((32,32));sources=[('direct','external','direct'),('Dual','external','Dual'),('AP','external','AP'),('DR','external','DR'),('Dual_full','external_continued','Dual')]
if (HERE/'SELECTION_SEAL.json').exists():raise FileExistsError('Do not overwrite diagnostic')
dump(HERE/'CONFIG.json',CONFIG);pools=[];selections=[];selected_arrays={}
for name,cohort,method in sources:
    out=HERE.parent/cohort/'RESULTS'/cid/method;r=json.loads((out/'RECEIPT.json').read_text())
    with np.load(out/'derived.npz')as z:A=z['candidates'].copy();final=z['m'].copy();raw=z['raw'].copy()
    h=r['history'];phi=np.array([x['Phi']for x in h],float)
    if A.shape!=(len(h),1024)or not np.isfinite(A).all()or not np.isfinite(phi).all():raise ValueError('Pool shape/finite failure')
    R=np.array([reg.evaluate(m)[0]for m in A]);d2=2*(phi-R)
    if d2.min() < -1e-10:raise ValueError('Recorded Phi incompatible with inheritedR')
    residual=np.sqrt(np.maximum(d2,0));eligible=np.flatnonzero(residual<=.12);chosen=int(eligible[np.argmin(R[eligible])])if len(eligible)else None
    original_matches=np.flatnonzero(np.all(A==final,axis=1));original_index=int(original_matches[0])if len(original_matches)else None
    if original_index is None:raise ValueError('Actual registeredfinal missing from saved pool; do not fabricate pair')
    original_phi=float(phi[original_index]);pool_min=int(np.argmin(phi))
    if abs(original_phi-r['original_Phi'])>1e-11:raise ValueError('Final objective mismatch')
    def item(i):
        if i is None:return None
        return dict(index=int(i),label=h[i]['label'],recorded_Phi=float(phi[i]),R=float(R[i]),relative_data_residual=float(residual[i]),eligible=bool(residual[i]<=.12),CPU_seconds=h[i]['CPU_seconds'])
    row=dict(method=name,cohort=cohort,actual_pool_count=len(A),duplicates_retained=True,counts_by_label={k:sum(x['label']==k for x in h)for k in sorted(set(x['label']for x in h))},pool_min_Phi=item(pool_min),original_final=item(original_index),minimum_R_eligible=item(chosen),eligible_count=len(eligible),coverage=dict(Phi=[float(phi.min()),float(phi.max())],R=[float(R.min()),float(R.max())],relative_data_residual=[float(residual.min()),float(residual.max())]),missing_pool_scope='Unselected Austria native-Dual stage arrays and all unaudited inner iterates are absent; only selected auxiliary plus recorded physical corrector candidates available'if method=='Dual'else'AP/DR physical audits only every50projectionsteps, plus recorded physical corrector candidates'if name in ['AP','DR']else'Recorded actual objective-call candidates; no invented intermediates')
    selections.append(row);selected_arrays[name+'_original']=final
    if chosen is not None:selected_arrays[name+'_minR']=A[chosen].copy()
    pools.append(dict(name=name,A=A,phi=phi,R=R,residual=residual,eligible=eligible,chosen=chosen,original=original_index,raw=raw))
    np.savez_compressed(HERE/(name+'_POOL_METRICS.npz'),Phi=phi,R=R,relative_data_residual=residual,eligible=residual<=.12)
np.savez_compressed(HERE/'SELECTED_MATERIALS.npz',**selected_arrays)
dump(HERE/'SELECTION_SEAL.json',dict(status='PUBLIC_PHI_R_NOISEBALL_SELECTIONS_FIXED_BEFORE_ORACLE_COMPUTATION',qualification='This task was already exposed to prior GT; seal is ordering of this diagnostic, not independent blinded confirmation',rows=selections))
truthpath=REPO.parent/'r8_worktree/research/research_max_struct_r8/data/paper_benchmarks_private'/(cid+'.npz');truth=np.load(truthpath)['m'];norm=np.linalg.norm(truth);oracle_arrays={};finalrows=[]
for pool,row in zip(pools,selections):
    errors=np.linalg.norm(pool['A']-truth,axis=1)/norm;oi=int(np.argmin(errors));ei=pool['eligible'];oe=int(ei[np.argmin(errors[ei])])if len(ei)else None
    def oracle(i):
        if i is None:return None
        return dict(index=i,NRMSE=float(errors[i]),Phi=float(pool['phi'][i]),R=float(pool['R'][i]),relative_data_residual=float(pool['residual'][i]),eligible=bool(pool['residual'][i]<=.12),scope='OFFLINE_ORACLE_ONLY_NOT_PREDICTION')
    row['original_final']['offline_NRMSE']=float(errors[pool['original']])
    if pool['chosen']is not None:row['minimum_R_eligible']['offline_NRMSE']=float(errors[pool['chosen']])
    row.update(offline_all_pool_oracle=oracle(oi),offline_eligible_pool_oracle=oracle(oe),saved_raw_offline_NRMSE=float(np.linalg.norm(pool['raw']-truth)/norm));finalrows.append(row);oracle_arrays[pool['name']+'_oracle']=pool['A'][oi]
    if oe is not None:oracle_arrays[pool['name']+'_eligible_oracle']=pool['A'][oe]
np.savez_compressed(HERE/'OFFLINE_ORACLE_MATERIALS.npz',**oracle_arrays)
dump(HERE/'RESULTS.json',dict(status='POSTHOC_EXPOSED_SELECTION_DIAGNOSTIC_COMPLETE',rows=finalrows,new_physical_solves=0,new_inversions=0,registered_endpoints_modified=False,CPU_seconds=time.process_time()-start,statistical_scope='Five dependent finite trajectory pools, duplicate evaluations retained. No success probability, confidence interval, general certificate or out-of-sample claim. Ball .12 is inherited discrepancy threshold on noisy observations, not known iidCN noise certificate. Feasible here means saved rawdata residual satisfies this scalar threshold, not material correctness/global inverse feasibility.'))
lines=['# Saved Austria pool: exposed posthoc selection diagnostic','', 'No new forward/inversion was run. Actual saved candidate arrays were paired one-to-one with recorded physical Phi values. Data residual is reconstructed only from the exact identity `r²=2(Phi−R)`, using inherited32² regularizer smooth=tv=1e-4,epsilon=.01. Eligibility is the unchanged rawdata relative ball r<=.12. The alternative rule minimizes inheritedR among eligible paid candidates, with first recorded candidate resolving ties. Registered final images remain unchanged.','', 'This is POSTHOC_EXPOSED_SELECTION_DIAGNOSTIC after prior GT scores, not a confirmatory algorithm or parameter search. The Phi/R/ball choices were saved before this script computes material oracles; prior exposure remains explicit. Offline all-pool and eligible-pool oracles only describe available material quality and are never predictions.','', '| Method | Saved / eligible | Original NRMSE | minR eligible NRMSE | Offline pool oracle NRMSE (residual) | Offline eligible oracle NRMSE |','|---|---:|---:|---:|---:|---:|']
for r in finalrows:
    mi=r['minimum_R_eligible'];oo=r['offline_all_pool_oracle'];eo=r['offline_eligible_pool_oracle'];a='NA'if mi is None else f"{mi['offline_NRMSE']:.6f}";b='NA'if eo is None else f"{eo['NRMSE']:.6f}"
    lines.append(f"|{r['method']}|{r['actual_pool_count']}/{r['eligible_count']}|{r['original_final']['offline_NRMSE']:.6f}|{a}|{oo['NRMSE']:.6f} ({oo['relative_data_residual']:.6f})|{b}|")
lines+=['', 'Coverage: every direct objective-call candidate; AP/DR physical audits every50steps plus physical-corrector objective calls; native Dual selected auxiliary and physical-corrector candidates. Unselected native-Dual auxiliary stage arrays were not saved by the caller and cannot be recovered from metadata. No missing stage/inner iterate was fabricated or replayed. Duplicate points are retained as recorded paid evaluations, not independent samples. Raw/final selected arrays are separately preserved; originals were verified present in each actual candidate pool.','', 'A minR result or oracle is limited to its saved finite pool. A pool may contain good material outside the .12 discrepancy ball while containing no equally good eligible candidate. Conversely, a large eligible-pool/original error gap demonstrates that registered Phi selection missed a better material in that particular paid pool, not that minR universally finds it. No generic basin, noise, feasibility, calibration or statistical certificate follows.']
(HERE/'REPORT.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(dict(status='COMPLETE',rows=[dict(method=r['method'],count=r['actual_pool_count'],eligible=r['eligible_count'],original=r['original_final']['offline_NRMSE'],minR=None if r['minimum_R_eligible']is None else r['minimum_R_eligible']['offline_NRMSE'],oracle=r['offline_all_pool_oracle'],eligible_oracle=r['offline_eligible_pool_oracle'])for r in finalrows]),indent=2))
