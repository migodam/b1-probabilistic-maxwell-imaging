#!/usr/bin/env python3
"""Independent read-only risk arithmetic after individual scoring Books close.

No model/scorer/solver imports. Independent complex-error arithmetic and paired
parent bootstrap; no source results are rewritten and no model is rerun.
"""
from pathlib import Path
import json,os,sys,time,csv,resource
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'tools')]
from p3r4.evaluation import checkpoint_identity,write_new_json
from r4_fresh_risk import component_check
S=ROOT/'research/p3l_r4_conditional_transport';F=S/'results/fresh';OUT=F/'independent_risk_statistical_audit_v1'
SEED=20261025099;DRAWS=4000;LEVEL=.9833333333333333
PRIMARY_REFS=['frozen_r2','strong_spectral','fusion_diag','fusion_permuted']
def read(p):return json.loads(Path(p).read_text())
def check(v):
 if checkpoint_identity(v['path'])!=v:raise ValueError('source identity changed: '+v['path'])
def csvwrite(path,rows):
 with path.open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def bootstrap(cand,ref):
 rng=np.random.default_rng(SEED);a=[];r=[]
 for start in range(0,DRAWS,100):
  idx=rng.integers(0,len(cand),size=(min(100,DRAWS-start),len(cand)))
  cm=cand[idx].mean(1);rm=ref[idx].mean(1);a.extend(rm-cm);r.extend((rm-cm)/rm)
 alpha=(1-LEVEL)/2
 return dict(n_parent=len(cand),candidate_mean=float(cand.mean()),reference_mean=float(ref.mean()),improvement_absolute=float((ref-cand).mean()),improvement_relative=float((ref.mean()-cand.mean())/ref.mean()),absolute_CI=np.quantile(a,[alpha,1-alpha]).tolist(),relative_CI=np.quantile(r,[alpha,1-alpha]).tolist(),CI_level=LEVEL,draws=DRAWS,seed=SEED,unit='independent parent; seed losses averaged before resampling')
def main():
 wall,cpu=time.monotonic(),time.process_time()
 if OUT.exists():raise FileExistsError('independent audit output already exists; no overwrite')
 release_path=S/'configs/FIRST_HELDOUT_GT_SCORE_RELEASE.json';release=read(release_path);releaseid=checkpoint_identity(release_path)
 if release.get('released') is not True or release.get('experiment')!='fresh_r4_isolated_scores_v1' or release.get('all_online_both_cohorts_closed_before_truth') is not True or release.get('new_models_after_truth') is not False:raise ValueError('explicit firstGT authorization missing')
 for v in release['closed_gates']:
  check(v)
  if v['mtime_ns']>releaseid['mtime_ns']:raise ValueError('prerequisite seal postdates authorization')
 comp_path=S/'AUDIT/risk/fresh_freeze_v3/RISK_COMPONENT_FROZEN.json';comp=read(comp_path);component_check(comp)
 gates={};sources=[releaseid,checkpoint_identity(comp_path)];failures=[]
 # Both scoring Books, complete source seals and chronological checks precede
 # the FIRST truth access by this audit, even though panels are now EXPOSED.
 for group,n,prefix in [('locked',512,'R4L1'),('confirmation',256,'R4L2')]:
  base=F/group/'risk45';d=S/'IV/jobs'/('fresh_risk_'+group+'_score_v1');closed=read(d/'CLOSED.json');started=read(d/'START.json');sc=read(base/'scores/RESULTS.json');g=read(base/'online/ONLINE_COMPLETE.json');bound=read(F/group/'BOUND_RISK.json');top=read(S/'configs'/('FROZEN_'+group.upper()+'.json'))
  parents=[f'{prefix}_{i:05d}' for i in range(n)]
  if closed['status']!='COMPLETE' or started['stage']!='IV' or started['kind']!='analysis' or closed['counters'].get('RHS',0)!=0 or int(started['start_unix']*1e9)<releaseid['mtime_ns']:raise ValueError('individual scoreBook/authorization gate failed')
  if sc['status']!='OFFLINE_SCORE_COMPLETE' or sc['scoring_pid']==os.getpid() or g['pid']==os.getpid() or sc['scoring_pid']==g['pid']:raise ValueError('independent scoring/inference PID gate failed')
  if not g['completed_time_ns']<releaseid['mtime_ns']<int(started['start_unix']*1e9) or (base/'scores/RESULTS.json').stat().st_mtime_ns>(d/'CLOSED.json').stat().st_mtime_ns:raise ValueError('GT chronology failed')
  if g['parents']!=parents or g['freeze']!=bound or bound['risk_rule']!=comp['risk_rule'] or top['risk_rule']!=comp['risk_rule'] or len(bound['methods'])!=45 or bound['adapt_after_L1'] is not False:raise ValueError('frozen roster/risk rule altered')
  if [(m['id'],m.get('checkpoint_identity')) for m in bound['methods']]!=[(m['id'],m.get('checkpoint_identity')) for m in comp['methods']]:raise ValueError('45 frozen methods changed')
  for identity in g['sealed_artifact_identities']+g['source_identities']:check(identity)
  if any(v['mtime_ns']<=bound['freeze_time_ns'] for v in read(base/'features/FEATURES_COMPLETE.json')['observation_sources']):raise ValueError('observation before binding')
  sources.extend([checkpoint_identity(d/'CLOSED.json'),checkpoint_identity(d/'START.json'),checkpoint_identity(base/'scores/RESULTS.json'),checkpoint_identity(base/'online/ONLINE_COMPLETE.json'),checkpoint_identity(F/group/'BOUND_RISK.json')]);gates[group]=(base,bound,g,sc,parents)
 OUT.mkdir(parents=True);truth_open_ns=time.time_ns();allloss={};family_rows=[];method_rows=[];numeric=[];truth_ids=[];batch_checks=[]
 for group in ['locked','confirmation']:
  base,bound,g,sc,parents=gates[group];roi=np.array(g['geometry']['roi'],int)
  truth=[]
  for parent in parents:
   p=S/'inputs/offline_only'/(parent+'.npz');truth_ids.append(checkpoint_identity(p))
   with np.load(p,allow_pickle=False) as z:chi=z['chi']
   if chi.shape!=(576,) or not np.iscomplexobj(chi) or not np.isfinite(chi).all():raise ValueError('invalid isolated full576 truth')
   truth.append(chi)
  truth=np.array(truth,dtype=np.complex128);loss={};flags={}
  for m in bound['methods']:
   mid=m['id'];path=base/'online'/(mid+'.npz')
   with np.load(path,allow_pickle=False) as z:
    if z['parent'].astype(str).tolist()!=parents:raise ValueError('canonical parent roster changed')
    arrays={action:z['prediction_'+action] for action in ['raw','clipped']};fallback=z['fallback_anchor'].copy()
   # Independently show aggregate is exactly the concatenation of frozen64
   # batches, rather than regenerate any network prediction through replay.
   count=0
   for start in range(0,len(parents),64):
    with np.load(base/'online/batches'/mid/f'batch_{start//64:03d}.npz',allow_pickle=False) as z:
     if z['parent'].astype(str).tolist()!=parents[start:start+64] or not np.array_equal(z['fallback_anchor'],fallback[start:start+64]):raise ValueError('batch parent/fallback mismatch')
     for action in ['raw','clipped']:
      if not np.array_equal(z['prediction_'+action],arrays[action][start:start+64]):raise ValueError('aggregate differs from canonical stored batch')
    count+=1
   batch_checks.append(dict(cohort=group,method=mid,batches=count,bitwise_aggregate_equals_batches=True))
   if len(fallback)!=len(parents) or int(fallback.sum())!=len(g['failures'][mid]) or int(fallback.sum())!=sc['methods'][mid]['fallback_count']:raise ValueError('fallback denominator/status mismatch')
   flags[mid]=fallback
   if fallback.any():failures.append(dict(cohort=group,method=mid,failed_parents=np.array(parents)[fallback].tolist(),fallback='anchor; included'))
   loss[mid]={}
   with np.load(base/'scores'/(mid+'_per_parent.npz'),allow_pickle=False) as saved:
    if saved['parent'].astype(str).tolist()!=parents:raise ValueError('scorer perparent order changed')
    for action,pred in arrays.items():
     if pred.shape!=truth.shape or not np.isfinite(pred).all():raise ValueError('invalid canonical material; cannot silently drop parents')
     # Independent real/imag arithmetic; no call to original material_losses.
     d=np.asarray(pred,dtype=np.complex128)-truth;sq=d.real*d.real+d.imag*d.imag;mag=np.sqrt(sq)
     metrics=dict(full_mse=sq.mean(1),full_mae=mag.mean(1),roi_mse=sq[:,roi].mean(1),roi_mae=mag[:,roi].mean(1));loss[mid][action]=metrics
     for metric,value in metrics.items():
      error=float(np.max(np.abs(value-saved[action+'__'+metric])));numeric.append(dict(cohort=group,method=mid,action=action,metric=metric,max_absolute_difference=error))
      if error>1e-12:raise ValueError('independent loss arithmetic disagrees with scorer')
     method_rows.append(dict(cohort=group,method=mid,family=m['family'],seed=m['seed'],action=action,n_parent=len(parents),fallbacks=int(fallback.sum()),**{k:float(v.mean()) for k,v in metrics.items()}))
  families={}
  for family in dict.fromkeys(m['family'] for m in bound['methods']):
   methods=[m for m in bound['methods'] if m['family']==family]
   if len(set(m['seed'] for m in methods))!=len(methods):raise ValueError('duplicate seed within family')
   families[family]={}
   for action in ['raw','clipped']:
    metrics={k:np.stack([loss[m['id']][action][k] for m in methods]).mean(0) for k in ['full_mse','full_mae','roi_mse','roi_mae']};families[family][action]=metrics
    family_rows.append(dict(cohort=group,family=family,action=action,n_parent=len(parents),seeds=len(methods),**{k:float(v.mean()) for k,v in metrics.items()}))
  allloss[group]=families
  np.savez_compressed(OUT/(group+'_seed_mean_parent_losses.npz'),parent=np.array(parents),**{family+'__'+action+'__'+metric:value for family,actions in families.items() for action,metrics in actions.items() for metric,value in metrics.items()})
 comparisons=[];compare_errors=[]
 for panel in ['locked','confirmation','pooled768_descriptive']:
  for action in ['raw','clipped']:
   for ref in PRIMARY_REFS:
    if panel.startswith('pooled'):
     cand=np.concatenate([allloss[g]['fusion_spline'][action]['full_mse'] for g in ['locked','confirmation']]);reference=np.concatenate([allloss[g][ref][action]['full_mse'] for g in ['locked','confirmation']])
    else:cand=allloss[panel]['fusion_spline'][action]['full_mse'];reference=allloss[panel][ref][action]['full_mse']
    row=dict(panel=panel,candidate_family='fusion_spline',reference_family=ref,action=action,scope='DESCRIPTIVE_POOLED_NOT_NEW_CONFIRMATORY_TEST' if panel.startswith('pooled') else 'FIXED_COHORT_PRIMARY' if ref=='frozen_r2' and action=='clipped' else 'SUPPORTING_POINTWISE_INTERVAL_NOT_EXTRA_FAMILYWISE_CLAIM',**bootstrap(cand,reference));comparisons.append(row)
    if not panel.startswith('pooled'):
     sc=gates[panel][3];matches=[v for v in sc['comparisons'] if v['candidate_family']=='fusion_spline' and v['reference_family']==ref and v['action']==action]
     if len(matches)!=1:raise ValueError('expected fixed comparison missing')
     old=matches[0];error=max(float(np.max(np.abs(np.asarray(row[k])-np.asarray(old[k])))) for k in ['candidate_mean','reference_mean','improvement_absolute','improvement_relative','absolute_CI','relative_CI']);compare_errors.append(dict(panel=panel,reference=ref,action=action,max_absolute_difference=error))
     if error>1e-12:raise ValueError('independent bootstrap disagrees with scorer')
 for family in ['fusion_spline',*PRIMARY_REFS]:
  for action in ['raw','clipped']:
   metrics={k:np.concatenate([allloss[g][family][action][k] for g in ['locked','confirmation']]) for k in ['full_mse','full_mae','roi_mse','roi_mae']}
   family_rows.append(dict(cohort='pooled768_descriptive',family=family,action=action,n_parent=768,seeds=1 if family=='frozen_r2' else 3,**{k:float(v.mean()) for k,v in metrics.items()}))
 csvwrite(OUT/'ALL45_METHOD_ACTIONS.csv',method_rows);csvwrite(OUT/'FAMILY_SEED_MEAN_TABLE.csv',family_rows)
 csvwrite(OUT/'PRIMARY_FIXED_COMPARISONS.csv',[{**{k:v for k,v in row.items() if k not in ['absolute_CI','relative_CI']},'absolute_CI_low':row['absolute_CI'][0],'absolute_CI_high':row['absolute_CI'][1],'relative_CI_low':row['relative_CI'][0],'relative_CI_high':row['relative_CI'][1]} for row in comparisons])
 result=dict(status='INDEPENDENT_READ_ONLY_STATISTICAL_AUDIT_COMPLETE',failures_first=failures,fallback_failure_count=sum(len(v['failed_parents']) for v in failures),source_or_chronology_mismatches=0,primary_family='fusion_spline',primary_selection='fixed original R2VAL256 means, not heldout/DEV; no switching for L2',risk_rule=comp['risk_rule'],heldout_label_status='EXPOSED_AFTER_CLOSED_SCORING; both panels now exposed',truth_first_access_ns=truth_open_ns,release_identity=releaseid,score_source_identities=sources,truth_identities=truth_ids,batch_aggregate_checks=batch_checks,numeric_score_reconciliation=numeric,bootstrap_reconciliation=compare_errors,comparisons=comparisons,methods=45,cohorts=dict(locked=512,confirmation=256),combined_scope='768 pooled descriptive only, no new confirmatory hypothesis',statistics=dict(draws=DRAWS,seed=SEED,level=LEVEL,unit='parent; arithmeticmean of three seed losses, never ensemble image'),multiplicity='only original3 broad predeclared primary comparisons retain their original claim scope; extra raw/control/pooled CIs are pointwise descriptive; no additional familywise claim',model_reruns=0,source_outputs_modified=False,new_physics_RHS=0,new_fits=0,wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,peak_RSS_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024),role='independent auditor of stored arrays/statistics, not scorer author or scientific judge')
 write_new_json(OUT/'AUDIT.json',result)
 lines=['# 独立只读风险统计审计','',f"**失败优先：预测/回退 {result['fallback_failure_count']}，源身份或时间链不一致0。** 两个原评分Book均已CLOSED COMPLETE后，本审计才首次打开隔离truth。现L1与确认集标签状态都是EXPOSED；没有重跑模型、评分wrapper或Maxwell，没有改原输出。",'', '| 面板 | 参照（均为clipped） | Spline MSE | 参照 MSE | 改善 | 98.333333%父样本bootstrap区间 |','|---|---|---:|---:|---:|---:|']
 for row in comparisons:
  if row['action']=='clipped':lines.append(f"| {row['panel']} | {row['reference_family']} | {row['candidate_mean']:.10f} | {row['reference_mean']:.10f} | {100*row['improvement_relative']:.3f}% | [{100*row['relative_CI'][0]:.3f}%, {100*row['relative_CI'][1]:.3f}%] |")
 lines+=['', '从全部45方法的canonical stored64批次独立复核aggregate逐项相同；逐父样本全576复MSE/MAE与ROI MSE/MAE用独立real/imag算式重算，逐项与原评分对账。先在每父样本内平均三个seed损失，再以512/256独立父样本进行4000次bootstrap（seed20261025099）。没有平均seed图像，也没有把seed当额外父样本。数值差异及bootstrap对账逐项记录AUDIT.json。','', '主族始终是验证集预选Spline，确认集没有更换；方法/checkpoint/provider/source/风险规则身份与观测前冻结一致，未发现测试后调参。对Diag的区间跨0仅是本次差异不确定，不能称等价或独特Flow优势；最终解释归主线程。','', '池化768仅作描述，不是新增确认性检验。原3项预声明大类主比较之外，本报告额外raw/控制/池化区间只按单项描述，不新增familywise覆盖或显著性声明。','',f"本审计wall {result['wall_seconds']:.3f}s、CPU {result['cpu_seconds']:.3f}s、RSS {result['peak_RSS_bytes']/1e6:.1f}MB，0新物理RHS/0训练。身份检查采用冻结路径/字节数/mtime，不读取模型张量或CFM/sample日志。",'', '全部方法及raw/clipped均值：ALL45_METHOD_ACTIONS.csv；seed损失均值表：FAMILY_SEED_MEAN_TABLE.csv；固定比较及区间：PRIMARY_FIXED_COMPARISONS.csv；原始独立逐父损失：locked/confirmation_seed_mean_parent_losses.npz。']
 (OUT/'REPORT_ZH.md').write_text('\n'.join(lines)+'\n');print(json.dumps(dict(status=result['status'],failures=result['fallback_failure_count'],wall=result['wall_seconds'],path=str(OUT))))
if __name__=='__main__':main()
