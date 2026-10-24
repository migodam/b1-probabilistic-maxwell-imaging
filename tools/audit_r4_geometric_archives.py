#!/usr/bin/env python3
"""Non-GT geometric archive build, then separately invoked isolated GT append.

No Full, model/backend construction, inverse, rescaling, training or posterior
weights. --score requires immutable choices and different process from --build.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json,csv,argparse,collections
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write,save
from p3r4.evaluation import checkpoint_identity
from p1_material.prior import MaterialPrior

PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(8));RADII=(.02,.04,.08)
CAMPAIGNS=('dynamic_geometry_first','dynamic_geometry_heun_r1')
FIRST_ARMS=('direct_bracket_forward','direct_bracket_reverse','frozen_weak','dynamic_weak','ordinary_field_curvature','ordinary_field_secant','circuit_forward','circuit_reverse','strong_compensated')
HEUN_ARMS=('direct_forward','direct_reverse','Heun_forward','Heun_reverse')
DEST=OUT/'results/III/geometric_candidate_archives_v2';OFFLINE=OUT/'results/offline_only/geometric_archive_GT_distance_v2'
def load(p):return json.loads(Path(p).read_text())
def lines(p):return [json.loads(s) for s in p.read_text().splitlines()] if p.exists() else []
def csv_write(p,data):
 if not data:return
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
def family(arm):return 'Euler' if arm.startswith('circuit') else 'Heun' if arm.startswith('Heun') else 'direct' if arm.startswith('direct') else arm

def build():
 if DEST.exists():raise FileExistsError('preserve first choices; no rebuild')
 quiet=OUT/'results/IV/interval_chain_locked_v1/COMPARISON.json';quiet_end=OUT/'results/IV/interval_chain_locked_pair/COMPLETE.json'
 if not quiet_end.exists() or load(quiet)['status']!='PASS_SAME_KERNEL_FOUR_ARM_NUMERICS':raise ValueError('quietlockedcomparisonmustclosebeforeanynumericalaudit')
 source_costs={}
 for campaign in CAMPAIGNS:
  receipt=OUT/'III/jobs'/campaign/'CLOSED.json';result=load(OUT/'results/III'/campaign/'RESULTS.json')
  if load(receipt)['status']!='COMPLETE' or result['status']!='COMPLETE' or result['no_GT'] is not True:raise ValueError('closed no-GT sources required')
  source_costs[campaign]=dict(receipt=checkpoint_identity(receipt),RHS=result['new_costs']['RHS'],configurations=result['configurations'])
 b=Book('III','geometric_archive_build_v2',max_wall=900,analysis_only=True);status='COMPLETE';error=None;seals=[]
 try:
  DEST.mkdir(parents=True);cfg=load(OUT/'configs/TARGET.json');sigma=cfg['noise']['sigma'];tau=load(OUT/'configs/REPRESENTATIVE_PROTOCOL.json')['Phi_threshold'];teacher=load(OUT/'configs/GEOMETRY_HEUN_PROTOCOL.json')['teacher_eligible']
  prior=MaterialPrior(cfg,load(ROOT/'research/p1_material_target_20261018/inputs/geometry.json')['points'])
  context_receipts=[]
  for p in sorted((OUT/'I/jobs').glob('contexts_dev*/CLOSED.json')):
   if load(p)['status']=='COMPLETE':
    source=p.parent/'SOURCE_REVISION.json'
    context_receipts.append(dict(receipt=checkpoint_identity(p),source=checkpoint_identity(source) if source.exists() else None,source_snapshot_status='PRESENT' if source.exists() else 'MISSING_ORIGINAL_AS_RUN_METADATA_NO_RECONSTRUCTION'))
  protocol=dict(status='FROZEN_NON_GT_ARCHIVE_RULE_BEFORE_CHOICES',role='GEOMETRIC_CANDIDATE',parents=PARENTS,
   candidate_scope='completed requestedradius task endpoints in originalfirst8DEV geometry+Heun,allcontrols;FD/smallh/intermediateprobes countedseparately,notnewtaskworlds',
   ordering='anchor FIRST; eligiblecandidates ascending originalonlineFullPhi thenUID',matched_radius_relative_error_max=.02,
   Phi_max=tau,Full_DeltaPhi_max=teacher['Full_DeltaPhi_max'],prior_delta_max=teacher['prior_energy_delta_max'],signed_progress_min_fraction=teacher['signed_material_progress_min_fraction'],
   backward_residual_max=teacher['backward_residual_max'],prior_legal=True,min_actual_full576_material_RMS=.02,
   separation='min(norm(chi-c)/24 for ALLselected)>=.02 STRICT,NOepsilon,NO.0196rescue,NOsnap/rescale',
   original_teacher_radius_relative_error=teacher['radius_relative_error_max'],source_rules_unchanged=True,weights=None,probabilities=None,candidate_only=True,
   posterior_sample=False,posterior_mean=False,mode_claim=False,GT_read=False,GT_append='onlyseparateprocessafterARCHIVE_SEALED',
   quiet_comparison_gate=checkpoint_identity(quiet),quiet_coordinator_gate=checkpoint_identity(quiet_end),context_owner_receipts=context_receipts)
  write(DEST/'FROZEN_ARCHIVE_PROTOCOL.json',protocol);allcases=[];summary=[];selected_metadata=[];owners=[];noncandidate_total=0
  for parent in PARENTS:
   cachepath=OUT/'results/contexts/cache'/f'{parent}.npz';obs=OUT/'inputs/online'/f'{parent}.npz';contextpath=OUT/'results/contexts/dev'/f'{parent}.json';context=load(contextpath)
   with np.load(cachepath,allow_pickle=False) as z:anchor={k:z[k].copy() for k in ('q','chi','prediction')}
   with np.load(obs,allow_pickle=False) as z:y=z['observation'].copy()
   phi0=float(np.vdot(anchor['prediction']-y,anchor['prediction']-y).real/sigma**2)
   if phi0>tau:raise ValueError('fixedanchoroutsideoriginalnoiseball; noanchorreplacement')
   anchor_uid='anchor:'+parent;candidates=[];parentcases=[];material_groups={};parentpaidkeys=set();taskkeys=set();owner={}
   for campaign in CAMPAIGNS:
    heun=campaign.endswith('heun_r1');arms=HEUN_ARMS if heun else FIRST_ARMS
    configs=source_costs[campaign]['configurations'];owner[campaign+'_parent_allconfig_RHS']=sum(c.get('actual_new_cost',c.get('actual_new_costs'))['RHS'] for c in configs if c['parent']==parent)
    for dimension in (64,66):
     for rank in (2,4,8):
      first=OUT/'results/III/dynamic_geometry_first'/parent/f'd{dimension}_r{rank}';folder=OUT/'results/III'/campaign/parent/f'd{dimension}_r{rank}'
      records=lines(folder/('TASK_ROWS.jsonl' if heun else 'ARM_ROWS.jsonl'));ep={r['endpoint']:r for r in lines(folder/'ENDPOINT_ROWS.jsonl')};parentpaidkeys.update((campaign,k) for k in ep)
      with np.load(first/'BASE_FRAME.npz',allow_pickle=False) as z:D=z['Dchi'].copy();X=z['X'].copy()
      with np.load(first/('BRACKET_h0.001.npz' if heun else 'BRACKET_h0.002.npz'),allow_pickle=False) as z:directions={k:z[k].copy() for k in ('bracket','curvature','secant')}
      for radius in RADII:
       for arm in arms:
        matches=[r for r in records if r['arm']==arm and r['radius']==radius]
        if len(matches)>1:raise ValueError('duplicate planned row')
        row=matches[0] if matches else None
        entry=dict(parent=parent,campaign=campaign,dimension=dimension,rank=rank,arm=arm,family=family(arm),target_RMS=radius,completed=row is not None,
         UID=f'{campaign}:{parent}:d{dimension}_r{rank}:{arm}:r{radius:g}',endpoint=None,actual_RMS=None,Phi=None,DeltaPhi=None,prior_delta=None,signed_material_progress=None,
         eligible_support=False,original_saved_radius_matched=None,original_saved_teacher_eligible=None,exact_material_duplicate_group_uid=None,
         archive_selected=False,archive_order=None,min_RMS_to_selected=None,nearest_selected_UID=None,separation_deficit=None,near_02_boundary_rejection=False,
         filter_reasons=None,source_path=None,path_RHS=None,endpoint_Full_RHS=8 if row else 0)
        reasons=[]
        if row is None:reasons.append('NOT_COMPLETED_OR_NOT_RUN_OR_CAP_SKIP')
        else:
         endpoint=row['endpoint'];uid=campaign+':'+endpoint;source=OUT/'III/jobs'/campaign/'endpoints'/f'{endpoint}.npz';meta=ep[endpoint];taskkeys.add((campaign,endpoint))
         with np.load(source,allow_pickle=False) as z:a={k:z[k].copy() for k in ('q','chi','prediction')}
         decoded,_=prior.decode(a['q'])
         if not np.array_equal(decoded,a['chi']):raise ValueError('cached material not original decoder')
         phi=float(np.vdot(a['prediction']-y,a['prediction']-y).real/sigma**2);actual=float(np.linalg.norm(a['chi']-anchor['chi'])/24)
         if not np.isclose(phi,row['Phi'],rtol=1e-12,atol=1e-10) or not np.isclose(actual,row['actual_RMS'],rtol=1e-10,atol=1e-12):raise ValueError('saved FullPhi/material mismatch')
         sign=-1 if arm.endswith('reverse') else 1
         direction=X if arm in ('frozen_weak','dynamic_weak','strong_compensated') else directions['curvature'] if arm=='ordinary_field_curvature' else directions['secant'] if arm=='ordinary_field_secant' else sign*directions['bracket']
         v=D@direction;ref=np.r_[v.real,v.imag]/24;delta=a['chi']-anchor['chi'];progress=float(np.r_[delta.real,delta.imag]@ref/24/max(np.linalg.norm(ref),1e-30))
         checks=dict(MATCHED_RMS_2PCT=abs(actual/radius-1)<=.02,ORIGINAL_NOISEBALL=phi<=tau,FULL_DELTA_PHI=phi-phi0<=teacher['Full_DeltaPhi_max'],
          PRIOR_DELTA=row['prior_delta']<=teacher['prior_energy_delta_max'],SIGNED_PROGRESS=progress>=teacher['signed_material_progress_min_fraction']*radius,
          PRIOR_LEGAL=meta['prior_legal'],BACKWARD_RESIDUAL=meta['backward_residual']<=teacher['backward_residual_max'])
         reasons.extend(k+'_FAIL' for k,value in checks.items() if not value)
         key=a['chi'].tobytes();material_groups.setdefault(key,[]).append(uid)
         entry.update(UID=uid,endpoint=endpoint,actual_RMS=actual,Phi=phi,DeltaPhi=phi-phi0,prior_delta=row['prior_delta'],signed_material_progress=progress,
          eligible_support=all(checks.values()),original_saved_radius_matched=row['radius_matched'],original_saved_teacher_eligible=row.get('teacher_eligible'),source_path=str(source),
          path_RHS=(row['actual_cost'] if heun else row['marginal_cost'])['RHS'])
         entry['_material_key']=key
         if entry['eligible_support']:candidates.append(dict(entry=entry,arrays=a))
        entry['filter_reasons']=';'.join(reasons);parentcases.append(entry)
   duplicate_uid={key:min(uids) for key,uids in material_groups.items() if len(uids)>1}
   for entry in parentcases:
    if '_material_key' in entry:entry['exact_material_duplicate_group_uid']=duplicate_uid.get(entry['_material_key'])
   selected=[dict(UID=anchor_uid,arrays=anchor,metadata=dict(UID=anchor_uid,role='GEOMETRIC_CANDIDATE',source_path=str(cachepath),parent=parent,Phi=phi0,archive_order=0,anchor=True,weights=None,probabilities=None,candidate_only=True))]
   for candidate in sorted(candidates,key=lambda c:(c['entry']['Phi'],c['entry']['UID'])):
    r=candidate['entry'];a=candidate['arrays'];distances=[float(np.linalg.norm(a['chi']-s['arrays']['chi'])/24) for s in selected];index=int(np.argmin(distances));distance=distances[index]
    r.update(min_RMS_to_selected=distance,nearest_selected_UID=selected[index]['UID'],separation_deficit=max(0.,.02-distance))
    if distance>=.02:
     r.update(archive_selected=True,archive_order=len(selected));m={k:v for k,v in r.items() if k!='_material_key'};m.update(role='GEOMETRIC_CANDIDATE',weights=None,probabilities=None,candidate_only=True,anchor=False)
     selected.append(dict(UID=r['UID'],arrays=a,metadata=m))
    else:r['filter_reasons']='GREEDY_MIN_FULL576_RMS_LT_0.02';r['near_02_boundary_rejection']=abs(distance-.02)<=1e-10
   p=DEST/parent/'ARCHIVE.npz';save(p,UID=np.array([s['UID'] for s in selected]),q=np.array([s['arrays']['q'] for s in selected]),chi=np.array([s['arrays']['chi'] for s in selected]),
    prediction=np.array([s['arrays']['prediction'] for s in selected]),Phi=np.array([s['metadata']['Phi'] for s in selected]),source_path=np.array([s['metadata']['source_path'] for s in selected]))
   selection=[s['metadata'] for s in selected];write(p.with_suffix('.json'),dict(parent=parent,role='GEOMETRIC_CANDIDATE',weights=None,probabilities=None,candidate_only=True,entries=selection))
   # Verify strict separation again on the exact arrays actually published.
   with np.load(p,allow_pickle=False) as z:published=z['chi'].copy()
   pairdist=[float(np.linalg.norm(published[i]-published[j])/24) for i in range(len(published)) for j in range(i)]
   if pairdist and min(pairdist)<.02:raise AssertionError('publishedarchiveviolatesstrictseparation')
   seals.extend([checkpoint_identity(p),checkpoint_identity(p.with_suffix('.json'))]);selected_metadata.extend(selection)
   eligible=[r for r in parentcases if r['eligible_support']]
   summary.append(dict(parent=parent,planned_task_attempts=len(parentcases),completed_task_endpoints=sum(r['completed'] for r in parentcases),support_eligible_endpoints=len(eligible),
    exact_unique_completed_materials=len(material_groups),exact_duplicate_completed_endpoints=sum(r['completed'] for r in parentcases)-len(material_groups),
    exact_unique_eligible_materials=len(set(r['_material_key'] for r in eligible)),eligible_dynamic_endpoints=sum(r['eligible_support'] and r['arm']=='dynamic_weak' for r in parentcases),
    strict_archive_worlds_including_anchor=len(selected),strict_archive_new_geometric_candidates=len(selected)-1,min_published_pair_RMS=min(pairdist) if pairdist else None,
    strict_separation_rejected=len(eligible)-(len(selected)-1),near_02_boundary_rejected=sum(r['near_02_boundary_rejection'] for r in parentcases),anchor_Phi=phi0,
    archive_role='GEOMETRIC_CANDIDATE',weights=None,probabilities=None,candidate_only=True))
   noncandidate=len(parentpaidkeys-taskkeys);noncandidate_total+=noncandidate
   owners.append(dict(parent=parent,context_source=str(contextpath),context_cache_identity=checkpoint_identity(cachepath),observation_identity=checkpoint_identity(obs),
    actual_shared_context_RHS=context['new_RHS'],context_native_cap=context['native_cap'],context_probit_cap=context['probit_cap'],
    **owner,noncandidate_paid_probe_endpoints=noncandidate,
    scope='originalsearchacquisitionowneronce,notperselectedworld;pathRHScontainsselectedFull8plusoriginalpathFullJ;FD480andHeungates1728alreadyincampaigncost,notrebilled'))
   for r in parentcases:r.pop('_material_key',None)
   allcases.extend(parentcases)
  csv_write(DEST/'ALL_PLANNED_CASES_AND_FILTER_REASONS.csv',allcases);csv_write(DEST/'PER_PARENT_WORLD_COUNTS.csv',summary)
  write(DEST/'ALL_PLANNED_CASES.json',dict(rows=allcases,GT_read=False));write(DEST/'COST_OWNERS.json',dict(campaigns={k:{f:v for f,v in s.items() if f!='configurations'} for k,s in source_costs.items()},
   entire_original_teacher_RHS=459712,first8_parent_campaign_RHS=sum(o['dynamic_geometry_first_parent_allconfig_RHS']+o['dynamic_geometry_heun_r1_parent_allconfig_RHS'] for o in owners),
   shared_context_RHS_first8=sum(o['actual_shared_context_RHS'] for o in owners),parent_owners=owners,context_owner_receipts=context_receipts,new_archive_Maxwell_RHS=0))
  reasons=collections.Counter(reason for r in allcases for reason in r['filter_reasons'].split(';') if reason)
  result=dict(status='ARCHIVE_CHOICES_COMPLETE_NON_GT',role='GEOMETRIC_CANDIDATE',parent_rows=summary,planned_task_attempts=len(allcases),completed_task_endpoints=sum(r['completed'] for r in allcases),
   support_eligible_endpoints=sum(r['eligible_support'] for r in allcases),eligible_dynamic_endpoints=sum(r['eligible_support'] and r['arm']=='dynamic_weak' for r in allcases),
   greedy_worlds_including_8anchors=sum(r['strict_archive_worlds_including_anchor'] for r in summary),filter_reason_counts=dict(reasons),noncandidate_paid_probe_endpoints=noncandidate_total,
   weights=None,probabilities=None,candidate_only=True,posterior_samples=False,posterior_mean=False,mode_count_claim=False,maximum_packing_claim=False,GT_read=False,new_Maxwell_RHS=0,new_LU=0,new_states=0,new_training=0)
  write(DEST/'RESULTS.json',result)
  report=['# 有限材料候选的严格分离档案：非后验','',
   '仅在quiet locked四臂比较完成后开始零物理分析。先固定非GT规则，再构造档案：原锚点始终第一；完成的有限task候选通过2%距离匹配、原Phi≤282.873792、DeltaPhi≤2、prior_delta≤2、signed progress≥.5target、原合法性/后向残差门禁，按online Full Phi再UID排序。然后逐一只在相对所有已选世界的 full576 chi RMS≥.02 时接纳。严格比较没有epsilon、.0196救回、snap或重缩放。','',
   '候选池包含原first8DEV两次封账几何研究的全部有限task控制，普通曲率/割线/strong控制也保留；不把FD、小h、Heun predictor或中间路径探针当新的task世界。所有计划跳过/未完成、重复材料、门禁失败和贪心距离拒绝均保留，不把83个dynamic合格端点直接称83个世界。共享阈值对原未定义通用教师资格的控制只作描述筛选，不改原1%教师flags。','',
   '|父对象|计划task|完成端点|support合格|其中dynamic端点|严格档案世界（含锚点）|新增候选|近.02边界严格拒绝|','|---|---:|---:|---:|---:|---:|---:|---:|']
  for r in summary:report.append(f"|{r['parent']}|{r['planned_task_attempts']}|{r['completed_task_endpoints']}|{r['support_eligible_endpoints']}|{r['eligible_dynamic_endpoints']}|{r['strict_archive_worlds_including_anchor']}|{r['strict_archive_new_geometric_candidates']}|{r['near_02_boundary_rejected']}|")
  report += ['', '每个ARCHIVE.npz保留完整chi/q/prediction和原付费端点来源；COST_OWNERS披露原459712RHS整研究、所用8对象几何取得、公共96+64上下文、路径/FD/门禁共享拥有，全部不重计。新Maxwell/LU/RHS/训练为0。相同chi但不同q仍是同材料世界，档案数是固定Phi优先贪心分离的数量，不是最大packing、独立样本数、模式数或后验复杂度。','',
   '角色始终GEOMETRIC_CANDIDATE，weights/probabilities为null，candidate_only=true；没有MMSE、后验均值/样本或几何优势结论。此构造完全不读取GT；choices封存并关闭后才允许另一个隔离进程追加GT距离，不能重排、删增或调权。']
  (OUT/'AUDIT/GEOMETRIC_CANDIDATE_ARCHIVE_REVIEW.md').write_text('\n'.join(report)+'\n')
  seals.extend(checkpoint_identity(DEST/name) for name in ('FROZEN_ARCHIVE_PROTOCOL.json','ALL_PLANNED_CASES_AND_FILTER_REASONS.csv','ALL_PLANNED_CASES.json','PER_PARENT_WORLD_COUNTS.csv','COST_OWNERS.json','RESULTS.json'))
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:
  b.close(status,error)
  if status=='COMPLETE':write(DEST/'ARCHIVE_SEALED.json',dict(status='CHOICES_SEALED_BEFORE_GT',build_pid=os.getpid(),GT_read=False,artifacts=seals,receipt=checkpoint_identity(b.jobdir/'CLOSED.json'),role='GEOMETRIC_CANDIDATE',weights=None,candidate_only=True))

def score():
 if OFFLINE.exists():raise FileExistsError('preserve isolated append')
 seal=load(DEST/'ARCHIVE_SEALED.json')
 if seal['status']!='CHOICES_SEALED_BEFORE_GT' or seal['GT_read'] is not False or seal['build_pid']==os.getpid():raise ValueError('isolatedprocessafterchoicesrequired')
 if load(seal['receipt']['path'])['status']!='COMPLETE' or checkpoint_identity(seal['receipt']['path'])!=seal['receipt']:raise ValueError('closedbuildreceiptchanged')
 for identity in seal['artifacts']:
  if checkpoint_identity(identity['path'])!=identity:raise ValueError('sealedchoiceschanged')
 b=Book('III','geometric_archive_isolated_GT_append_v2',max_wall=600,analysis_only=True);status='COMPLETE';error=None
 try:
  OFFLINE.mkdir(parents=True);data=[]
  for parent in PARENTS:
   # First truth read is confined to this separate post-seal branch/process.
   with np.load(OUT/'inputs/offline_only'/f'{parent}.npz',allow_pickle=False) as z:truth=z['chi'].copy()
   with np.load(DEST/parent/'ARCHIVE.npz',allow_pickle=False) as z:chi=z['chi'].copy();uids=z['UID'].astype(str);phi=z['Phi'].copy()
   for i,(uid,world) in enumerate(zip(uids,chi)):data.append(dict(parent=parent,UID=uid,archive_order=i,anchor=i==0,Phi=float(phi[i]),material_RMS_to_GT=float(np.linalg.norm(world-truth)/24),
    role='GEOMETRIC_CANDIDATE',weights=None,candidate_only=True,GT_use='appendonlydescriptive;choicesunchanged',source_archive=str(DEST/parent/'ARCHIVE.npz')))
  csv_write(OFFLINE/'SEALED_WORLD_GT_DISTANCES.csv',data);write(OFFLINE/'RESULTS.json',dict(status='ISOLATED_GT_DISTANCE_APPEND_COMPLETE',rows=data,
   choices_seal=checkpoint_identity(DEST/'ARCHIVE_SEALED.json'),choices_modified=False,GT_used_for_archive=False,new_Maxwell_RHS=0,new_training=0))
  with (OUT/'AUDIT/GEOMETRIC_CANDIDATE_ARCHIVE_REVIEW.md').open('a') as f:f.write('\n隔离GT距离已在choices封存后由不同PID追加至 results/offline_only/geometric_archive_GT_distance_v2/SEALED_WORLD_GT_DISTANCES.csv；原档案身份复核未变，GT没有参与顺序、数量、角色、权重或候选选择。\n')
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('build','score'));a=p.parse_args()
 if a.action=='build':build()
 else:score()
if __name__=='__main__':main()
