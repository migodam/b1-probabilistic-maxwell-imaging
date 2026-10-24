#!/usr/bin/env python3
"""Postboundary saved-array eligibility audit; no geometry rerun/rescaling."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,csv
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write

PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(8));RADII=(.02,.04,.08)
FIRST_ARMS=('direct_bracket_forward','direct_bracket_reverse','frozen_weak','dynamic_weak','circuit_forward','circuit_reverse')
HEUN_ARMS=('direct_forward','direct_reverse','Heun_forward','Heun_reverse')
def load(p):return json.loads(Path(p).read_text())
def rows(p):return [json.loads(s) for s in p.read_text().splitlines()] if p.exists() else []
def csv_write(p,data):
 if not data:return
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
def family(arm):return 'Euler' if arm.startswith('circuit') else 'Heun' if arm.startswith('Heun') else 'direct' if arm.startswith('direct') else 'frozen' if arm=='frozen_weak' else 'dynamic'

def main():
 if not (OUT/'results/II/interval_chain_boundary_v1/COORDINATOR_COMPLETE.json').exists():raise ValueError('parent requested this audit only afterboundaryclosed')
 dest=OUT/'results/offline_only/geometry_matched_gate_audit_v1'
 if dest.exists():raise FileExistsError('preserve first gateaudit')
 for name in ('dynamic_geometry_first','dynamic_geometry_heun_r1'):
  if load(OUT/'III/jobs'/name/'CLOSED.json')['status']!='COMPLETE':raise ValueError('closedgeometryonly')
 b=Book('III','geometry_matched_gate_audit_v1',max_wall=600,analysis_only=True);status='COMPLETE';error=None
 try:
  dest.mkdir(parents=True);first_protocol=load(OUT/'configs/GEOMETRY_PROTOCOL.json');heun_protocol=load(OUT/'configs/GEOMETRY_HEUN_PROTOCOL.json');teacher=heun_protocol['teacher_eligible']
  white_radius=first_protocol['literal_noiseball_white_radius'];noise_phi=white_radius**2/2;contract_phi=load(OUT/'configs/REPRESENTATIVE_PROTOCOL.json')['Phi_threshold'];allrows=[]
  write(dest/'FROZEN_ANALYSIS_PROTOCOL.json',dict(distance_gate='abs(actualmaterialRMS/requested-1)<=.02',distance_gate_source='explicit parent analysis instruction afterclosedgeometry;notalgorithmretuning',
   original_teacher_distance_tolerance=teacher['radius_relative_error_max'],saved_Heun_teacher_thresholds=teacher,saved_geometry_noiseball_white_radius=white_radius,derived_geometry_noiseball_Phi=noise_phi,original_contract_Phi=contract_phi,
   noiseball_scope='historical diagnostic,notposteriorregion',threshold_transfer='Heun saved numerical/progress criteria applied DESCRIPTIVELY tofirstpilot classes;firstpilot didnot registergenericteachereligible rule',
   planned_independent_parents=8,planned_attempts='all8parents×dims64/66×ranks2/4/8×radii.02/.04/.08×declaredarms;missing/partialcountasnotcompleted',
   no_rescale=True,no_geometry_source_change=True,GT_gates=False,GT_read_after_eligibility_saved=True,new_RHS=0))
  for campaign,arms in (('dynamic_geometry_first',FIRST_ARMS),('dynamic_geometry_heun_r1',HEUN_ARMS)):
   isheun=campaign.endswith('heun_r1')
   for parent in PARENTS:
    with np.load(OUT/'results/contexts/cache'/f'{parent}.npz',allow_pickle=False) as z:anchor=z['chi'].copy()
    for dimension in (64,66):
     for rank in (2,4,8):
      config=OUT/'results/III'/campaign/parent/f'd{dimension}_r{rank}';records=rows(config/('TASK_ROWS.jsonl' if isheun else 'ARM_ROWS.jsonl'));endpoint_meta={r['endpoint']:r for r in rows(config/'ENDPOINT_ROWS.jsonl')}
      first=OUT/'results/III/dynamic_geometry_first'/parent/f'd{dimension}_r{rank}'
      with np.load(first/'BASE_FRAME.npz',allow_pickle=False) as z:D=z['Dchi'].copy();X=z['X'].copy()
      with np.load(first/('BRACKET_h0.001.npz' if isheun else 'BRACKET_h0.002.npz'),allow_pickle=False) as z:B=z['bracket'].copy()
      for radius in RADII:
       for arm in arms:
        selected=[r for r in records if r['arm']==arm and r['radius']==radius]
        if len(selected)>1:raise ValueError('duplicateplannedattempt')
        row=selected[0] if selected else None;f=family(arm);r=dict(campaign=campaign,parent=parent,dimension=dimension,rank=rank,radius=radius,arm=arm,family=f,completed=row is not None,
         endpoint=None,actual_RMS=None,relative_radius_error=None,matched_2pct=False,original_saved_radius_matched=None,Phi=None,DeltaPhi=None,prior_delta=None,signed_material_progress=None,
         endpoint_noiseball=False,endpoint_contract_noiseball=False,sampled_path_noiseball=False,all_paid_stages_noiseball=False,likelihood_delta_gate=False,progress_gate=False,prior_delta_gate=False,
         backward_gate=False,prior_legal=False,all_support_and_match_2pct=False,original_saved_teacher_eligible=None,actual_cost_RHS=None)
        if row is not None:
         endpoint=row['endpoint'];meta=endpoint_meta[endpoint]
         with np.load(OUT/'III/jobs'/campaign/'endpoints'/f'{endpoint}.npz',allow_pickle=False) as z:chi=z['chi'].copy()
         actual=float(np.linalg.norm(chi-anchor)/24)
         saved=row['actual_RMS']
         if not np.isclose(actual,saved,rtol=1e-10,atol=1e-12):raise ValueError('savedvsarrayRMSmismatch')
         sign=-1 if arm.endswith('reverse') else 1;direction=X if f in ('frozen','dynamic') else sign*B
         material_direction=D@direction;ref=np.r_[material_direction.real,material_direction.imag]/24
         delta=chi-anchor;v=np.r_[delta.real,delta.imag]/24;progress=float(v@ref/max(np.linalg.norm(ref),1e-30))
         if isheun and not np.isclose(progress,row['signed_material_progress'],rtol=1e-10,atol=1e-12):raise ValueError('savedsignedprogressmismatch')
         if f=='Heun':paid=row['stage_endpoint_keys'];sampled=paid[1::2]
         elif f=='Euler':paid=row['intermediate_endpoints'];sampled=paid
         elif f=='dynamic':paid=[row['midpoint']['endpoint'],endpoint];sampled=paid
         else:paid=[endpoint];sampled=paid
         completed_noise=all(endpoint_meta[e]['Phi']<=noise_phi for e in sampled);paid_noise=all(endpoint_meta[e]['Phi']<=noise_phi for e in paid)
         r.update(endpoint=endpoint,actual_RMS=actual,relative_radius_error=abs(actual/radius-1),matched_2pct=abs(actual/radius-1)<=.02,
          original_saved_radius_matched=row['radius_matched'],Phi=row['Phi'],DeltaPhi=row['DeltaPhi'],prior_delta=row['prior_delta'],signed_material_progress=progress,
          endpoint_noiseball=row['Phi']<=noise_phi,endpoint_contract_noiseball=row['Phi']<=contract_phi,sampled_path_noiseball=completed_noise,all_paid_stages_noiseball=paid_noise,
          likelihood_delta_gate=row['DeltaPhi']<=teacher['Full_DeltaPhi_max'],progress_gate=progress>=teacher['signed_material_progress_min_fraction']*radius,
          prior_delta_gate=row['prior_delta']<=teacher['prior_energy_delta_max'],backward_gate=meta['backward_residual']<=teacher['backward_residual_max'],prior_legal=meta['prior_legal'],
          original_saved_teacher_eligible=row.get('teacher_eligible'),actual_cost_RHS=(row['actual_cost'] if isheun else row['marginal_cost'])['RHS'])
         r['all_support_and_match_2pct']=all(r[k] for k in ('matched_2pct','endpoint_noiseball','likelihood_delta_gate','progress_gate','prior_delta_gate','backward_gate','prior_legal'))
        allrows.append(r)
  csv_write(dest/'ALL_PLANNED_ATTEMPT_GATES_NON_GT.csv',allrows);write(dest/'NON_GT_ELIGIBILITY_SAVED.json',dict(rows=allrows,GT_read=False,allplanned=len(allrows)))
  groups=[]
  for key in sorted(set((r['campaign'],r['dimension'],r['rank'],r['radius'],r['family']) for r in allrows)):
   selected=[r for r in allrows if (r['campaign'],r['dimension'],r['rank'],r['radius'],r['family'])==key];n=len(selected)
   g=dict(campaign=key[0],dimension=key[1],rank=key[2],radius=key[3],family=key[4],planned_attempts=n,independent_parents=8,completed=sum(r['completed'] for r in selected))
   for field in ('matched_2pct','endpoint_noiseball','sampled_path_noiseball','likelihood_delta_gate','progress_gate','prior_delta_gate','all_support_and_match_2pct'):
    g[field+'_count']=sum(r[field] for r in selected);g[field+'_fraction_allplanned']=g[field+'_count']/n
   g['all_support_fraction_among_matched']=g['all_support_and_match_2pct_count']/g['matched_2pct_count'] if g['matched_2pct_count'] else None;groups.append(g)
  csv_write(dest/'GATE_FRACTIONS_BY_RANK_RADIUS.csv',groups)
  index={(r['campaign'],r['parent'],r['dimension'],r['rank'],r['radius'],r['arm']):r for r in allrows};pairs=[]
  for r in allrows:
   if r['family']=='direct':continue
   refarm=('direct_' if r['campaign'].endswith('heun_r1') else 'direct_bracket_')+('reverse' if r['arm'].endswith('reverse') else 'forward')
   ref=index[(r['campaign'],r['parent'],r['dimension'],r['rank'],r['radius'],refarm)]
   pairs.append(dict(campaign=r['campaign'],parent=r['parent'],dimension=r['dimension'],rank=r['rank'],radius=r['radius'],arm=r['arm'],family=r['family'],reference_arm=refarm,
    candidate_endpoint=r['endpoint'],reference_endpoint=ref['endpoint'],both_completed=r['completed'] and ref['completed'],both_RMS_matched_2pct=r['matched_2pct'] and ref['matched_2pct'],
    both_full_support_gates=r['all_support_and_match_2pct'] and ref['all_support_and_match_2pct']))
  csv_write(dest/'PAIRED_ELIGIBILITY_NON_GT.csv',pairs)
  # First GT-derived material score access, after all eligibility artifacts.
  source=OUT/'results/offline_only/geometry_material_paths_v1/ALL_PAID_ENDPOINT_MATERIAL_SCORES.csv'
  with source.open() as f:scores={(r['campaign'],r['parent'],int(r['dimension']),int(r['rank']),r['endpoint']):float(r['GT_RMS']) for r in csv.DictReader(f)}
  eligible=[]
  for r in pairs:
   for gate in ('both_RMS_matched_2pct','both_full_support_gates'):
    if not r[gate]:continue
    key=(r['campaign'],r['parent'],r['dimension'],r['rank']);a=scores[key+(r['candidate_endpoint'],)];v=scores[key+(r['reference_endpoint'],)]
    eligible.append(dict(**r,eligibility_gate=gate,candidate_GT_RMS=a,reference_GT_RMS=v,GT_RMS_difference=a-v,scope='exposedposthocconditionalcomparison,notGTfiltered;allplanneddenominatorsretained'))
  csv_write(dest/'ELIGIBLE_ONLY_PAIRED_GT_SCORES.csv',eligible)
  totals=[]
  for key in sorted(set((r['campaign'],r['family']) for r in allrows)):
   selected=[r for r in allrows if (r['campaign'],r['family'])==key];n=len(selected);m=sum(r['matched_2pct'] for r in selected)
   totals.append(dict(campaign=key[0],family=key[1],planned=n,completed=sum(r['completed'] for r in selected),matched_2pct=m,matched_fraction_allplanned=m/n,
    endpoint_noiseball_and_matched=sum(r['endpoint_noiseball'] and r['matched_2pct'] for r in selected),likelihood_and_matched=sum(r['likelihood_delta_gate'] and r['matched_2pct'] for r in selected),
    progress_and_matched=sum(r['progress_gate'] and r['matched_2pct'] for r in selected),all_support_and_matched=sum(r['all_support_and_match_2pct'] for r in selected)))
  write(dest/'RESULTS.json',dict(status='NON_GT_GATE_AUDIT_AND_ELIGIBLE_OFFLINE_SCORE_COMPLETE',totals=totals,groups=groups,allplanned=len(allrows),completed=sum(r['completed'] for r in allrows),
   eligible_GT_comparison_rows=len(eligible),fixed_distance_gate=.02,original_teacher_gate=teacher['radius_relative_error_max'],new_RHS=0,new_states=0,new_LU=0,new_training=0,
   geometry_source_or_tolerances_changed=False,unmatched_GT_average_ranking_allowed=False,selection_uses_GT=False))
  report=['# 几何有限距离门禁：保留全部尝试的离线复核','',
   '按父线程明确指定的 abs(actual RMS/target−1)≤.02 做距离适用门禁；这不是把原算法/教师的1%条件改成2%，原始 radius_matched 与 teacher_eligible 全部原样保存。统计仅原前8DEV；四Euler已知自身场缺陷依旧不作为交换子理论反证。','',
   '|研究/类|计划尝试|完成|2%距离合格|合格率|合格且noiseball|合格且DeltaPhi≤2|合格且进展≥.5r|所有支持门禁合格|','|---|---:|---:|---:|---:|---:|---:|---:|---:|']
  for r in totals:report.append(f"|{r['campaign']}/{r['family']}|{r['planned']}|{r['completed']}|{r['matched_2pct']}|{r['matched_fraction_allplanned']:.3%}|{r['endpoint_noiseball_and_matched']}|{r['likelihood_and_matched']}|{r['progress_and_matched']}|{r['all_support_and_matched']}|")
  report += ['', '未完成/部分/跳过都保留在原计划分母。分维数、秩、半径的完整分数见 GATE_FRACTIONS_BY_RANK_RADIUS.csv；重复配置不是新增独立父对象。noiseball直接用原几何白化半径，合同Tau另列；DeltaPhi≤2、进展≥.5r、prior_delta≤2、残差≤1e-10来自原Heun教师阈值。对首轮没有通用教师资格规则的类仅作同阈值描述，不追认原先已通过教师门禁。采样路径和所有付费Heun predictor的noiseball标志分开，均无连续证书。','',
   '先保存不读真值的全部逐尝试门禁及固定直接控制配对资格，再读取原隔离离线GT评分。ELIGIBLE_ONLY_PAIRED_GT_SCORES.csv只含双方均满足2%距离、或更严格全部支持条件的同父对象配对。选择不使用GT，缺失/不合格分母仍在计划表；这是已暴露开发数据的事后条件描述，不能推广成盲验证。','',
   '原离线表中的Heun均值.0229/.0231与直接法.0200不构成同材料距离比较；严禁用那些不匹配平均GT误差排名算法。未满足距离门禁本身就是该算法在所请求半径的失败证据，不能通过离线缩放/挑近点救回。原始GT评分和失败结果都保留，没有新物理或rescale。']
  (OUT/'AUDIT/GEOMETRY_MATCHED_DISTANCE_GATE_REVIEW.md').write_text('\n'.join(report)+'\n')
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)

if __name__=='__main__':main()
