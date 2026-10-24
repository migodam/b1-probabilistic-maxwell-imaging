#!/usr/bin/env python3
"""Isolated postclose GT scoring of already-paid geometry arrays; no Full/fit."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,csv
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write

def load(p):return json.loads(Path(p).read_text())
def lines(p):return [json.loads(z) for z in p.read_text().splitlines()] if p.exists() else []
def csv_write(p,rows):
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
 dest=OUT/'results/offline_only/geometry_material_paths_v1'
 if dest.exists():raise FileExistsError('fresh offline score only')
 campaigns=('dynamic_geometry_first','dynamic_geometry_heun_r1');gates={}
 # Validate every online close before the first GT access.
 for campaign in campaigns:
  result=load(OUT/'results/III'/campaign/'RESULTS.json');receipt=load(OUT/'III/jobs'/campaign/'CLOSED.json');start=load(OUT/'III/jobs'/campaign/'START.json')
  if result['status']!='COMPLETE' or receipt['status']!='COMPLETE' or result['no_GT'] is not True or start['pid']==os.getpid():raise ValueError('separate isolated process and sealed observation-only online work required')
  gates[campaign]=dict(result_status=result['status'],receipt=str(OUT/'III/jobs'/campaign/'CLOSED.json'),physical_pid=start['pid'],new_costs=result['new_costs'])
 b=Book('III','geometry_material_paths_offline_v1',max_wall=900,analysis_only=True);status='COMPLETE';error=None
 try:
  dest.mkdir(parents=True);write(dest/'OFFLINE_GATE.json',dict(status='CLOSED_ONLINE_GATES_VERIFIED_BEFORE_GT',sources=gates,scoring_pid=os.getpid(),GT_directory=str(OUT/'inputs/offline_only'),new_physics_RHS=0))
  endpoint_rows=[];path_rows=[];missing=[]
  for parent_i in range(8):
   parent=f'R4DEV_{parent_i:05d}';truth_path=OUT/'inputs/offline_only'/f'{parent}.npz'
   if not truth_path.exists():missing.append(dict(parent=parent,reason='MISSING_GT_NO_SUBSTITUTION'));continue
   with np.load(truth_path,allow_pickle=False) as z:truth=z['chi'].copy()
   if truth.shape!=(576,) or not np.iscomplexobj(truth) or not np.isfinite(truth).all():raise ValueError('invalid isolated material truth')
   with np.load(OUT/'results/contexts/cache'/f'{parent}.npz',allow_pickle=False) as z:anchor=z['chi'].copy()
   anchor_gt=float(np.linalg.norm(anchor-truth)/24);context=load(OUT/'results/contexts/dev'/f'{parent}.json')
   for campaign in campaigns:
    root=OUT/'results/III'/campaign/parent;job=OUT/'III/jobs'/campaign;heun=campaign.endswith('heun_r1')
    for config in sorted(root.glob('d*_r*')):
     parts=config.name.split('_');dimension=int(parts[0][1:]);rank=int(parts[1][1:]);metrics={}
     for row in lines(config/'ENDPOINT_ROWS.jsonl'):
      endpoint=row['endpoint'];source=job/'endpoints'/f'{endpoint}.npz'
      if not source.exists():missing.append(dict(parent=parent,campaign=campaign,endpoint=endpoint,reason='MISSING_PAID_ARRAY'));continue
      with np.load(source,allow_pickle=False) as z:chi=z['chi'].copy()
      gt_rms=float(np.linalg.norm(chi-truth)/24);actual=float(np.linalg.norm(chi-anchor)/24)
      item=dict(campaign=campaign,parent=parent,dimension=dimension,rank=rank,arm_or_phase=row.get('arm',row.get('phase')),label=row['label'],endpoint=endpoint,
       GT_RMS=gt_rms,anchor_GT_RMS=anchor_gt,GT_RMS_change=gt_rms-anchor_gt,actual_material_RMS_from_anchor=actual,Phi=row['Phi'],prior_delta=row['prior_delta'],source=str(source),offline_oracle_only=True)
      metrics[endpoint]=item;endpoint_rows.append(item)
     arms=lines(config/('TASK_ROWS.jsonl' if heun else 'ARM_ROWS.jsonl'))
     referenced=set()
     for row in arms:
      arm=row['arm'];end=row['endpoint']
      if end not in metrics:continue
      if heun and arm.startswith('Heun'):
       sampled=row['stage_endpoint_keys'][1::2];paid=row['stage_endpoint_keys'];pathkind='four corrected Heun endpoints; predictors separately audited'
      elif arm.startswith('circuit'):
       sampled=row['intermediate_endpoints'];paid=sampled;pathkind='four Euler segment endpoints; leading self-defect retained'
      elif 'midpoint' in row:
       sampled=[row['midpoint']['endpoint'],end];paid=sampled;pathkind='two sampled dynamic weak endpoints'
      else:sampled=[end];paid=sampled;pathkind='one finite endpoint, no continuous path evaluated'
      available=[e for e in sampled if e in metrics];allavailable=[e for e in paid if e in metrics];referenced.update(allavailable)
      candidates=[('anchor',anchor_gt)]+[(e,metrics[e]['GT_RMS']) for e in available];best=min(candidates,key=lambda x:x[1]);allbest=min([('anchor',anchor_gt)]+[(e,metrics[e]['GT_RMS']) for e in allavailable],key=lambda x:x[1])
      cost=row['actual_cost'] if heun else row['marginal_cost'];sharedFD=480;gatecost=1728 if heun else 0
      path_rows.append(dict(campaign=campaign,parent=parent,dimension=dimension,rank=rank,arm=arm,radius=row['radius'],actual_final_material_RMS=metrics[end]['actual_material_RMS_from_anchor'],
       radius_matched=row['radius_matched'],anchor_GT_RMS=anchor_gt,final_GT_RMS=metrics[end]['GT_RMS'],final_GT_RMS_change=metrics[end]['GT_RMS_change'],
       closest_sampled_path_GT_RMS=best[1],closest_sampled_path_endpoint=best[0],closest_any_paid_stage_GT_RMS=allbest[1],closest_any_paid_stage_endpoint=allbest[0],
       sampled_points=len(available),all_paid_stage_points=len(allavailable),path_kind=pathkind,path_RHS=cost['RHS'],path_states=cost['states'],path_LU=cost['LU'],
       preceding_shared_threeh_FD_RHS=sharedFD,preceding_Heun_small_h_gate_RHS=gatecost,shared_online_context_RHS=context['new_RHS'],
       shared_setup_inclusive_research_path_RHS=context['new_RHS']+sharedFD+gatecost+cost['RHS'],
       setup_scope='actual original shared research setup,NOT marginal standalone-method necessity; repeated acrosspaths,do not sum',
       source_row_file=str(config/('TASK_ROWS.jsonl' if heun else 'ARM_ROWS.jsonl')),source_endpoint=metrics[end]['source'],oracle_only=True,continuous_path_certificate=False))
     # Partial/nonterminal probes stay in the endpoint table, with explicit
     # unmatched count; never discard them or turn them into completed paths.
   print(parent,flush=True)
  csv_write(dest/'ALL_PAID_ENDPOINT_MATERIAL_SCORES.csv',endpoint_rows);csv_write(dest/'COMPLETED_FINITE_PATH_SCORES.csv',path_rows)
  grouped=[]
  keys=sorted(set((r['campaign'],r['dimension'],r['rank'],r['arm'],r['radius']) for r in path_rows))
  for key in keys:
   rows=[r for r in path_rows if (r['campaign'],r['dimension'],r['rank'],r['arm'],r['radius'])==key]
   grouped.append(dict(campaign=key[0],dimension=key[1],rank=key[2],arm=key[3],radius=key[4],parents=len(rows),
    mean_actual_final_material_RMS=float(np.mean([r['actual_final_material_RMS'] for r in rows])),matched_radius_count=sum(r['radius_matched'] for r in rows),
    mean_anchor_GT_RMS=float(np.mean([r['anchor_GT_RMS'] for r in rows])),mean_final_GT_RMS=float(np.mean([r['final_GT_RMS'] for r in rows])),
    mean_final_GT_RMS_change=float(np.mean([r['final_GT_RMS_change'] for r in rows])),mean_closest_sampled_path_GT_RMS=float(np.mean([r['closest_sampled_path_GT_RMS'] for r in rows])),
    mean_path_RHS=float(np.mean([r['path_RHS'] for r in rows])),mean_shared_setup_inclusive_research_path_RHS=float(np.mean([r['shared_setup_inclusive_research_path_RHS'] for r in rows]))))
  csv_write(dest/'GROUPED_DESCRIPTIVE.csv',grouped)
  result=dict(status='ISOLATED_OFFLINE_SCORE_COMPLETE',parents_requested=8,paid_endpoint_scores=len(endpoint_rows),completed_finite_paths=len(path_rows),missing=missing,
   groups=grouped,source_gate=gates,source_total_teacher_RHS=459712,all_scored_materials_already_paid=True,new_Full_RHS=0,new_LU=0,new_states=0,neural_fits=0,
   GT_read=True,GT_use='offlineonly descriptive endpoint/path error; no online selection/proposal/teacher changes',continuous_path_certification=False,
   closest_distance='posthocoracle among originalanchor and savedpathpoints only; Heunpredictors reported separately',Euler_scope='knownfourEuler selfdefect;notexactcommutator/theoremfalsification',parent_final_review=True)
  write(dest/'RESULTS.json',result)
  report=['# 已封账几何材料路径：隔离离线评分','',f'先验证两次在线物理作业 COMPLETE/no_GT 和不同 PID，再在本独立 BookIII 分析进程读取原前8DEV真值。已评分 {len(endpoint_rows)} 个已付费端点、{len(path_rows)} 个完成的有限路径；缺失 {len(missing)} 项。没有 Full、LU、RHS 或训练。','',
   '每条路径列出原锚点、终点和已采样路径点到真值的材料 RMS=||chi−chiGT||/24，以及实际有限材料移动与匹配标记。最接近真值的点仅为离线事后 oracle 描述，不反向挑选候选、修正在线图像或改变教师资格。Heun 四个校正端点与四个 predictor 付费场探针分开，predictor 的较近距离不冒充物理路径终点；Euler 仍保留已知 O(h²) 自身场缺陷，不能当成精确交换子或理论反证。','',
   'COMPLETED_FINITE_PATH_SCORES.csv 列原路径所有中间 Full-J/最终 Full 成本。原三档FD共享480RHS，Heun小步门禁原1728RHS，原在线上下文实际费用也分列；setup-inclusive列是原共享研究取得成本加路径，不是宣称每方法都必需同样冷成本，各路径共享项不能相加。ALL_PAID_ENDPOINT_MATERIAL_SCORES.csv 保留包括 FD 和未完成末段在内的全部端点，未把部分路径当完成。','',
   'GROUPED_DESCRIPTIVE.csv 按原算法/维数/秩/半径汇总，不以重复配置增加独立父对象数；不同实际移动不能假称同材料距离比较。输入都是已暴露开发对象，没有后验参考或连续路径证书。完整机器证据位于 results/offline_only/geometry_material_paths_v1/RESULTS.json；严禁把这里的真值误差或 oracle 索引作为在线特征。']
  (OUT/'AUDIT/GEOMETRY_MATERIAL_PATHS_OFFLINE.md').write_text('\n'.join(report)+'\n')
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)

if __name__=='__main__':main()
