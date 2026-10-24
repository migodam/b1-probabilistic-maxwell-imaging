#!/usr/bin/env python3
"""Descriptive, physics-free Heun revision aggregation after its job closes."""
from pathlib import Path
from collections import Counter,defaultdict
import json,numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'
RESULT=OUT/'results/III/dynamic_geometry_heun_r1';FIRST=OUT/'results/III/dynamic_geometry_first'


def stats(values):
 a=np.asarray(values,float)
 return dict(n=len(a),minimum=float(a.min()),median=float(np.median(a)),q25=float(np.quantile(a,.25)),q75=float(np.quantile(a,.75)),maximum=float(a.max())) if len(a) else dict(n=0)


def main():
 receipt=json.loads((OUT/'III/jobs/dynamic_geometry_heun_r1/CLOSED.json').read_text());result=json.loads((RESULT/'RESULTS.json').read_text())
 if receipt['status']!='COMPLETE' or result['status']!='COMPLETE':raise ValueError('completed Heun job required')
 small=[];task=[];skips=[]
 for path in RESULT.glob('R4DEV_*/d*_r*/SMALL_H_ROWS.jsonl'):
  parent=path.parent.parent.name;parts=path.parent.name.split('_');d=int(parts[0][1:]);rank=int(parts[1][1:])
  for line in path.read_text().splitlines():small.append(dict(parent=parent,dimension=d,rank=rank,**json.loads(line)))
 for path in RESULT.glob('R4DEV_*/d*_r*/TASK_ROWS.jsonl'):
  parent=path.parent.parent.name;parts=path.parent.name.split('_');d=int(parts[0][1:]);rank=int(parts[1][1:])
  for line in path.read_text().splitlines():task.append(dict(parent=parent,dimension=d,rank=rank,**json.loads(line)))
 for path in RESULT.glob('R4DEV_*/d*_r*/SKIPS.jsonl'):
  for line in path.read_text().splitlines():skips.append(dict(configuration=str(path.parent.relative_to(RESULT)),**json.loads(line)))
 small_groups=[]
 for d in (64,66):
  for rank in (2,4,8):
   for h in (.001,.002,.004):
    rows=[z for z in small if z['dimension']==d and z['rank']==rank and z['step']==h]
    small_groups.append(dict(dimension=d,rank=rank,h=h,count=len(rows),relative_material_error=stats([z['relative_material_bracket_error'] for z in rows]),relative_latent_error=stats([z['relative_latent_bracket_error'] for z in rows]),seam_flags=sum(z['any_seam_crossing'] for z in rows)))
 groups=defaultdict(list)
 for z in task:groups[(z['dimension'],z['rank'],z['arm'],z['radius'])].append(z)
 task_groups=[]
 for (d,rank,arm,radius),rows in sorted(groups.items()):
  task_groups.append(dict(dimension=d,rank=rank,arm=arm,radius=radius,n=len(rows),
   actual_RMS_ratio=stats([z['actual_RMS']/radius for z in rows]),matched_within_1pct=sum(abs(z['actual_RMS']/radius-1)<=.01 for z in rows),
   DeltaPhi=stats([z['DeltaPhi'] for z in rows]),prior_delta=stats([z['prior_delta'] for z in rows]),
   actual_normal_RMS=stats([z['actual_endpoint_normal']['actual_normal_RMS'] for z in rows]),eligible_teachers=sum(z['teacher_eligible'] for z in rows),
   marginal_RHS=stats([z['actual_cost']['RHS'] for z in rows]),wall_seconds=stats([z['wall_seconds'] for z in rows])))
 index={(z['parent'],z['dimension'],z['rank'],z['arm'],z['radius']):z for z in task}
 paired=[]
 for z in task:
  if not z['arm'].startswith('Heun'):continue
  orientation=z['arm'].split('_')[-1];direct=index.get((z['parent'],z['dimension'],z['rank'],'direct_'+orientation,z['radius']))
  if direct:
   pair=dict(parent=z['parent'],dimension=z['dimension'],rank=z['rank'],radius=z['radius'],orientation=orientation,
    Heun_actual_RMS_ratio=z['actual_RMS']/z['radius'],direct_actual_RMS_ratio=direct['actual_RMS']/z['radius'],
    Heun_DeltaPhi=z['DeltaPhi'],direct_DeltaPhi=direct['DeltaPhi'],Heun_prior_delta=z['prior_delta'],direct_prior_delta=direct['prior_delta'],
    marginal_RHS_ratio=z['actual_cost']['RHS']/direct['actual_cost']['RHS'],
    shared_FD_and_small_h_inclusive_RHS_ratio=(480+1728+z['actual_cost']['RHS'])/(480+1728+direct['actual_cost']['RHS']),
    no_same_material_precision_claim=True,continuous_path_certificate=False)
   oldpath=FIRST/z['parent']/f'd{z["dimension"]}_r{z["rank"]}'/'ARM_ROWS.jsonl'
   old=[json.loads(line) for line in oldpath.read_text().splitlines()]
   old=[r for r in old if r['arm']=='circuit_'+orientation and r['radius']==z['radius']]
   if old:pair.update(first_Euler_RMS_ratio=old[0]['actual_RMS']/z['radius'],first_Euler_DeltaPhi=old[0]['DeltaPhi'],Euler_reference_h=.002,Heun_reference_h=.001,causal_same_step_integrator_comparison=False)
   paired.append(pair)
 heun=[z for z in task if z['arm'].startswith('Heun')];direct=[z for z in task if z['arm'].startswith('direct')]
 summary=dict(status='DESCRIPTIVE_HEUN_REPAIR_COMPLETED_PARENT_REVIEW_REQUIRED',independent_fixed_parents=8,configurations=len(result['configurations']),
  passed_gates=result['passed_gates'],failed_gates=result['failed_gates'],gate_failure_reasons=dict(Counter(reason for c in result['configurations'] for reason in c['gate']['reasons'])),
  small_h_completed=len(small),fine_h_material_relative_error=stats([z['relative_material_bracket_error'] for z in small if z['step']==.001]),
  fine_h_latent_relative_error=stats([z['relative_latent_bracket_error'] for z in small if z['step']==.001]),
  small_h_seam_flags=sum(z['any_seam_crossing'] for z in small),rank_changes=sum(z['any_rank_change'] for z in small),near_crossings=sum(z['any_near_crossing'] for z in small),
  Heun_task_endpoints=len(heun),direct_task_endpoints=len(direct),
  Heun_task_actual_RMS_ratio=stats([z['actual_RMS']/z['radius'] for z in heun]),Heun_within_1pct=sum(abs(z['actual_RMS']/z['radius']-1)<=.01 for z in heun),
  direct_within_1pct=sum(abs(z['actual_RMS']/z['radius']-1)<=.01 for z in direct),
  eligible_Heun_teachers=sum(z['teacher_eligible'] for z in heun),eligible_direct_teachers=sum(z['teacher_eligible'] for z in direct),
  skip_count=len(skips),skip_reasons=dict(Counter(z['reason'] for z in skips)),new_costs=result['new_costs'],
  previously_paid_selected_shared_FD_RHS=result['previously_paid_selected_shared_FD_RHS'],
  controlled_wall_seconds=receipt['controlled_wall_seconds'],profile=receipt['profile'],small_groups=small_groups,task_groups=task_groups,
  paired_Heun_direct=len(paired),marginal_RHS_ratio=stats([z['marginal_RHS_ratio'] for z in paired]),
  shared_setup_inclusive_RHS_ratio=stats([z['shared_FD_and_small_h_inclusive_RHS_ratio'] for z in paired]),
  no_GT=True,continuous_path_certificate=False,posterior_MCSE='NOT_RUN',exactly_reversible_map=False,
  qualification='Heun removes leading Euler self-defect in tested regime; empirical convergence only; no global smooth O(h^3) at C1 seams; no theorem falsification or posterior success claim')
 for filename,data in [('OBSERVED_METRICS.json',summary),('PAIRED_HEUN_DIRECT.json',paired),('ALL_SKIPS.json',skips)]:
  (RESULT/filename).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
 print(json.dumps({k:v for k,v in summary.items() if k not in ('small_groups','task_groups')},indent=2))

if __name__=='__main__':main()
