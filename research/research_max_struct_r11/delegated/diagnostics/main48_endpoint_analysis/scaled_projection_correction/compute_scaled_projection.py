"""Pure saved-public-m/g arithmetic; no physics, GT, optimizer or runtime probe."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
import csv,json,time
import numpy as np
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;ROOT=HERE.parents[3]
METHODS=('cie_single','cie_dual','global_dual','fullG19','freezeA','contractedT')
def read(p):return json.loads(Path(p).read_text())
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def write(p,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with Path(p).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)
def public_path(p):
 p=Path(p)
 if not p.is_absolute() or not p.is_relative_to(ROOT/'RESULTS') or {'sealed_truth','data_private','private'} & set(p.parts):raise ValueError('Only explicit public RESULTS files allowed')
 return p
def main():
 wall=time.perf_counter();source=ROOT/'delegated/twod/independent48/endpoint_bookkeeping_corrected_path/ENDPOINT_BOOKKEEPING.json';end=read(source)
 with (BASE/'JOINED_FINAL288.csv').open(newline='') as f:old=list(csv.DictReader(f))
 lookup={(r['case_id'],r['method']):r for r in old}
 if len(end)!=288 or len(lookup)!=288 or {(e['case_id'],e['method']) for e in end}!=set(lookup):raise ValueError('288 exact public endpoint rows required')
 rows=[];merged=[]
 for e in end:
  if not e['available']:raise RuntimeError('No exact saved gradient, do not use a replacement')
  d=public_path(e['output_dir']);gp=public_path(e['gradient_source']);capfile=d/'R11_CAP_PREDICTION.npz'
  with np.load(capfile,allow_pickle=False) as z:m=z['m'].copy()
  with np.load(gp,allow_pickle=False) as z:mg=z['m'].copy();g=z['gradient'].copy()
  if not np.array_equal(m,mg) or m.shape!=(576,) or g.shape!=m.shape or np.iscomplexobj(g) or not np.isfinite(g).all():raise ValueError('Requires elementwise exact cap/saved gradient identity')
  u=float(np.max(np.abs(m-np.clip(m-g,.02,20.))));s=float(np.max(np.abs(m-np.clip(m-1e6*g,.02,20.))));unscaled_pass=u<=1e-8;scaled_pass=s<=.01
  prior=lookup[e['case_id'],e['method']];matches={Path(c['path']).name:c.get('material_exact') for c in e['candidate_match_records']}
  row=dict(task_id=e['task_id'],case_id=e['case_id'],method=e['method'],sealed_incumbent_unscaled_projected_mapping_inf=u,unscaled_mapping_le_1e_minus8=unscaled_pass,sealed_incumbent_scaled_projected_mapping_inf=s,scaled_mapping_le_solver_gtol_point01=scaled_pass,criteria_disagree=unscaled_pass!=scaled_pass,scaled_pass_but_unscaled_exceeds=scaled_pass and not unscaled_pass,unscaled_pass_but_scaled_exceeds=unscaled_pass and not scaled_pass,raw_gradient_inf=float(np.max(np.abs(g))),original_corrector_stop=prior['corrector_stop'],PGTOL_stop_string='NORM OF PROJECTED GRADIENT <= PGTOL' in prior['corrector_stop'],PGTOL_stop_but_sealed_incumbent_scaled_exceeds=('NORM OF PROJECTED GRADIENT <= PGTOL' in prior['corrector_stop']) and not scaled_pass,sealed_cap_matches_ACTUAL=matches.get('R11_ACTUAL_PAID_INCUMBENT.npz'),sealed_cap_matches_LAST_PAID_STATE=matches.get('R11_LAST_PAID_STATE.npz'),optimizer_terminal_point_identity_established=False,optimizer_terminal_point_note='No result.x material exported by refine_observed; LAST_PAID_STATE is only last paid evaluation. No terminal-point equality inferred.',gradient_source=str(gp),cap_material_source=str(capfile),unscaled_difference_from_prior_bookkeeping=u-e['projected_gradient_inf'],scope='Original solver-type scaled mapping evaluated at sealed minPhi incumbent; not a reconstructed terminal stopping event or new success gate')
  rows.append(row);r=dict(prior);r.update(row);r['stationarity_scope']='SUPERSEDES_V1_LABEL: 1e-8 mapping is unscaled descriptor only; actual scaled mapping<=.01 is a separate sealed-incumbent quantity, not proof of optimizer-terminal identity';merged.append(r)
 write(HERE/'SCALED_UNSCALED_ALL288.csv',rows);dump(HERE/'SCALED_UNSCALED_ALL288.json',rows);write(HERE/'JOINED_FINAL288_CORRECTED_PROJECTION_LABELS.csv',merged)
 disagree=[r for r in rows if r['criteria_disagree']];write(HERE/'CRITERION_DISAGREEMENT_CASES.csv',disagree);dump(HERE/'CRITERION_DISAGREEMENT_CASES.json',disagree)
 summaries=[]
 for method in METHODS:
  rr=[r for r in rows if r['method']==method];mm=[r for r in rr if r['criteria_disagree']];qualified=[r for r in rr if r['scaled_mapping_le_solver_gtol_point01']]
  summaries.append(dict(method=method,objects=48,unscaled_descriptor_le_1e_minus8=sum(r['unscaled_mapping_le_1e_minus8'] for r in rr),scaled_mapping_le_point01_at_sealed_incumbent=sum(r['scaled_mapping_le_solver_gtol_point01'] for r in rr),criteria_disagreements=len(mm),scaled_pass_unscaled_fail=sum(r['scaled_pass_but_unscaled_exceeds'] for r in rr),unscaled_pass_scaled_fail=sum(r['unscaled_pass_but_scaled_exceeds'] for r in rr),max_unscaled_mapping_all=max(r['sealed_incumbent_unscaled_projected_mapping_inf'] for r in rr),max_unscaled_mapping_while_scaled_passes=max((r['sealed_incumbent_unscaled_projected_mapping_inf'] for r in qualified),default=None),max_unscaled_mapping_among_disagreements=max((r['sealed_incumbent_unscaled_projected_mapping_inf'] for r in mm),default=None),PGTOL_stop_strings=sum(r['PGTOL_stop_string'] for r in rr),PGTOL_strings_but_incumbent_scaled_fails=sum(r['PGTOL_stop_but_sealed_incumbent_scaled_exceeds'] for r in rr),sealed_matches_last_paid=sum(r['sealed_cap_matches_LAST_PAID_STATE'] is True for r in rr),optimizer_terminal_identity_established=0,original_stop_counts=dict(Counter(r['original_corrector_stop'] or 'UNRECORDED' for r in rr))))
 write(HERE/'METHOD_SUMMARY_CORRECTED.csv',summaries);dump(HERE/'METHOD_SUMMARY_CORRECTED.json',summaries)
 receipt=dict(status='COMPLETE_ZERO_FW_SCALED_PROJECTION_CORRECTION',created_utc=datetime.now(timezone.utc).isoformat(),rows=288,criteria_disagreements=len(disagree),all_original_v1_unscaled_numbers_reproduced=all(r['unscaled_difference_from_prior_bookkeeping']==0 for r in rows),maximum_unscaled_among_scaled_passes=max(r['sealed_incumbent_unscaled_projected_mapping_inf'] for r in rows if r['scaled_mapping_le_solver_gtol_point01']),optimizer_terminal_identity_established=0,GT_arrays_read=0,new_FW=0,new_adjoint=0,new_solver=0,new_scoring=0,new_bootstrap=0,runtime_version_probe=False,SciPy_runtime_version=None,CPU_seconds_near_write=time.process_time(),wall_seconds=time.perf_counter()-wall)
 dump(HERE/'RECEIPT.json',receipt);print(json.dumps(summaries));print(json.dumps(receipt))
if __name__=='__main__':main()
