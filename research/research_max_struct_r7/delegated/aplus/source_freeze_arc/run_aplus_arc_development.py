"""Final structural test: state stationary proposals with unchanged prior MM acceptance."""
from pathlib import Path
from datetime import datetime,timezone
import json,sys,subprocess
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'research/research_max_struct_r7/delegated/aplus'

def main():
 manifest=OUT/'ARC_DEVELOPMENT_MANIFEST.json'
 if not manifest.is_file():
  old=json.loads((OUT/'DEVELOPMENT_MANIFEST.json').read_text());cfg=json.loads((OUT/'mobius_CONFIG.json').read_text());cfg.update(run_id='r7_aplus_development_mobius_arc_candidates',classification='FINAL_STRUCTURAL_QUADRATIC_ARC_PROPOSAL_COMPARISON_NO_GT');cfg['aplus_options']['aplus_inherited_gtol']=True;cfg['aplus_options']['aplus_arc_candidates']=True
  config=OUT/'mobius_arc_candidates_CONFIG.json';config.write_text(json.dumps(cfg,indent=2)+'\n');tasks=[]
  for original in old['tasks']:
   if original['development_mode']!='mobius':continue
   case=original['case_id'];tasks.append(dict(original,task_id=case+'__mobius_arc_candidates',config=str(config),output_dir=str(OUT/'trajectories'/(case+'__mobius_arc_candidates')),development_mode='mobius_arc_candidates'))
  data=dict(status='FROZEN_BEFORE_REFINEMENT_TRAJECTORIES',freeze_utc=datetime.now(timezone.utc).isoformat(),truth_access='NONE',justification='Replace regularized stationary quintic roots by state stationary quadratic proposals, endpoints and old coordinate; accept lowest complete prior-MM surrogate. This is an approximate regularized conditional proposal. Unchanged exact objective/prior/two-sweep maximum/inherited gtol.',selection_rule='Complete all four predetermined disclosed cases under original budgets, with the two old bars first. Interpret full original Phi, public observation/noise residual contract, order-of-magnitude changes and paid CPU together; no arbitrary relative-Phi threshold or sample omission. No R7 primary results or GT. No further hyperparameter variants.',tasks=tasks)
  manifest.write_text(json.dumps(data,indent=2)+'\n')
 data=json.loads(manifest.read_text());rows=[]
 for task in data['tasks']:
  out=Path(task['output_dir'])
  if not out.exists():
   with (OUT/(task['task_id']+'.log')).open('w') as stream:
    result=subprocess.run([sys.executable,str(ROOT/'research/research_max_struct_r7/experiments/run_case.py'),'--manifest',str(manifest),'--task-id',task['task_id']],cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
   if result.returncode:raise RuntimeError('Failed charged refinement task '+task['task_id'])
  receipt=json.loads((out/'RECEIPT.json').read_text())
  if receipt['status']!='SEALED_PUBLIC_PREDICTION':raise RuntimeError('Incomplete refinement task '+task['task_id'])
  row=dict(task_id=task['task_id'],case_id=task['case_id'],mode=task['development_mode'],status=receipt['status'],cpu_seconds=receipt['cpu_seconds'],wall_seconds=receipt['wall_seconds'],auxiliary_cpu_seconds=receipt['auxiliary_cpu_seconds'],original_Phi=receipt['physical_objective'],physical_relative_residual=receipt['physical_relative_residual'],physical_cost=receipt['physical_model_cost'],auxiliary_cost=receipt['auxiliary_execution_receipt'],output_dir=str(out));rows.append(row)
  (OUT/'ARC_DEVELOPMENT_RESULTS.json').write_text(json.dumps(dict(status='COMPLETE' if len(rows)==len(data['tasks']) else 'PARTIAL',truth_access='NONE',rows=rows),indent=2)+'\n')
  print(json.dumps({k:row[k] for k in ('task_id','cpu_seconds','auxiliary_cpu_seconds','original_Phi','physical_relative_residual')}),flush=True)
if __name__=='__main__':main()
