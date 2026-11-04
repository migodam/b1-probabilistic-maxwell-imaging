"""One declared Möbius refinement: inherited projected-gradient tolerance."""
from pathlib import Path
from datetime import datetime,timezone
import json,sys,subprocess
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'research/research_max_struct_r7/delegated/aplus'

def main():
 manifest=OUT/'GUARDED_DEVELOPMENT_MANIFEST.json'
 if not manifest.is_file():
  old=json.loads((OUT/'DEVELOPMENT_MANIFEST.json').read_text());cfg=json.loads((OUT/'mobius_CONFIG.json').read_text());cfg.update(run_id='r7_aplus_development_mobius_guarded',classification='ONE_DECLARED_ENGINEERING_REFINEMENT_INHERITED_GTOL_NO_GT');cfg['aplus_options']['aplus_inherited_gtol']=True
  config=OUT/'mobius_guarded_CONFIG.json';config.write_text(json.dumps(cfg,indent=2)+'\n');tasks=[]
  for original in old['tasks']:
   if original['development_mode']!='mobius':continue
   case=original['case_id'];tasks.append(dict(original,task_id=case+'__mobius_guarded',config=str(config),output_dir=str(OUT/'trajectories'/(case+'__mobius_guarded')),development_mode='mobius_guarded'))
  data=dict(status='FROZEN_BEFORE_REFINEMENT_TRAJECTORIES',freeze_utc=datetime.now(timezone.utc).isoformat(),truth_access='NONE',justification='Earlier candidates continue work below inherited L-BFGS projected-gradient tolerance. Honor exactly material_inner_gtol=1e-8 before roots and after each complete sweep; unchanged two-sweep maximum.',selection_rule='Use original public Phi and paid complete CPU only. Retention summary uses per-object Phi ratios; do not infer GT recovery. Prefer guarded Mobius over unguarded only if old-bar public Phi retained and aggregate charged CPU reduced, reporting any adverse control-object difference.',tasks=tasks)
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
  (OUT/'GUARDED_DEVELOPMENT_RESULTS.json').write_text(json.dumps(dict(status='COMPLETE' if len(rows)==len(data['tasks']) else 'PARTIAL',truth_access='NONE',rows=rows),indent=2)+'\n')
  print(json.dumps({k:row[k] for k in ('task_id','cpu_seconds','auxiliary_cpu_seconds','original_Phi','physical_relative_residual')}),flush=True)
if __name__=='__main__':main()
