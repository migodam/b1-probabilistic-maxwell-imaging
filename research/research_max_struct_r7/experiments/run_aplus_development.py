"""Frozen disclosed four-object / three-readout development, one CPU process."""
from pathlib import Path
import json,sys,subprocess,time
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'research/research_max_struct_r7/delegated/aplus'
CASES=('obs_7b8d5df517f1497d','obs_b20d5a55677d4eca','obs_40ac5b79e8a74ac3','obs_00f185ff5525406b')

def main():
 OUT.mkdir(parents=True,exist_ok=True)
 manifest=OUT/'DEVELOPMENT_MANIFEST.json'
 if not manifest.is_file():
  source=json.loads((ROOT/'research/research_max_struct_r6/configs/primary_final_manifest.json').read_text())
  base=json.loads((ROOT/'research/research_max_struct_r7/configs/primary.json').read_text())
  tasks=[]
  for mode in ('original_a','mobius','projected'):
   cfg=dict(base,run_id='r7_aplus_development_'+mode,classification='DISCLOSED_DEVELOPMENT_NO_GT_SELECTION',aplus_options=dict(aplus_mode=mode,aplus_sweeps=2,aplus_projected_steps=8,aplus_step_tol=1e-8))
   config=OUT/(mode+'_CONFIG.json');config.write_text(json.dumps(cfg,indent=2)+'\n')
   for case in CASES:
    original=next(t for t in source['tasks'] if t['case_id']==case and t['arm']=='a_material')
    task=dict(original,task_id=case+'__'+mode,arm='a_material' if mode=='original_a' else 'a_plus',config=str(config),output_dir=str(OUT/'trajectories'/(case+'__'+mode)),development_mode=mode)
    tasks.append(task)
  manifest.write_text(json.dumps(dict(status='FROZEN_BEFORE_FULL_TRAJECTORIES',truth_access='NONE',selection_rule='Report all three variants. Choose A+ by disclosed public original-Phi retention and paid CPU only; never GT.',cases=list(CASES),tasks=tasks),indent=2)+'\n')
 data=json.loads(manifest.read_text());receipts=[]
 for task in data['tasks']:
  out=Path(task['output_dir'])
  if out.exists():
   receipt=json.loads((out/'RECEIPT.json').read_text())
   if receipt['status']!='SEALED_PUBLIC_PREDICTION':raise RuntimeError('Existing incomplete task; preserve and audit before retry '+str(out))
  else:
   log=OUT/(task['task_id']+'.log')
   with log.open('w') as stream:
    result=subprocess.run([sys.executable,str(ROOT/'research/research_max_struct_r7/experiments/run_case.py'),'--manifest',str(manifest),'--task-id',task['task_id']],cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
   if result.returncode:raise RuntimeError('Failed charged development task '+str(log))
   receipt=json.loads((out/'RECEIPT.json').read_text())
  row=dict(task_id=task['task_id'],case_id=task['case_id'],mode=task['development_mode'],status=receipt['status'],cpu_seconds=receipt['cpu_seconds'],wall_seconds=receipt['wall_seconds'],auxiliary_cpu_seconds=receipt['auxiliary_cpu_seconds'],original_Phi=receipt['physical_objective'],physical_relative_residual=receipt['physical_relative_residual'],physical_cost=receipt['physical_model_cost'],auxiliary_cost=receipt['auxiliary_execution_receipt'],output_dir=str(out))
  receipts.append(row)
  (OUT/'DEVELOPMENT_RESULTS.json').write_text(json.dumps(dict(status='COMPLETE' if len(receipts)==len(data['tasks']) else 'PARTIAL',truth_access='NONE',rows=receipts),indent=2)+'\n')
  print(json.dumps({k:row[k] for k in ('task_id','cpu_seconds','auxiliary_cpu_seconds','original_Phi','physical_relative_residual')}),flush=True)
if __name__=='__main__':main()
