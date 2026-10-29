"""Release holdout GT only after all63 prediction snapshots and source checks."""
import json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
 checks=[];predictions=[]
 for split,expected in [('holdout',27),('budget_filled',27),('data_adjoint',9)]:
  p=ROOT/f'audit/EXPERIMENT_PREDICTIONS_FREEZE_{split.upper()}.json'
  if not p.exists():raise RuntimeError('Missing '+str(p))
  receipt=json.loads(p.read_text())
  if receipt['predictions']!=expected:raise RuntimeError('Wrong prediction count')
  for f in receipt['files']:
   equal=(ROOT/f['source']).read_bytes()==(ROOT/f['snapshot']).read_bytes();predictions.append(dict(split=split,observation_id=f['observation_id'],method=f['method'],live_equals_snapshot=equal))
 primary=json.loads((ROOT/'audit/HOLDOUT_FREEZE.json').read_text())
 for f in primary['files']:checks.append(dict(group='primary16',path=f['path'],byte_equal=(ROOT/f['path']).read_bytes()==(ROOT/f['snapshot']).read_bytes()))
 for name in ['BUDGET_FILL','DATA_ADJOINT_CONTROL']:
  folder=ROOT/f'audit/{name}_SOURCE_SNAPSHOT';receipt=json.loads((folder/'FREEZE.json').read_text())
  for f in receipt['files']:checks.append(dict(group=name,path=f['path'],byte_equal=(ROOT/f['path']).read_bytes()==(folder/f['path']).read_bytes()))
 if not all(c['byte_equal'] for c in checks) or not all(p['live_equals_snapshot'] for p in predictions):raise RuntimeError('Source or prediction drift; do not read GT')
 release=dict(released_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),primary_prediction_count=27,budget_filled_prediction_count=27,ordinary_control_prediction_count=9,total_prediction_snapshots=63,all_predictions_snapshot_equal=True,all_sources_snapshot_equal=True,truth_read_outside_generator_before_all_predictions=False,primary_frozen_before_Y=True,supplements_frozen_after_Y_before_GT=True,source_checks=checks,prediction_checks=predictions)
 (ROOT/'audit/EXPERIMENT_ALL_PREDICTIONS_FROZEN_GT_RELEASE.json').write_text(json.dumps(release,indent=2)+'\n')
 env=os.environ.copy()
 for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:env[key]='1'
 commands=[['experiments/score_predictions.py','score','holdout'],['experiments/score_budget_filled.py','score'],['experiments/score_data_adjoint.py','score'],['experiments/aggregate_delivery.py'],['experiments/render_attribution_comparison.py']]
 for i,args in enumerate(commands):
  r=subprocess.run([sys.executable,*args],cwd=ROOT,env=env,capture_output=True,text=True);(ROOT/f'audit/EXPERIMENT_POSTSCORE_{i}.log').write_text(r.stdout+'\nSTDERR\n'+r.stderr)
  print(json.dumps(dict(command=args,exit_code=r.returncode,stdout=r.stdout[-5000:],stderr=r.stderr[-2000:])),flush=True)
  if r.returncode:raise RuntimeError('Postscore failure; see retained log')
if __name__=='__main__':main()
