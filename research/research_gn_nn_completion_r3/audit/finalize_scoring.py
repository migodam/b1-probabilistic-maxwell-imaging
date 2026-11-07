"""Seal every predeclared prediction before first TEST truth scoring."""
from pathlib import Path
import json,os,subprocess,sys,time
S=Path(__file__).resolve().parents[1];R=S.parents[1]
sys.path.insert(0,str(R));os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['VECLIB_MAXIMUM_THREADS']='1'
while not ((S/'audit/CONFIRMATION_PHYSICS_DONE.json').exists()and(S/'audit/SECONDARY_DONE.json').exists()):time.sleep(30)
from research.research_gn_nn_completion_r3.src.score import seal_predictions,score_predictions
plan=json.loads((S/'audit/CONFIRMATION_PLAN_FROZEN.json').read_text())
runs=[p['run']for p in plan['plans']]+plan['controls']
core_seal=seal_predictions(S,S/'data/manifests/evaluation.json',runs,splits=['test'],output=S/'evaluation/CONFIRMATION_SEALED.json')
secondary=json.loads((S/'audit/SECONDARY_DONE.json').read_text())
second_seal=None
if secondary['status']=='COMPLETED':
    manifest=json.loads((S/'data/manifests/evaluation.json').read_text());ids=set(secondary['obs_ids'])
    manifest['records']=[r for r in manifest['records']if r['obs_id']in ids]
    subset=S/'data/manifests/secondary_evaluation.json';subset.write_text(json.dumps(manifest,indent=2))
    second_runs=[p['run']for p in secondary['plans']]+['confirm1024_v_physics_seed7','confirm1024_j_null_seed7']+plan['controls']
    second_seal=seal_predictions(S,subset,second_runs,splits=['test'],output=S/'evaluation/SECONDARY_SEALED.json')
# Both cohorts are now closed before any TEST material or current is opened.
score_predictions(core_seal,output=S/'evaluation/CONFIRMATION_SCORES.json',include_oracle=True)
if second_seal:score_predictions(second_seal,output=S/'evaluation/SECONDARY_SCORES.json',include_oracle=False)
for source,dest in [('CONFIRMATION_SCORES.json','results'),('SECONDARY_SCORES.json','results/secondary')]:
    if (S/'evaluation'/source).exists():
        subprocess.run([sys.executable,str(S/'src/summarize_results.py'),'--scores',str(S/'evaluation'/source),'--output',str(S/dest)],cwd=S,check=True)
subprocess.run([sys.executable,str(S/'src/cost_ledger.py')],cwd=S,check=True)
(S/'audit/FINAL_SCORING_DONE.json').write_text(json.dumps(dict(unix=time.time(),test_truth='OPENED_ONLY_AFTER_ALL_PLANNED_COHORTS_SEALED',core=str(core_seal),secondary=None if second_seal is None else str(second_seal)),indent=2))
