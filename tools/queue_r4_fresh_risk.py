#!/usr/bin/env python3
"""Explicit released zero-physics queue after parent generation/context closure.

Does not bind/generate observations, acquire contexts, train, or run physics.
'plan' is metadata only. 'run' requires a distinct immutable parent release.
"""
from pathlib import Path
import argparse,json,subprocess,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.evaluation import checkpoint_identity,write_new_json
from r4_fresh_risk import bound
S=ROOT/'research/p3l_r4_conditional_transport'
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['plan','run']);p.add_argument('--cohort',required=True,choices=['locked','confirmation'])
 for k in ['component','top-freeze','freeze','release','out']:p.add_argument('--'+k,type=Path)
 p.add_argument('--stage',default='IV');a=p.parse_args();c,top=bound(a);parents=c['cohorts'][a.cohort]['parents'];dest=a.out
 base=[sys.executable,str(ROOT/'tools/r4_fresh_risk.py')];common=['--cohort',a.cohort,'--component',str(a.component),'--top-freeze',str(a.top_freeze),'--freeze',str(a.freeze)]
 steps=[('features',base+['features',*common,'--observations',str(S/'inputs/online'),'--out',str(dest/'features')]),('infer',base+['infer',*common,'--contexts',str(S/'results/contexts'/a.cohort),'--native',str(S/'results/contexts/details'),'--features',str(dest/'features'),'--out',str(dest/'online')]),('export',base+['export',*common,'--online',str(dest/'online'),'--observations',str(S/'inputs/online'),'--out',str(dest/'exports')]),('score',base+['score',*common,'--online',str(dest/'online'),'--truth',str(S/'inputs/offline_only'),'--out',str(dest/'scores')])]
 if a.command=='plan':
  print(json.dumps(dict(status='PLAN_ONLY_NO_EXECUTION',steps=steps,required_parents=len(parents),bind_before_generation=True,no_auto_physics=True)));return
 release=json.loads(a.release.read_text())
 if release.get('released') is not True or release.get('experiment')!='fresh_risk_panels_v1' or release.get('cohort')!=a.cohort or release.get('component_identity')!=checkpoint_identity(a.component) or release.get('top_freeze_identity')!=checkpoint_identity(a.top_freeze) or release.get('allow_features_infer') is not True:raise ValueError('separate parent zero-physics release required')
 if not a.freeze.exists():raise ValueError('parent must bind risk before observation generation')
 for parent in parents:
  if not (S/'inputs/online'/(parent+'.npz')).exists() or not (S/'results/contexts'/a.cohort/(parent+'.npz')).exists():raise ValueError('all fixed online observations/contexts must be present; never generate or replace missing parents')
 dest.mkdir(parents=True,exist_ok=True);write_new_json(dest/'QUEUE_FROZEN.json',dict(component_identity=checkpoint_identity(a.component),top_freeze_identity=checkpoint_identity(a.top_freeze),release_identity=checkpoint_identity(a.release),steps=steps,parents=parents,new_physics=0,new_fits=0))
 for label,cmd in steps:
  if label=='score' and release.get('allow_offline_score') is not True:continue
  wrapped=[sys.executable,str(ROOT/'tools/run_r4_job.py'),'--stage',a.stage,'--job',f'fresh_risk_{a.cohort}_{label}_v1','--max-wall','1200','--rss-limit-mib','3072','--',*cmd]
  subprocess.run(wrapped,check=True,cwd=ROOT)
 write_new_json(dest/'QUEUE_COMPLETE.json',dict(status='COMPLETE',cohort=a.cohort,scored=release.get('allow_offline_score') is True,new_physics=0,new_fits=0))
if __name__=='__main__':main()
