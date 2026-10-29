"""Separate postfreeze GENERAL6 root-set scorer. Never import in solver code.
Run only after parent/engineer confirms every prediction is frozen.
"""
import argparse,json,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
P=Path(__file__).resolve().parents[1]
def main():
 a=argparse.ArgumentParser();a.add_argument('--prediction-freeze-receipt',required=True);args=a.parse_args();receipt=Path(args.prediction_freeze_receipt).resolve();output=P/'results/general6_scores.json'
 if not receipt.exists():raise RuntimeError('Need externally confirmed prediction freeze receipt')
 if output.exists():raise RuntimeError('GENERAL6 already scored; preserve original postfreeze score')
 gt=np.load(P/'private/general6/sealed_targets.npz');rows=[]
 for f in sorted((P/'results').glob('general6_*.json')):
  r=json.loads(f.read_text());case=r.get('case')
  if case not in gt.files:continue
  roots=[]
  for root in r.get('independently_validated',[]):
   v=root.get('validation',{})
   if root.get('status')=='QUALIFIED' or v.get('qualified') is True or v.get('response_fit_pass') is True:roots.append(root)
  errors=[float(np.linalg.norm(np.asarray(root['m'],dtype=float)-gt[case])/np.sqrt(6)) for root in roots]
  nearest=min(errors) if errors else None;rows.append({'file':f.name,'case':case,'method':r.get('method'),'solver_status':r.get('status'),'qualified_root_count':len(roots),'all_qualified_root_RMS_to_hidden_truth':errors,'minimum_root_RMS':nearest,'contains_GT_RMS1e5':bool(nearest is not None and nearest<=1e-5),'closest_truth_role':'evaluation of frozen root-set coverage only; not a solver-selected image','all_validation':[{k:v for k,v in root.get('validation',{}).items() if k!='m'} for root in roots],'method_wall_seconds':r.get('method_wall_seconds'),'method_CPU_seconds':r.get('method_CPU_seconds'),'accounting':r.get('accounting',{}),'GT_disclosed':False})
 if not rows:raise RuntimeError('No frozen GENERAL6 predictions found; no score emitted')
 result={'schema':'Z-ALG-R2-GENERAL6-SCORE-v1','status':'SCORED_AFTER_EXTERNAL_PREDICTION_FREEZE_CONFIRMATION','UTC':datetime.now(timezone.utc).isoformat(),'freeze_receipt':str(receipt),'rows':rows,'GT_values_disclosed':False,'GT_seed_disclosed':False,'scope':'Independent new6-cell3-source support;finite numerical root-set coverage,not rigorous completeness; noGT-driven tuning or root selection','generation_actual_complex_RHS':9};output.write_text(json.dumps(result,indent=2)+'\n');print(f'GENERAL6 scored {len(rows)} frozen root-set outputs;GT/seed suppressed.')
if __name__=='__main__':main()
