"""Post-freeze root-set evaluation only; no selection feedback to solvers."""
import argparse,json,time
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parents[2]
def main():
 a=argparse.ArgumentParser();a.add_argument('--freeze-receipt',required=True);args=a.parse_args();receipt=Path(args.freeze_receipt)
 if not receipt.exists():raise RuntimeError('external algorithm+prediction freeze required')
 gt=np.load(P/'private/high4/sealed_targets.npz');rows=[]
 for f in sorted((P/'results').glob('high4_*.json')):
  r=json.loads(f.read_text());case=r.get('case')
  if case not in gt.files:continue
  roots=[x for x in r.get('independently_validated',[]) if x.get('status')=='QUALIFIED' or x.get('validation',{}).get('qualified')]
  errors=[float(np.linalg.norm(np.asarray(x['m'])-gt[case])/2) for x in roots]
  rows.append({'file':f.name,'case':case,'method':r.get('method'),'qualified_roots':len(roots),'all_root_material_RMS':errors,'minimum_root_RMS':min(errors) if errors else None,'contains_GT_within_RMS_1e5':bool(errors and min(errors)<=1e-5),'root_selection':'none; nearest truth is evaluation of coverage only, not solver output selection','GT_disclosed':False})
 out={'status':'SCORED_AFTER_FREEZE','freeze_receipt':str(receipt),'scored_unix':time.time(),'rows':rows,'scope':'finite numerical root-set coverage; not proof of completeness; no GT-selected image claim'}
 (P/'results/high4_scores.json').write_text(json.dumps(out,indent=2)+'\n');print(f'Scored {len(rows)} high4 frozen root-set outputs; truth values suppressed.')
if __name__=='__main__':main()
