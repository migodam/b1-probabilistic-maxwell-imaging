"""Isolated post-freeze truth scorer. Does not expose GT arrays."""
import argparse,json,time
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parents[1]
def main():
 a=argparse.ArgumentParser();a.add_argument('--freeze-receipt',required=True);args=a.parse_args();receipt=Path(args.freeze_receipt)
 if not receipt.exists():raise RuntimeError('Need external solver/run freeze receipt before truth scoring')
 gt=np.load(P/'private/sealed_targets.npz');rows=[]
 for path in sorted((P/'results').glob('blind_*.json')):
  if path.name=='blind_scores.json':continue
  r=json.loads(path.read_text())
  if r.get('case') not in gt.files or 'm' not in r:continue
  m=np.asarray(r['m']);truth=gt[r['case']];rms=float(np.linalg.norm(m-truth)/np.sqrt(len(m)))
  rows.append({'file':path.name,'case':r['case'],'method':r['method'],'material_RMS':rms,'material_MAE':float(np.mean(abs(m-truth))),'material_max_absolute_error':float(np.max(abs(m-truth))),'imaging_pass':bool(rms<=1e-5 and r['validation']['response_fit_pass']),'validation':r['validation'],'accounting':r.get('accounting',{})})
 result={'status':'SCORED_AFTER_EXTERNAL_FREEZE','freeze_receipt':str(receipt),'scored_unix':time.time(),'rows':rows,'GT_disclosed':False}
 (P/'results/blind_scores.json').write_text(json.dumps(result,indent=2)+'\n');print(f'Scored {len(rows)} frozen results; GT values not disclosed.')
if __name__=='__main__':main()
