"""Publish local-reproduction training readiness from public anchor receipts."""
import json,sys,time
from pathlib import Path
S=Path(__file__).resolve().parents[1]
records=json.loads((S/'data/manifests/online.json').read_text())['records']
counts={'train':0,'val':0};missing=[];nonfull=[];wrong=[]
for r in records:
    if r['split']not in counts:continue
    p=S/'anchors'/f"{r['obs_id']}.json"
    if not p.exists()or not p.with_suffix('.npz').exists():missing.append(r['obs_id']);continue
    a=json.loads(p.read_text())
    if a['config'].get('discrepancy_stop',False):wrong.append(r['obs_id']);continue
    counts[r['split']]+=1
    if a['final_stage']!=24 or not a['all_unknowns_released']:nonfull.append(r['obs_id'])
ready=counts=={'train':1008,'val':156}and not missing and not nonfull and not wrong
out={'status':'TRAIN_VAL_STRICT_ANCHORS_READY'if ready else'HOLD_INCOMPLETE_PUBLIC_ANCHORS',
     **counts,'all_primary_full24':not nonfull,'all_primary_discrepancy_stop':False,
     'missing':missing,'nonfull':nonfull,'wrong_primary':wrong,'truth_access':'NONE',
     'time_unix':time.time(),'source':'local metadata-only release; no numerical physics revalidation'}
(S/'audit/DATA_READY.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
sys.exit(0 if ready else 1)
