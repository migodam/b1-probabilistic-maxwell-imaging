"""Independent stored-evidence checks, no Maxwell solves or model selection."""
from core import *
import csv,datetime

def main():
 idx=load_setup()['ROI'];issues=[];checks=[];lastcharge=0.;rawseals=[]
 stop=json.loads((P/'audit/PHYSICS_STOP.json').read_text())['utc']
 for f in sorted((P/'results').glob('*/*/*/PREDICTION_SEALED.json')):
  m=json.loads(f.read_text());z=dict(np.load(f.parent/'PREDICTION.npz'));y=np.load(P/'data/online'/f"{m['parent']}.npz")['y'];a=z['occupancy'];chi=z['material'];delta=complex(z['delta']);expected=np.full(576,HOST);expected[idx]+=delta*a[idx]
  good=dict(parent=m['parent'],method=m['method'],array_shapes=bool(a.shape==(576,) and chi.shape==(576,) and z['F'].shape==(32,8)),finite=bool(all(np.isfinite(z[k]).all() for k in ['x','occupancy','material','F'])),occupancy_legal=bool(a.min()>=0 and a.max()<=1),property_legal=bool(.8<=delta.real<=1.8 and .03<=delta.imag<=.18),same_material_map=float(np.max(abs(expected-chi))),passivity=bool(chi.imag.min()>=HOST.imag-1e-14),Phi_recomputed=phi(z['F'],y),Phi_receipt=m['Phi'],prediction_blake2b_matches=blake(f.parent/'PREDICTION.npz')==m['prediction_blake2b'],GT_used=m['GT_used'])
  if not all(good[k] for k in ['array_shapes','finite','occupancy_legal','property_legal','passivity','prediction_blake2b_matches']) or good['GT_used'] or abs(good['Phi_recomputed']-m['Phi'])>1e-8:issues.append(good)
  checks.append(good);rawseals.append(f.stat().st_mtime)
 for line in (P/'COST_LEDGER.jsonl').read_text().splitlines():
  q=json.loads(line);lastcharge=max(lastcharge,float(q['utc'] or 0))
 score=json.loads((P/'EXPERIMENT_RESULTS.json').read_text())['GT_scoring_utc']
 records=dict(status='STORED_EVIDENCE_PASS_WITH_RESOURCE_PROTOCOL_FAILURE' if not issues else 'BLOCKED',new_physics_RHS=0,new_LU=0,predictions=len(checks),checks=checks,issues=issues,last_paid_event_utc=lastcharge,STOP_utc=stop,no_paid_event_after_STOP=lastcharge<=stop,all_prediction_receipts_before_GT_score=max(rawseals)<score,all_blind_records_before_score=all(json.loads((P/'audit'/n).read_text())['recorded_utc']<score for n in ['BLIND_READ_N.json','BLIND_READ_P.json']),resource_protocol_status='FAIL_CAP_EXCEEDED')
 write(P/'audit/FINAL_EVIDENCE_QA.json',records)
 print(json.dumps({k:v for k,v in records.items() if k!='checks'},indent=2))
if __name__=='__main__':main()
