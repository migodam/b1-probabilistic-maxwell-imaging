import sys,json,time
from pathlib import Path
OUT=Path(__file__).resolve().parent;sys.path.insert(0,str(OUT))
import numpy as np,torch
from inference import Action,OLD
from p3r2.residual_risk import read_anchor
torch.set_num_threads(1);d=dict(np.load(OUT/'prepared_validation.npz'));a=read_anchor(OLD/'research/p3l_r4_conditional_transport/inputs/strong_context_1024_256/anchor_validation.npz');records=[]
for method in ['deterministic_matched','distilled_w128_n32768']:
 for seed in [2026102601,2026102602,2026102603]:
  act=Action(method,seed);pred=act.predict(d['y'],a,batch_size=64);folder=OUT/'experiments'/f'deterministic_matched_s{seed}' if method=='deterministic_matched' else OUT/'learned_context'/f'distilled_w128_n32768_s{seed}';ref=np.load(folder/'validation.npz')['prediction'];records.append(dict(method=method,seed=seed,parent_count=256,max_abs_prediction_difference=float(np.max(abs(pred-ref))),absolute_mean_risk_difference=abs(float(np.mean(abs(pred-d['truth'])**2))-float(np.mean(abs(ref-d['truth'])**2)))))
r=dict(status='PASS' if max(x['max_abs_prediction_difference'] for x in records)<3e-6 else 'FAIL',tolerance=3e-6,scope='CPUinferenceentry vsstoredMPSvalidationpredictions allfinalists3seeds; existingVALnotnewtruth',records=records,physics_RHS=0);(OUT/'INFERENCE_REPLAY.json').write_text(json.dumps(r,indent=2)+'\n');print(r['status'],max(x['max_abs_prediction_difference'] for x in records))
