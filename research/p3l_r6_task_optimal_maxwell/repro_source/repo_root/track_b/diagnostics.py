"""Finite trained-model swap and feature controls; historical VAL only."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
import sys,json,time,resource
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from inference import Action
from run import OUT,HIST,OLD,SEEDS,ARMS,write
import numpy as np,torch

def main():
 torch.set_num_threads(1);start=time.monotonic();va=dict(np.load(HIST/'prepared_validation.npz'));out={};controls={};replay={}
 for arm in ['original_feature32768']+ARMS:
  if arm!='original_feature32768' and not all((OUT/'models'/f'{arm}_s{s}'/'RESULT.json').exists() for s in SEEDS):continue
  matrix=[];preds=[];swapped=[]
  for i,pi in enumerate(SEEDS):
   row=[];row_predictions=[]
   for j,hj in enumerate(SEEDS):
    action=Action(arm,hj,provider_seed=pi);mom=action.moments(va['y']);prediction=action.predict_prepared(va['image'],va['anchor'],mom);risk=np.mean(abs(prediction-va['truth'])**2,1);row.append(float(risk.mean()));row_predictions.append(prediction);np.savez_compressed(OUT/'diagnostics'/f'{arm}_p{pi}_h{hj}.npz',parent=va['parent'],prediction=prediction,risk=risk)
    if i==j:
     preds.append(prediction);cc={};mm=action.hp['moment_mean'];base=mm+np.zeros_like(mom)
     for control in ['mean_only','variance_only','constant','parent_shuffle']:
      mc=mom.copy()
      if control=='mean_only':mc[:,2:]=base[:,2:]
      elif control=='variance_only':mc[:,:2]=base[:,:2]
      elif control=='constant':mc=base
      else:mc=mc[np.random.default_rng(2026102701).permutation(len(mc))]
      cp=action.predict_prepared(va['image'],va['anchor'],mc);cc[control]=float(np.mean(abs(cp-va['truth'])**2));np.savez_compressed(OUT/'diagnostics'/f'{arm}_s{pi}_{control}.npz',parent=va['parent'],risk=np.mean(abs(cp-va['truth'])**2,1))
     controls[f'{arm}_s{pi}']=cc
   matrix.append(row);swapped.append(row_predictions)
  P9=np.asarray(swapped);grand=P9.mean((0,1));pe=P9.mean(1)-grand;he=P9.mean(0)-grand;interaction=P9-P9.mean(1)[:,None]-P9.mean(0)[None]+grand;totalrisk=float(np.mean(abs(P9-va["truth"][None,None])**2));grandrisk=float(np.mean(abs(grand-va["truth"])**2));pv=float(np.mean(abs(pe)**2));hv=float(np.mean(abs(he)**2));iv=float(np.mean(abs(interaction)**2))
  M=np.asarray(matrix);P=np.asarray(preds);ensemble=P.mean(0);variance=np.mean(abs(P-ensemble)**2);out[arm]=dict(provider_seed_by_head_seed_risk=matrix,diagonal_mean=float(np.trace(M)/3),offdiagonal_mean=float((M.sum()-np.trace(M))/6),row_range=float(np.ptp(M.mean(1))),column_range=float(np.ptp(M.mean(0))),ensemble_risk=float(np.mean(abs(ensemble-va['truth'])**2)),mean_seed_risk=float(np.mean(abs(P-va['truth'][None])**2)),finite_ensemble_prediction_variance=float(variance),finite_bias_variance_identity_gap=float(np.mean(abs(P-va['truth'][None])**2)-np.mean(abs(ensemble-va['truth'])**2)-variance),finite_crossed_decomposition=dict(all9_mean_risk=totalrisk,grand_prediction_risk=grandrisk,provider_prediction_variation=pv,head_prediction_variation=hv,interaction_prediction_variation=iv,identity_gap=totalrisk-grandrisk-pv-hv-iv),interpretation='finite trained-model identity; no population bias/variance inference')
 # Actual Action loading with original full VAL anchors, all candidates and seeds.
 with np.load(OLD/'research/p3l_r4_conditional_transport/inputs/strong_context_1024_256/anchor_validation.npz') as z:
  online={k:z[k] for k in ['observation','residual','q0','chi0','J']};parents=z['parent']
 assert np.array_equal(parents,va['parent'])
 for arm in out:
  for seed in SEEDS:
   action=Action(arm,seed);pr=action.predict(va['y'],online);folder=HIST/'learned_context'/f'distilled_w128_n32768_s{seed}' if arm=='original_feature32768' else OUT/'models'/f'{arm}_s{seed}';stored=np.load(folder/'validation.npz')['prediction'];replay[f'{arm}_s{seed}']=dict(max_abs_prediction_difference=float(np.max(abs(pr-stored))),risk=float(np.mean(abs(pr-va['truth'])**2)))
 write(OUT/'SEED_COUPLING.json',dict(status='COMPLETE',seeds=SEEDS,data='historical selection VAL256 only',records=out,physics_RHS=0));write(OUT/'CHANNEL_CONTROLS.json',dict(status='COMPLETE',controls=controls,note='inference ablations with matched head frozen; no retraining; parent_shuffle destroys alignment across parents',physics_RHS=0));write(OUT/'ACTION_REPLAY.json',dict(status='COMPLETE',records=replay,scope='CPU Action versus stored MPS predictions; drift quantified; no invented acceptance threshold'))
 write(OUT/'jobs/diagnostics.json',dict(status='COMPLETE',device='cpu',CPU_threads=1,wall_seconds=time.monotonic()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,physics_RHS=0))
if __name__=='__main__':
 (OUT/'diagnostics').mkdir(exist_ok=True);main()
