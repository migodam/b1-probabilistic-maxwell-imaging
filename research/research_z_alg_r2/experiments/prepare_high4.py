"""Freeze a NEW positive higher-contrast E2 support, then generate blind inputs.
Generator-only process; do not import this from solver code.
"""
import os,sys,json,time,secrets
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.stats import qmc
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))
from physics import Model

def main():
 contract_path=P/'experiments/HIGH4_CONTRACT.json'
 if contract_path.exists():raise RuntimeError('Frozen high4 contract exists')
 old=json.loads((P/'experiments/contract.json').read_text());g=old['presets']['E1_LOW6'].copy();g['centers']=g['centers'][:4];g.update(radius=.04,k=2*np.pi)
 midpoint=15.25;starts=[np.full(4,midpoint)]
 for x in qmc.Sobol(4,scramble=False).random_base2(4):
  m=.5+29.5*x
  if all(np.linalg.norm(m-y)>1e-12 for y in starts):starts.append(m)
  if len(starts)==8:break
 c={'schema':'Z-ALG-R2-E2-HIGH4-v1','frozen_utc':datetime.now(timezone.utc).isoformat(),'scope':'NEW model support and higher positive contrast regime; not a rescue/replacement of old [0.5,4] evidence; not Full66','geometry':g,'material_box':[.5,30.],'hidden_sampling':'3 independent uniform product targets, private cryptographic seed','noise':'none','initialization':{'midpoint':midpoint,'FullJ8_starts':[x.tolist() for x in starts],'policy':'fixed midpoint plus first7 distinct unscrambled Sobol4 points; no GT use'},'methods':['six_path_bilinear_dual','FullJ8_exact_TRF'],'homotopy':{'seed':2026102904,'tolerance':2e-11,'max_steps':4000,'initial_step':.01,'min_step':1e-9,'max_step':.04,'max_newton':12,'wall_seconds':60},'baseline':{'max_nfev_per_start':120,'ftol':1e-12,'xtol':1e-12,'gtol':1e-12,'loss':'linear','x_scale':1.,'residual_scaling':'fixed observed target norm','bounds':[.5,30.]},'validation':{'fresh_numpy_solve':True,'response_relative':1e-8,'state_relative':1e-10,'resolvent_sigma_min':1e-6,'merge_material_RMS':1e-7,'material_imaging_RMS':1e-5},'blindness':'computational input separation; shared filesystem is not an OS sandbox','failure_policy':'save all paths,roots,starts,rejections; no truth scores until implementation+run freeze; no tuning after scores','budget':{'method_case_complex_RHS':15000,'method_case_wall_seconds':120,'single_BLAS_thread':True}}
 contract_path.write_text(json.dumps(c,indent=2)+'\n')
 seed=secrets.randbits(128);rng=np.random.default_rng(seed);model=Model(None,**g);targets={};cases=[];t=time.perf_counter();cpu=time.process_time()
 for i in range(3):
  case=f'E2_HIGH4_{i:02d}';m=rng.uniform(.5,30.,4);targets[case]=m;state=model.forward(m,independent=True)
  np.savez(P/f'experiments/high4_public/{case}.npz',y=state['y']);cases.append({'case':case,'observation_file':f'{case}.npz'})
 np.savez(P/'private/high4/sealed_targets.npz',**targets)
 secret=P/'private/high4/generator_secret.json';secret.write_text(json.dumps({'seed':str(seed),'generation_wall_seconds':time.perf_counter()-t,'generation_CPU_seconds':time.process_time()-cpu,'source_rhs':6},indent=2)+'\n');os.chmod(secret,0o600);os.chmod(P/'private/high4/sealed_targets.npz',0o600)
 (P/'experiments/high4_public/manifest.json').write_text(json.dumps({'schema':c['schema'],'cases':cases},indent=2)+'\n');print('New high4 protocol frozen;3 observation-only cases generated;private values suppressed.')
if __name__=='__main__':main()
