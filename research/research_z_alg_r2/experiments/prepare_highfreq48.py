"""Separate private image generator; not imported by solver processes."""
import json,sys,secrets,os,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))
from physics import Model

def main():
 cp=P/'experiments/HIGHFREQ48_CONTRACT.json'
 if cp.exists():raise RuntimeError('Frozen highfreq protocol exists')
 audit=json.loads((P/'experiments/HIGHFREQ48_GEOMETRY_AUDIT.json').read_text());g=audit['geometry']
 if not audit['condition_under_1e6']:raise RuntimeError('public geometry gate failed')
 c={'schema':'Z-ALG-R2-HIGHFREQ48-v1','frozen_utc':datetime.now(timezone.utc).isoformat(),'scope':'NEW high-frequency 48-cell fixed pulse support; all48 real cells unknown; not old low-frequency support, not Full66','geometry':g,'material_box':[.5,4.],'private_case_families':['smooth_two_blobs','piecewise_two_regions','independent_uniform'],'cases_per_family':1,'noise_modes':['clean','relative_complex_1e-3'],'noise_definition':'draw isotropic independent complex Gaussian direction, normalize its l2 norm to1e-3*norm(clean complex y); one fixed noisy observation shared by all arms','arms':{'exactSVD':'full fixed receiver SVD rank cutoff1e-12; raw shared-real constitutive estimate retained even if illegal; no optimization','FullJ_midpoint':'all48 material coordinates free, fixed midpoint2.25, bounded exact full-J TRF','SVD_FullJ':'same TRF, initialize raw exactSVD shared-real estimate clipped to[.5+1e-10,4-1e-10] without GT; no other tuning'},'SVD_parameters':{'rank_rtol':1e-12,'node_threshold':1e-20,'constitutive_rtol':1e-8},'optimizer':{'method':'trf','max_nfev':120,'ftol':1e-12,'xtol':1e-12,'gtol':1e-12,'x_scale':1.,'loss':'linear','residual':'all labeled observed channels divided by fixed observed target norm','J':'exact all-cell tangent solve48 columns vs96 receiver adjoints; cheaper chosen from public dimensions','trial_state':'LU plus source solves; no pertrial resolvent SVD','wall_seconds_per_arm_case':120},'validation':{'independent':'fresh numpy solve and full resolvent singular values only at final endpoints','clean_response_relative':1e-8,'state_relative':1e-10,'legal_tolerance':1e-10,'sigma_floor':1e-6,'clean_imaging_RMS':1e-5,'noise_evaluation':'report all material/response errors; no GT-tuned noise threshold'},'accounting':{'setup':'one shared geometry Green/SVD setup for all cases/arms, separately timed and added once to method-arm totals','incremental':'per-case solve+finalvalidation excludes shared geometry','RHS':'8 source RHS per forward;48 tangent RHS per J;8 fresh source RHS pervalidation','hardware':'single BLAS thread CPU','per_arm_case_complex_RHS':15000},'blindness':'computational observation-only process contract; no OS sandbox; separate generator/scorer/plotter','freeze_policy':'freeze implementations and all18predictions before any GT score; preserve all failed candidates; no tuning from GT'}
 cp.write_text(json.dumps(c,indent=2)+'\n');seed=secrets.randbits(128);rng=np.random.default_rng(seed);xy=np.asarray(g['centers']);x,y=xy.T;targets={};rows=[];model=Model(None,**g);t=time.perf_counter();cpu=time.process_time()
 for i,family in enumerate(c['private_case_families']):
  case=f'HF48_{i:02d}'
  if i==0:
   v=np.zeros(48)
   for _ in range(2):
    bx,by=rng.uniform([-.25,-.17],[.25,.17]);sx,sy=rng.uniform(.08,.18,2);v+=rng.uniform(.6,1.2)*np.exp(-((x-bx)**2/(2*sx*sx)+(y-by)**2/(2*sy*sy)))
   m=.7+2.9*v/max(float(v.max()),1e-12)
  elif i==1:
   m=np.full(48,rng.uniform(.6,1.0))
   for _ in range(2):
    bx,by=rng.uniform([-.2,-.15],[.2,.15]);wx,wy=rng.uniform([.12,.10],[.25,.20]);mask=(abs(x-bx)<wx)&(abs(y-by)<wy);m[mask]=rng.uniform(2.,3.8)
  else:m=rng.uniform(.5,4.,48)
  targets[case]=m;y0=model.forward(m,independent=True)['y'];z=rng.normal(size=y0.shape)+1j*rng.normal(size=y0.shape);noise=z/np.linalg.norm(z)*(.001*np.linalg.norm(y0))
  for mode,obs in [('clean',y0),('relative_complex_1e-3',y0+noise)]:
   name=f'{case}_{mode}';np.savez(P/f'experiments/highfreq48_public/{name}.npz',y=obs);rows.append({'case':case,'family':family,'noise_mode':mode,'input_id':name,'observation_file':f'{name}.npz','relative_noise_l2':0. if mode=='clean' else .001})
 np.savez(P/'private/highfreq48/sealed_targets.npz',**targets)
 secret=P/'private/highfreq48/generator_secret.json';secret.write_text(json.dumps({'seed':str(seed),'wall_seconds':time.perf_counter()-t,'CPU_seconds':time.process_time()-cpu,'source_rhs':24,'scope':'private generator only'},indent=2)+'\n');os.chmod(secret,0o600);os.chmod(P/'private/highfreq48/sealed_targets.npz',0o600)
 (P/'experiments/highfreq48_public/manifest.json').write_text(json.dumps({'schema':c['schema'],'cases':rows},indent=2)+'\n');print('Highfreq48 contract frozen;3 images6 observation-only inputs generated;truth/seed suppressed.')
if __name__=='__main__':main()
