import sys,json,time,resource
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];OLD=ROOT.parent/'b1-p3l-r4-conditional-transport';sys.path.insert(0,str(OLD/'src'))
import numpy as np,torch
from material_moments import *
from p3r4.density import load_density
from p3r4.field_moments import load_field_moments
from p3r4.probability_eval import decode_gaussianized
from p3r4.risk_fusion import parent_seed,FusionRisk
from p3r4.fusion_provider_compat import load_fusion_compatible
from p3r4.risk import RiskActionPredictor,feature_image,roi_components
from p3r2.residual_risk import read_anchor
S=OLD/'research/p3l_r4_conditional_transport';torch.set_num_threads(1);d=dict(np.load(OUT/'prepared_validation.npz'));decoder=json.loads((OUT/'ASSET_INVENTORY.json').read_text())['decoder'];model=load_density(S/'models/density/extended_diag_gaussian_s2026102501/BEST.pt');fm=load_field_moments(S/'models/field/extended_mlp_small_moments_s2026102501/BEST.pt');anchor=read_anchor(S/'inputs/strong_context_1024_256/anchor_validation.npz');rows=[]
# Numerical speed audit predeclared tolerance1e-6 maxabsolute in all4channels versus128nodes.
qref=[]
for begin in range(0,256,32):
 m,C=raw_gaussian(model,d['y'][begin:begin+32]);ref=diagonal_image(m,C,decoder,128)[0]
 for q in [16,32,64]:qref.append(dict(nodes=q,max_abs=float(np.max(abs(diagonal_image(m,C,decoder,q)[0]-ref)))))
qaudit={str(q):max(r['max_abs'] for r in qref if r['nodes']==q) for q in [16,32,64]};qpass=[q for q in [16,32,64] if qaudit[str(q)]<=1e-6];chosen=min(qpass)
heads={}
for arm in ['deterministic_matched','direct_moments_matched']:
 p=torch.load(OUT/'experiments'/f'{arm}_s2026102601/BEST.pt',map_location='cpu',weights_only=False);net=FusionRisk(p['geometry']);net.load_state_dict(p['state_dict']);net.eval();heads[arm]=(p,net)
original=load_fusion_compatible(S/'models/risk/riskfusion_diag_gaussian_s2026102501/BEST.pt');strong=RiskActionPredictor(S/'models/risk/strongrisk_spectral_s2026102501/BEST.pt')

def fourdirect(y):
 r=fm.predict_moments(y,batch_size=64);mean=r['material_mean'];cv=r['psd_pointwise_covariance'];return np.stack((mean.real,mean.imag,cv[:,:,0,0],cv[:,:,1,1]),1).reshape(-1,4,24,24)

def provider(y,parents,label):
 if label=='original_K32_exact_parent_seed':
  ims=[]
  for yy,parent in zip(y,parents):
   x=model.sample(yy[None],32,parent_seed(parent,20261025121),batch_size=8)[0];chi,_,_=decode_gaussianized(x,decoder);ims.append(np.stack((chi.real.mean(0),chi.imag.mean(0),chi.real.var(0),chi.imag.var(0))).reshape(4,24,24))
  return np.array(ims,dtype='float32')
 if label.startswith('det_'):
  m,C=raw_gaussian(model,y);return diagonal_image(m,C,decoder,int(label.split('_')[1]))[0]
 return fourdirect(y).astype('float32')

for batch in [1,64]:
 y=d['y'][:batch];par=d['parent'][:batch];anc={k:v[:batch] for k,v in anchor.items() if k in ['observation','residual','q0','chi0','J']}
 for label in ['original_K32_exact_parent_seed','det_16','det_32','det_64','direct']:
  if label=='original_K32_exact_parent_seed':p,net=original.p,original.model
  else:p,net=heads['direct_moments_matched' if label=='direct' else 'deterministic_matched']
  payload=p['feature_payload'] if 'feature_payload' in p else p;g=p['geometry']
  def call():
   moments=provider(y,par,label);image=feature_image(y,g,anc,payload['feature_normalization']);base=roi_components(anc['chi0'],g['roi']);moments=(moments-p['moment_mean'])/p['moment_scale']
   with torch.inference_mode():return net(torch.tensor(image),torch.tensor(moments),torch.tensor(base)).numpy()
  begin=time.monotonic();call();cold=time.monotonic()-begin;times=[]
  for rep in range(12):begin=time.monotonic();call();times.append(time.monotonic()-begin)
  rows.append(dict(arm=label,batch=batch,first_inference_seconds=cold,warm_median_seconds=float(np.median(times)),warm_per_parent_seconds=float(np.median(times)/batch)))
 # Common strong reference completehead+alloriginal online featureimage assembly.
 def basecall():return strong.predict(y,anc,batch_size=64)
 basecall();times=[]
 for rep in range(12):begin=time.monotonic();basecall();times.append(time.monotonic()-begin)
 rows.append(dict(arm='strong_reference',batch=batch,warm_median_seconds=float(np.median(times)),warm_per_parent_seconds=float(np.median(times)/batch)))
 for arm in ['raw_y','decoder_only','physical_residual','same_spectrum_random']:
  fit=dict(np.load(OUT/'low_rank_correction'/f'{arm}_s2026102601/BEST.npz'))
  def rrrcall():
   base=strong.predict(y,anc,batch_size=64)
   if arm in ['raw_y','same_spectrum_random']:raw=y.astype(float)
   elif arm=='decoder_only':raw=np.concatenate((roi_components(anc['chi0'],strong.geometry['roi']).reshape(len(y),-1),anc['q0']),1).astype(float)
   else:
    im=feature_image(y,strong.geometry,anc,strong.p['feature_normalization']);raw=np.concatenate((im[:,2:6].reshape(len(y),-1),anc['q0']),1).astype(float)
   if arm=='same_spectrum_random':X=(((raw-fit['raw_mean'])/fit['raw_scale']@fit['random_projection'])-fit['projection_mean'])@fit['spectrum_factor']
   else:X=(raw-fit['raw_mean'])/fit['raw_scale']@fit['projection']
   corr=((X-fit['input_mean'])@fit['input_factor']@fit['output_factor']+fit['output_mean'])*24
   return base,corr
  rrrcall();times=[]
  for rep in range(12):begin=time.monotonic();rrrcall();times.append(time.monotonic()-begin)
  rows.append(dict(arm='RRR_'+arm,batch=batch,warm_median_seconds=float(np.median(times)),warm_per_parent_seconds=float(np.median(times)/batch),stored_float_scalar_count=int(sum(v.size for v in fit.values() if v.dtype.kind=='f')),baseline_parameters=5842))
res=dict(status='COMPLETE',device='CPU1',scope='ALL providerencoder+momentgeneration+originalyonlyfeatureassembly+6418head; RRRstrongbaseline+fullinputtransform+RRR; excludes paidanchor acquisition common toallarms and interpreter imports',rows=rows,quadrature_numerical_audit=dict(reference_nodes=128,max_abs_tolerance=1e-6,all256VAL_4channels_errors=qaudit,fastest_numerically_qualified_nodes=chosen,feature_only_no_truth_read=True),physics_RHS=0,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
(OUT/'FULL_COST_BENCHMARK.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res))
