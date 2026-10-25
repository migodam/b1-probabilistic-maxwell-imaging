import sys,json,time,resource
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];OLD=ROOT.parent/'b1-p3l-r4-conditional-transport';sys.path.insert(0,str(OLD/'src'))
import numpy as np,torch
from material_moments import *
from p3r4.density import load_density
from p3r4.field_moments import load_field_moments
from p3r4.probability_eval import decode_gaussianized
S=OLD/'research/p3l_r4_conditional_transport';torch.set_num_threads(1);data=dict(np.load(OUT/'prepared_validation.npz'));inventory=json.loads((OUT/'ASSET_INVENTORY.json').read_text());decoder=inventory['decoder'];rows=[]
for provider in ['gaussian_K32','gaussian_deterministic','no_density_direct']:
 start=time.monotonic()
 model=(load_field_moments(S/'models/field/extended_mlp_small_moments_s2026102501/BEST.pt') if provider=='no_density_direct' else load_density(S/'models/density/extended_diag_gaussian_s2026102501/BEST.pt'))
 load=time.monotonic()-start
 def call(y):
  if provider=='gaussian_K32':
   x=model.sample(y,32,2026102621,batch_size=2048);chi,_,_=decode_gaussianized(x.reshape(-1,66),decoder);chi=chi.reshape(len(y),32,576);return np.stack((chi.real.mean(1),chi.imag.mean(1),chi.real.var(1),chi.imag.var(1)),1)
  if provider=='gaussian_deterministic':
   m,C=raw_gaussian(model,y);return diagonal_image(m,C,decoder,64)[0]
  return model.predict_moments(y,batch_size=64)
 for batch in [1,64]:
  y=data['y'][:batch];start=time.monotonic();call(y);first=time.monotonic()-start;times=[]
  for rep in range(30):
   start=time.monotonic();call(y);times.append(time.monotonic()-start)
  rows.append(dict(provider=provider,batch=batch,cold_load_seconds=load,cold_first_inference_seconds=first,warm_median_seconds=float(np.median(times)),warm_per_parent_seconds=float(np.median(times)/batch),warm_p10_seconds=float(np.quantile(times,.1)),warm_p90_seconds=float(np.quantile(times,.9))))
 # numerical-order validation across all256 independentexistingVAL parents
 if provider=='gaussian_deterministic':
  errors=[]
  for b in range(0,len(data['y']),32):
   m,C=raw_gaussian(model,data['y'][b:b+32]);a=diagonal_image(m,C,decoder,64)[0];ref=diagonal_image(m,C,decoder,128)[0];errors.append(np.max(abs(a-ref)))
  quadrature=dict(nodes=64,reference_nodes=128,parents=256,max_abs_float32_difference=float(max(errors)))
 if provider=='no_density_direct':direct_count=sum(v.numel() for v in model.model.parameters())
res=dict(status='COMPLETE',device='CPU1',rows=rows,production_quadrature_check=quadrature,direct_provider_parameter_count=direct_count,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,physics_RHS=0,cold_scope='checkpointload+firstinference; excludes common interpretertorchimport',warm_scope='providerencoder+sampling/integration+fourchannels; excludes shared paid anchor acquisition and riskhead')
(OUT/'COST_BENCHMARK.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res))
