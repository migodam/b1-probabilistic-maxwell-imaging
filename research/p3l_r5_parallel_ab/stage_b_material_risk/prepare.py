import sys,json,time,resource
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];OLD=ROOT.parent/'b1-p3l-r4-conditional-transport';sys.path.insert(0,str(OLD/'src'))
import numpy as np,torch
from material_moments import *
from p3r4.density import load_density
from p3r4.probability_eval import decoder_from_geometry,decode_gaussianized
from p3r4.risk_fusion import read_moments
from p3r4.fusion_provider_compat import load_fusion_compatible
from p3r4.risk import feature_image,roi_components,RiskActionPredictor
from p3r2.residual_risk import read_anchor
S=OLD/'research/p3l_r4_conditional_transport';torch.set_num_threads(1);start=time.monotonic()
def readj(path):
 with np.load(path,allow_pickle=False) as z:return {k:z[k] for k in ['parent','y','material']}
model=load_density(S/'models/density/extended_diag_gaussian_s2026102501/BEST.pt');cfg=json.loads((S/'configs/TARGET.json').read_text());points=json.loads((OLD/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'];decoder=decoder_from_geometry(cfg,points)
head=load_fusion_compatible(S/'models/risk/riskfusion_diag_gaussian_s2026102501/BEST.pt');g=head.geometry
stats=model.stats;F=np.asarray(stats['frame_inverse']);recs=[];meta=dict(provider_parameters=sum(z.numel() for z in model.model.parameters()),frame_cross_block_max=float(max(np.max(abs(F[:64,64:])),np.max(abs(F[64:,:64])))),frame_inverse_rule='x=(normalized_x*x_scale+x_mean)@frame_inverse',noise=cfg['noise'],variance_convention='R4 K32 empirical ddof0; deterministic population q variance, no31/32 correction',device='cpu',MPS_available=torch.backends.mps.is_available(),new_physics_RHS=0,source_root=str(S),decoder=decoder)
for fold,filename in [('train','train32768.npz'),('validation','validation.npz')]:
 joint=readj(S/'inputs/joint'/filename);a=read_anchor(S/f'inputs/strong_context_1024_256/anchor_{fold}.npz');lookup={x:i for i,x in enumerate(joint['parent'])};ids=np.array([lookup[x] for x in a['parent']]);y=joint['y'][ids];truth=joint['material'][ids];parents=joint['parent'][ids];rows=[];cross=0;wall=time.monotonic()
 for begin in range(0,len(y),32):
  m,C=raw_gaussian(model,y[begin:begin+32]);im,c=diagonal_image(m,C,decoder);rows.append(im);cross=max(cross,c)
 det=np.concatenate(rows);elapsed=time.monotonic()-wall
 old=read_moments(S/f'inputs/fusion_diag_gaussian/{fold}',parents);raw_image=feature_image(y,g,a,head.p['feature_normalization']);anchor=roi_components(a['chi0'],g['roi']);target=roi_components(truth,g['roi']);extras={}
 for seed in [2026102501,2026102502,2026102503]:
  h=load_fusion_compatible(S/f'models/risk/riskfusion_diag_gaussian_s{seed}/BEST.pt');base=RiskActionPredictor(S/f'models/risk/strongrisk_spectral_s{seed}/BEST.pt');extras[f'strong_s{seed}']=roi_components(base.predict(y,a,batch_size=64),g['roi'])
  for label,mom in [('K32',old),('det_frozen',det)]:extras[f'{label}_s{seed}']=roi_components(h.predict(y,a,mom,parents,batch_size=64,moment_provider='diag_gaussian',moment_K=32),g['roi'])
 np.savez_compressed(OUT/f'prepared_{fold}.npz',parent=parents,y=y,truth=truth,target=target,anchor=anchor,image=raw_image,det=det,K32=old,q0=a['q0'],**extras)
 recs.append(dict(fold=fold,parents=len(y),deterministic_seconds=elapsed,seconds_per_parent=elapsed/len(y),field_property_cov_max=cross,moment_K32_RMSE=float(np.sqrt(np.mean((old-det)**2)))))
 # Scalar high-K diagnostic, joint material covariance and exact source law.
 if fold=='validation':
  refs=[]
  for i in range(8):
   sums=np.zeros((4,576));n=0
   for chunk in range(32):
    sam=model.sample(y[i:i+1],4096,2026102690+i*100+chunk,batch_size=2048)[0];chi,_,_=decode_gaussianized(sam,decoder);sums+=np.stack([chi.real.sum(0),chi.imag.sum(0),(chi.real**2).sum(0),(chi.imag**2).sum(0)]);n+=len(chi)
   ref=sums/n;ref[2:]-=ref[:2]**2;ref=ref.reshape(4,24,24);refs.append(dict(parent=str(parents[i]),K=n,max_abs=float(np.max(abs(ref-det[i]))),RMSE=float(np.sqrt(np.mean((ref-det[i])**2)))))
  meta['provider_MC_reference']=refs
meta.update(records=recs,wall_seconds=time.monotonic()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
(OUT/'ASSET_INVENTORY.json').write_text(json.dumps(meta,indent=2)+'\n');torch.save(dict(geometry=g,reference_head_payload=head.p),OUT/'geometry.pt');print(json.dumps({k:v for k,v in meta.items() if k!='decoder'}))
