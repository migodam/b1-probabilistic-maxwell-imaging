"""Offline-only fixed digit8 generator; all shape templates remain outside online inputs."""
from pathlib import Path
import importlib.util,time,json
import numpy as np
from digit_stress import ROOT,EX,TX,TY,FREQ,SIGMA,GridEval,dump
p=EX/'digit_patterns.py';spec=importlib.util.spec_from_file_location('digit_patterns',p);patterns=importlib.util.module_from_spec(spec);spec.loader.exec_module(patterns)
def main():
 private=EX/'private';private.mkdir(exist_ok=True);online=EX/'online';online.mkdir(exist_ok=True)
 if (private/'truth.npz').exists():raise RuntimeError('No truth regeneration allowed')
 r=np.random.default_rng(33106008);mask=np.flipud(patterns.get_pattern(8,padded=True)).astype(bool);we=r.uniform(2,4.5,6);wl=r.uniform(.02,.12,6);tc=np.where(mask,r.uniform(2.5,4.0,(7,7)),0.);tl=np.where(mask,r.uniform(.02,.12,(7,7)),0.);z=np.r_[we,wl,tc.ravel(),tl.ravel()];model=GridEval();clean,diag=model.run(z)
 np.savez_compressed(private/'truth.npz',z=z,mask=mask,clean=clean,target_contrast=tc,target_loss=tl,target_x_edges=TX,target_y_edges=TY)
 r=np.random.default_rng(33106108);y=clean+SIGMA*(r.normal(size=clean.shape)+1j*r.normal(size=clean.shape));np.savez_compressed(online/'observation.npz',y=y)
 dump(online/'config.json',dict(model='independent piecewise constant target material, all49cells free; unknownwall12',target_x_edges=TX,target_y_edges=TY,frequencies=FREQ,TX_count=4,RX_count=16,sigma_real_imag=SIGMA,target_chi_box=[0,8],target_loss_box=[0,.25],no_shape_or_digit_input=True))
 dump(EX/'GENERATION_RECEIPT.json',dict(time_unix=time.time(),generated_glyph_ids=list(range(1,10)),PDE_generated_ids=[8],reconstructed_ids=[],same_discrete_generator_inverse=True,independent_finer_validation=False,truth_private=True,public_observation_keys=['y'],diagnostics=diag))
if __name__=='__main__':main()
