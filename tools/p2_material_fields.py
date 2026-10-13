"""Posterior location/material uncertainty fields and properly weighted PPC."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.b1_joint.material import geometry
from src.b1_joint.response import Response
from src.b1_joint.package import predictive_check
from src.b1_joint.accounting import dump
def main():
    out=ROOT/'research/joint_p2/results';r=np.linspace(0,1,257);names=[];probs=[];means=[];sds=[];response=Response(out/'COMMON_RESPONSE_BANK.npz');rng=np.random.default_rng(261013219)
    for p in sorted((out/'packages').glob('*.json')):
        d=json.loads(p.read_text());a=np.load(p.with_name(p.stem+'_samples.npz'));x=a['latent'];w=a['weights'];g=geometry(x);v=(r[None,:]<g[:,0,None])|((g[:,1,None]>0)&(abs(r[None,:]-.65)<g[:,1,None]));field=np.where(v,1.,x[:,2,None]+1j*x[:,3,None]);m=w@field;vr=w@((field.real-m.real)**2);vi=w@((field.imag-m.imag)**2);names.append(p.stem);probs.append(w@v);means.append(m);sds.append(np.array([np.sqrt(vr),np.sqrt(vi)]))
        old=d['posterior_predictive_consistency'];d['posterior_predictive_consistency']=predictive_check(response,x,np.array(d['measurement']),np.array(d['noise_sigma']),rng,weights=w);d['posterior_predictive_consistency']['old_uniform_sample_list_diagnostic']=old;d['posterior_predictive_consistency']['sampling_correction']='use supplied whole-posterior weights; appended rare-event particles must not be over-weighted';d['material_location_probability_field']=dict(path='../MATERIAL_FIELDS.npz',case=p.stem,field='probability of exact vacuum at radial position; epsilon mean/real&imag SD',measure='whole modeled event masses + finite proxy event-conditional weighted particles',full_field_certificate=False);dump(p,d)
    np.savez_compressed(out/'MATERIAL_FIELDS.npz',case_ids=np.array(names),radius=r,vacuum_probability=np.array(probs),epsilon_mean=np.array(means),epsilon_SD_real_imag=np.array(sds))
    sensitivity=out/'PRIOR_NUISANCE_SENSITIVITY_REFINED.json';d=json.loads(sensitivity.read_text())
    for case in d['cases'].values():
        a=case['radius_uniform_gain015'];a['old_log_evidence_coefficient_nodes']=a['log_evidence'];a['log_evidence']=a['independent_coordinate_adaptive']['log_evidence']
    dump(sensitivity,d)
if __name__=='__main__':main()
