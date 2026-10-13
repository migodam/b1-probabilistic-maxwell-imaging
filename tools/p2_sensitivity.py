"""Different prior/nuisance identities, same paid fields; no algorithm ranking."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.b1_joint.material import RMAX,SCALES,geometry
from src.b1_joint.response import Response
from src.b1_joint.reference import Reference
from src.b1_joint.likelihood import log_likelihood
from src.b1_joint.accounting import dump
def main():
    out=ROOT/'research/joint_p2/results';response=Response(out/'COMMON_RESPONSE_BANK.npz');ref=Reference(response);cases=json.loads((out/'DEVELOPMENT_CASES.json').read_text())['cases'];rows={}
    for item in cases:
        y=np.array(item['y'])[:32];sigma=np.array(item['sigma'])[:32];r=ref.solve(y,sigma,order=17);bank=r['bank'];x=bank['x'];g=geometry(x);ratio=np.ones(len(x))
        for j in range(2):
            active=g[:,j]>0;v=g[active,j];original=v/SCALES[j]**2*(1-(v/SCALES[j])**2)**-3;ratio[active]*=(.5/RMAX[j])/original
        variants={}
        for name,gs,lr in [('coefficient_uniform_gain015',.15,0.),('radius_uniform_gain015',.15,np.log(ratio)),('coefficient_uniform_gain010',.1,0.),('coefficient_uniform_gain020',.2,0.)]:
            lw=log_likelihood(bank['F1'],y,sigma,gs)+np.log(bank['prior'])+lr;z=logsumexp(lw);w=np.exp(lw-z);p=np.bincount(bank['events'],weights=w,minlength=12)
            variants[name]=dict(target_identity='joint-p2-sensitivity-v1:'+name,posterior_event_mass=p,log_evidence=float(z),prior_change=name.startswith('radius'),noise_unchanged=True,gain_prior=dict(sd=gs,bounds=[.7,1.3],truncation_normalizer_retained=True),numerical_status='same continuous proxy fields; prior change kept in separate identity')
        rows[item['id']]=variants
    dump(out/'PRIOR_NUISANCE_SENSITIVITY.json',dict(cases=rows,new_physics=0,method_comparison=False,not_noise_tuning=True,prior_measure='radius-uniform variant retains .5 zero-radius atom per geometry coordinate; positive radius is .5/Rmax, eps unchanged',numerical_limit='integrated by17th-order coefficient-node rule; reweighting near r=0 has integrable cusp, independent refinement required for strong radius-prior claims'))
    comp=json.loads((out/'COMPENSATION.json').read_text());a=comp['local_projection'][:3];rad=np.array([v['radius'] for v in a]);raw=np.array([v['raw'] for v in a]);projected=np.array([v['projected'] for v in a]);dump(out/'LOCAL_MECHANISM_STATUS.json',dict(fit_interval=[.02,.08],raw_log_slope=float(np.polyfit(np.log(rad),np.log(raw),1)[0]),projected_log_slope=float(np.polyfit(np.log(rad),np.log(projected),1)[0]),projection_scope='white unconstrained epsilonReal/epsilonImag/shared-gain Jacobian at fixedsolid anchor, both frequencies',leading_term_cancellation_supported=False,meaning='projection reduces amplitude but observed small-radius response retains approximately cubic scaling; no universal r^5 mechanism claim',new_physics=0))
if __name__=='__main__':main()
