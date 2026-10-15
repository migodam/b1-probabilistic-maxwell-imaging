"""One full-validated pair, shared across the whole registered action menu."""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
REPO=Path(__file__).resolve().parents[3];sys.path.insert(0,str(REPO))
from cost import Cost,ROOT,dump
from prepare_e1 import forward
from src.b1_joint.material import from_geometry,void_volume
from src.b1_joint.response import Response
from sensing import field_l1

def run(cost):
    bank=np.load(ROOT/'results/E1_NOISE_ANCHORS.npz');ks=list(bank['frequencies']);anchor=bank['F'];sd=bank['sigma']
    cal_sd=sd[ks.index(1.6)]/np.linalg.norm(anchor[ks.index(1.6)]);response=Response(ROOT/'results/E1_RESPONSE_BANK.npz')
    ref=from_geometry(0,0)
    def residual(p):
        f=response.predict(from_geometry(.08,0,p[0],p[1]),ks).reshape(6,32)
        return np.r_[((p[2]*f-anchor)/sd[:,None]).ravel(),(p[2]-1)/cal_sd]
    opt=least_squares(residual,[2.5,.1,1.],bounds=([2.2,.04,.7],[2.8,.16,1.3]),max_nfev=100,ftol=1e-12,xtol=1e-12,gtol=1e-12)
    x=from_geometry(.08,0,*opt.x[:2]);paid=forward(cost);before=cost.count()
    full=np.array([paid(x,k,kind='SMALL_CAVITY_ONE_WORLD_ALL_MENU')[0] for k in ks])
    assert cost.count()-before<=6
    d=np.linalg.norm((opt.x[2]*full-anchor)/sd[:,None],axis=1);cal=abs(opt.x[2]-1)/cal_sd
    context=abs(float(void_volume(x)-void_volume(ref)))/.01
    d0=d[ks.index(1.6)];maxd=max(float(d.max()),cal,context);D2=d0*d0+3*maxd*maxd
    result=dict(status='ONE_PAIR_FULL_VALIDATED_ALL_REGISTERED_ACTIONS',reference=ref,candidate=x,shared_gain=float(opt.x[2]),
        frequency_white_distances={f'{k:g}':float(v) for k,v in zip(ks,d)},calibration_white_distance=cal,external_geometry_white_distance=context,
        full_response=full,reference_response=anchor,sigma=sd,material_component_L1=float(field_l1(ref,x)[0,0]),
        optimizer_success=bool(opt.success),optimizer_nfev=opt.nfev,new_attempts=cost.count()-before,
        max_additional_blocks=3,gaussian_chain_KL_upper=.5*D2,two_world_equal_design_prior_error_lower=.5*(1-min(1,np.sqrt(D2)/2)),
        bound_scope='Numerical full-Mie fixed two-world experiment, known independent thermal noise and lawful persistent per-world gains; chain-rule KL/Pinsker covers adaptive menu choice, not a full-space posterior or hardware certificate',
        not_a_global_minimum=True,unknown_world_mass=None)
    dump(ROOT/'results/SMALL_CAVITY_HORIZON_WITNESS.json',result)
    np.savez_compressed(ROOT/'results/SMALL_CAVITY_HORIZON_WITNESS.npz',reference_latent=ref,candidate_latent=x,F=full,F_reference=anchor,gain=opt.x[2],sigma=sd)
    print('all-menu small-cavity witness',d,'error lower',result['two_world_equal_design_prior_error_lower'],'new',cost.count()-before)

if __name__=='__main__':
    with Cost().job('e2-smallcavity-horizon-audit-v1','E2') as cost:run(cost)
