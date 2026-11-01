"""One audit of the newly wired joint derivative and profile implementation."""
from pathlib import Path
import json,sys,time
import numpy as np
from scipy.linalg import lu_solve
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
from src.physics import MaxwellModel,make_grid
from src.current_normal import spatial_modes,normal_projection
from src.joint_current_model import JointCurrentModel
from src.profile_schur import profile_schur,release_candidate


def main():
    started=time.perf_counter()
    model=MaxwellModel(make_grid(4))
    rng=np.random.default_rng(202611010531)
    m=np.linspace(.6,1.4,model.n)
    anchor=model.state(m)
    V=rng.normal(size=(model.n,model.ns))+1j*rng.normal(size=(model.n,model.ns))
    V*=.025
    j=anchor['j']+V
    Y=anchor['Y']+.01*(rng.normal(size=anchor['Y'].shape)+1j*rng.normal(size=anchor['Y'].shape))
    modes=spatial_modes(np.linalg.qr(rng.normal(size=(model.n,2))+1j*rng.normal(size=(model.n,2)))[0],model.ns)
    cfg=dict(chart='cie',beta=.2,rho=1.,data_scale=np.linalg.norm(Y),state_weight=1/np.linalg.norm(model.E))
    joint=JointCurrentModel(model,Y,cfg)
    out=joint.linearize(m,j,modes)
    h=rng.normal(size=model.n)*.1
    c=rng.normal(size=modes.shape[2])*.1
    current=np.einsum('npd,d->np',modes,c)
    t=1e-6
    fd=(joint.residual(m+t*h,j+t*current)['b']-joint.residual(m-t*h,j-t*current)['b'])/(2*t)
    exact=out['B']@h+out['A']@c
    derivative=float(np.linalg.norm(fd-exact)/np.linalg.norm(exact))
    # An interior frozen model; no clipping, so the stated step identities apply.
    prof=profile_schur(out['b'],out['B'],out['A'][:,:12],Lambda=.03*np.eye(model.n))
    AC=out['A'][:,12:]
    hard=release_candidate(prof,AC,'hard')
    direct=profile_schur(out['b'],out['B'],out['A'],Lambda=.03*np.eye(model.n))
    dressing=AC+prof['A']@rng.normal(size=(12,12))
    dressed=release_candidate(prof,dressing,'hard')
    transform=np.diag(np.linspace(.3,1.7,12))
    soft=release_candidate(prof,AC,'soft',eta_target=.5)
    softmix=release_candidate(prof,AC@transform+prof['A']@rng.normal(size=(12,12)),'soft',eta_target=.5)
    def error(a,b):return float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(a),np.linalg.norm(b)))
    # Real-material normal annihilates the genuine all-source Maxwell tangent.
    T=model.tangent(m,h)[0]
    LV=T-m[:,None]*(model.G@T)
    normal,_=normal_projection(anchor['Q'],LV)
    R=Y-anchor['Y'];ref=anchor
    grad={}
    for name,beta in [('data',None),('feedback',.2),('globaldual',0.)]:
        value,g=model.data_value_gradient(m,Y,beta,ref)
        fp=model.data_value_gradient(m+t*h,Y,beta,ref)[0]
        fm=model.data_value_gradient(m-t*h,Y,beta,ref)[0]
        dd=(fp-fm)/(2*t)
        grad[name]=dict(relative_error=float(abs(dd-g@h)/max(abs(dd),abs(g@h),1e-12)))
    results=dict(scope='new production wiring; supplied algebra checks separately replayed once',
        joint_CIE_direction_relative_error=derivative,
        physical_tangent_normal_relative=float(np.linalg.norm(normal)/np.linalg.norm(LV)),
        hard_H_error=error(hard['H'],direct['H']),hard_g_error=error(hard['g'],direct['g']),
        hard_step_error=error(hard['p'],direct['p']),
        dressing_step_error=error(hard['p'],dressed['p']),
        soft_coordinate_step_error=error(soft['p'],softmix['p']),
        scalar_objective_gradients=grad,wall=time.perf_counter()-started,cost=dict(model.cost))
    checks=[derivative<1e-6,results['physical_tangent_normal_relative']<1e-9,
            results['hard_H_error']<1e-9,results['hard_g_error']<1e-9,
            results['hard_step_error']<1e-9,results['dressing_step_error']<1e-9,
            results['soft_coordinate_step_error']<1e-9]+[r['relative_error']<1e-6 for r in grad.values()]
    results['status']='PASS' if all(checks) else 'FAIL'
    (HERE/'audit/E0_PRODUCTION_WIRING.json').write_text(json.dumps(results,indent=2))
    print(json.dumps(results,indent=2))


if __name__=='__main__':main()
