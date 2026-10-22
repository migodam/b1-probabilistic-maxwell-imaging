"""Pure wrong-force MH controls and saved-array physical-proxy derivatives.

Reads only paid DEV030 anchor arrays, online y, geometry/config. Never
constructs Target or Book, calls a Maxwell backend, trains, or reads GT.
"""
from pathlib import Path
import importlib.util
import ast
import json
import sys
import numpy as np
from scipy.linalg import solve_triangular,lu_factor

ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'src'))
from p3l.hmc import PreconditionedHMC
from p3l.kernels import ChainState
from p3l.core import to_probit
from p1_material.prior import MaterialPrior

OUT=ROOT/'research/p3l_local_posterior_imaging'
spec=importlib.util.spec_from_file_location('C_surrogate_controls',OUT/'C/surrogate_code.py')
proxy=importlib.util.module_from_spec(spec);spec.loader.exec_module(proxy)


def flipped(q,p,C,h,steps,force):
    q=q.copy();p=p.copy();p+=.5*h*force(q)
    for i in range(steps):
        q+=h*(C@p);p+=(.5 if i==steps-1 else 1)*h*force(q)
    return q,-p


def fd(fun,x,h=1e-5):
    E=np.eye(len(x))*h
    return np.asarray([(fun(x+e)-fun(x-e))/(2*h) for e in E])


def main():
    runner=ast.parse((ROOT/'tools/run_p3l_surrogate_hmc.py').read_text())
    for node in ast.walk(runner):
        if not isinstance(node,ast.Call) or not isinstance(node.func,ast.Attribute):continue
        assert node.func.attr!='phi_gradient'
        assert not (node.func.attr=='gradient' and isinstance(node.func.value,ast.Name) and node.func.value.id=='post')
    rng=np.random.default_rng(2026102349)
    A=rng.normal(size=(9,66))*.04;H=np.eye(66)+A.T@A
    D=rng.normal(size=(66,66))*.006;C=np.eye(66)*.7+D@D.T
    skew=np.zeros((66,66));skew[0,1]=.04;skew[1,0]=-.04
    offset=rng.normal(size=66)*.03
    # Deliberately wrong, even nonconservative, deterministic fixed force.
    wrong=lambda q:-.8*q+.05*np.sin(q)+skew@q+offset
    target=lambda q:-.5*float(q@H@q)
    q=rng.normal(size=66)*.3;p=rng.normal(size=66)
    q1,p1=flipped(q,p,C,.025,8,wrong)
    q2,p2=flipped(q1,p1,C,.025,8,wrong)
    reversal=max(float(np.max(np.abs(q2-q))),float(np.max(np.abs(p2-p))))
    joint=lambda a,b:target(a)-.5*float(b@C@b)
    ratio=joint(q1,p1)-joint(q,p)
    flux_error=abs((joint(q,p)+min(0.,ratio))-(joint(q1,p1)+min(0.,-ratio)))
    position=np.r_[q,p]
    def mapping(v):return np.r_[flipped(v[:66],v[66:],C,.025,4,wrong)]
    eps=1e-5;E=np.eye(132)*eps
    jac=np.column_stack([(mapping(position+e)-mapping(position-e))/(2*eps) for e in E])
    LU,pivots=lu_factor(jac);diagonal=np.diag(LU)
    logdet=float(np.log(np.abs(diagonal)).sum())
    sign=float(np.prod(np.sign(diagonal))*(-1)**np.sum(pivots!=np.arange(132)))
    seed=2026102350;oracle_rng=np.random.default_rng(seed)
    momentum=solve_triangular(np.linalg.cholesky(C).T,oracle_rng.normal(size=66),lower=False)
    expected_q,expected_p=flipped(q,momentum,C,.025,8,wrong)
    expected_ratio=joint(expected_q,expected_p)-joint(q,momentum)
    accept=np.log(max(oracle_rng.random(),np.finfo(float).tiny))<min(0.,expected_ratio)
    calls=dict(full_target=0,surrogate_force=0)
    def counted_target(x):calls['full_target']+=1;return target(x)
    def counted_force(x):calls['surrogate_force']+=1;return wrong(x)
    transition=PreconditionedHMC(counted_target,counted_force,.025,C,8).step(
        ChainState(q,target(q),wrong(q)),np.random.default_rng(seed))
    endpoint_error=float(np.max(np.abs(transition.state.q-(expected_q if accept else q))))
    energy_error=float(abs(transition.metadata['log_mh_ratio']-expected_ratio))
    # Physical proxy checks use only arrays from the previously paid anchor.
    anchor_path=OUT/'A/anchors/DEV_030.npz'
    with np.load(anchor_path,allow_pickle=False) as z:
        a={k:z[k].copy() for k in ('q','chi','prediction','field','receiver_T')}
    if a['receiver_T'].shape!=(32,576):raise ValueError('paid receiver transfer array required')
    cfg=json.loads((OUT/'configs/TARGET.json').read_text())
    points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
    prior=MaterialPrior(cfg,points);sigma=cfg['noise']['sigma']
    with np.load(OUT/'inputs/online/DEV_030.npz',allow_pickle=False) as z:y=z['observation'].copy()
    x0=to_probit(a['q'],initialization=True)
    Cchi=np.einsum('qp,ps->qsp',a['receiver_T'],a['field']).reshape(256,576)
    material=proxy.MaterialResponseForce(prior,x0,a['chi'],a['prediction'],Cchi,y,sigma)
    delta=rng.normal(size=576)*.01+1j*rng.normal(size=576)*.01
    direct=np.einsum('rp,ps,p->rs',a['receiver_T'],a['field'],delta)
    packing_error=float(np.max(np.abs(direct-(Cchi@delta).reshape(32,8))))
    center_q=proxy.native(x0);complex_J=Cchi@prior.derivative(center_q)
    Jx=np.column_stack([proxy.white(complex_J[:,j].reshape(32,8),sigma) for j in range(66)])
    Jx*=proxy.probit_derivative(x0)[None,:]
    # If inward initialization adjusted a boundary, use the material proxy's
    # arithmetic Fhat(x0) as the linear control center; no new Full is claimed.
    linear=proxy.LinearResponseForce(x0,material.prediction(x0),Jx,y,sigma)
    x=x0+rng.normal(size=66)*.015
    gradient_results={}
    for name,force in [('linear',linear),('material',material)]:
        analytic=force.gradient(x);finite=-fd(force.energy,x)
        gradient_results[name]=dict(max_abs_error=float(np.max(np.abs(analytic-finite))),
            relative_norm_error=float(np.linalg.norm(analytic-finite)/max(1.,np.linalg.norm(analytic))))
    local_q=proxy.native(x);local_complex_J=Cchi@prior.derivative(local_q)
    local_Jx=np.column_stack([proxy.white(local_complex_J[:,j].reshape(32,8),sigma) for j in range(66)])
    local_Jx*=proxy.probit_derivative(x)[None,:]
    oracle=-local_Jx.T@proxy.white(material.prediction(x)-y,sigma)-x
    packing_gradient_error=float(np.linalg.norm(material.gradient(x)-oracle)/max(1.,np.linalg.norm(oracle)))
    prior_only=proxy.MaterialResponseForce(prior,x0,a['chi'],a['prediction'],np.zeros_like(Cchi),y,sigma)
    unchanged_prior_error=float(np.max(np.abs(prior_only.gradient(x)+x)))
    whiten_roundtrip=float(np.max(np.abs(proxy.unwhite(proxy.white(y,sigma),sigma)-y)))
    assert reversal<1e-12 and flux_error<1e-12
    assert abs(logdet)<2e-7 and sign==1
    assert endpoint_error<1e-12 and energy_error<1e-12
    assert calls==dict(full_target=1,surrogate_force=8)
    assert packing_error<1e-10 and packing_gradient_error<1e-10 and unchanged_prior_error==0
    assert all(v['relative_norm_error']<5e-6 for v in gradient_results.values())
    assert whiten_roundtrip<1e-12
    result=dict(status='PASS_PURE_SURROGATE_FORCE_AND_SAVED_ARRAY_CONTROLS',
        wrong_force_is_nonconservative=True,reversibility_max_error=reversal,
        detailed_balance_joint_log_flux_error=float(flux_error),
        finite_difference_phase_space_log_abs_det=logdet,phase_space_determinant_sign=sign,
        exact_HMC_endpoint_oracle_error=endpoint_error,exact_Full_Hamiltonian_ratio_error=energy_error,
        invocation_counts=calls,physical_proxy_gradient_finite_difference=gradient_results,
        runner_AST_no_Full_gradient_calls=True,
        material_response_source_packing_error=packing_error,
        material_gradient_white_J_oracle_relative_error=packing_gradient_error,
        full_Gaussian66_prior_error=unchanged_prior_error,white_roundtrip_error=whiten_roundtrip,
        paid_array_source=anchor_path.relative_to(OUT).as_posix(),
        actual_Maxwell_executed=False,new_physical_RHS=0,Target_constructed=False,Book_constructed=False,
        GT_read=False,training=False,novelty_claim=False,mixing_certificate=False)
    Path(__file__).with_name('RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
