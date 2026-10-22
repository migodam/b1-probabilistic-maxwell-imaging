"""Synthetic runner contracts only; never construct Book, Target, or a trainer."""
from pathlib import Path
from types import SimpleNamespace
import ast
import importlib.util
import json
import sys
import tempfile
import numpy as np

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'src'))
from p3l.learned_proposal import GaussianWeakDensity,LearnedFullProposal
from p3l.core import pack,from_probit


def main():
    source=ROOT/'tools/run_p3l_learned_correction.py'
    ast.parse(source.read_text())
    spec=importlib.util.spec_from_file_location('correction_runner',source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    rng=np.random.default_rng(2026102283)
    q0=np.r_[np.zeros(64),[.4,.6]]
    J=rng.normal(size=(512,66))*.01;residual=rng.normal(size=512)*.02
    model=GaussianWeakDensity(np.zeros(10),np.eye(10)*.7,rng=rng)
    proposal=LearnedFullProposal(model,np.zeros(1154),np.zeros((4,24,24)),
        np.eye(64)[:,:8],np.eye(64)[:,8:],q0,J,residual,defensive_weight=.1)
    R=rng.normal(size=(17,66))*.01
    H=np.eye(66)+R.T@R
    logpi=lambda x:-.5*float(x@H@x)
    errors=[];wrong_errors=[]
    for _ in range(12):
        x,y=proposal.sample(2,rng=rng)
        qx,qy=float(proposal.log_prob(x)),float(proposal.log_prob(y))
        ratio=module.independent_log_ratio(logpi(y),logpi(x),qy,qx)
        reverse=module.independent_log_ratio(logpi(x),logpi(y),qx,qy)
        flux_forward=logpi(x)+qy+min(0.,ratio)
        flux_reverse=logpi(y)+qx+min(0.,reverse)
        errors.append(abs(flux_forward-flux_reverse))
        wrong_errors.append(abs((logpi(x)+qy+min(0.,logpi(y)-logpi(x)))-
            (logpi(y)+qx+min(0.,logpi(x)-logpi(y)))))
    # The frozen context feature API exactly reproduces the context runner.
    prior=SimpleNamespace(A=rng.normal(size=(576,64))*.01)
    derivative=rng.normal(size=(576,66))*.01+1j*rng.normal(size=(576,66))*.01
    prior.derivative=lambda q:derivative
    y=rng.normal(size=(32,8))+1j*rng.normal(size=(32,8))
    state=SimpleNamespace(data=y+.01,chi=np.ones(576)*(.2+.1j))
    point=SimpleNamespace(q=q0,state=state)
    target=SimpleNamespace(prior=prior,jacobian=lambda p:J,white=lambda a:pack(a)*np.sqrt(2)/.1)
    scalar,image,gotJ,gotres=module.context_arrays(target,point,y)
    singular=np.linalg.svd(J[:,:64],compute_uv=False)
    expected=np.r_[target.white(y),target.white(state.data-y),q0,np.log1p(singular*singular)]
    feature_error=float(np.max(np.abs(scalar-expected)))
    expected_image=np.stack((state.chi.real,state.chi.imag,prior.A@q0[:64],
        np.sqrt(np.sum(np.abs(derivative[:,:64])**2,axis=1)))).reshape(4,24,24)
    image_error=float(np.max(np.abs(image-expected_image)))
    with tempfile.TemporaryDirectory() as temporary:
        original=module.OUT;module.OUT=Path(temporary)
        boundary=np.r_[np.zeros(64),[0.,1.]]
        x,metadata=module.finite_start('DEV_000',boundary)
        fallback_error=float(np.max(np.abs(from_probit(x)[64:]-[1e-4,1-1e-4])))
        modefile=module.OUT/'A/jobs/reference_DEV_000_mala_mode_c0/PROBIT_MODE_INITIALIZATION.json'
        modefile.parent.mkdir(parents=True)
        mode=np.arange(66,dtype=float)*.001
        modefile.write_text(json.dumps(dict(mode_endpoint='mode',mode_history=[dict(endpoint='mode',x=mode.tolist())])))
        got,metadata=module.finite_start('DEV_000',boundary)
        mode_error=float(np.max(np.abs(got-mode)))
        module.OUT=original
    result=dict(status='PASS_PURE_NORMALIZED_MH_AND_CONTEXT_CONTRACT_CONTROLS',
        detailed_balance_log_flux_max_error=float(max(errors)),
        deliberately_omitted_proposal_ratio_flux_max_error=float(max(wrong_errors)),
        scalar_feature_error=feature_error,image_feature_error=image_error,
        inward_anchor_fallback_error=fallback_error,saved_mode_preference_error=mode_error,
        fixed_parents=list(module.PARENTS),fixed_families=list(module.FAMILIES),
        fixed_training_seed=module.SEED,transitions=module.TRANSITIONS,
        Maxwell_executed=False,Book_constructed=False,Target_constructed=False,
        training=False,real_checkpoints_loaded=False,GT_read=False)
    assert max(errors)<1e-10 and max(wrong_errors)>1e-3
    assert feature_error==0 and image_error==0 and fallback_error<1e-14 and mode_error==0
    Path(__file__).with_name('RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
