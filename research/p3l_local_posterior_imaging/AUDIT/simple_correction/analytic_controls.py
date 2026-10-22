"""Synthetic simple-density/MH contracts; no paid models, context, or physics."""
from pathlib import Path
import ast
import json
import sys
import tempfile
import numpy as np

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tools'))
import run_p3l_simple_correction as runner
import run_p3l_learned_correction as common
from p3l.learned_proposal import LearnedFullProposal


def main():
    rng=np.random.default_rng(2026102299)
    scalar=rng.normal(size=(1,1154));image=np.zeros((1,4,24,24))
    arrays=dict(scalar_mean=np.zeros(1154),scalar_scale=np.ones(1154),
        label_mean=np.zeros(10),label_scale=np.linspace(.8,1.2,10),
        projection=rng.normal(size=(1154,3))*.01,weights=rng.normal(size=(3,10))*.1,
        covariance=np.eye(10)*.7)
    with tempfile.TemporaryDirectory() as temporary:
        path=Path(temporary)/'synthetic_fit.npz';np.savez(path,**arrays)
        m1=runner.load_simple(path,3,rng=np.random.default_rng(42))
        m2=runner.load_simple(path,3,rng=np.random.default_rng(42))
        sample1=m1.sample(5,scalar,image);sample2=m2.sample(5,scalar,image)
        reproducibility_error=float(np.max(np.abs(sample1-sample2)))
        mean=(scalar@arrays['projection'])@arrays['weights']
        standard=(sample1-arrays['label_mean'])/arrays['label_scale']-mean[:,None,:]
        oracle=-.5*(np.sum(standard**2/.7,axis=-1)+10*np.log(2*np.pi*.7))-np.log(arrays['label_scale']).sum()
        normalized_density_error=float(np.max(np.abs(m1.log_prob(sample1,scalar,image)-oracle)))
        frame=np.eye(64);q0=np.r_[np.zeros(64),[.5,.5]]
        J=rng.normal(size=(512,66))*.005;residual=rng.normal(size=512)*.01
        proposal=LearnedFullProposal(m1,scalar,image,frame[:,:8],frame[:,8:],q0,J,residual,defensive_weight=.1)
        flux=[]
        for _ in range(8):
            x,y=proposal.sample(2,rng=rng)
            px,py=-.5*x@x,-.5*y@y;qx,qy=proposal.log_prob(x),proposal.log_prob(y)
            ratio=common.independent_log_ratio(py,px,qy,qx)
            reverse=common.independent_log_ratio(px,py,qx,qy)
            flux.append(abs((px+qy+min(0.,ratio))-(py+qx+min(0.,reverse))))
    source=ROOT/'tools/run_p3l_simple_correction.py';tree=ast.parse(source.read_text())
    calls=[n for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='reconstruct']
    assert len(calls)==1
    options={k.arg:ast.literal_eval(k.value) for k in calls[0].keywords}
    assert options==dict(method='coarse_TRF',max_evaluations=16,initial_point=None,polish=False)
    assert common.neural_selection_status(dict(best_update=0))=='NOT_RUN_VALIDATION_SELECTED_INITIALIZATION_NO_TRAINED_GAIN'
    assert common.neural_selection_status(dict(best_update=100))=='AVAILABLE_TRAINED_VALIDATION_SELECTION'
    assert runner.run_pair is common.run_pair
    result=dict(status='PASS_SYNTHETIC_FITTED_DENSITY_AND_FULL_MH_CONTRACT_CONTROLS',
        deterministic_numpy_model_RNG_error=reproducibility_error,
        normalized_affine_Gaussian_log_density_error=normalized_density_error,
        normalized_defensive_Full_MH_log_flux_error=float(max(flux)),
        exact_shared_run_pair=True,context_options=options,
        update0_neural_status='NOT_RUN_VALIDATION_SELECTED_INITIALIZATION_NO_TRAINED_GAIN',
        fixed_families=list(runner.FAMILIES),job=runner.JOB,total_RHS_upper=19008,
        real_fitted_checkpoints_loaded=False,Target_constructed=False,Book_constructed=False,
        physics_executed=False,training_or_fitting=False,GT_read=False)
    assert reproducibility_error==0 and normalized_density_error<1e-12 and max(flux)<1e-10
    Path(__file__).with_name('RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
