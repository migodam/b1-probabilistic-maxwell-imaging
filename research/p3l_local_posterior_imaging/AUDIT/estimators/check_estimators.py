"""Standalone supplied-sample controls; no physics/GT/backend imports."""
import json
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT/'src'))
from p3l.estimators import (posterior_mean_variance, complex_geometric_median,
                          component_quantiles, material_diagnostics,
                          project_legal_bayes_action)
from p1_material.prior import MaterialPrior
OUT = Path(__file__).resolve().parent
results = {}
def check(name, condition, **details):
    results[name] = dict(pass_=bool(condition), **details)
    if not condition:
        raise AssertionError(name)
def constant_pixels(x):
    return np.repeat(np.asarray(x, complex)[:,None],576,axis=1)

# An antithetic finite ensemble has exactly the Gaussian mean and symmetry.
rng=np.random.default_rng(20261022)
center=.4+.7j
noise=rng.normal(size=(32,576))+1j*rng.normal(size=(32,576))
samples=np.r_[center+noise,center-noise]
a=posterior_mean_variance(samples)
m=complex_geometric_median(samples)
check('symmetric_gaussian_ensemble',np.max(abs(a['mean']-center))<2e-15 and np.max(abs(m.value-center))<2e-9 and np.all(m.converged),mean_error=float(np.max(abs(a['mean']-center))),median_error=float(np.max(abs(m.value-center))))
# Complex variance must use modulus, including imaginary directions.
a=posterior_mean_variance(constant_pixels([-1j,1j]))
check('modulus_variance',np.allclose(a['variance'],1) and np.allclose(a['mean'],0))
# Coordinatewise medians are zero, but the Euclidean complex median is not.
t=(3-np.sqrt(3))/6
m=complex_geometric_median(constant_pixels([0,1,1j]))
check('triangle_median_coincident_nonoptimum',np.max(abs(m.value-(t+1j*t)))<2e-8 and np.all(m.converged),maximum_error=float(np.max(abs(m.value-(t+1j*t)))),subgradient_residual=float(np.max(m.normalized_subgradient_residual)))
for name,points,expected in [('repeated_majority',[0,0,0,1+2j],0),('all_identical',[2+3j]*5,2+3j),('two_point_nonunique',[0,2],1)]:
    m=complex_geometric_median(constant_pixels(points))
    check(name,np.allclose(m.value,expected) and np.all(m.converged))
q=component_quantiles(constant_pixels([0,2+4j]),[0,.5,1])
check('component_quantile_semantics',np.allclose(q['real'][:,0],[0,1,2]) and np.allclose(q['imag'][:,0],[0,2,4]))
d=material_diagnostics(constant_pixels([-1j,1j]))
check('spread_without_truth',d['mean_material_variance']==1 and d['material_spread_rms']==1)

class LinearDecoder:
    dimension=2
    gaussian_dimension=2
    def decode(self,q):return np.full(576,q[0]+1j*q[1])
    def derivative(self,q):return np.tile([1,1j],(576,1))
a=project_legal_bayes_action(np.full(576,center),LinearDecoder(),[0,0],posterior_latent_samples=[[1,1],[-1,1],[1,1]],maximum_starts=4)
check('linear_decoder_mean_is_legal',a['projection_mse']<1e-20,projection_mse=a['projection_mse'],starts=len(a['trials']))

class ParabolaDecoder:
    dimension=1
    gaussian_dimension=1
    def decode(self,q):return np.tile([q[0],q[0]**2],288).astype(complex)
    def derivative(self,q):return np.tile([1,2*q[0]],288)[:,None].astype(complex)
p=ParabolaDecoder()
mean=(p.decode([-1])+p.decode([1]))/2
projected=project_legal_bayes_action(mean,p,[0],posterior_latent_samples=[[-1],[1]],maximum_starts=3)
check('nonlinear_mean_outside_decoder_image',abs(projected['projection_mse']-.375)<1e-8,projection_mse=projected['projection_mse'],q=projected['q'].tolist(),status=projected['status'])
check('multiple_starts_escape_stationary_ridge',projected['trials'][0]['objective']>.49 and projected['projection_mse']<.38)

# Supplied P1 decoder: different global property phases yield a mean with
# spatially varying phase, which cannot obey its shared increment law.
cfg=json.loads((ROOT/'research/p1_material_target_20261018/configs/ONLINE_TARGET.json').read_text())
points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
prior=MaterialPrior(cfg,points)
q1=np.zeros(66);q2=np.zeros(66)
q1[0]=-.3;q1[3]=.5;q1[64:]=[.1,.9]
q2[0]=.5;q2[1]=-.2;q2[64:]=[.9,.1]
mean=(prior.decode(q1)[0]+prior.decode(q2)[0])/2
inc=mean-prior.host
mask=prior.mask & (inc.real>1e-8)
ratio=inc.imag[mask]/inc.real[mask]
phase_range=float(np.ptp(ratio))
projection=project_legal_bayes_action(mean,prior,q1,posterior_latent_samples=np.array([q1,q2]),maximum_starts=2,max_nfev=40)
check('p1_mean_not_single_global_increment_law',phase_range>1e-3,phase_ratio_range=phase_range,projection_mse=projection['projection_mse'],optimization_success=projection['trials'][projection['selected_start']]['success'],status=projection['status'])
check('p1_projection_bound_and_no_prior_penalty',prior.legal(projection['q']) and abs(np.mean(abs(projection['material']-mean)**2)-projection['projection_mse'])<1e-14)
# Intentionally capped projection must retain an explicit unsuccessful status.
a=project_legal_bayes_action(np.full(576,center),LinearDecoder(),[4,4],max_nfev=1)
check('cap_not_claimed_optimum',not a['trials'][0]['success'] and a['status']=='LOCAL_NUMERICAL_PROJECTION_NOT_GLOBAL_ARGMIN')
for name,call in [('nonfinite',lambda:posterior_mean_variance(np.full((1,576),np.nan))),('wrong_shape',lambda:component_quantiles(np.zeros((2,2)))),('illegal_latent',lambda:project_legal_bayes_action(mean,prior,np.r_[np.zeros(64),2,0]))]:
    try:call()
    except ValueError:check(name,True)
    else:check(name,False)
report=dict(status='PASS_SUPPLIED_SAMPLE_ZERO_PHYSICS_CONTROLS_ONLY',seed=20261022,new_maxwell=0,new_rhs=0,truth_access=False,checks=results)
(OUT/'results.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
