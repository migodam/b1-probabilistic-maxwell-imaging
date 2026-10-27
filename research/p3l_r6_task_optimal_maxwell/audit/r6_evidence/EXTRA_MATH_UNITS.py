"""Zero-Maxwell tests; imports only numerical libraries and copied static prior."""
from pathlib import Path
import json,importlib.util,time
import numpy as np
from scipy.special import ndtr
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parent

def main():
 start=time.perf_counter();src=ROOT/'decoder_replay/decoder_sources'
 spec=importlib.util.spec_from_file_location('static_decoder',src/'prior.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 cfg=json.loads((src/'ONLINE_TARGET.json').read_text());geom=json.loads((src/'geometry.json').read_text());prior=mod.MaterialPrior(cfg,geom['points'])
 # Legal 3D frozen submodel: field offset plus both property coordinates.
 qfrac=.95;delta=complex(np.mean(prior.real_interval),np.mean(prior.imag_interval));scaled=delta/qfrac
 def props(d):return [(d.real-prior.real_interval[0])/np.diff(prior.real_interval)[0],(d.imag-prior.imag_interval[0])/np.diff(prior.imag_interval)[0]]
 x0=brentq(lambda x:float(mod.saturation(1.3*x))-qfrac,-10/13,10/13)
 pure=np.r_[2.,np.zeros(63),props(delta)];trans=np.r_[x0,np.zeros(63),props(scaled)]
 cp,ip=prior.decode(pure);ct,it=prior.decode(trans);gerr=float(np.max(np.abs(cp-ct)));assert prior.legal(pure) and prior.legal(trans) and gerr<1e-13
 gauge=dict(status='PASS_DECODER_ONLY',pure_offset=2.,transition_offset=x0,fraction=qfrac,pure_property=[delta.real,delta.imag],transition_property=[scaled.real,scaled.imag],maximum_material_difference=gerr,scope='new 3D fixed-offset+two-property submodel; not Full66 positive-mass overlap')
 # Exact omitted-mass identity and a conservative task-error bound with signed numerators.
 rng=np.random.default_rng(202610270611);bound_cases=[]
 for k in range(100):
  prior_mass=rng.dirichlet(np.ones(9));likelihood=np.exp(-rng.uniform(0,8,9));phi=rng.normal(size=9);covered=np.arange(9)<5
  ev=prior_mass*likelihood;Z=ev.sum();mu=float(ev@phi/Z);Zc=ev[covered].sum();mc=float(ev[covered]@phi[covered]/Zc);eps=float(ev[~covered].sum()/Z);mm=float(ev[~covered]@phi[~covered]/ev[~covered].sum());identity=abs((mu-mc)-eps*(mm-mc))
  action=float(rng.normal());N=ev[covered]@phi[covered];M=ev[~covered]*1.01;b=np.abs(phi[~covered]-action);lower=.98*Zc
  bound=(abs(N-action*Zc)+b@M)/lower;error=abs(mu-action);assert identity<2e-15 and error<=bound+1e-14
  bound_cases.append(dict(identity_error=identity,error=error,bound=float(bound)))
 missing=dict(status='PASS_FINITE_ALGEBRA',cases=100,max_identity_error=max(x['identity_error'] for x in bound_cases),minimum_bound_slack=min(x['bound']-x['error'] for x in bound_cases),qualification='known finite source partition; no unvisited Maxwell domain certified')
 # Compensation: integrate z in [-2,2] analytically for f=z+h^3.
 n,w=leggauss(256);h=(n+1)/2;wh=w/2;comp=[]
 for sigma in [.1,.05,.02]:
  integral=sigma*np.sqrt(2*np.pi)*(ndtr((2+h**3)/sigma)-ndtr((-2+h**3)/sigma));weights=wh*integral;weights/=weights.sum();mean=float(weights@h);var=float(weights@(h-mean)**2);mass=float(weights[h>.5].sum());assert abs(mean-.5)<1e-12 and abs(var-1/12)<1e-12 and mass>.49
  comp.append(dict(sigma=sigma,mean_h=mean,var_h=var,posterior_mass_h_gt_half=mass))
 # Split at x=.5 for event M>.125. Atom x<=0 has exact prior mass .5.
 xx=[];ww=[]
 for lo,hi in [(0.,.5),(.5,1.)]:xx.extend((lo+(hi-lo)*(n+1)/2).tolist());ww.extend((w*(hi-lo)/4).tolist())
 material=np.r_[0.,np.array(xx)**3];base=np.r_[.5,np.array(ww)];tasks=np.vstack([material,material**2,(material>.125).astype(float)])
 sigma=.2;y=.13
 def moments(yy):
  weights=base*np.exp(-(material-yy)**2/(2*sigma*sigma));weights/=weights.sum();return weights,tasks@weights
 weights,means=moments(y);cov=(tasks*material)@weights-means*(material@weights);deriv=cov/sigma**2;rows=[]
 for step in [1e-3,1e-4,1e-5]:
  fd=(moments(y+step)[1]-moments(y-step)[1])/(2*step);err=float(np.max(np.abs(fd-deriv)));rows.append(dict(step=step,max_abs_error=err,finite_difference=fd.tolist()))
 assert rows[-1]['max_abs_error']<1e-7 and rows[-1]['max_abs_error']<rows[0]['max_abs_error']
 # A deliberately y-independent q has zero actual derivative but positive covariance.
 uq=np.ones(3)/3;mq=np.array([0.,.5,1.]);wrong_cov=float((uq@mq**2-(uq@mq)**2)/sigma**2);assert wrong_cov>1
 derivative=dict(status='PASS_FINITE_QUADRATURE',sigma=sigma,y=y,tasks=['material_mean','material_second_moment','material_event_gt_0.125'],quadrature_nodes=len(material),analytic_covariance_derivative=deriv.tolist(),finite_difference_checks=rows,wrong_q_counterexample=dict(actual_derivative=0.,cov_q_over_sigma2=wrong_cov),qualification='toy likelihood-tilt posterior incl atom; finite quadrature consistency only, no directed rounding')
 result=dict(status='PASS_ZERO_PHYSICS_EXTRA_MATH',same_material_gauge=gauge,missing_mass=missing,compensated_scalar=dict(status='PASS_NO_NEW_NORMAL_SHRINKAGE',model='z Uniform[-2,2], h Uniform[0,1], f=z+h^3, observation0',rows=comp,reason='z=-h^3 remains interior for every h; c lies in range A'),derivative_identity=derivative,cost=dict(new_Maxwell_states=0,new_Maxwell_RHS=0,training_updates=0,wall_seconds=time.perf_counter()-start))
 (ROOT/'EXTRA_MATH_RESULTS.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
