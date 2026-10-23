"""Deterministic synthetic algebra checks only; no data/truth/physics access."""
from pathlib import Path
import json
import time
import numpy as np
from scipy.special import logsumexp
from wr.posterior_math import (log_marginal_likelihood,finite_grid_posterior,
    finite_grid_summary,weighted_cdf,weighted_crps,brier_score)


def direct(f,y,n,shape,sigma=.01,return_weights=False):
    # Independent complex-residual evaluation and independent quadrature nodes.
    # Sum all observations before integrating the one shared nuisance.
    f,y=np.broadcast_arrays(np.asarray(f,dtype=np.clongdouble),np.asarray(y,dtype=np.clongdouble))
    z,w=np.polynomial.legendre.leggauss(n)
    a=np.asarray(1+.01*z,dtype=np.longdouble)
    phi=np.asarray(.01*z,dtype=np.longdouble)
    values=[];aa=[];pp=[]
    count=int(np.prod(shape));constant=-2*count*np.log(np.longdouble(sigma)*np.sqrt(np.longdouble(2)*np.pi))
    for j in range(n):
        for k in range(n):
            gamma=a[j]*(np.cos(phi[k])+1j*np.sin(phi[k]))
            residual=gamma*f-y
            square=np.sum(residual.real**2+residual.imag**2,axis=(-3,-2,-1))
            values.append(constant-square/(2*np.longdouble(sigma)**2)+np.log(np.longdouble(w[j]*w[k]/4)))
            aa.append(a[j]);pp.append(phi[k])
    lognodes=np.stack(values,axis=-1);evidence=logsumexp(lognodes,axis=-1)
    weights=np.exp(lognodes-evidence[...,None])
    if return_weights:return evidence,weights,np.asarray(aa),np.asarray(pp)
    return evidence


def scalar(z):return float(np.asarray(z))
def maxerr(a,b):return float(np.max(np.abs(np.asarray(a)-np.asarray(b))))


def main():
    start=time.process_time();rows={}
    zeros=np.zeros((5,4,12),dtype=complex)
    out=log_marginal_likelihood(zeros,zeros,n=11,return_node_weights=True)
    expected=-480*np.log(.01*np.sqrt(2*np.pi))
    assert abs(float(out['logZ'])-expected)<1e-11
    mm=out['posterior_nuisance_moments']
    assert abs(float(mm['amplitude_mean'])-1)<1e-14
    assert abs(float(mm['phase_mean']))<1e-16
    assert abs(float(mm['amplitude_variance'])-(.02**2/12))<1e-16
    assert abs(float(mm['phase_variance'])-(.02**2/12))<1e-16
    huge_y=log_marginal_likelihood(zeros,np.full_like(zeros,1e6),n=11)
    assert abs(float(huge_y['posterior_nuisance_moments']['amplitude_variance'])-.02**2/12)<1e-16
    rows['absolute_normalization']=dict(large_common_residual_preserves_uniform_nuisance=True,shape=[5,4,12],complex_count=240,real_count=480,
        sigma_per_real=.01,expected_log_constant=expected,computed_logZ=scalar(out['logZ']),
        uniform_amplitude_variance=scalar(mm['amplitude_variance']),
        uniform_phase_variance=scalar(mm['phase_variance']),status='ALGEBRA_CHECK_OK')
    shape=(2,4,3);grid=np.arange(np.prod(shape)).reshape(shape)
    base=.2*np.exp(1j*(.13*grid+.2))
    f=np.stack((base,base+.001*np.cos(grid),base*.997),axis=0)
    y=1.002*np.exp(.003j)*base+.002*np.sin(.3*grid)+.001j*np.cos(.2*grid)
    r11=log_marginal_likelihood(f,y,n=11,observation_shape=shape,return_node_weights=True)
    r21=log_marginal_likelihood(f,y,n=21,observation_shape=shape)
    exact,pw,a,phi=direct(f,y,11,shape,return_weights=True)
    e=maxerr(r11['logZ'],exact)
    assert e<1e-11
    am=np.sum(pw*a,axis=-1);pm=np.sum(pw*phi,axis=-1)
    moment_errors=dict(amplitude_mean=maxerr(am,r11['posterior_nuisance_moments']['amplitude_mean']),
        phase_mean=maxerr(pm,r11['posterior_nuisance_moments']['phase_mean']),
        amplitude_variance=maxerr(np.sum(pw*(a-am[...,None])**2,axis=-1),r11['posterior_nuisance_moments']['amplitude_variance']))
    assert max(moment_errors.values())<1e-12
    rows['independent_complex_residual']=dict(observation_shape=shape,batch_size=3,
        maximum_logZ_error=e,moment_errors=moment_errors,
        accumulation_dtype=r11['accumulation_dtype'],accumulation_epsilon=float(r11['accumulation_epsilon']),
        residual_square_roundoff_scale=float(r11['maximum_residual_square_roundoff_scale']),
        quadrature11_vs21_max_logZ_difference=maxerr(r11['logZ'],r21['logZ']),status='ALGEBRA_CHECK_OK')
    blocked=log_marginal_likelihood(f,y,n=11,observation_shape=shape,source_blocks=True)
    paired=log_marginal_likelihood(f,y,n=11,observation_shape=shape,source_blocks=((0,2),(1,3)))
    err=max(maxerr(blocked['logZ'],r11['logZ']),maxerr(paired['logZ'],r11['logZ']))
    assert err<1e-11
    independent_source_nuisance=sum(direct(f[..., :,j:j+1,:],y[:,j:j+1,:],11,(2,1,3)) for j in range(4))
    difference=maxerr(independent_source_nuisance,r11['logZ'])
    assert difference>1e-6
    try:log_marginal_likelihood(f,y,observation_shape=shape,source_blocks=((0,1),(1,2,3)))
    except ValueError:duplicate_rejected=True
    else:raise AssertionError('Duplicate likelihood source accepted')
    rows['source_partition']=dict(maximum_shared_nuisance_block_error=err,
        wrong_independent_per_source_nuisance_logZ_difference=difference,
        duplicate_source_block_rejected=duplicate_rejected,status='ALGEBRA_CHECK_OK')
    fb=f[:2,None];yb=np.stack([y,y+.001j,y-.002],axis=0)
    batch=log_marginal_likelihood(fb,yb,n=11,observation_shape=shape)
    batchdirect=direct(fb,yb,11,shape)
    assert np.shape(batch['logZ'])==(2,3) and maxerr(batch['logZ'],batchdirect)<1e-11
    rows['broadcasting']=dict(f_batch=[2,1],y_batch=[3],output_batch=list(np.shape(batch['logZ'])),
        maximum_logZ_error=maxerr(batch['logZ'],batchdirect),status='ALGEBRA_CHECK_OK')
    # Stress case is intentionally sharply concentrated. Its11-vs21 discrepancy
    # is disclosed, not relabeled as physical failure or universal convergence.
    sf=10*np.exp(1j*(.03*grid));sy=1.003*np.exp(.004j)*sf
    s11=log_marginal_likelihood(sf,sy,n=11,observation_shape=shape)
    s21=log_marginal_likelihood(sf,sy,n=21,observation_shape=shape)
    d11=direct(sf,sy,11,shape);d21=direct(sf,sy,21,shape)
    assert max(maxerr(s11['logZ'],d11),maxerr(s21['logZ'],d21))<1e-8
    rows['quadrature_stress']=dict(logZ11=scalar(s11['logZ']),logZ21=scalar(s21['logZ']),
        difference11_vs21=maxerr(s11['logZ'],s21['logZ']),
        maximum_independent_direct_error=max(maxerr(s11['logZ'],d11),maxerr(s21['logZ'],d21)),
        quadrature_convergence='NOT_CERTIFIED;11/21 can be inadequate for sharply concentrated nuisance',
        status='DISCREPANCY_RETAINED')
    residual=(1.003*np.exp(.004j)*base-y)
    complex_phi=.5*np.sum(abs(residual)**2)/.01**2
    real_phi=.5*np.sum(np.r_[residual.real.ravel(),residual.imag.ravel()]**2)/.01**2
    assert abs(complex_phi-real_phi)<1e-12
    rows['noise_whitening']=dict(complex_phi=float(complex_phi),realified_phi=float(real_phi),
        divided_by_channel_count=False,status='ALGEBRA_CHECK_OK')
    ll=np.array([0.,np.log(2),-1.]);prior=np.array([2.,.5,1.]);vol=np.array([.1,.6,.3])
    p=finite_grid_posterior(ll,prior,vol);manual=np.exp(ll)*prior*vol;manual/=manual.sum()
    assert maxerr(p['posterior_weights'],manual)<1e-15
    assert abs(float(p['prior_quadrature_mass'])-.8)<1e-15
    assert abs(float(p['log_evidence_normalized_finite_prior']-p['log_evidence_raw_quadrature'])+np.log(.8))<1e-15
    x=np.array([-1.,.25,2.]);weights=p['posterior_weights']
    summary=finite_grid_summary(weights,branch_ids=['left','left','right'],
        event_masks={'x_positive':x>0},target_positions=np.c_[x,x*x])
    assert maxerr(summary['target_position_mean'],np.sum(manual[:,None]*np.c_[x,x*x],axis=0))<1e-15
    cdf=weighted_cdf(x,weights,[-2.,-1.,.3,3.])
    assert maxerr(cdf['values'],[0,manual[0],manual[0]+manual[1],1])<1e-15
    outcome=.1;crps=weighted_crps(x,weights,outcome)
    brute=np.sum(manual*abs(x-outcome))-.5*np.sum(manual[:,None]*manual[None,:]*abs(x[:,None]-x[None,:]))
    assert abs(float(crps)-brute)<1e-15
    brier=brier_score(summary['event_weights']['x_positive'],1)
    assert abs(float(brier)-(1-manual[1]-manual[2])**2)<1e-15
    rows['finite_quadrature_posterior']=dict(prior_density=prior.tolist(),cell_weights=vol.tolist(),
        posterior_weights=np.asarray(weights,dtype=float).tolist(),prior_quadrature_mass=scalar(p['prior_quadrature_mass']),
        raw_log_evidence=scalar(p['log_evidence_raw_quadrature']),
        normalized_finite_prior_log_evidence=scalar(p['log_evidence_normalized_finite_prior']),
        branch_weights=[dict(branch=z['branch'],weight=scalar(z['weight'])) for z in summary['branch_weights']],
        target_mean=np.asarray(summary['target_position_mean'],dtype=float).tolist(),
        target_variance=np.asarray(summary['target_position_variance'],dtype=float).tolist(),
        finite_atom_cdf=np.asarray(cdf['values'],dtype=float).tolist(),crps=scalar(crps),brier=scalar(brier),
        interpretation=cdf['interpretation'],status='ALGEBRA_CHECK_OK')
    # Extreme log likelihood verifies that normalization does not underflow.
    stable=finite_grid_posterior([-1e6,-1e6-1000],[1.,1.],[.25,.75])
    assert abs(float(np.sum(stable['posterior_weights']))-1)<1e-12
    rows['logsumexp_stability']=dict(log_evidence=scalar(stable['log_evidence_raw_quadrature']),
        posterior_weight_sum=scalar(np.sum(stable['posterior_weights'])),status='ALGEBRA_CHECK_OK')
    result=dict(status='SYNTHETIC_PROBABILITY_ALGEBRA_CHECKS_COMPLETE',P1_gate='NOT_EVALUATED',
        Maxwell_status='NOT_RUN',physical_solves=0,RHS=0,LU=0,data_generation=False,
        truth_access=False,sampling=False,synthetic_algebra_arrays_created=True,checks=rows,cpu_seconds=time.process_time()-start,
        scope='deterministic small synthetic arrays only; not Maxwell result, posterior result or qualification')
    path=Path(__file__).with_name('CHECK_RESULTS.json');path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','P1_gate','Maxwell_status','cpu_seconds')}))

if __name__=='__main__':main()
