"""Material-event handoff, context update and explicitly bounded prototypes.

An event mass belongs to a measurable region, never its displayed medoid.
This module doesn't import or construct a solver. A metered forward is passed
only for the selected representative and best-found MAP endpoints.
"""
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import truncnorm,chi2
from .material import BOUNDS,WIDTH,geometry,void_volume,event_ids,prior_event_masses,material_distance,decode,from_geometry,RMAX
from .likelihood import log_likelihood,conditional_gain,context_loglik
from .accounting import dump

def update_context(latent,weights,log_context_likelihood):
    """Full-support weighted-sample update; callers declare conditional independence.

    The callback consumes materials only. It must not reuse scattering data.
    Finite sample support and mixing/importance errors remain caller obligations.
    """
    x=np.asarray(latent);w=np.asarray(weights,float)
    if w.shape!=(len(x),) or np.any(w<0) or abs(w.sum()-1)>1e-8:raise ValueError('normalized whole modeled posterior required')
    ll=np.asarray(log_context_likelihood(x),float)
    if ll.shape!=w.shape or not np.all(np.isfinite(ll)):raise ValueError('finite independent context log likelihood required')
    lw=np.where(w>0,np.log(np.maximum(w,np.finfo(float).tiny))+ll,-np.inf);logz=logsumexp(lw);a=np.exp(lw-logz)
    return dict(weights=a,event_mass=np.bincount(event_ids(x),weights=a,minlength=12),context_log_evidence=float(logz),ESS=float(1/np.sum(a*a)),outside_model_support_mass=None)

def medoid(x,rng):
    """Finite material-distance medoid on a 257-cell radial volume grid."""
    x=np.asarray(x);c=x[rng.choice(len(x),min(96,len(x)),replace=False)];b=x[rng.choice(len(x),min(512,len(x)),replace=False)]
    edges=np.linspace(0,1,258);mid=(edges[:-1]+edges[1:])/2;vol=np.diff(edges**3)
    def field(v):
        g=geometry(v);void=(mid[None,:]<g[:,0,None])|((g[:,1,None]>0)&(abs(mid[None,:]-.65)<g[:,1,None]))
        eps=v[:,2,None]+1j*v[:,3,None];return np.where(void,1.,eps)
    a=field(c);bb=field(b);dist=np.sqrt(np.maximum(((abs(a)**2*vol).sum(1)[:,None]+(abs(bb)**2*vol).sum(1)[None,:]-2*np.real((a*vol)@bb.conj().T)),0))/abs(1.5+.1j)
    return c[np.argmin(dist.mean(1))]

def endpoint(forward,x,y,sigma):
    F=np.concatenate([forward.predict(x,k,kind='SELECTED_MATERIAL_HYPOTHESIS_ENDPOINT')[0] for k in [1.6,2.4]])
    m,v=conditional_gain(F[:32],y[:32],sigma[:32]);a,b=(.7-m)/np.sqrt(v),(1.3-m)/np.sqrt(v)
    gain,variance=truncnorm.stats(a,b,loc=m,scale=np.sqrt(v),moments='mv')
    return dict(unscaled_field=F.tolist(),packing='each frequency real32: Re16 then Im16',frequencies=[1.6,2.4],posterior_gain_mean=float(gain),posterior_gain_SD=float(np.sqrt(variance)),gain_support=[.7,1.3],gain_used_for_both_frequencies=True,mean_prediction=(gain*F).tolist(),whitened_residual_C=float(np.linalg.norm((gain*F[:32]-y[:32])/sigma[:32])),marginal_log_likelihood_C=float(log_likelihood(F[:32],y[:32],sigma[:32])),likelihood_status='full vector-Mie n14; proper bounded-gain marginal',new_or_paid_cache_cost='sole PHYSICS_LEDGER receipt')

def map_local(response,y,sigma,rng):
    starts=[from_geometry(0,0),from_geometry(.08,0),from_geometry(.25,0),from_geometry(0,.08),from_geometry(.2,.06),*rng.uniform(BOUNDS[:,0],BOUNDS[:,1],size=(3,4))];runs=[]
    for x in starts:
        r=minimize(lambda v:-float(log_likelihood(response.predict(v)[0],y[:32],sigma[:32])),x,method='L-BFGS-B',bounds=BOUNDS,options={'maxiter':100,'ftol':1e-10})
        runs.append(dict(latent=r.x,log_likelihood=-r.fun,success=bool(r.success),proxy_evaluations=r.nfev))
    best=max(runs,key=lambda a:a['log_likelihood'])
    return best,runs

def predictive_check(response,x,y,sigma,rng,weights=None):
    # Physical approximation is explicitly the checked proxy. No extra solve.
    xx=x[rng.choice(len(x),min(len(x),1024),replace=weights is not None,p=weights)];F=response.predict(xx);m,v=conditional_gain(F,y[:32],sigma[:32]);g=truncnorm.ppf(rng.random(len(xx)),(.7-m)/np.sqrt(v),(1.3-m)/np.sqrt(v),loc=m,scale=np.sqrt(v))
    t=((g[:,None]*F-y[:32])/sigma[:32])**2
    return dict(posterior_predictive_tail_probability=float(chi2.sf(t.sum(1),32).mean()),statistic='sum of fixed-noise whitened squared residuals conditional on sampled material/gain',samples=len(xx),numerical_target='proxy; empirically full-field checked separately',not_model_adequacy_certificate=True)

def build_packages(results,response,reference,forward):
    import json
    results=Path(results);(results/'packages').mkdir(exist_ok=True);items=json.loads((results/'DEVELOPMENT_CASES.json').read_text())['cases'];pilot=json.loads((results/'PILOT_REFERENCE.json').read_text())['cases'];corr=json.loads((results/'FULL_MIE_CORRECTED_REFERENCE.json').read_text())['cases'];rng=np.random.default_rng(261013215);packages={}
    for item in items:
        name=item['id'];y=np.array(item['y']);sigma=np.array(item['sigma']);raw=np.load(results/f'{name}_REFERENCE_PARTICLES.npz');x=raw['latent'];eid=event_ids(x)
        proxy=pilot[name]['C']['orders'][-1];p=np.array(corr[name]['event_mass'] if name in corr else proxy['event_mass']);prior=prior_event_masses();order=np.argsort(-p,kind='stable');k=1+np.searchsorted(np.cumsum(p[order]),.95);selected=order[:min(12,k)];hyp=[]
        for event in selected:
            event=int(event);xx=x[eid==event];rep=medoid(xx,rng);c,s=divmod(event,3)
            pred=endpoint(forward,rep,y,sigma);conditional=(corr[name]['conditional'] if name in corr else proxy['conditional'])[event]
            mass_error=dict(MC_residual_SE=float(corr[name]['event_mass_MC_se'][event]) if name in corr else None,MC_interval95=corr[name]['event_mass_interval95'][event] if name in corr else None,quadrature_refinement_max=pilot[name]['C']['refinement']['event_mass'],interval_scope='residual Monte Carlo only, cubature and untested full-field discrepancy separate')
            hyp.append(dict(hypothesis_id=f'{name}-E{event}',event_id=event,event_definition=dict(central_radius_bin=c,annular_halfwidth_bin=s,central_edges=[0,.1,.2,float(RMAX[0])],annular_edges=[0,.06,float(RMAX[1])],zero_is_atom=True,disjoint=True),representative_latent=rep.tolist(),representative_material={kk:vv.tolist() if hasattr(vv,'tolist') else vv for kk,vv in decode(rep).items() if kk!='eps_layers'},representative_epsilon=dict(real=float(rep[2]),imag=float(rep[3])),representative_role='finite event-medoid in volume-weighted complex-material distance; not a point-probability',posterior_mass=float(p[event]),mass_error=mass_error,prior_mass=float(prior[event]),posterior_prior_odds_ratio=float((p[event]/(1-p[event]))/(prior[event]/(1-prior[event]))),log_event_evidence=float((corr[name]['log_evidence'] if name in corr else proxy['log_evidence'])+np.log(p[event])-np.log(prior[event])),material_uncertainty=dict(coordinates=['central_radius','annular_halfwidth','epsilon_real','epsilon_imag'],mean=conditional['mean'],SD=conditional['SD'],covariance=conditional['covariance'],continuous_quantile_reference=pilot[name]['C']['continuous'],quantile_target='proxy continuous density; full-Mie correction verified for moments/event masses, not full-CDF uniform certification'),forward_prediction=pred,compensation_diagnostics=dict(raw_material_parameters=rep.tolist(),gain_is_instrument_nuisance=True,minimum_residual_is_not_evidence=True),indistinguishable_from=[]))
        for i,a in enumerate(hyp):
            for b in hyp[i+1:]:
                fa=np.array(a['forward_prediction']['unscaled_field']);fb=np.array(b['forward_prediction']['unscaled_field']);ga=a['forward_prediction']['posterior_gain_mean'];gb=b['forward_prediction']['posterior_gain_mean'];d=float(np.linalg.norm((ga*fa[:32]-gb*fb[:32])/sigma[:32]));raw_d=float(np.linalg.norm((fa[:32]-fb[:32])/sigma[:32]));md=material_distance(a['representative_latent'],b['representative_latent'])
                row=dict(raw_noise_distance=raw_d,conditional_gain_adjusted_noise_distance=d,material_distance=md,frequencies=[1.6],noise_scale='fixed sigma, joint real-channel whitened norm',status='representative posterior worlds; not global compensation minimum or exact equivalence')
                if d<=2:
                    a['indistinguishable_from'].append(dict(hypothesis_id=b['hypothesis_id'],**row));b['indistinguishable_from'].append(dict(hypothesis_id=a['hypothesis_id'],**row))
        best,runs=map_local(response,y,sigma,rng);best['full_endpoint']=endpoint(forward,best['latent'],y,sigma);best['status']='best of eight bounded proxy multistart fits; selected full endpoint; not global MAP certification';best['event_id']=int(event_ids(best['latent']))
        # Whole posterior sample handoff has exact event masses plus finite
        # within-event proxy quadrature resampling. Full-Mie arrays separate.
        # Restore a finite conditional approximation for tiny events that the
        # 8192 quadrature resamples missed, without assigning zero mass to them.
        adaptive=json.loads((results/'ADAPTIVE_PHYSICAL_REFERENCE.json').read_text())['cases'][name]
        from .correction import EventProposal
        extra={}
        for event in range(12):
            if not np.any(eid==event):
                q=EventProposal(adaptive['leaves'],event);xx=q.sample(rng,32);lw=log_likelihood(response.predict(xx),y[:32],sigma[:32])-q.logpdf(xx)-np.log(WIDTH).sum();extra[event]=np.exp(lw-logsumexp(lw));x=np.concatenate((x,xx));eid=np.concatenate((eid,np.full(32,event)))
        w=np.zeros(len(x))
        for event in range(12):
            mask=eid==event
            if mask.any():w[mask]=p[event]*(extra[event] if event in extra else np.full(mask.sum(),1/mask.sum()))
        missing=1-w.sum();w/=w.sum()
        np.savez_compressed(results/'packages'/f'{name}_samples.npz',latent=x,weights=w,event_ids=eid,event_mass=p,unrepresented_sample_event_mass=missing)
        package=dict(observation_id=name,case_id=name,target_identity='joint-p2-radial-v1:C:full-vector-Mie-n14' if name in corr else 'joint-p2-radial-v1:C:tensor-Chebyshev-five-per-physical-axis-radial-v1',measurement=y[:32].tolist(),noise_sigma=sigma[:32].tolist(),hypotheses=hyp,represented_mass=float(p[selected].sum()),unshown_modeled_mass=float(1-p[selected].sum()),all_event_mass=p.tolist(),all_prior_mass=prior.tolist(),unexplored_mass_status='declared radial support integrated; no guarantee for omitted modes outside tested support',outside_model_support_mass=None,outside_model_support_mass_status='unquantified',numerical_status=corr[name]['numerical_status'] if name in corr else 'EMPIRICAL_PROXY_REFERENCE; 13/17 and independent adaptive integration',model_adequacy_status='radial concentric 3D vector model only; finite noise/frequency/calibration assumptions, no real-world certification',posterior_samples=dict(path=f'{name}_samples.npz',measure='latent dx / whole twelve-event support',weights_status='corrected event masses with proxy event-conditional quadrature samples; not exact full-Mie posterior draws',sample_mass_numerical_gap=float(missing),full_corrected_samples=f'../{name}_FULL_CORRECTION.npz' if name in corr else None),single_MAP=best,multiple_start_fits=runs,simple_residual_ranking='marginal likelihood best-fit endpoints saved; residual ranking is not evidence and no mass assigned',posterior_predictive_consistency=predictive_check(response,x,y,sigma,rng,weights=w),downstream_context_update=dict(interface='b1_joint.package.update_context(latent,weights,log_context_likelihood)',condition='independent geometry sensor given material; do not reuse y',scope='finite weighted approximation; omitted/model-outside mass unquantified'),truth_latent=item['latent'],truth_usage='simulation/evaluation only, never inference feature')
        dump(results/'packages'/f'{name}.json',package);packages[name]=package
    return packages
