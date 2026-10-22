"""Full-target, audited MCMC; transform and preconditioner never change prior."""
from collections import OrderedDict
import json,time
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from .core import OUT,save,write,from_probit,to_probit,probit_derivative
from .kernels import ChainState,PreconditionedMALA,PCN,EllipticalSlice,PropertyReflection

class FrozenMaterialMobility:
    """Changing decoder geometry, one explicitly paid frozen Maxwell factor.

    The mobility is only a deterministic proposal rule. It never replaces the
    full likelihood or source-adjoint gradient, and its two endpoint Gaussian
    determinants must enter MH. Identity precision belongs to probit Gaussian66.
    """
    def __init__(self,target,anchor):
        self.target=target;target.jacobian(anchor)
        self.Cchi=np.einsum('qp,ps->qsp',anchor.T,anchor.state.field).reshape(256,576)
        self.anchor_key=anchor.key;self.cache=OrderedDict()
    def __call__(self,x):
        key=np.asarray(x).tobytes()
        if key not in self.cache:
            q=from_probit(x);D=self.target.prior.derivative(q)
            Z=(self.Cchi@D).reshape(32,8,66)
            J=np.concatenate((Z.real.transpose(1,0,2),Z.imag.transpose(1,0,2)),axis=1).reshape(512,66)*np.sqrt(2)/self.target.sigma
            J*=probit_derivative(x)[None,:]
            H=J.T@J+np.eye(66);H=(H+H.T)/2
            self.cache[key]=H
            while len(self.cache)>4:self.cache.popitem(last=False)
            self.target.book.counts['cached_material_decoder_metric_builds']+=1
            self.target.book.counts['cached_material_decoder_metric_complex_columns']+=66
        return self.cache[key]

class Posterior:
    def __init__(self,target,y,*,probit):
        self.target,self.y,self.probit=target,y,probit;self.cache=OrderedDict();self.current=None
    def point(self,x):
        x=np.asarray(x);k=x.tobytes()
        if k not in self.cache:
            q=from_probit(x) if self.probit else x
            p=self.target.evaluate(q,'FULL_POSTERIOR')
            self.cache[k]=p
            while len(self.cache)>5:
                old=next(v for v in self.cache if v!=self.current and v!=k);del self.cache[old]
        else:self.target.book.counts['identical_target_point_reuses']+=1
        return self.cache[k]
    def __call__(self,x):
        if not np.all(np.isfinite(x)):return -np.inf
        if not self.probit and not self.target.prior.legal(x):return -np.inf
        p=self.point(x);z=x if self.probit else x[:64]
        return -self.target.phi(p,self.y)-.5*float(z@z)
    def gradient(self,x):
        p=self.point(x);g=self.target.phi_gradient(p,self.y)
        if self.probit:return -g*probit_derivative(x)-x
        return -g-np.r_[np.asarray(x)[:64],np.zeros(2)]
    def pin(self,x):self.current=np.asarray(x).tobytes()

def metric(target,p,x):
    J=target.jacobian(p)*probit_derivative(x)[None,:]
    H=J.T@J+np.eye(66)
    C=cho_solve(cho_factor(H,lower=True),np.eye(66));return (C+C.T)/2

def probit_mode_initialization(target,y,anchor,chain,rng,*,max_evaluations=48):
    """Initialization only: optimize the correctly transformed target.

    A native-box MAP can be on u=1, whose probit coordinate is infinite.
    This finite source of bad initial energies is distinct from production
    mixing. Optimize Phi(Tx)+||x||²/2 with its exact chain rule; this is not
    another native-q MAP image and does not alter the posterior target.
    """
    from scipy.optimize import least_squares
    post=Posterior(target,y,probit=True)
    x0=to_probit(anchor.q,initialization=True)
    x0[64:]=to_probit(np.r_[np.zeros(64),np.array([[.2,.2],[.4,.8],[.65,.35],[.85,.65]])[chain%4]])[64:]
    history=[];points={};best=[np.inf,None,None]
    def point(x):
        k=np.asarray(x).tobytes()
        if k not in points:
            p=post.point(x);points[k]=p
            value=target.phi(p,y)+.5*float(x@x)
            history.append(dict(endpoint=p.key,objective=value,x=np.asarray(x).copy()))
            if value<best[0]:best[:]=[value,np.asarray(x).copy(),p]
        return points[k]
    def fun(x):return np.r_[target.white(point(x).state.data-y),x]
    def jac(x):return np.vstack((target.jacobian(point(x))*probit_derivative(x)[None,:],np.eye(66)))
    result=least_squares(fun,x0,jac=jac,method='trf',tr_solver='exact',x_scale='jac',max_nfev=max_evaluations,ftol=1e-8,xtol=1e-8,gtol=1e-6)
    mode_value,mode_x,mode_p=best;C=metric(target,mode_p,mode_x);L=np.linalg.cholesky(C)
    perturb=[.25,.5,.75,1.][chain%4]*(L@rng.normal(size=66));trials=[]
    # Finite burn-in starts; not rejection samples and not posterior draws.
    for shrink in range(7):
        x=mode_x+perturb*(.5**shrink);p=point(x);value=target.phi(p,y)+.5*float(x@x)
        trials.append(dict(shrink=shrink,objective=value,endpoint=p.key))
        if value<=mode_value+40:break
    C=metric(target,p,x)
    write(target.book.jobdir/'PROBIT_MODE_INITIALIZATION.json',dict(kind='transformed-target optimization and finite energy-dispersed starts, burn-in only',native_MAP_unchanged=True,chain=chain,GT_used=False,mode_objective=mode_value,mode_endpoint=mode_p.key,scipy_optimality=float(result.optimality),scipy_status=int(result.status),mode_history=history,perturbation_trials=trials,max_optimization_evaluations=max_evaluations,not_posterior_initial_draws=True))
    return x,p,C

def run_reference(target,parent,chain,*,method='mala',warmup=1024,production=4096,seed=2026102230):
    book=target.book;y=target.observation(parent)
    anchor=target.load_cache(OUT/f'A/anchors/{parent}.npz',owner='P3L A deterministic inversion')
    rng=np.random.default_rng(np.random.SeedSequence([seed,int(parent.split('_')[-1]),chain,0 if method=='mala' else 1]))
    transformed=method.startswith('mala');dynamic=method in ['mala_decoder','mala_decoder_mode'];mode_init=method.endswith('_mode');retracted=method in ['mala_retracted','mala_decoder'] or mode_init;post=Posterior(target,y,probit=transformed)
    if mode_init:
        init_x,pinit,C=probit_mode_initialization(target,y,anchor,chain,rng)
        x=init_x if transformed else from_probit(init_x)
        post.cache[x.tobytes()]=pinit
    elif transformed:
        if retracted:
            # Initialization, not a posterior transition. Explore weak prior
            # coordinates/property interiors, then repair strong energy only.
            from .conditional import make_frame,conditional_reference
            frame=make_frame(anchor.q,target.jacobian(anchor),8,kind='posterior_material_risk',material_derivative=target.prior.derivative(anchor.q))
            a=frame['gaussian_prior_weak_center']+.5*rng.normal(size=8)
            u=np.array([[.2,.2],[.4,.8],[.65,.35],[.85,.65]])[chain%4]
            init=conditional_reference(target,y,frame,a,u,max_evaluations=8)
            pinit=init['point'];x=to_probit(pinit.q);C=metric(target,pinit,x)
            write(book.jobdir/'INITIALIZATION.json',dict(kind='weak-coordinate conditional repair, not posterior draw',weak_coordinates=a,properties=u,conditional_initial_Phi=target.phi(pinit,y),conditional_gradient_inf=init['strong_gradient_inf'],GT_used=False,initial_endpoint=pinit.key,init_cost=init['counts']))
            post.cache[x.tobytes()]=pinit
        else:
            center=to_probit(anchor.q,initialization=True)
            C=metric(target,anchor,center);L=np.linalg.cholesky(C)
            x=center+[.5,1,1.5,2][chain%4]*(L@rng.normal(size=66))
    else:
        x=anchor.q.copy();x[:64]+=.03*rng.normal(size=64);C=None
    logp=post(x);g=post.gradient(x) if transformed else None;current=ChainState(x,logp,g);post.pin(x)
    mobility=FrozenMaterialMobility(target,anchor) if dynamic else None
    h=.05 if retracted else .25;beta=.02;prop_sd=.01
    # Standard dual averaging; reset only in the preregistered warmup windows.
    da_t=0;da_mu=np.log(10*h);da_hbar=0.;da_logbar=np.log(h)
    qs=[];logps=[];phis=[];accepted=[];evals=[];keys=[];xs=[];metadata=[]
    total=warmup+production;started=time.perf_counter()
    def checkpoint(i):
        arrays=dict(q=np.asarray(qs),x=np.asarray(xs),log_density=np.asarray(logps),Phi=np.asarray(phis),accepted=np.asarray(accepted),target_evaluations=np.asarray(evals),endpoint_key=np.asarray(keys),warmup=np.array(warmup),target_family=np.array('FULL'),log_mh_ratio=np.array([v.get('log_mh_ratio',np.nan) for v in metadata]),log_proposal_ratio=np.array([v.get('log_proposal_ratio',np.nan) for v in metadata]),gradient_evaluations=np.array([v.get('gradient_evaluations',0) for v in metadata]))
        save(book.jobdir/'CHAIN.npz',**arrays)
        write(book.jobdir/'CHECKPOINT.json',dict(parent=parent,method=method,chain=chain,iteration=i,total=total,seed=seed,RNG_state=rng.bit_generator.state,step_size=h,beta=beta,property_sd=prop_sd,warmup=warmup,production=production,target_family='FULL',all_rejections_retained=True,elapsed=time.perf_counter()-started,metric_frozen=i>=warmup))
        p=post.point(current.q);target.persist_cache(p,book.jobdir/'cache/CURRENT.npz')
        print(json.dumps(dict(parent=parent,chain=chain,method=method,iteration=i,total=total,acceptance=float(np.mean(accepted[-256:])),Phi=phis[-1],step_size=h)),flush=True)
    kernel=None
    for i in range(total):
        if transformed:
            if i in (256,512,768) and i<warmup and not dynamic:
                C=metric(target,post.point(current.q),current.q)
                save(book.jobdir/f'METRIC_WARMUP_{i}.npz',covariance=C,q=from_probit(current.q))
                if retracted:
                    h=.05;da_t=0;da_mu=np.log(10*h);da_hbar=0.;da_logbar=np.log(h)
            if kernel is None or i<=warmup:
                if dynamic:
                    from .kernels import PositionDependentMALA
                    kernel=PositionDependentMALA(post,post.gradient,h,mobility,dimension=66)
                else:kernel=PreconditionedMALA(post,post.gradient,h,C,dimension=66)
            tr=kernel.step(current,rng)
        elif method in ['pcn','pcn_mode']:
            tr=PCN(post,beta=beta,property_sd=prop_sd).step(current,rng)
        elif method in ['ess','ess_mode']:
            field=EllipticalSlice(post,max_evaluations=None).step(current,rng)
            prop=PropertyReflection(post,sd=prop_sd).step(field.state,rng)
            tr=prop
            tr.metadata['target_evaluations']+=field.metadata['target_evaluations']
            tr.metadata['field_rejected_proposals']=field.metadata['rejected_proposals']
            tr.metadata['field_accepted']=field.metadata['accepted']
            tr.metadata['property_accepted']=prop.metadata['accepted']
        else:raise ValueError(method)
        current=tr.state;post.pin(current.q)
        a=bool(tr.metadata['accepted'])
        if i<warmup:
            rate=.04/(1+i/100)**.6
            if transformed and retracted:
                da_t+=1;alpha=float(np.exp(min(0.,tr.metadata['log_mh_ratio'])));eta=1/(da_t+10)
                da_hbar=(1-eta)*da_hbar+eta*(.574-alpha)
                logh=da_mu-np.sqrt(da_t)/.05*da_hbar;weight=da_t**(-.75)
                da_logbar=weight*logh+(1-weight)*da_logbar
                h=float(np.clip(np.exp(np.clip(da_logbar if i==warmup-1 else logh,-20,2)),1e-6,2.))
            elif transformed:h=float(np.clip(h*np.exp(rate*(a-.574)),1e-5,2.))
            else:
                beta=float(np.clip(beta*np.exp(rate*(a-.234)),1e-5,.4))
                prop_sd=float(np.clip(prop_sd*np.exp(rate*(a-.234)),1e-5,.2))
        p=post.point(current.q);q=from_probit(current.q) if transformed else current.q
        qs.append(q.copy());xs.append(current.q.copy());logps.append(current.log_density);phis.append(target.phi(p,y));accepted.append(a);evals.append(tr.metadata['target_evaluations']);keys.append(p.key)
        metadata.append(tr.metadata)
        if (i+1)%256==0 or i+1==total:checkpoint(i+1)
    save(book.jobdir/'METRIC.npz',covariance=C if transformed else np.zeros((0,0)),is_production_covariance=np.array(not dynamic),description=np.array('initial covariance only; production uses saved decoder mobility rule' if dynamic else 'last warmup frozen covariance'))
    if dynamic:save(book.jobdir/'DECODER_MOBILITY.npz',frozen_material_operator=mobility.Cchi,final_precision=mobility(current.q),anchor_key=np.array(mobility.anchor_key))
    write(book.jobdir/'TRANSITION_SUMMARY.json',dict(parent=parent,chain=chain,method=method,warmup_acceptance=float(np.mean(accepted[:warmup])),production_acceptance=float(np.mean(accepted[warmup:])),total_target_evaluations=int(sum(evals)),warmup_adaptation='step only; fixed deterministic decoder-mobility rule' if dynamic else 'step and GN metric only',production_adaptation=False,position_dependent_mobility=dynamic,anchor_material_factor=anchor.key if dynamic else None,target_family='FULL',full_likelihood=True,posterior_normalizer_not_required=True,probability_status='EXACT_INVARIANT_KERNEL_NOT_A_MIXING_CERTIFICATE',coordinate='x=(z,eta);u=normal_CDF(eta)' if transformed else 'q=(z,u)',cost_scope='new physical actions in campaign ledger; same paid y anchor shared as explicitly disclosed setup'))
    return book.jobdir/'CHAIN.npz'
