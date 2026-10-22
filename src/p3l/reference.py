"""Full-target, audited MCMC; transform and preconditioner never change prior."""
from collections import OrderedDict
import json,time
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from .core import OUT,save,write,from_probit,to_probit,probit_derivative
from .kernels import ChainState,PreconditionedMALA,PCN,EllipticalSlice,PropertyReflection

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

def run_reference(target,parent,chain,*,method='mala',warmup=1024,production=4096,seed=2026102230):
    book=target.book;y=target.observation(parent)
    anchor=target.load_cache(OUT/f'A/anchors/{parent}.npz',owner='P3L A deterministic inversion')
    rng=np.random.default_rng(np.random.SeedSequence([seed,int(parent.split('_')[-1]),chain,0 if method=='mala' else 1]))
    transformed=method=='mala';post=Posterior(target,y,probit=transformed)
    if transformed:
        center=to_probit(anchor.q,initialization=True)
        C=metric(target,anchor,center);L=np.linalg.cholesky(C)
        x=center+[.5,1,1.5,2][chain%4]*(L@rng.normal(size=66))
    else:
        x=anchor.q.copy();x[:64]+=.03*rng.normal(size=64);C=None
    logp=post(x);g=post.gradient(x) if transformed else None;current=ChainState(x,logp,g);post.pin(x)
    h=.25;beta=.02;prop_sd=.01
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
            if i in (256,512,768) and i<warmup:
                C=metric(target,post.point(current.q),current.q)
                save(book.jobdir/f'METRIC_WARMUP_{i}.npz',covariance=C,q=from_probit(current.q))
            if kernel is None or i<=warmup:kernel=PreconditionedMALA(post,post.gradient,h,C,dimension=66)
            tr=kernel.step(current,rng)
        elif method=='pcn':
            tr=PCN(post,beta=beta,property_sd=prop_sd).step(current,rng)
        elif method=='ess':
            field=EllipticalSlice(post,max_evaluations=None).step(current,rng)
            prop=PropertyReflection(post,sd=prop_sd).step(field.state,rng)
            tr=prop
            tr.metadata['target_evaluations']+=field.metadata['target_evaluations']
            tr.metadata['field_rejected_proposals']=field.metadata['rejected_proposals']
        else:raise ValueError(method)
        current=tr.state;post.pin(current.q)
        a=bool(tr.metadata['accepted'])
        if i<warmup:
            rate=.04/(1+i/100)**.6
            if transformed:h=float(np.clip(h*np.exp(rate*(a-.574)),1e-5,2.))
            else:
                beta=float(np.clip(beta*np.exp(rate*(a-.234)),1e-5,.4))
                prop_sd=float(np.clip(prop_sd*np.exp(rate*(a-.234)),1e-5,.2))
        p=post.point(current.q);q=from_probit(current.q) if transformed else current.q
        qs.append(q.copy());xs.append(current.q.copy());logps.append(current.log_density);phis.append(target.phi(p,y));accepted.append(a);evals.append(tr.metadata['target_evaluations']);keys.append(p.key)
        metadata.append(tr.metadata)
        if (i+1)%256==0 or i+1==total:checkpoint(i+1)
    save(book.jobdir/'METRIC.npz',covariance=C if transformed else np.zeros((0,0)))
    write(book.jobdir/'TRANSITION_SUMMARY.json',dict(parent=parent,chain=chain,method=method,warmup_acceptance=float(np.mean(accepted[:warmup])),production_acceptance=float(np.mean(accepted[warmup:])),total_target_evaluations=int(sum(evals)),warmup_adaptation='step and GN metric only',production_adaptation=False,target_family='FULL',full_likelihood=True,posterior_normalizer_not_required=True,probability_status='EXACT_INVARIANT_KERNEL_NOT_A_MIXING_CERTIFICATE',coordinate='x=(z,eta);u=normal_CDF(eta)' if transformed else 'q=(z,u)',cost_scope='new physical actions in campaign ledger; same paid y anchor shared as explicitly disclosed setup'))
    return book.jobdir/'CHAIN.npz'
