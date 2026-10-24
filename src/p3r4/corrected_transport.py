"""Explicit-q independence MH and constant mixtures with prior pCN; pure algebra.

No physical backend, GT, adaptation, sampler filtering, or automatic expansion.
"""
import numpy as np

KERNELS=('pcn','independence_q','mixture_half')
BETA=.002
MIXTURE_Q_PROBABILITY=.5


def component_for(kernel,selector):
    if kernel not in KERNELS or not 0<=selector<1:raise ValueError('fixedkernel/independent uniformselector required')
    return 'q' if kernel=='independence_q' or (kernel=='mixture_half' and selector<.5) else 'pcn'


def log_acceptance(component,x,current_phi,proposed,proposed_phi,logq_current=None,logq_proposed=None):
    x,proposed=np.asarray(x,float),np.asarray(proposed,float)
    if x.shape!=(66,) or proposed.shape!=(66,) or not all(np.isfinite(a).all() for a in (x,proposed)):raise ValueError('finite full66x')
    if not np.isfinite(current_phi) or not np.isfinite(proposed_phi):raise ValueError('finite paidFullpotentials')
    ratio=float(current_phi-proposed_phi)
    if component=='q':
        if not np.isfinite(logq_current) or not np.isfinite(logq_proposed):raise ValueError('both explicit rawdx66 density values required')
        ratio+=.5*float(x@x-proposed@proposed)+float(logq_current-logq_proposed)
    elif component!='pcn':raise ValueError('q/pcn component required; never use an incorrect mixture-density ratio')
    return min(0.,ratio)


def seed_for(seed,parent_index,chain,step,role):
    return int(np.random.SeedSequence([int(seed),parent_index,chain,step,role]).generate_state(1)[0])


def step_randomness(seed,parent_index,chain,step):
    # Separate state-independent streams let common pCN innovations/U and q
    # proposal seeds agree across kernels whenever that component is used.
    innovation=np.random.default_rng(seed_for(seed,parent_index,chain,step,1)).normal(size=66)
    logu=float(np.log(np.random.default_rng(seed_for(seed,parent_index,chain,step,2)).uniform()))
    selector=float(np.random.default_rng(seed_for(seed,parent_index,chain,step,3)).uniform())
    qseed=seed_for(seed,parent_index,chain,step,4)
    return innovation,logu,selector,qseed


def pcn_proposal(x,innovation):return np.sqrt(1-BETA*BETA)*np.asarray(x)+BETA*np.asarray(innovation)


def multichain_diagnostics(values):
    """One parent/kernel only; preserve every poststep rejection repeat.

    Conventional split-Rhat and initial-positive ACF diagnostics are descriptive
    finite-run numbers. Constant/stuck chains are never assigned ESS=N/MCSE0.
    """
    a=np.asarray(values,float)
    if a.ndim!=3 or a.shape[0]!=4 or a.shape[1]<2 or not np.isfinite(a).all():raise ValueError('4chains x frozenproductionsteps x functions')
    chains,n,functions=a.shape;half=n//2
    if half>=2:
        split=np.concatenate((a[:,:half],a[:,-half:]),axis=0);split_variance=split.var(axis=1,ddof=1)
        split_variance[np.all(split==split[:,0:1,:],axis=1)]=0.
        within=split_variance.mean(0);between=half*split.mean(1).var(0,ddof=1)
        estimate=(half-1)/half*within+between/half
        rhat=np.sqrt(np.divide(estimate,within,out=np.full_like(within,np.nan),where=within>0));rhat[(within==0)&(between>0)]=np.inf
    else:rhat=np.full(functions,np.nan)
    ess=np.full((chains,functions),np.nan);mcse=np.full_like(ess,np.nan)
    for chain in range(chains):
        for f in range(functions):
            x=a[chain,:,f]
            if np.all(x==x[0]):continue
            variance=float(x.var(ddof=1))
            if variance<=0:continue
            centered=x-x.mean();cov=np.correlate(centered,centered,mode='full')[n-1:]/np.arange(n,0,-1);acf=cov/cov[0];total=0.
            for lag in range(1,min(n-1,512),2):
                pair=acf[lag]+acf[lag+1]
                if pair<=0:break
                total+=pair
            ess[chain,f]=min(n,n/max(1.,1+2*total));mcse[chain,f]=np.sqrt(variance/ess[chain,f])
    return dict(status='FINITE_RUN_UNQUALIFIED',production_steps=n,chains=chains,split_Rhat=rhat.tolist(),ESS_per_chain=ess.tolist(),MCSE_mean_per_chain=mcse.tolist(),
      constant_or_stuck_chain_function_count=int(np.isnan(ess).sum()),split_Rhat_infinite_function_indices=np.flatnonzero(np.isinf(rhat)).tolist(),split_Rhat_unavailable_function_indices=np.flatnonzero(np.isnan(rhat)).tolist(),Rhat_method='conventional2halvesofeach4chains, same parent/target only; not ranknormalized',
      ESS_method='initial-positive paired ACF perchain, repeats retained',common_Full_reference='NONE',posterior_qualification='NOT_ESTABLISHED')
