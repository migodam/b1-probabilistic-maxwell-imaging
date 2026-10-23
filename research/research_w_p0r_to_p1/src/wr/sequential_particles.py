"""Generic four-source finite-quadrature SMC; no data/truth/PDE/Book/I/O.

Node ravel order is(branch,iu,ic,ia,iphi). The supplied proper continuous prior
is represented by caller-supplied finite quadrature masses. Results concern
that finite approximation, not arbitrary-wall or continuous posterior proofs.
Local torus moves wrap GRID INDICES only; no physical periodicity is assumed.
No function runs at import; tiny_algebra_checks() is an optional tiny audit.
"""
from __future__ import annotations
import math
import numpy as np

STAGES=4


def _lse(x):
    x=np.asarray(x,dtype=float);maximum=np.max(x)
    if not np.isfinite(maximum):return float(maximum)
    return float(maximum+np.log(np.sum(np.exp(x-maximum))))


def _shape(shape,n):
    shape=tuple(shape)
    if len(shape)!=5 or shape[0]!=2 or any(isinstance(v,(bool,np.bool_)) or int(v)!=v or v<1 for v in shape):
        raise ValueError('grid_shape must be(2,nu,nc,ng,ng) with positive integer axes')
    shape=tuple(int(v) for v in shape)
    if shape[3]!=shape[4] or math.prod(shape)!=n:
        raise ValueError('grid_shape must match node count and use equal amplitude/phase quadrature counts')
    return shape


def _prior(prior,n=None):
    p=np.asarray(prior)
    if np.iscomplexobj(p) or p.ndim!=1 or (n is not None and len(p)!=n) or not len(p):
        raise ValueError('prior_weights must be a real nonempty node vector')
    p=np.asarray(p,dtype=float)
    if not np.isfinite(p).all() or np.any(p<0) or not np.any(p>0):
        raise ValueError('prior weights must be finite nonnegative and have positive mass')
    mass=float(np.sum(p))
    if not np.isfinite(mass) or mass<=0:raise ValueError('Prior mass must be finite positive')
    return p/mass,mass


def _likelihood(blocks,n=None):
    x=np.asarray(blocks)
    if np.iscomplexobj(x) or x.ndim!=2 or x.shape[1]!=4 or (n is not None and len(x)!=n):
        raise ValueError('loglike_blocks must have real shape(Nnodes,4)')
    x=np.asarray(x,dtype=float)
    if np.isnan(x).any() or np.isposinf(x).any():raise ValueError('Log likelihood may contain-inf, but not NaN/+inf')
    return x


def beta_schedule(schedule):
    """Exactly four fixed stages. Both schedules terminate at all-source beta1."""
    if schedule=='direct':return np.repeat((np.arange(1,5,dtype=float)/4)[:,None],4,axis=1)
    if schedule=='progressive':return np.tril(np.ones((4,4),dtype=float))
    raise ValueError('schedule must be direct or progressive')


def _at_beta(blocks,beta):
    # Explicitly skip beta0 so(-inf)*0 never creates NaN.
    result=np.zeros(len(blocks),dtype=float)
    for source,coefficient in enumerate(beta):
        if coefficient:result+=coefficient*blocks[:,source]
    return result


def _draw_prior(p,uniforms):
    cdf=np.cumsum(p);cdf[-1]=1.
    return np.searchsorted(cdf,uniforms,side='right').astype(np.int64)


def systematic_resampling(weights,offset_uniform):
    """One U[0,1) offset; exactly P systematic offspring indexes."""
    w=np.asarray(weights,dtype=float)
    if w.ndim!=1 or not len(w) or not np.isfinite(w).all() or np.any(w<0) or not np.isclose(w.sum(),1,rtol=0,atol=1e-12):
        raise ValueError('Systematic weights must be finite normalized probabilities')
    if not np.isscalar(offset_uniform) or not 0<=offset_uniform<1:
        raise ValueError('offset_uniform must lie in[0,1)')
    points=(float(offset_uniform)+np.arange(len(w)))/len(w)
    cdf=np.cumsum(w);cdf[-1]=1.
    return np.searchsorted(cdf,points,side='right').astype(np.int64)


def make_proposal_library(prior_weights,grid_shape,*,rngseed=0,nparticles=2048,sweeps=4,max_array_bytes=128_000_000):
    """Fixed common random numbers, drawn UP FRONT for both schedules.

    Includes unused proposals/acceptance uniforms, so rejection or resampling
    never changes later RNG consumption. Global proposals use proper prior
    masses. Local coordinate/sign choices are independent of current state.
    """
    p,mass=_prior(prior_weights);shape=_shape(grid_shape,len(p))
    if isinstance(nparticles,bool) or int(nparticles)!=nparticles or nparticles<1:
        raise ValueError('nparticles must be a positive integer')
    if isinstance(sweeps,bool) or int(sweeps)!=sweeps or sweeps<0:
        raise ValueError('sweeps must be a nonnegative integer')
    P,S=int(nparticles),int(sweeps)
    estimate=16*len(p)+8*P+4*S*P*(8+8+1+1+8)+32
    if not np.isfinite(max_array_bytes) or max_array_bytes<=0 or estimate>max_array_bytes:
        raise MemoryError('Proposal library array budget exceeded')
    rng=np.random.default_rng(rngseed);size=(4,S,P)
    initial=_draw_prior(p,rng.random(P))
    global_nodes=_draw_prior(p,rng.random(size))
    library=dict(grid_shape=shape,prior_weights=p.copy(),prior_input_mass=mass,
        nparticles=P,sweeps=S,rngseed=rngseed,initial_nodes=initial,estimated_array_bytes=int(estimate),
        global_nodes=global_nodes,global_uniforms=rng.random(size),
        coordinates=rng.integers(0,5,size=size,dtype=np.int8),
        signs=(2*rng.integers(0,2,size=size,dtype=np.int8)-1),
        acceptance_uniforms=rng.random(size),resampling_uniforms=rng.random(4))
    # A library is reusable without mutation across direct/progressive runs.
    for value in library.values():
        if isinstance(value,np.ndarray):value.flags.writeable=False
    return library


def _library_check(library,p,shape,P,S):
    if tuple(library['grid_shape'])!=shape or library['nparticles']!=P or library['sweeps']!=S or not np.array_equal(library['prior_weights'],p):
        raise ValueError('Proposal library differs from grid/prior/particle/sweep contract')
    size=(4,S,P)
    for key in ('global_nodes','global_uniforms','coordinates','signs','acceptance_uniforms'):
        if np.shape(library[key])!=size:raise ValueError('Malformed proposal library:'+key)
    if np.shape(library['initial_nodes'])!=(P,) or np.shape(library['resampling_uniforms'])!=(4,):
        raise ValueError('Malformed initial/resampling library')
    for key in ('initial_nodes','global_nodes'):
        v=np.asarray(library[key])
        if v.dtype.kind not in 'iu' or np.any(v<0) or np.any(v>=len(p)) or np.any(p[v]<=0):
            raise ValueError('Prior proposal library includes invalid/zero-prior nodes')
    for key in ('global_uniforms','acceptance_uniforms','resampling_uniforms'):
        v=np.asarray(library[key])
        if not np.isfinite(v).all() or np.any((v<0)|(v>=1)):raise ValueError('Invalid library uniform')
    if np.asarray(library['coordinates']).dtype.kind not in 'iu' or np.asarray(library['signs']).dtype.kind not in 'iu':
        raise ValueError('Local coordinate/sign arrays must be integer typed')
    if np.any((library['coordinates']<0)|(library['coordinates']>4)) or np.any((library['signs']!=1)&(library['signs']!=-1)):
        raise ValueError('Invalid local coordinate/sign library')


def _neighbors(nodes,coordinates,signs,shape):
    result=nodes.copy()
    for axis,length in enumerate(shape):
        use=coordinates==axis
        stride=math.prod(shape[axis+1:])
        coordinate=(nodes[use]//stride)%length
        neighbor=(coordinate+signs[use])%length
        result[use]+=(neighbor-coordinate)*stride
    return result


def _log_ratio(new,old):
    result=np.full(np.shape(new),-np.inf,dtype=float)
    both=np.isfinite(new)&np.isfinite(old)
    result[both]=new[both]-old[both]
    result[np.isfinite(new)&np.isneginf(old)]=np.inf
    # Both zero-target states => reject; a zero-target particle has zero weight.
    return result


def exact_finite_reference(loglike_blocks,prior_weights):
    blocks=_likelihood(loglike_blocks);p,mass=_prior(prior_weights,len(blocks))
    with np.errstate(divide='ignore'):lp=np.log(p)
    total=np.sum(blocks,axis=1);valid=(p>0)&np.isfinite(total)
    if not np.any(valid):raise ValueError('Finite quadrature target has zero evidence')
    anchor=float(np.max(total[valid]));score=lp+(total-anchor);relative=_lse(score);logz=anchor+relative
    if not np.isfinite(logz):raise ValueError('Finite quadrature target has zero evidence')
    weights=np.exp(score-relative);weights/=weights.sum()
    return dict(log_evidence=logz,log_evidence_raw_prior_quadrature=logz+np.log(mass),
        posterior_weights=weights,prior_input_mass=mass,
        prior_handling='normalized finite quadrature prior; raw prior mass exposed',
        interpretation='exact sum on supplied finite quadrature, not exact continuous/arbitrary-wall posterior',
        P1_gate='NOT_EVALUATED_ALGEBRA_ONLY')


def smc(loglike_blocks,prior_weights,grid_shape,*,schedule='direct',rngseed=0,
        nparticles=2048,sweeps=4,global_prob=.25,proposal_library=None,
        max_array_bytes=128_000_000):
    """Four fixed incremental likelihood stages, weighted final SMC particles.

    ESS<.5P triggers systematic resampling. Mutation kernels leave each stage's
    finite target invariant: prior-independence proposals accept delta(logL),
    symmetric coordinate torus moves accept delta(logL)+delta(logprior).
    Rejections/self-proposals are recorded as repeated states. Resampling is
    recorded even when the identity map is used. Do not pool trajectories from
    changing beta targets as if they were draws from the final posterior.
    """
    blocks=_likelihood(loglike_blocks);p,mass=_prior(prior_weights,len(blocks));shape=_shape(grid_shape,len(p))
    if isinstance(nparticles,bool) or int(nparticles)!=nparticles or nparticles<1:raise ValueError('Positive integer nparticles required')
    if isinstance(sweeps,bool) or int(sweeps)!=sweeps or sweeps<0:raise ValueError('Nonnegative integer sweeps required')
    if not np.isscalar(global_prob) or not np.isfinite(global_prob) or not 0<=global_prob<=1:raise ValueError('global_prob must lie in[0,1]')
    P,S=int(nparticles),int(sweeps);T=1+4*(2+S)
    # Conservative working-array estimate, not a claim about total process RSS.
    estimated=8*len(p)*12+(4*S*P)*(8+8+8+1+1)+T*P*16+4*P*8+P*128
    if not np.isfinite(max_array_bytes) or max_array_bytes<=0 or estimated>max_array_bytes:
        raise MemoryError('SMC array budget exceeded before proposal/trajectory allocation')
    library=proposal_library if proposal_library is not None else make_proposal_library(
        prior_weights,shape,rngseed=rngseed,nparticles=P,sweeps=S,max_array_bytes=max_array_bytes)
    _library_check(library,p,shape,P,S)
    beta=beta_schedule(schedule);previous_beta=np.zeros(4)
    states=np.array(library['initial_nodes'],copy=True);weights=np.full(P,1/P)
    ancestors=np.arange(P,dtype=np.int64)
    trajectories=np.empty((T,P),dtype=np.int64);trajectory_weights=np.empty((T,P),dtype=float)
    stages=np.empty(T,dtype=np.int8);kinds=[];sweep_ids=np.empty(T,dtype=np.int16);cursor=0
    resample_indexes=np.empty((4,P),dtype=np.int64);resampled=np.zeros(4,dtype=bool)
    accepted=np.zeros((4,S,P),dtype=bool);global_flags=np.asarray(library['global_uniforms'])<global_prob
    stage_rows=[];log_evidence=0.
    with np.errstate(divide='ignore'):logprior=np.log(p)
    def record(stage,kind,sweep=-1):
        nonlocal cursor
        trajectories[cursor]=states;trajectory_weights[cursor]=weights
        stages[cursor]=stage;sweep_ids[cursor]=sweep;kinds.append(kind);cursor+=1
    record(0,'initial_prior')
    for stage in range(4):
        delta=beta[stage]-previous_beta
        increment=_at_beta(blocks,delta)[states]
        with np.errstate(divide='ignore'):logw=np.log(weights)
        valid=(weights>0)&np.isfinite(increment)
        if not np.any(valid):raise ValueError('All particles have zero incremental likelihood mass')
        anchor=float(np.max(increment[valid]));relative_increment=increment-anchor
        relative_logz=_lse(logw+relative_increment);increment_logz=anchor+relative_logz
        if not np.isfinite(increment_logz):raise ValueError('All particles have zero incremental likelihood mass')
        log_evidence+=increment_logz
        weights=np.exp(logw+relative_increment-relative_logz);weights/=weights.sum()
        ess=float(1/np.sum(weights**2));record(stage+1,'incremental_weight')
        if ess<.5*P:
            indexes=systematic_resampling(weights,library['resampling_uniforms'][stage])
            states=states[indexes];ancestors=ancestors[indexes];weights=np.full(P,1/P);resampled[stage]=True
        else:indexes=np.arange(P,dtype=np.int64)
        resample_indexes[stage]=indexes;record(stage+1,'resampling_or_identity')
        ll=_at_beta(blocks,beta[stage]);counts={'global':0,'local':0,'global_proposals':0,'local_proposals':0,'self_proposals':0,'accepted_state_changes':0}
        for sweep in range(S):
            global_move=global_flags[stage,sweep]
            proposals=_neighbors(states,library['coordinates'][stage,sweep],library['signs'][stage,sweep],shape)
            proposals[global_move]=library['global_nodes'][stage,sweep,global_move]
            old_like,new_like=ll[states],ll[proposals]
            ratio=_log_ratio(new_like,old_like)
            local=~global_move
            valid_local=local&(p[proposals]>0)
            ratio[valid_local]+=logprior[proposals[valid_local]]-logprior[states[valid_local]]
            # A zero-prior local node is always rejected, including old zeroL.
            ratio[p[proposals]<=0]=-np.inf
            with np.errstate(divide='ignore'):log_uniform=np.log(library['acceptance_uniforms'][stage,sweep])
            accept=log_uniform<np.minimum(0.,ratio)
            accepted[stage,sweep]=accept
            counts['global']+=int(np.count_nonzero(accept&global_move));counts['local']+=int(np.count_nonzero(accept&local))
            counts['global_proposals']+=int(np.count_nonzero(global_move));counts['local_proposals']+=int(np.count_nonzero(local))
            counts['self_proposals']+=int(np.count_nonzero(proposals==states))
            counts['accepted_state_changes']+=int(np.count_nonzero(accept&(proposals!=states)))
            states[accept]=proposals[accept];record(stage+1,'mutation_with_rejections',sweep)
        stage_rows.append(dict(stage=stage+1,beta=beta[stage].tolist(),delta_beta=delta.tolist(),
            incremental_log_evidence=increment_logz,ESS_before_resampling=ess,
            ESS_after_resampling=float(1/np.sum(weights**2)),resampled=bool(resampled[stage]),
            surviving_initial_ancestors=int(len(np.unique(ancestors))),unique_current_nodes=int(len(np.unique(states))),
            accept_counts=counts))
        previous_beta=beta[stage]
    assert cursor==T
    return dict(particle_indexes=states.copy(),weights=weights.copy(),log_evidence_estimate=log_evidence,
        log_evidence_raw_prior_quadrature_estimate=log_evidence+np.log(mass),
        evidence_note='SMC normalizer estimate; its log is not an unbiased log-evidence estimator',
        ancestor_indexes=ancestors.copy(),surviving_initial_ancestors=int(len(np.unique(ancestors))),
        stages=stage_rows,trajectory_indexes=trajectories,trajectory_weights=trajectory_weights,
        trajectory_stage=stages,trajectory_kind=kinds,trajectory_sweep=sweep_ids,
        mutation_accepted=accepted,mutation_global=global_flags,resampling_indexes=resample_indexes,
        resampling_performed=resampled,final_beta=beta[-1].copy(),schedule=schedule,
        proposal_library_seed=library['rngseed'],proposal_library_reusable=True,
        prior_input_mass=mass,prior_handling='normalized proper finite quadrature masses',
        estimated_working_array_bytes=int(estimated),max_array_bytes=int(max_array_bytes),
        interpretation='finite quadrature approximation to the declared controlled prior only; no arbitrary-wall/continuous-posterior certification',
        P1_gate='NOT_EVALUATED_ALGEBRA_ONLY',physics='NOT_RUN')


def tiny_algebra_checks():
    """Enumerable balance/ratio checks and a tiny structural path audit only.

    No convergence/accuracy conclusion is drawn from the16-particle toy paths.
    This does not qualify any Maxwell model, likelihood or P1 posterior.
    """
    shape=(2,2,1,1,1);n=math.prod(shape);prior=np.array([.1,.2,.3,.4])
    blocks=np.array([[0.,-.2,-.1,-.3],[-.4,0.,-.2,-.1],[-.1,-.7,0.,-.2],[-.3,-.2,-.5,0.]])
    qlocal=np.zeros((n,n))
    # Independent, fully enumerated neighbor proposal table.
    for i in range(n):
        coordinate=list(np.unravel_index(i,shape))
        for axis in range(5):
            for sign in (-1,1):
                target=coordinate.copy();target[axis]=(target[axis]+sign)%shape[axis]
                j=int(np.ravel_multi_index(tuple(target),shape));qlocal[i,j]+=.1
    assert np.max(abs(qlocal-qlocal.T))<1e-15
    worst_balance=0.;worst_ratio=0.;worst_local_ratio=0.
    for name in ('direct','progressive'):
        betas=beta_schedule(name);increments=np.diff(np.vstack([np.zeros(4),betas]),axis=0)
        assert np.array_equal(betas[-1],np.ones(4))
        assert np.max(abs(np.sum(blocks@increments.T,axis=1)-np.sum(blocks,axis=1)))<1e-15
        for beta in betas:
            ll=blocks@beta;target=prior*np.exp(ll);target/=target.sum();transition=np.zeros((n,n))
            for i in range(n):
                for j in range(n):
                    global_alpha=min(1.,np.exp(ll[j]-ll[i]))
                    local_alpha=min(1.,np.exp(ll[j]-ll[i])*prior[j]/prior[i])
                    proposal=.25*prior[j]*global_alpha+.75*qlocal[i,j]*local_alpha
                    transition[i,j]+=proposal
                    transition[i,i]+=.25*prior[j]*(1-global_alpha)+.75*qlocal[i,j]*(1-local_alpha)
                    independent_ratio=(target[j]*prior[i])/(target[i]*prior[j])
                    worst_ratio=max(worst_ratio,abs(np.log(independent_ratio)-(ll[j]-ll[i])))
                    if qlocal[i,j]>0:
                        worst_local_ratio=max(worst_local_ratio,abs(np.log(target[j]/target[i])-
                            (ll[j]-ll[i]+np.log(prior[j])-np.log(prior[i]))))
            assert np.max(abs(transition.sum(axis=1)-1))<1e-14
            balance=target[:,None]*transition
            worst_balance=max(worst_balance,float(np.max(abs(balance-balance.T))))
    assert worst_balance<1e-14 and worst_ratio<1e-14 and worst_local_ratio<1e-14
    assert np.array_equal(systematic_resampling(prior,.25),[0,2,2,3])
    library=make_proposal_library(prior,shape,rngseed=17,nparticles=16,sweeps=1)
    paths={name:smc(blocks,prior,shape,schedule=name,nparticles=16,sweeps=1,proposal_library=library)
           for name in ('direct','progressive')}
    for result in paths.values():
        assert result['trajectory_indexes'].shape==(13,16)
        assert abs(result['weights'].sum()-1)<1e-14
        assert np.array_equal(result['final_beta'],np.ones(4))
        assert len(result['trajectory_kind'])==13
        for stage in range(4):
            weight_frame=1+stage*3;resampling_frame=weight_frame+1;mutation_frame=weight_frame+2
            assert np.array_equal(result['trajectory_indexes'][resampling_frame],
                result['trajectory_indexes'][weight_frame][result['resampling_indexes'][stage]])
            rejected=~result['mutation_accepted'][stage,0]
            assert np.array_equal(result['trajectory_indexes'][mutation_frame,rejected],
                result['trajectory_indexes'][resampling_frame,rejected])
    exact=exact_finite_reference(blocks,prior)
    manual=np.sum(prior*np.exp(np.sum(blocks,axis=1)))
    assert abs(exact['log_evidence']-np.log(manual))<1e-14
    common=exact_finite_reference(np.full((4,4),-1e20),prior)
    assert np.max(abs(common['posterior_weights']-prior))<1e-15
    return dict(status='TINY_FINITE_ALGEBRA_CHECKS_COMPLETE',maximum_detailed_balance_error=worst_balance,
        maximum_global_forward_reverse_log_ratio_error=float(worst_ratio),
        maximum_local_forward_reverse_log_ratio_error=float(worst_local_ratio),
        large_common_loglikelihood_preserves_prior=True,
        incremental_likelihood_counted_once=True,direct_progressive_same_final_target=True,
        proposal_table_symmetric=True,systematic_example_indexes=[0,2,2,3],
        tiny_path_particles=16,tiny_path_sweeps=1,tiny_path_schedule_count=2,
        final_target_exact_log_evidence=exact['log_evidence'],
        rejected_and_resampled_state_records_verified=True,tiny_MC_accuracy_claim=False,P1_gate='NOT_EVALUATED',Maxwell='NOT_RUN',I_O='NONE')
