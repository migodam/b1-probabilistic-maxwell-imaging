"""Explicit finite-prior inference. No continuous-integration certificate.

All probabilities refer to FINITE_PRIOR_ENUM4_V1 physical states. Material
quadrature nodes define a discrete prior, not a certified continuous integral.
"""
import numpy as np
from .geometry import enumerate_supports

PRIOR_ID = 'FINITE_PRIOR_ENUM4_V1'
STALE_LABEL = 'STALE_WALL_MARGINAL'


def _normalize(weights):
    w = np.asarray(weights,dtype=float)
    if w.ndim != 1 or not w.size or np.any(~np.isfinite(w)) or np.any(w<0) or w.sum()<=0:
        raise ValueError('nonnegative finite weights with positive total required')
    return w/w.sum()


def _lognormalize(logweights):
    a = np.asarray(logweights,dtype=float)
    if a.ndim != 1 or a.size == 0 or np.any(np.isnan(a)) or np.any(np.isposinf(a)):
        raise ValueError('one-dimensional logweights must be finite or negative infinity')
    m = a.max()
    if not np.isfinite(m):
        raise ValueError('empty posterior support')
    return np.exp(a-m)/np.exp(a-m).sum()


def _grid(eps,ell,material_weight,pose_x,pose_y,pose_weight,bank=False):
    patterns = enumerate_supports(2,.65)
    legal = np.stack([s for s,_ in patterns])
    support_prior = np.array([p for _,p in patterns])
    rows = []
    # Wall-group IDs are stable across poses, enabling forward basis reuse.
    walls = []
    for sid,s in enumerate(legal):
        if not s.any():
            walls.append((sid,3.75,.125,float(support_prior[sid])))
        else:
            for i,e in enumerate(eps):
                for j,l in enumerate(ell):
                    walls.append((sid,e,l,float(support_prior[sid]*material_weight[i]*material_weight[j])))
    for iy,y in enumerate(pose_y):
        for ix,x in enumerate(pose_x):
            for wid,(sid,e,l,mass) in enumerate(walls):
                rows.append((sid,e,l,x,y,mass*pose_weight[ix]*pose_weight[iy],wid))
    a = np.asarray(rows,dtype=float)
    ids = a[:,0].astype(int)
    return dict(label='PROPOSAL_BANK_ONLY' if bank else PRIOR_ID,
                prior_id=PRIOR_ID,discrete_prior=True,
                continuous_uniform_integration_certificate=False,
                posterior_claims_allowed=not bank,
                support_id=ids,support=legal[ids].reshape(len(a),4),
                epsinf=a[:,1],ell=a[:,2],pose=a[:,3:5],
                prior_weight=_normalize(a[:,5]),wall_group_id=a[:,6].astype(int),
                legal_supports=legal,support_prior=support_prior,
                material_nodes=dict(epsinf=np.asarray(eps),ell=np.asarray(ell)),
                state_count=len(a),wall_group_count=len(walls),pose_outermost=True)


def build_grid():
    nodes = np.sort(np.cos(np.pi*np.arange(5)/4))
    x = nodes
    y = 2.1+.5*nodes
    cc5 = np.array([1/30,4/15,2/5,4/15,1/30])
    return _grid([1.5,3.75,6],[0,.125,.25],[1/6,4/6,1/6],x,y,cc5)


def build_coarse_grid():
    """A proposal bank; its weights do not establish posterior claims."""
    return _grid([2.5,5],[.0625,.1875],[.5,.5],[-.7,0,.7],[1.75,2.1,2.45],[1/3]*3,bank=True)


def joint_coordinates(grid):
    z = np.asarray(grid['support'],dtype=float).reshape(-1,4)
    contrast = z*(np.asarray(grid['epsinf'])[:,None]-1)
    loss = z*np.asarray(grid['ell'])[:,None]
    return np.column_stack((contrast,loss,np.asarray(grid['pose'])))


def _distribution(values,w,quantiles=(.025,.5,.975)):
    values = np.asarray(values)
    unique,inverse = np.unique(values,return_inverse=True)
    mass = np.bincount(inverse,weights=w,minlength=len(unique))
    cdf = np.cumsum(mass)
    q = np.asarray(quantiles,dtype=float)
    if np.any((q<0)|(q>1)):
        raise ValueError('quantiles must lie in [0,1]')
    indices = np.minimum(np.searchsorted(cdf,q,side='left'),len(unique)-1)
    return dict(values=unique,mass=mass,cdf=cdf,quantile_levels=q,quantiles=unique[indices])


def summarize(weights,grid,quantiles=(.025,.5,.975)):
    w = _normalize(weights)
    if len(w)!=len(grid['prior_weight']):
        raise ValueError('grid and weight lengths differ')
    z = np.asarray(grid['support'],dtype=float).reshape(-1,4)
    pose = np.asarray(grid['pose'],dtype=float)
    coords = joint_coordinates(grid)
    mean = w@coords
    centered = coords-mean
    cov = (centered*w[:,None]).T@centered
    mode_mass = np.bincount(grid['support_id'],weights=w,minlength=len(grid['legal_supports']))
    conditionals = []
    for region in range(4):
        occupied = z[:,region]>0
        probability = float(w[occupied].sum())
        if probability == 0:
            conditionals.append(dict(occupied_probability=0.,defined=False,epsinf=None,ell=None,joint=None))
            continue
        conditional_w = w[occupied]/probability
        pair = np.column_stack((grid['epsinf'][occupied],grid['ell'][occupied]))
        unique,inverse = np.unique(pair,axis=0,return_inverse=True)
        mass = np.bincount(inverse,weights=conditional_w,minlength=len(unique))
        conditionals.append(dict(occupied_probability=probability,defined=True,
            epsinf=_distribution(grid['epsinf'][occupied],conditional_w,quantiles),
            ell=_distribution(grid['ell'][occupied],conditional_w,quantiles),
            joint=dict(values=unique,mass=mass)))
    return dict(label=grid['label'],occupancy=w@z,conditional_material=conditionals,
                full_wall_mean=mean[:8],full_wall_cov=cov[:8,:8],
                joint_mean=mean,joint_cov=cov,pose_mean=mean[8:],pose_cov=cov[8:,8:],
                mode_probability=mode_mass,event_px_ge_0=float(w[pose[:,0]>=0].sum()),
                coordinate_distributions=[_distribution(coords[:,j],w,quantiles) for j in range(10)],
                coordinate_names=['z0_eps_contrast','z1_eps_contrast','z2_eps_contrast','z3_eps_contrast','z0_ell','z1_ell','z2_ell','z3_ell','px','py'])


def update_logweights(prior,uniqueLLlist):
    """Normalize prior × each observation likelihood exactly once.

    Prefer a mapping observation_id→loglikelihood, or (id, array) pairs; duplicate
    IDs are rejected. A bare array list is supported when caller guarantees IDs.
    Returns normalized logweights, retaining impossible states as negative infinity.
    """
    p = _normalize(prior)
    with np.errstate(divide='ignore'):
        total = np.log(p)
    items = uniqueLLlist.items() if isinstance(uniqueLLlist,dict) else uniqueLLlist
    seen = set()
    for item in items:
        if isinstance(item,tuple) and len(item)==2 and isinstance(item[0],(str,int)):
            identifier,ll = item
            if identifier in seen:
                raise ValueError('duplicate observation likelihood ID')
            seen.add(identifier)
        else:
            ll = item
        ll = np.asarray(ll,dtype=float)
        if ll.shape!=p.shape or np.any(np.isnan(ll)) or np.any(np.isposinf(ll)):
            raise ValueError('likelihood log grid must match prior and contain no NaN/+infinity')
        total = total+ll
    _lognormalize(total)  # Validate nonempty support without losing finite log tails.
    maximum = np.max(total)
    return total-(maximum+np.log(np.exp(total-maximum).sum()))


def conditional_target_only(oldw,newLL,wall_group):
    """Return stale-wall approximation: retain q(wall), update q(p|wall)."""
    w = _normalize(oldw)
    ll = np.asarray(newLL,dtype=float)
    groups = np.asarray(wall_group)
    if ll.shape!=w.shape or groups.shape!=w.shape or np.any(np.isnan(ll)) or np.any(np.isposinf(ll)):
        raise ValueError('matching weight, loglikelihood, and wall-group arrays required')
    result = np.zeros_like(w)
    for group in np.unique(groups):
        selected = groups==group
        mass = w[selected].sum()
        if mass>0:
            with np.errstate(divide='ignore'):
                local = np.log(w[selected])+ll[selected]
            result[selected] = mass*_lognormalize(local)
    return dict(label=STALE_LABEL,weights=result,exact_new_target_posterior=False,
                interpretation='old wall marginal retained; target conditional updated within each wall group')


def posterior_metrics(approx,ref,grid):
    """Finite-grid discrepancies; primary squared-loss action is target pose only."""
    a,r = _normalize(approx),_normalize(ref)
    if a.shape!=r.shape or len(a)!=len(grid['prior_weight']):
        raise ValueError('grid and posterior shapes differ')
    coords = joint_coordinates(grid)
    cdf_sup=[]
    for j in range(coords.shape[1]):
        _,inverse = np.unique(coords[:,j],return_inverse=True)
        difference = np.bincount(inverse,weights=a-r)
        cdf_sup.append(float(np.abs(np.cumsum(difference)).max()))
    ids = grid['support_id']
    am = np.bincount(ids,weights=a,minlength=len(grid['legal_supports']))
    rm = np.bincount(ids,weights=r,minlength=len(grid['legal_supports']))
    selected = rm>=.05
    pose = np.asarray(grid['pose'], dtype=float)
    action,optimal = a@pose,r@pose
    ref_risk = float(r@np.sum((pose-optimal)**2,axis=1))
    approx_action_ref_risk = float(r@np.sum((pose-action)**2,axis=1))
    occupancy_error = np.abs(a@grid['support']-r@grid['support'])
    return dict(label='FINITE_GRID_DISCREPANCY_NO_GLOBAL_CERTIFICATE',coordinate_cdf_sup=np.asarray(cdf_sup),
                event_px_ge_0_error=float(abs(a[grid['pose'][:,0]>=0].sum()-r[grid['pose'][:,0]>=0].sum())),
                mode_mass_max_error_ref_ge_005=float(np.abs(am[selected]-rm[selected]).max()) if np.any(selected) else 0.,
                evaluated_mode_ids=np.flatnonzero(selected),occupancy_error=occupancy_error,
                occupancy_max_error=float(occupancy_error.max()),action=action,reference_action=optimal,
                loss='target pose Euclidean squared loss: (px-action_x)^2 + (py-action_y)^2',
                reference_bayes_risk=ref_risk,approx_action_reference_risk=approx_action_ref_risk,
                action_regret=max(0.,approx_action_ref_risk-ref_risk),
                approximate_bayes_risk=float(a@np.sum((pose-action)**2,axis=1)))


def importance_sample(grid,log_likelihood,n_particles=1000,seed=0):
    """Self-normalized prior-proposal importance sampling on the finite prior."""
    p = _normalize(grid['prior_weight'])
    ll = _check_sampling_inputs(p,log_likelihood,n_particles)
    rng = np.random.default_rng(seed)
    indices = rng.choice(len(p),size=n_particles,p=p)
    weights = _lognormalize(ll[indices])
    return dict(label='FINITE_PRIOR_IMPORTANCE_MONTE_CARLO',indices=indices,particle_weights=weights,
                weights=np.bincount(indices,weights=weights,minlength=len(p)),
                effective_sample_size=float(1/np.sum(weights**2)),seed=seed,
                global_posterior_certificate=False)


def _check_sampling_inputs(p,log_likelihood,n_particles):
    ll = np.asarray(log_likelihood,dtype=float)
    if ll.shape!=p.shape or np.any(~np.isfinite(ll)):
        raise ValueError('sampling requires matching finite loglikelihood values')
    if not isinstance(n_particles,int) or isinstance(n_particles,bool) or n_particles<=0:
        raise ValueError('positive integer particle count required')
    return ll


def tempered_smc(grid,log_likelihood,n_particles=1000,temperatures=(0.,.25,.5,.75,1.),seed=0,mh_steps=2,resample_ess_fraction=.5):
    """Tempered finite-prior SMC with invariant independence-MH moves.

    Proposal is the exact finite prior. At inverse temperature beta, prior terms
    cancel in MH and acceptance is min(1, exp(beta*(LL_new-LL_old))). Resampling
    only replaces the particle representation and adds no observation factors.
    """
    p = _normalize(grid['prior_weight'])
    ll = _check_sampling_inputs(p,log_likelihood,n_particles)
    beta = np.asarray(temperatures,dtype=float)
    if beta.ndim!=1 or len(beta)<2 or beta[0]!=0 or beta[-1]!=1 or np.any(~np.isfinite(beta)) or np.any(np.diff(beta)<=0):
        raise ValueError('temperatures must strictly increase from 0 to 1')
    if not isinstance(mh_steps,int) or isinstance(mh_steps,bool) or mh_steps<0 or not 0<resample_ess_fraction<=1:
        raise ValueError('nonnegative integer MH steps and ESS fraction in (0,1] required')
    rng = np.random.default_rng(seed)
    indices = rng.choice(len(p),size=n_particles,p=p)
    weights = np.full(n_particles,1/n_particles)
    history=[]
    for previous,current in zip(beta[:-1],beta[1:]):
        with np.errstate(divide='ignore'):
            weights = _lognormalize(np.log(weights)+(current-previous)*ll[indices])
        ess = float(1/np.sum(weights**2))
        resampled = ess < n_particles*resample_ess_fraction
        if resampled:
            indices = indices[rng.choice(n_particles,size=n_particles,p=weights)]
            weights.fill(1/n_particles)
        accepted = 0
        for _ in range(mh_steps):
            proposal = rng.choice(len(p),size=n_particles,p=p)
            log_ratio = current*(ll[proposal]-ll[indices])
            accept = np.log(rng.random(n_particles)) < np.minimum(0,log_ratio)
            indices[accept] = proposal[accept]
            accepted += int(accept.sum())
        history.append(dict(beta=float(current),effective_sample_size=ess,resampled=resampled,
                            mh_accepted=accepted,mh_proposals=mh_steps*n_particles))
    return dict(label='FINITE_PRIOR_TEMPERED_SMC_MONTE_CARLO',indices=indices,particle_weights=weights,
                weights=np.bincount(indices,weights=weights,minlength=len(p)),history=history,
                seed=seed,global_posterior_certificate=False,
                observation_likelihood_insertions=1,final_temperature=1.)
