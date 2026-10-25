"""Measurement-only bounded joint MAP proposals and reversible support helpers.

No truth I/O, posterior-mass claim, or experiment driver. Born ranking is only
an uncharged-in-RHS analytic proposal calculation; full-wave acceptance is the
caller's responsibility. Algorithmic TRF damping is never a statistical prior.
"""
from functools import lru_cache
import time

import numpy as np
from scipy.optimize import least_squares, minimize
from scipy.special import hankel1

from .geometry import ROI, enumerate_supports, validate_geometry
from .world import World

LOWER = np.array([1.5, 0., -1., 1.6, -.01, -.01])
UPPER = np.array([6., .25, 1., 2.6, .01, .01])
X_SCALE = np.array([4.5, .25, 2., 1., .02, .02])
GAIN_PRIOR_SIGMA = .005
NOISE_SIGMA = .01
RANKING_LABEL = 'BOUNDED_JOINT_MAP_CANDIDATE_NOT_POSTERIOR_MASS'


def _pack(data):
    data = np.asarray(data, complex).ravel()
    return np.r_[data.real, data.imag] / NOISE_SIGMA


def _packj(data):
    data = np.asarray(data, complex)
    return np.vstack((data.real, data.imag)) / NOISE_SIGMA


def _data(y):
    y = np.asarray(y, complex)
    if y.size != 32 or not np.isfinite(y).all():
        raise ValueError('Expected finite base acquisition: 2 Tx by 16 Rx')
    return y.ravel()


def _support(support):
    s = np.asarray(support)
    if s.shape not in ((2, 2), (3, 3)) or not validate_geometry(s)['legal']:
        raise ValueError('Expected legal public 2x2 or 3x3 occupancy')
    return s.astype(bool)


def support_id(support):
    """Stable row-major binary ID within the explicitly named grid shape."""
    s = _support(support)
    return f'{s.shape[0]}x{s.shape[1]}:' + ''.join(str(int(x)) for x in s.ravel())


@lru_cache(maxsize=2)
def _support_masses(size):
    return {s.tobytes(): p for s, p in enumerate_supports(size, .65)}


def support_mass(support):
    s = _support(support)
    return float(_support_masses(s.shape[0])[s.tobytes()])


def fit_joint(y, support, theta0, book=None, *, evaluator=None, stage='A0', max_nfev=6):
    """Bounded TRF in [epsinf,ell,target_x,target_y,gain_offset,phase].

    The optional evaluator allows dependency injection; its book charges all
    actual solver actions. max_nfev is at most six; callers may allocate three
    per stage/start to compose an equal-budget strategy. Initial J counts as
    one paid Jacobian evaluation. Returned arrays include every residual trial.
    Uniform material coordinates remain latent even for an empty wall.
    """
    y, support = _data(y), _support(support)
    theta0 = np.asarray(theta0, float)
    if theta0.shape != (6,) or not np.isfinite(theta0).all() or np.any(theta0 < LOWER) or np.any(theta0 > UPPER):
        raise ValueError('Initial theta must already satisfy physical bounds')
    if not isinstance(max_nfev, int) or isinstance(max_nfev, bool) or not 1 <= max_nfev <= 6:
        raise ValueError('max_nfev must be an integer from one through six')
    if evaluator is None:
        if book is None:
            raise ValueError('A shared accounting book is required')
        from .forward import Evaluator
        evaluator = Evaluator(book, stage=stage)
    book = evaluator.book
    if getattr(evaluator, 'stage', stage) != stage:
        raise ValueError('Evaluator accounting stage mismatch')
    before = book.totals(stage)
    began, cpu_began = time.perf_counter(), time.process_time()
    cache, traces = {}, []
    physical_forward_calls = physical_jacobian_calls = 0
    evaluator_wall_seconds = evaluator_cpu_seconds = 0.

    def get(theta, want_j=False):
        nonlocal physical_forward_calls, physical_jacobian_calls
        nonlocal evaluator_wall_seconds, evaluator_cpu_seconds
        # Gain changes are algebraic: share physical data/J across them.
        key = tuple(float(x) for x in theta[:4])
        item = cache.get(key)
        needs_f = item is None
        if needs_f or want_j and item['J'] is None:
            world = World(support.copy(), theta[0], theta[1], theta[2:4].copy(),
                          radius=.12, contrast=3., target_loss=.02)
            action_wall, action_cpu = time.perf_counter(), time.process_time()
            ev = evaluator.evaluate(world, jacobian=want_j)
            evaluator_wall_seconds += time.perf_counter()-action_wall
            evaluator_cpu_seconds += time.process_time()-action_cpu
            f = _data(ev.data)
            if item is None:
                item = dict(F=f.copy(), J=None, metadata=dict(ev.metadata))
                cache[key] = item
            if want_j:
                physical_jacobian_calls += 1
                item['J'] = np.asarray(ev.jacobian, complex)[:, :4].copy()
                if item['J'].shape != (32, 4):
                    raise ValueError('Expected physical Jacobian with four leading columns')
                item['metadata'] = dict(ev.metadata)
            # This count describes evaluator calls; the Book is authoritative
            # if a revisited F cache needs an actual field rebuild for its J.
            physical_forward_calls += int(needs_f)
        return item

    def fun(theta):
        item = get(theta)
        gamma = (1+theta[4])*np.exp(1j*theta[5])
        raw = _pack(gamma*item['F']-y)
        residual = np.r_[raw, theta[4:6]/GAIN_PRIOR_SIGMA]
        traces.append(dict(kind='residual_trial', theta=theta.copy(),
                           raw_train_white=raw.copy(), objective=float(.5*residual@residual),
                           pred=(gamma*item['F']).copy(), physical_key=theta[:4].copy()))
        return residual

    def jac(theta):
        item = get(theta, True)
        phase = np.exp(1j*theta[5])
        gamma = (1+theta[4])*phase
        complex_j = np.column_stack((gamma*item['J'], phase*item['F'], 1j*gamma*item['F']))
        prior_j = np.zeros((2, 6))
        prior_j[:, 4:6] = np.eye(2)/GAIN_PRIOR_SIGMA
        result = np.vstack((_packj(complex_j), prior_j))
        traces.append(dict(kind='jacobian_paid_or_reused', theta=theta.copy(),
                           physical_jacobian=item['J'].copy(), residual_jacobian=result.copy(),
                           metadata=dict(item['metadata'])))
        return result

    solution = least_squares(fun, theta0, jac=jac, method='trf', bounds=(LOWER, UPPER),
                             x_scale=X_SCALE, max_nfev=max_nfev)
    # No clipping, no extra physics, and no additional residual call at exit.
    theta = solution.x.copy()
    item = cache[tuple(float(x) for x in theta[:4])]
    pred = (1+theta[4])*np.exp(1j*theta[5])*item['F']
    raw = _pack(pred-y)
    objective = float(.5*raw@raw + .5*np.sum((theta[4:6]/GAIN_PRIOR_SIGMA)**2))
    mass = support_mass(support)
    after = book.totals(stage)
    return dict(theta=theta, support=support.copy(), support_id=support_id(support),
                pred=pred, objective=objective, support_negative_log_prior=float(-np.log(mass)),
                ranking_objective=objective-float(np.log(mass)), support_prior_mass=mass,
                rawtrainwhite=raw, rawtrainwhite_norm=float(np.linalg.norm(raw)),
                success=bool(solution.success), status=int(solution.status), message=solution.message,
                nfev=int(solution.nfev), njev=int(solution.njev or 0), traces=traces,
                saved_evaluations=[dict(physical_theta=np.array(k), **v) for k,v in cache.items()],
                stage_cost_diffs={k:after[k]-before[k] for k in after},
                evaluator_forward_requests=physical_forward_calls,
                evaluator_jacobian_requests=physical_jacobian_calls,
                total_wall_seconds=time.perf_counter()-began,
                total_process_cpu_seconds=time.process_time()-cpu_began,
                evaluator_wall_seconds=evaluator_wall_seconds,
                evaluator_process_cpu_seconds=evaluator_cpu_seconds,
                optimizer_and_bookkeeping_cpu_seconds=max(0.,time.process_time()-cpu_began-evaluator_cpu_seconds),
                label=RANKING_LABEL, proper_posterior_claim=False,
                physical_cost_authority='SHARED_STAGE_BUDGET_BOOK',
                prior='Uniform bounded physical coordinates; truncated independent Gaussian gain/phase sigma .005')


def profile_gain_phase(y, prediction):
    """Cheap bounded 2D gain/phase MAP profiling with analytic derivatives."""
    y, f = _data(y), _data(prediction)
    def objective(z):
        phase = np.exp(1j*z[1])
        gamma = (1+z[0])*phase
        delta = gamma*f-y
        value = .5*np.vdot(delta, delta).real/NOISE_SIGMA**2 + .5*np.sum((z/GAIN_PRIOR_SIGMA)**2)
        j = np.column_stack((phase*f, 1j*gamma*f))
        gradient = np.real(j.conj().T@delta)/NOISE_SIGMA**2 + z/GAIN_PRIOR_SIGMA**2
        return float(value), gradient
    result = minimize(objective, np.zeros(2), jac=True, method='L-BFGS-B', bounds=[(-.01,.01)]*2)
    value, gradient = objective(result.x)
    return dict(gain_offset=float(result.x[0]), phase=float(result.x[1]), objective=value,
                success=bool(result.success), nfev=int(result.nfev), gradient=gradient,
                physics_rhs=0, proper_posterior_claim=False)


def bank_rank(y, bankdata, bankgrid):
    """Profile every row; return globally sorted rows and best per support.

    All legal supports must occur in the bank. No Top-K support restriction is
    returned or imposed; best_by_support retains every legal support (14 at 2x2).
    Finite bank quadrature weights are not posterior masses or continuous MAP
    prior densities; ranking uses the conditional Bernoulli support mass only.
    """
    y = _data(y)
    bankdata = np.asarray(bankdata, complex)
    n = len(bankgrid['epsinf'])
    if n == 0:
        raise ValueError('Nonempty proposal bank required')
    if bankdata.shape[0] != n or bankdata.size != n*32:
        raise ValueError('Bank data/grid shape mismatch')
    bankdata = bankdata.reshape(n,32)
    rows = []
    for index in range(n):
        raw_support = np.asarray(bankgrid['support'][index])
        size = int(round(np.sqrt(raw_support.size)))
        s = _support(raw_support.reshape(size,size))
        profiled = profile_gain_phase(y, bankdata[index])
        theta = np.r_[bankgrid['epsinf'][index],bankgrid['ell'][index],bankgrid['pose'][index],
                      profiled['gain_offset'],profiled['phase']]
        mass = support_mass(s)
        rows.append(dict(row_id=index, support_id=support_id(s), support=s.copy(), theta=theta,
                         objective=profiled['objective'], support_prior_mass=mass,
                         ranking_objective=profiled['objective']-np.log(mass), profile=profiled,
                         label=RANKING_LABEL, proper_posterior_claim=False))
    rows.sort(key=lambda row:(row['ranking_objective'],row['support_id'],row['row_id']))
    best = {}
    for row in rows:
        best.setdefault(row['support_id'],row)
    expected = {support_id(s) for s,_ in enumerate_supports(size,.65)}
    if set(best) != expected:
        raise ValueError('Bank must retain all legal supports, including empty wall')
    return dict(rows=rows, best_by_support=[best[k] for k in sorted(best)],
                best=rows[0], all_supports_retained=True, label=RANKING_LABEL,
                proper_posterior_claim=False, physics_rhs=0)


def coarse_support_to3x3(support):
    """Explicit center-sampled representation change, not prior preservation.

    The returned mapping may be illegal; callers must validate before physics.
    """
    s = _support(support)
    if s.shape != (2,2):
        raise ValueError('Expected coarse 2x2 occupancy')
    index = np.floor((np.arange(3)+.5)*2/3).astype(int)
    return dict(support=s[np.ix_(index,index)].copy(),
                representation_change='2x2_TO_3x3_REGION_CENTER', prior_preserved=False)


def legal_support_neighbors(support, *, include_noop=True):
    """All legal single births/deaths and occupied-air exchanges, stable IDs."""
    s = _support(support)
    flat = s.ravel()
    candidates = {support_id(s):s.copy()} if include_noop else {}
    edits = [(i,) for i in range(flat.size)]
    edits += [(i,j) for i in np.flatnonzero(flat) for j in np.flatnonzero(~flat)]
    for indices in edits:
        candidate = flat.copy()
        for i in indices:
            candidate[i] = ~candidate[i]
        candidate = candidate.reshape(s.shape)
        if validate_geometry(candidate)['legal']:
            candidates[support_id(candidate)] = candidate
    return [dict(support_id=k, support=candidates[k], is_noop=np.array_equal(candidates[k],s))
            for k in sorted(candidates)]


def born_region_screen(y, theta, supports, *, quadrature_degree=4):
    """Free-space first Born proposal ranking over public 3x3 rectangles.

    Analytic unit line sources and target disk area quadrature. No numerical
    Maxwell solve; not a truth model, acceptance test, or qualified surrogate.
    Uses only supplied current material/target hypotheses and observed y.
    """
    y = _data(y)
    theta = np.asarray(theta,float)
    if theta.shape != (6,) or not np.isfinite(theta).all() or np.any(theta<LOWER) or np.any(theta>UPPER):
        raise ValueError('Expected current bounded six-parameter hypothesis')
    if not isinstance(quadrature_degree,int) or not 2 <= quadrature_degree <= 12:
        raise ValueError('Use an explicit quadrature degree from 2 through 12')
    tx = np.array([[-1.5,0.],[1.5,0.]])
    rx = np.column_stack((np.linspace(-2,2,16),np.zeros(16)))
    k = 2*np.pi
    def green(points,sources):
        distance = np.linalg.norm(points[:,None,:]-sources[None,:,:],axis=-1)
        if np.any(distance==0):
            raise ValueError('Singular free-space quadrature geometry')
        return .25j*hankel1(0,k*distance)
    def integral(points,weights):
        return k*k*np.einsum('qt,qr,q->tr',green(points,tx),green(points,rx),weights)
    nodes,weights = np.polynomial.legendre.leggauss(quadrature_degree)
    xmin,xmax,ymin,ymax = ROI
    blocks = []
    for row in range(3):
        for col in range(3):
            dx,dy = (xmax-xmin)/3,(ymax-ymin)/3
            x = xmin+(col+.5)*dx+nodes*dx/2
            yy = ymin+(row+.5)*dy+nodes*dy/2
            xx,zz = np.meshgrid(x,yy)
            points = np.column_stack((xx.ravel(),zz.ravel()))
            blocks.append(integral(points,np.outer(weights,weights).ravel()*dx*dy/4))
    nr = quadrature_degree
    radius = (nodes+1)*.12/2
    angle = np.arange(4*nr)*2*np.pi/(4*nr)
    xx,aa = np.meshgrid(radius,angle,indexing='ij')
    points = theta[2:4]+np.column_stack((xx.ravel()*np.cos(aa.ravel()),xx.ravel()*np.sin(aa.ravel())))
    disk_weights = np.repeat(weights*.12/2*radius,4*nr)*2*np.pi/(4*nr)
    base = green(rx,tx).T+(3.+.02j)*integral(points,disk_weights)
    rows = []
    for support in supports:
        s = _support(support)
        if s.shape != (3,3):
            raise ValueError('Born screen uses public 3x3 regions')
        f = base+(theta[0]-1+1j*theta[1])*np.tensordot(s.ravel().astype(float),np.array(blocks),axes=1)
        pred = (1+theta[4])*np.exp(1j*theta[5])*f.ravel()
        residual = _pack(pred-y)
        obj = .5*residual@residual+.5*np.sum((theta[4:6]/GAIN_PRIOR_SIGMA)**2)-np.log(support_mass(s))
        rows.append(dict(support=s.copy(),support_id=support_id(s),ranking_objective=float(obj),
                         prediction=pred,proposal_only=True,physics_rhs=0,
                         label='ANALYTIC_FREE_SPACE_BORN_SCREEN_NOT_MAXWELL_TRUTH'))
    return sorted(rows,key=lambda row:(row['ranking_objective'],row['support_id']))
