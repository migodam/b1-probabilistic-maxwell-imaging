"""Orthogonally split frozen compact illuminations between two physical bounds."""
import json

from flint import acb, arb, ctx
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

from directional import compare_bases
from finite_states import pure_shift_score, spectral_enclosure
from pure_shift import center_ball
from run import (HERE, adj, cap_coefficients, clock, elapsed, finite_geometry, herm,
                 relax_receiver, tail, write_json)
from verify import (FiniteExact, ball, blocks, dual_evaluate, exact_number, eye,
                    gram_interval, hi, lo, norm, num, primal_evaluate, star)


TANGENT = .2
CAPS = (1., 1.25, 1.5)
IDS = ('loss_tm_01__fixed_t0.2__pair2pct', 'loss_tm_02__fixed_t0.2__pair2pct')
BASES = ('right_singular', 'best_common_shift')


def numerical_state(geo, y, cap, profile):
    mu, optical = profile['mu'], profile['weight']
    beta, alpha = cap_coefficients(cap,TANGENT)
    k = adj(geo.s)@geo.s/mu-geo.hk+optical*geo.dk+(beta+optical*alpha)*np.eye(geo.n)
    l = adj(geo.s)@y/mu+(1.+1.j*optical)/2*geo.b
    kval, vec = eigh(herm(k))
    nu = int((kval<0).sum())
    if min(abs(kval)) < 1.e-11:
        raise ArithmeticError('Numerical state spectral gap is too small.')
    j = vec @ ((adj(vec)@l)/kval[:,None])
    return herm(adj(y)@y-mu*adj(l)@j),nu


def verify_state(geo, exact, y, cap, profile):
    start = clock()
    beta = 1/((1+exact_number(TANGENT)**2)*exact_number(cap))
    alpha = exact_number(TANGENT)*beta
    mu,optical = arb(float(profile['mu'])),arb(float(profile['weight']))
    if not (mu>0 and optical>=0):
        raise ArithmeticError('Whole-shift weights must be physically eligible.')
    y = ball(y)
    k = (star(exact.s)*exact.s/mu - (exact.g+star(exact.g))/2
         + optical*exact.r + (beta+optical*alpha)*eye(geo.n))
    l = star(exact.s)*y/mu + acb(1,optical)*exact.b/2
    lows, highs, state_error = spectral_enclosure(k)
    if any(not (a>0 or b<0) for a,b in zip(lows,highs,strict=True)):
        raise ArithmeticError('Exact state inertia unresolved.')
    nu = sum(b<0 for b in highs)
    gap = min((-b if b<0 else a).lower() for a,b in zip(lows,highs,strict=True))
    vals, vec = eigh(herm(num(k)))
    trial = ball(vec@((adj(vec)@num(l))/vals[:,None]))
    kj = k*trial
    error = norm(kj-l)
    h = star(y)*y+mu*(star(trial)*kj-star(trial)*l-star(l)*trial)
    displacement = (mu*error**2/gap).upper()
    all_sector = tail_interval(h,nu,displacement)
    return h,nu,displacement,dict(inertia=nu,gap_lower=lo(gap),state_spectral_error_upper=state_error,
        stationarity_error_upper=hi(error),source_error_upper=hi(displacement),
        total_sector_distance_interval=distance_interval(*all_sector),**elapsed(start))


def tail_interval(h, nu, error):
    retained = h.nrows()-nu
    if retained<=0:
        return arb(0),arb(0)
    lower, upper, _ = spectral_enclosure(h)
    return (sum((max(arb(0),(x-error).lower()) for x in lower[:retained]),arb(0)).lower(),
            sum((max(arb(0),(x+error).upper()) for x in upper[:retained]),arb(0)).upper())


def distance_interval(lower, upper):
    return [max(0.,lo(lower.sqrt())),hi(upper.sqrt())]


def mask_indices(mask, n):
    rec = tuple(i for i in range(n) if mask & (1<<i))
    state = tuple(i for i in range(n) if not mask & (1<<i))
    return rec,state


def screen(y, rs, cb, fc, basis, prices, h, nu):
    start = clock()
    d = (y-fc)@basis
    source_quadratic = np.diag(adj(basis)@cb@basis).real
    energy = np.array([
        tau*(np.vdot(d[:,i],np.linalg.solve(rs+tau*np.eye(rs.shape[0]),d[:,i])).real-source_quadratic[i])
        if tau>0 else 0. for i,tau in enumerate(prices)])
    energy = np.maximum(energy,0.)
    h_rotated = herm(adj(basis)@h@basis)
    values = []
    for mask in range(1,2**len(prices)-1):
        receiver_ids, state_ids = mask_indices(mask,len(prices))
        sector = tail(h_rotated[np.ix_(state_ids,state_ids)],nu)**2
        receiver = float(energy[list(receiver_ids)].sum())
        values.append(dict(mask=mask,receiver_columns=list(receiver_ids),sector_columns=list(state_ids),
                           receiver_squared=receiver,sector_squared=sector,distance=float(np.sqrt(receiver+sector))))
    best_gen = max(values,key=lambda r:r['distance'])
    both = [v for v in values if min(v['receiver_squared'],v['sector_squared'])>1.e-20]
    return dict(max_nontrivial=best_gen,max_two_positive=max(both,key=lambda r:r['distance']) if both else None,
                candidate_masks=len(values),positive_receiver_prices=int((energy>0).sum()),
                all_directional_receiver_numeric=float(np.sqrt(energy.sum())),**elapsed(start))


def interval_receivers(geo, y, basis, prices, trial, error, polar_e):
    """Bound the attained fixed-price expression, with both sides of the exact Gram."""
    rslo,cblo,fc = blocks(geo,trial,arb(0))
    rshi,cbhi,_ = blocks(geo,trial,error)
    D = ball(y)-center_ball(fc,error)
    trial_receiver_floor = spectral_enclosure(rslo)[0][0].lower()
    positive_prices = [float(t) for t in prices if t>0]
    resolvent_floor = trial_receiver_floor+arb(min(positive_prices)) if positive_prices else None
    if resolvent_floor is not None and not resolvent_floor>0:
        raise ArithmeticError('Lower receiver Gram plus selected price is not verified positive.')
    polar_distance = polar_e/(1+(1-polar_e).sqrt())
    radius_factor = (1+(1+polar_e).sqrt())*polar_distance
    bounds = []
    for i,tau in enumerate(prices):
        if tau<=0:
            bounds.append((arb(0),arb(0)))
            continue
        t,v = arb(float(tau)),ball(basis[:,[i]])
        dv = D*v
        negative = (star(dv)*(rshi+t*eye(geo.p)).solve(dv) - star(v)*cbhi*v)[0,0].real*t
        positive = (star(dv)*(rslo+t*eye(geo.p)).solve(dv) - star(v)*cblo*v)[0,0].real*t
        drift = radius_factor*(norm(D)**2+t*norm(cbhi))
        bounds.append((max(arb(0),(negative-drift).lower()),max(arb(0),(positive+drift).upper())))
    return bounds,polar_distance,resolvent_floor


def verify_split(geo,y,h,nu,source_error,trial,error,basis,prices,proposals):
    start = clock()
    vb = ball(basis)
    polar_e = norm(star(vb)*vb-eye(geo.inputs))
    if not polar_e<1:
        raise ArithmeticError('Polar unitary reconstruction not enclosed.')
    prices_bound,polar_drift,resolvent_floor = interval_receivers(geo,y,basis,prices,trial,error,polar_e)
    s_error = (source_error + (1+(1+polar_e).sqrt())*polar_drift*norm(h)).upper()
    results = {}
    for name in ('max_nontrivial','max_two_positive'):
        proposed = proposals[name]
        if proposed is None:
            continue
        inds,sids = proposed['receiver_columns'],proposed['sector_columns']
        q = ball(basis[:,sids])
        slo,shi = tail_interval(star(q)*h*q,nu,s_error)
        plo = sum((prices_bound[i][0] for i in inds),arb(0)).lower()
        phi = sum((prices_bound[i][1] for i in inds),arb(0)).upper()
        results[name] = dict(distance_interval=distance_interval(plo+slo,phi+shi),
            receiver_squared_interval=[max(0.,lo(plo)),hi(phi)],
            sector_squared_interval=[max(0.,lo(slo)),hi(shi)])
    return dict(polar_defect_upper=hi(polar_e),sector_polar_and_stationarity_error_upper=hi(s_error),
                shifted_lower_receiver_gram_spectral_floor=None if resolvent_floor is None else lo(resolvent_floor),
                selections=results,**elapsed(start))


def polish_sector(geo,y,cap,profile):
    start = clock()
    initial = np.array([np.log(profile['mu']),np.log1p(profile['weight'])])

    def objective(q):
        settings = dict(mu=float(np.exp(q[0])),weight=float(np.expm1(q[1])))
        try:
            h,nu = numerical_state(geo,y,cap,settings)
        except ArithmeticError:
            return 0.
        return -tail(h,nu)

    res = minimize(objective,initial,method='Nelder-Mead',bounds=[(-20,1),(0,8)],
                   options=dict(maxiter=110,maxfev=180,xatol=5.e-5,fatol=1.e-10))
    return dict(mu=float(np.exp(res.x[0])),weight=float(np.expm1(res.x[1])),
                objective_calls=int(res.nfev),numerical_sector_distance=-float(res.fun),
                solver_message=str(res.message),**elapsed(start))


def comparison(geo,exact,y,cap,profile,trial,error,rs,cb,fc,bases,all_receivers,radius,arrays,key):
    h,nu,stationary, state = verify_state(geo,exact,y,cap,profile)
    hnum,_ = numerical_state(geo,y,cap,profile)
    variants = {}
    for name in BASES:
        v,prices,details = bases[name]
        selected = screen(y,rs,cb,fc,v,prices,hnum,nu)
        checked = verify_split(geo,y,h,nu,stationary,trial,error,v,prices,selected)
        arrays[key+'__'+name+'_basis'] = v
        arrays[key+'__'+name+'_prices'] = prices
        variants[name] = dict(proposals=selected,outward=checked,price_cost=details)
    both = [ (v['outward']['selections']['max_two_positive']['distance_interval'][0],name)
             for name,v in variants.items() if 'max_two_positive' in v['outward']['selections'] ]
    if not both:
        return dict(profile=profile,all_sector=state,bases=variants,status='no_two_positive_terms')
    found, name = max(both)
    state_lo,state_hi = state['total_sector_distance_interval']
    recv_lo,recv_hi = all_receivers
    proof = variants[name]['outward']['selections']['max_two_positive']
    return dict(profile=profile,all_sector=state,bases=variants,selected=name,
                numerical_selection=variants[name]['proposals']['max_two_positive'],outward_selection=proof,
                strict_gain_over_matched_all_sector=(found>state_hi),
                strict_gain_over_full_receiver_ellipsoid=(found>recv_hi),
                new_strict_exclusion=(found>radius and max(state_hi,recv_hi)<=radius),
                fixed_sector_lower_over_radius=state_lo/radius,
                split_lower_over_radius=found/radius)


def main():
    start = clock()
    geo = finite_geometry(.6)
    exact = FiniteExact(geo,.6)
    original = json.loads((HERE/'finite_0_state_verified.json').read_text())['rows']
    runs = []
    gram_setup = []
    arrays = {}
    with np.load(HERE/'finite_0_inputs.npz',allow_pickle=False) as data:
        for cap in CAPS:
            trial,error,_,setup = gram_interval(geo,exact,exact_number(cap),exact_number(TANGENT))
            gram_setup.append(dict(cap=cap,**setup))
            rs,cb,fc = map(num,blocks(geo,trial,arb(0)))
            for ident in IDS:
                y = data[ident]
                row = next(r for r in original if r['id']==ident and r['cap']==cap and r['tangent']==TANGENT)
                r = row['rho']
                common = pure_shift_score(geo,y,trial)
                bases,base_cost = compare_bases(y,rs,cb,fc,common['tau'])
                numeric,shadow = relax_receiver(y,rs,cb,fc)
                dual_proof = dual_evaluate(geo,y,shadow,trial,error)
                primal_proof = primal_evaluate(geo,y,trial,error)
                recv = [dual_proof['dual_lower'],primal_proof['primal_upper']]
                key=ident+'__u'+str(cap)
                frozen = comparison(geo,exact,y,cap,row['state']['shift']['profile'],trial,error,
                                    rs,cb,fc,bases,recv,r,arrays,key+'__frozen')
                result=dict(id=ident,cap=cap,tangent=TANGENT,radius=r,old_passive_sector_distance_lower=row['state']['old']['distance_lower'],
                            frozen_selected=frozen,all_receiver_interval=recv,
                            all_receiver_proof=dict(numerical=numeric,dual=dual_proof,primal=primal_proof),
                            independently_priced_basis_setup=base_cost,common_shift_basis_seed=common)
                if ident==IDS[0] and cap in (1.,1.25):
                    polish = polish_sector(geo,y,cap,row['state']['shift']['profile'])
                    result['sector_polish_sensitivity'] = comparison(geo,exact,y,cap,polish,trial,error,
                                    rs,cb,fc,bases,recv,r,arrays,key+'__sector_polish')
                runs.append(result)
                print(ident,cap,'receiver',*[v/r for v in recv],
                      'frozen sector',frozen['fixed_sector_lower_over_radius'],
                      'split',frozen['split_lower_over_radius'],
                      'new exclusion',frozen['new_strict_exclusion'],flush=True)
                if 'sector_polish_sensitivity' in result:
                    sensitivity = result['sector_polish_sensitivity']
                    print('  polished sector',sensitivity['fixed_sector_lower_over_radius'],
                          'split',sensitivity['split_lower_over_radius'],flush=True)
    np.savez_compressed(HERE/'source_partition.npz',**arrays)
    write_json(HERE/'source_partition.json',dict(scope='12² full TM disk collocation, 0.6 m, 400 MHz, known-loss-angle fixed observation, 10 illuminations with fixed fluence',
              receipts='finite_0_inputs.npz and finite_0_state_verified.json, no regenerated observation or material; exact represented observations, analytic operators',
              selection='Data-only right singular basis or common-scalar Hermitian eigenbasis. All nonempty proper subsets as receiver, complement as sector. Best each basis by squared distance; also require both contributions positive for the primary split.',
              strict_claim_scope='Matched attained endpoints, including full receiver ellipsoid upper and fixed selected all-sector quadratic. Numeric search across subsets and optional below-cap sensitivity; no maximum over every sector weight or source basis.',
              precision_bits=ctx.prec,setup_by_cap=gram_setup,runs=runs,**elapsed(start)))


if __name__=='__main__':
    ctx.prec,ctx.threads=128,1
    with threadpool_limits(1):
        main()
