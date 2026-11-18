"""Good deflated Gauss–Newton (Bloor Riley, Webb & Baker 2025, alg. 3.2), adapted to a box."""
from argparse import ArgumentParser
import csv
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize

from high_endpoints import R12, bounds, historical_path, load_endpoint, noise_radius, physics
from high_score import material_scores, truth_for

CASES = ('ring', 'shell', 'disconnected')
MAX_FORWARD, MAX_RHS = 150, 4800
DEFLATION_EPSILON = .01
INNER_DAMPING = 1e-10


class BudgetReached(Exception):
    pass


def norm(x):
    return float(np.linalg.norm(x))


def rms(x):
    return norm(x)/np.sqrt(x.size)


class PhysicalLeastSquares:
    def __init__(self, public, max_forward, max_rhs):
        self.public = public
        self.model, self.reg = physics(public)
        self.target = public['Y']
        self.data_norm = norm(self.target)
        self.difference = self.reg.D.toarray()
        self.low, self.high = bounds(public)
        self.smooth = np.sqrt(2*self.reg.smooth/self.reg.edges)
        self.tv = np.sqrt(2*self.reg.tv/self.reg.edges)
        self.delta = noise_radius(public)
        self.max_forward, self.max_rhs = max_forward, max_rhs

    def guard(self, fw, rhs):
        counts = self.model.counts
        if counts['FULL_FORWARD']+fw > self.max_forward or counts['actual_complex_RHS']+rhs > self.max_rhs:
            raise BudgetReached

    def state(self, m):
        if self.model.state is None or not np.array_equal(self.model.state['m'],m):
            self.guard(1,self.model.ns)
            self.model.forward(m)
        return self.model.state

    def evaluate(self, m, jacobian=False):
        state = self.state(m)
        d = self.reg.D@m
        soft_norm = np.sqrt(d*d+self.reg.epsilon**2)
        data = self.model.pack(state['y']-self.target)
        residual = np.r_[data/self.data_norm, self.smooth*d,
                         self.tv*d/np.sqrt(soft_norm+self.reg.epsilon)]
        values = dict(phi=.5*float(residual@residual), residual_over_noise=norm(data)/self.delta)
        if not jacobian:
            return values
        if 'J' not in state:
            self.guard(0,self.model.nr)
        derivative = self.model.jacobian(state)/self.data_norm
        tv_slope = .5*self.tv*np.sqrt(soft_norm+self.reg.epsilon)/soft_norm
        jac = np.vstack((derivative, self.smooth*self.difference, tv_slope[:,None]*self.difference))
        return values, residual, jac


def step_subproblem(jac, residual, lo, hi):
    """Bound-constrained damped linear least squares via its convex quadratic."""
    H = jac.T@jac
    H.flat[::len(lo)+1] += INNER_DAMPING
    g = jac.T@residual
    diagonal = np.sqrt(H.diagonal())
    p0 = np.clip(cho_solve(cho_factor(H,check_finite=False), -g,check_finite=False), lo, hi)
    projected = np.clip(-g/diagonal**2,lo,hi)
    curvature = float(projected @ H @ projected)
    cauchy = projected * min(1.,-float(g@projected)/curvature) if curvature>0 else projected
    p0 = min((p0,cauchy,np.zeros_like(p0)),key=lambda p:float((g+.5*(H@p))@p))
    scale = 1e6

    def surrogate(z):
        p = z/diagonal
        hp = H@p
        return scale*float((g+.5*hp)@p), scale*(g+hp)/diagonal

    solution = minimize(surrogate, diagonal*p0, method='L-BFGS-B', jac=True,
                        bounds=list(zip(diagonal*lo,diagonal*hi)),
                        options=dict(maxiter=240,maxls=30,maxcor=30,ftol=2e-14,gtol=5e-6))
    p = min((p0,solution.x/diagonal),key=lambda p:float((g+.5*(H@p))@p))
    true_gradient = g+H@p
    violation = np.max(abs(p-np.clip(p-true_gradient,lo,hi)))
    return p, g, dict(inner_evaluations=int(solution.nfev), inner_success=bool(solution.success),
                      inner_projected_gradient=float(violation), inner_message=str(solution.message))


def log_mu_dot(m, center, step):
    # Shifted inverse-square Euclidean deflation in normalized m/sqrt(N) units.
    d = m-center
    d2 = float(d@d)
    if d2 == 0:
        raise ZeroDivisionError('Cannot initialize exactly at the deflated endpoint')
    return -2*float(d@step)/(d2*(1+d2/m.size))


def trial_run(public, material, center, mode, max_forward, max_rhs):
    t0,c0 = perf_counter(), process_time()
    fit = PhysicalLeastSquares(public,max_forward,max_rhs)
    m = material.copy()
    best = (float('inf'),m.copy(),None)
    values = fit.evaluate(m)
    raw = []
    accepted_steps = 0
    visited = [m.copy()]
    outcome = 'evaluation ceiling'

    def log(x, info):
        nonlocal best
        entry = dict(phi=info['phi'],residual_over_noise=info['residual_over_noise'],
                     distance_to_deflated_rms=rms(x-center))
        if info['phi'] < best[0]:
            best = (info['phi'],x.copy(),entry)
        return entry

    initial = log(m,values)
    while True:
        try:
            values, residual, jac = fit.evaluate(m,jacobian=True)
            step, gradient, inner = step_subproblem(jac,residual,fit.low-m,fit.high-m)
            row = dict(iteration=len(raw),base=log(m,values),GN_step_rms=rms(step),**inner)
            raw.append(row)
            if rms(step) < 1e-9:
                outcome = 'bounded GN step tolerance'
                break
            slope = log_mu_dot(m,center,step) if mode=='good_deflated_gn' else 0.
            beta = 1.-slope
            row.update(log_mu_directional=slope,beta=beta,trials=[])
            if mode=='good_deflated_gn' and slope > DEFLATION_EPSILON and abs(beta) >= 1e-8:
                # Paper zone 2/3: original scalar correction, no merit line search.
                # Project the resulting full trial into this application's box.
                step_scaled = step/beta
                x = np.clip(m+step_scaled,fit.low,fit.high)
                row.update(update='deflated',unprojected_step_rms=rms(step_scaled),
                           projection_rms=rms(x-(m+step_scaled)))
                if rms(x-m) < 1e-12:
                    outcome = 'deflated step clipped to current box point'
                    break
                values = fit.evaluate(x)
                row['trials'].append(log(x,values))
                m = x
            else:
                row['update'] = 'ordinary_line_search'
                derivative = float(gradient@step)
                if derivative >= 0:
                    outcome = 'bounded GN non-descent at numerical precision'
                    break
                for j in range(12):
                    a = 2.**(-j)
                    x = np.clip(m+a*step,fit.low,fit.high)
                    candidate = fit.evaluate(x)
                    row['trials'].append(log(x,candidate))
                    if candidate['phi'] <= values['phi']+1e-4*a*derivative:
                        m,values = x,candidate
                        break
                else:
                    outcome = 'original-objective line search exhausted'
                    break
            accepted_steps += 1
            visited.append(m.copy())
            if mode=='good_deflated_gn' and np.array_equal(m,center):
                outcome = 'box projection returned exactly to deflated point'
                break
        except BudgetReached:
            break
    work = dict(wall_seconds=perf_counter()-t0,cpu_seconds=process_time()-c0,counts=dict(fit.model.counts))
    return best[1], m, np.stack(visited), dict(start=initial,best=best[2],last=log(m,values),
              status=outcome,accepted_steps=accepted_steps,trace=raw,cost=work)


def case_run(label):
    spec, public, start = load_endpoint(label)
    lo,hi = bounds(public)
    center = np.clip(start,lo,hi)
    path = R12/'results/high_fold'/f'{label}.npz'
    with np.load(path,allow_pickle=False) as a:
        # Reuse ordinary multistart points exactly, selected without any material truth.
        seeds = [np.clip(a[f'multistart_seed{i}'],lo,hi) for i in range(3)]
    if any(np.array_equal(m,center) for m in seeds):
        raise ValueError('Shared start equals a deflated endpoint')
    t0,c0 = perf_counter(), process_time()
    incumbent_fit = PhysicalLeastSquares(public,1,public['E0'].shape[1])
    incumbent = incumbent_fit.evaluate(center)
    baseline_cost = dict(wall_seconds=perf_counter()-t0,cpu_seconds=process_time()-c0,counts=dict(incumbent_fit.model.counts))
    data = dict(deflated_endpoint=center,**{f'start{i}':m for i,m in enumerate(seeds)})
    arms = {}
    for method in ('ordinary_gn','good_deflated_gn'):
        candidates = []
        used_fw,used_rhs = 1,incumbent_fit.model.ns
        chosen,chosen_score = center,incumbent['phi']
        chosen_name = 'archived_incumbent'
        best_distinct = None
        for i,seed in enumerate(seeds):
            best,last,trajectory,detail = trial_run(public,seed,center,method,
              (MAX_FORWARD-used_fw)//(len(seeds)-i),(MAX_RHS-used_rhs)//(len(seeds)-i))
            counts = detail['cost']['counts']
            used_fw += counts['FULL_FORWARD']
            used_rhs += counts['actual_complex_RHS']
            candidates.append(detail)
            data[f'{method}_trial{i}_best'] = best
            data[f'{method}_trial{i}_last'] = last
            data[f'{method}_trial{i}_trajectory'] = trajectory
            if detail['best']['phi'] < chosen_score:
                chosen_name = f'trial{i}_best'
                chosen,chosen_score = best,detail['best']['phi']
            if detail['best']['distance_to_deflated_rms'] >= .05 and (
                best_distinct is None or detail['best']['phi'] < best_distinct['phi']):
                best_distinct = dict(trial=i,**detail['best'])
        data[method] = chosen
        arms[method] = dict(selection=chosen_name,selection_phi=chosen_score,
            selection_distance_to_archived_rms=rms(chosen-center),
            selected_residual_over_noise=incumbent['residual_over_noise'] if chosen_name=='archived_incumbent' else
                next(d['best']['residual_over_noise'] for i,d in enumerate(candidates) if chosen_name==f'trial{i}_best'),
            best_distinct_best_trial=best_distinct, trials=candidates,
            charged_full_forward=used_fw,charged_complex_RHS=used_rhs,
            wall_seconds=baseline_cost['wall_seconds']+sum(d['cost']['wall_seconds'] for d in candidates),
            cpu_seconds=baseline_cost['cpu_seconds']+sum(d['cost']['cpu_seconds'] for d in candidates),
            deflated_jumps=sum(t.get('update')=='deflated' and bool(t['trials']) for d in candidates for t in d['trace']),
            uphill_deflated_jumps=sum(t.get('update')=='deflated' and bool(t['trials']) and t['trials'][0]['phi']>t['base']['phi']
                                      for d in candidates for t in d['trace']))
    metadata = dict(case=label,case_id=spec.case_id,archived_method=spec.method,
                    source_observations=historical_path(spec.input_path),start_provenance=f'{path.name}:multistart_seed0..2',
                    archived_incumbent=incumbent,shared_incumbent_cost=baseline_cost,arms=arms)
    return data, metadata


def main(output):
    arrays, rows = {},[]
    for label in CASES:
        prediction, metadata = case_run(label)
        arrays.update({f'{label}_{key}':value for key,value in prediction.items()})
        rows.append(metadata)
        print(label,[(name,arm['selected_residual_over_noise'],arm['deflated_jumps'],arm['charged_full_forward'],arm['charged_complex_RHS'])
                     for name,arm in metadata['arms'].items()],flush=True)
    settings = dict(paper='10.1137/24M1695191, algorithm 3.2',
        reference='https://github.com/AlbanBloorRiley/DeflatedGaussNewton/blob/615c185f3d61f044eab7cfd27f97f698ff341d8f/Methods/Good_Deflated_GaussNewton.m',
        starts='exact ordinary multi-start seeds from high_fold',
        box='same archived physical box; solve ordinary GN linear model with box step bound; componentwise clip p/beta for legal deflated trials',
        theta=2.,shift=1.,epsilon=DEFLATION_EPSILON,deflation_distance='Euclidean in m/sqrt(N) coordinates',
        regularization='original quadratic + smooth TV carried as exact least-squares residual and Jacobian',
        inner_solver='box-constrained full J linearized least squares with 1e-10 damping; diagonally scaled convex quadratic L-BFGS-B with a feasible Cauchy initialization candidate',
        merit='original Phi; backtracking Armijo for ordinary steps; deflated steps taken without a merit check',
        budget_full_forward=MAX_FORWARD,budget_complex_RHS=MAX_RHS,
        incumbent='archived bad endpoint available to both arms; initial physical evaluation charged to both')
    # Persist every new material before the separate scoring pass opens GT.
    np.savez_compressed(output.with_suffix('.npz'),**arrays)
    output.with_suffix('.json').write_text(json.dumps(dict(config=settings,cases=rows),indent=2,allow_nan=False)+'\n')
    scores=[]
    with np.load(output.with_suffix('.npz'),allow_pickle=False) as predictions:
        for row in rows:
            label = row['case']
            truth,support,amplitude,source = truth_for(row['case_id'])
            prefix = label+'_'
            for method,arm in row['arms'].items():
                material = predictions[prefix+method]
                scores.append(dict(case=label,method=method,choice=arm['selection'],
                  residual_over_noise=arm['selected_residual_over_noise'],phi=arm['selection_phi'],
                  **material_scores(material,truth,support,amplitude),truth_source=source,
                  physical_forward=arm['charged_full_forward'],complex_RHS=arm['charged_complex_RHS'],
                  cpu_seconds=arm['cpu_seconds'],wall_seconds=arm['wall_seconds'],
                  deflated_jumps=arm['deflated_jumps'],uphill_deflated_jumps=arm['uphill_deflated_jumps']))
    with output.with_suffix('.csv').open('w') as f:
        writer = csv.DictWriter(f,fieldnames=list(scores[0]))
        writer.writeheader()
        writer.writerows(scores)


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=R12/'results/high_deflated_gn')
    main(parser.parse_args().output)
