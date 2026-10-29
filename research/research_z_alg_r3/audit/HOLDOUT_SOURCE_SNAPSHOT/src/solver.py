"""All-unknown R3 solvers. Input is restricted to public operators and data.

No private target, material truth, generator secret, or grading module is read.
Quotient roots are candidates on exact coordinate slices; no full-dimensional
root completeness or global optimality is asserted.
"""
import itertools
import json
import time
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, minimize
from scipy.linalg import null_space, lstsq
from scipy.stats import qmc
from threadpoolctl import threadpool_limits, threadpool_info
from .model import PublicModel, subset_kernel
from .physics import ReducedKernel
from .quotient_solver import enumerate_quotient
from .branch_enumerator import filter_bilinear_dual_roots
from .general_quotient_solver import enumerate_general_quotient, filter_general_roots


class BudgetExhausted(RuntimeError):
    pass


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, complex):
        return dict(real=value.real, imag=value.imag)
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def load_public(source):
    """Read one explicitly supplied public NPZ, or copy its public field names."""
    allowed = ('y', 'G', 'E0', 'S', 'centers', 'sources', 'receivers', 'radius', 'k',
               'geometry', 'box', 'case_id')
    if isinstance(source, (str, Path)):
        with np.load(source, allow_pickle=False) as data:
            return {key: data[key] for key in allowed if key in data.files}
    return {key: source[key] for key in allowed if key in source}


DEFAULTS = dict(wall_seconds=60., seed=2026103001, fullj_nstarts=32,
                start_policy='strong_constants_sobol',
                fullj_max_nfev=200, algebraic_initial_starts=1,
                initial_max_nfev=80, correction_max_nfev=20,
                correction_per_block=1, beam_width=3, rounds=8,
                blocks_per_round=None, active_size=4, taus=[1., .5, .25, 0.],
                merge_rms=1e-5, fit_ftol=1e-11, fit_xtol=1e-11, fit_gtol=1e-11,
                full_rhs_limit=None, full_lu_limit=None, stop_relative_residual=1e-10,
                coordinate_sweeps=8,
                coordinate_correction_max_nfev=0, coordinate_pre_sweeps=0,
                beam_diversity_rms=0., parameterization='material',
                globaldual_nstarts=1024, globaldual_maxiter=200,
                globaldual_maxfun=1000, globaldual_correction_max_nfev=120,
                globaldual_ftol=1e-12, globaldual_gtol=1e-10)


class Run:
    def __init__(self, public, method, config):
        self.started = time.perf_counter()
        self.config = dict(DEFAULTS)
        self.config.update(config or {})
        self.method = method
        self.model = PublicModel(public)
        self.y = np.asarray(public['y'], dtype=complex)
        if self.y.shape != (self.model.nr, self.model.ns) or not np.all(np.isfinite(self.y)):
            raise ValueError('y must be finite complex[nr,ns]')
        box = np.asarray(public.get('box', [.5, 30.]), dtype=float)
        if box.shape == (2,):
            self.lower = np.full(self.model.n, box[0])
            self.upper = np.full(self.model.n, box[1])
        elif box.shape == (2, self.model.n):
            self.lower, self.upper = box
        else:
            raise ValueError('box shape must be(2,) or(2,n)')
        if not np.all(self.upper > self.lower):
            raise ValueError('Every material upper bound must exceed lower bound')
        self.width, self.midpoint = self.upper-self.lower, (self.lower+self.upper)/2
        self.scale = max(np.linalg.norm(self.y), 1e-300)
        self.target = self.model.pack(self.y)
        self.archive, self.events, self.best_trace, self.rejections = [], [], [], Counter()
        self.best_seen = np.inf
        self.best_seen_material = None
        self.kernel = None
        self.fits = []
        self.root_calls = 0
        self.corrected_candidates = set()
        self.status, self.stop_reason = 'COMPLETED_NUMERICAL', None
        original_operation = self.model.operation

        @contextmanager
        def checked_operation(kind, rhs=0, **cost):
            self.check_budget(rhs=rhs, full_lu=cost.get('LU', 0))
            with original_operation(kind, rhs=rhs, **cost) as action:
                yield action
        self.model.operation = checked_operation

    def elapsed(self):
        return time.perf_counter()-self.started

    def check_budget(self, rhs=0, full_lu=0):
        cfg = self.config
        if cfg['wall_seconds'] is not None and self.elapsed() >= cfg['wall_seconds']:
            raise BudgetExhausted('WALL_SECONDS')
        if cfg['full_rhs_limit'] is not None and self.model.counts['actual_complex_RHS']+rhs > cfg['full_rhs_limit']:
            raise BudgetExhausted('FULL_COMPLEX_RHS')
        if cfg['full_lu_limit'] is not None and self.model.counts['LU']+full_lu > cfg['full_lu_limit']:
            raise BudgetExhausted('FULL_FACTORIZATIONS')

    def discover(self, state, origin):
        residual = float(np.linalg.norm(state['y']-self.y)/self.scale)
        if residual < self.best_seen:
            self.best_seen = residual
            self.best_seen_material = state['m'].copy()
            self.best_trace.append(dict(seconds=self.elapsed(), relative_residual=residual,
                                        m=state['m'].copy(),
                                        full_complex_rhs=self.model.counts['actual_complex_RHS'],
                                        full_lu=self.model.counts['LU'],
                                        reduced_complex_rhs=self.model.counts['reduced_complex_RHS'],
                                        reduced_lu=self.model.counts['reduced_LU'],
                                        algebraic_calls=self.root_calls, origin=origin))
        return residual

    def add(self, m, origin, state=None, metadata=None):
        self.check_budget()
        m = np.asarray(m, dtype=float)
        if not np.all(np.isfinite(m)) or np.any(m < self.lower-1e-8) or np.any(m > self.upper+1e-8):
            self.rejections['INVALID_OR_OUTSIDE_FULL_BOX'] += 1
            return None
        for index, old in enumerate(self.archive):
            if np.sqrt(np.mean(((m-old['m'])/self.width)**2)) <= self.config['merge_rms']:
                self.rejections['DUPLICATE_ARCHIVE_MATERIAL'] += 1
                return index
        if state is None:
            state = self.model.forward(m, label=origin)
        residual = self.discover(state, origin)
        row = dict(m=m.copy(), relative_residual=residual,
                   residual_norm=float(np.linalg.norm(state['y']-self.y)),
                   discovery_seconds=self.elapsed(),
                   discovery_full_complex_rhs=self.model.counts['actual_complex_RHS'],
                   discovery_full_lu=self.model.counts['LU'],
                   discovery_reduced_complex_rhs=self.model.counts['reduced_complex_RHS'],
                   discovery_reduced_lu=self.model.counts['reduced_LU'],
                   origin=origin, metadata=metadata or {}, _y=state['y'].copy())
        self.archive.append(row)
        return len(self.archive)-1

    def fit(self, start, max_nfev, origin):
        cache = {}
        reciprocal = self.config['parameterization'] == 'reciprocal'
        if self.config['parameterization'] not in ('material','reciprocal'):
            raise ValueError('parameterization must be material or reciprocal')
        if reciprocal and np.any(self.lower <= 0):
            raise ValueError('Reciprocal parameterization requires strictly positive box')

        def material(x):
            return 1/x if reciprocal else x

        def state(x):
            self.check_budget()
            if 'x' not in cache or not np.array_equal(cache['x'], x):
                current = self.model.forward(material(x), label=origin)
                self.discover(current, origin)
                cache.update(x=x.copy(), state=current)
            return cache['state']

        start = np.clip(start, self.lower+1e-12*self.width, self.upper-1e-12*self.width)
        fit_start = 1/start if reciprocal else start
        fit_bounds = (1/self.upper,1/self.lower) if reciprocal else (self.lower,self.upper)
        before = dict(self.model.counts)
        begun = self.elapsed()
        try:
            fit = least_squares(lambda x: (self.model.pack(state(x)['y'])-self.target)/self.scale,
                                fit_start, jac=lambda x: self.model.jacobian(state(x))*(-material(x)**2 if reciprocal else 1)/self.scale,
                                bounds=fit_bounds, method='trf', max_nfev=int(max_nfev),
                                ftol=self.config['fit_ftol'], xtol=self.config['fit_xtol'],
                                gtol=self.config['fit_gtol'])
            endpoint = state(fit.x)
            index = self.add(endpoint['m'], origin, state=endpoint)
            self.fits.append(dict(origin=origin, nfev=fit.nfev, njev=fit.njev,
                                  status=int(fit.status), message=fit.message,
                                  started_seconds=begun, wall_seconds=self.elapsed()-begun,
                                  candidate=index, counts={key:value-before.get(key,0) for key,value in self.model.counts.items()}))
            return index
        except BudgetExhausted:
            if 'state' in cache:
                # Budget stop preserves the last already-paid state without new work.
                saved = cache['state']
                self.archive.append(dict(m=saved['m'].copy(), relative_residual=float(np.linalg.norm(saved['y']-self.y)/self.scale),
                                         residual_norm=float(np.linalg.norm(saved['y']-self.y)),
                                         discovery_seconds=self.elapsed(), origin=origin+'_budget_endpoint',
                                         discovery_full_complex_rhs=self.model.counts['actual_complex_RHS'],
                                         discovery_full_lu=self.model.counts['LU'], metadata={}, _y=saved['y'].copy()))
            self.fits.append(dict(origin=origin, status='BUDGET_STOP', started_seconds=begun,
                                  wall_seconds=self.elapsed()-begun))
            raise

    def starts(self, count):
        if count < 1:
            return []
        # Same deterministic start prefix is shared by both methods.
        policy = self.config['start_policy']
        if policy == 'midpoint_sobol':
            prefix = [self.midpoint.copy()]
        elif policy == 'strong_constants_sobol':
            prefix = [self.lower+fraction*self.width for fraction in (.01,.1,.5,.99)]
        else:
            raise ValueError('Unknown public start policy: '+str(policy))
        engine = qmc.Sobol(self.model.n, scramble=True, seed=int(self.config['seed']))
        power = int(np.ceil(np.log2(max(count-len(prefix),1))))
        sample = engine.random_base2(power)
        return (prefix+[self.lower+u*self.width for u in sample[:max(count-len(prefix),0)]])[:count]

    def fullj(self):
        for index, start in enumerate(self.starts(int(self.config['fullj_nstarts']))):
            self.fit(start, self.config['fullj_max_nfev'], 'fullj_start_'+str(index))
            if self.best_seen <= self.config['stop_relative_residual']:
                self.stop_reason = 'DATA_RESIDUAL_GOAL'
                break

    def beam(self):
        ranked = sorted(range(len(self.archive)), key=lambda i: self.archive[i]['relative_residual'])
        chosen = []
        distance = float(self.config['beam_diversity_rms'])
        for index in ranked:
            if all(np.sqrt(np.mean((self.archive[index]['m']-self.archive[old]['m'])**2)) >= distance for old in chosen):
                chosen.append(index)
                if len(chosen) >= int(self.config['beam_width']):
                    break
        return chosen

    def algebraic(self):
        if self.model.ns < 2:
            raise ValueError('Active4 Q6 requires at least2 sources')
        for index, start in enumerate(self.starts(int(self.config['algebraic_initial_starts']))):
            if self.config['initial_max_nfev']:
                self.fit(start, self.config['initial_max_nfev'], 'algebraic_initial_'+str(index))
            else:
                self.add(start, 'algebraic_initial_'+str(index))
        if self.best_seen <= self.config['stop_relative_residual']:
            self.stop_reason = 'DATA_RESIDUAL_GOAL_BEFORE_ALGEBRAIC_SETUP'
            return
        self.check_budget()
        self.kernel = ReducedKernel(self.model, self.midpoint, active=np.arange(self.model.n))
        if self.config['coordinate_pre_sweeps']:
            self.coordinate(sweeps=int(self.config['coordinate_pre_sweeps']), initialize=False)
        rng = np.random.default_rng(int(self.config['seed']))
        order = rng.permutation(self.model.n)
        size = int(self.config['active_size'])
        if size not in (4,6) or size > self.model.n:
            raise ValueError('Active size must be4 or6')
        blocks = int(self.config['blocks_per_round'] or np.ceil(self.model.n/size))
        source_sets = list(itertools.combinations(range(self.model.ns), size//2))
        if not source_sets:
            raise ValueError('Not enough sources for chosen active size')
        for round_index in range(int(self.config['rounds'])):
            if self.best_seen <= self.config['stop_relative_residual']:
                self.stop_reason = 'DATA_RESIDUAL_GOAL'
                break
            round_beam = self.beam()
            for block_index in range(blocks):
                # Shifted cyclic covers activate every cell across successive blocks.
                offset = (round_index+block_index*size)%self.model.n
                active = order[(offset+np.arange(size))%self.model.n]
                receiver = (round_index*blocks+block_index)%self.model.nr
                sources = source_sets[(round_index+block_index)%len(source_sets)]
                for beam_index in round_beam:
                    self.slice_step(beam_index, active, receiver, sources, round_index, block_index)
            # Refresh the material ordering without any observation-independent oracle.
            order = np.roll(order, 1)

    def slice_step(self, incumbent_index, active, receiver, sources, round_index, block_index):
        self.check_budget()
        incumbent = self.archive[incumbent_index]
        block = self.kernel.restrict(active, incumbent['m'])
        view = subset_kernel(block, receiver, sources)
        observed = self.y[np.ix_([receiver], sources)]
        current = incumbent['_y'][np.ix_([receiver], sources)]
        corrections = []
        for tau in self.config['taus']:
            self.check_budget()
            target = current+float(tau)*(observed-current)
            event = dict(round=round_index, block=block_index, incumbent=incumbent_index,
                         active=active.copy(), receiver=int(receiver), sources=list(sources), tau=float(tau),
                         started_seconds=self.elapsed(), expected_degree=6 if len(active)==4 else 20)
            self.root_calls += 1
            self.model.counts['algebraic_expected_roots'] += event['expected_degree']
            with self.model.operation('R3_ALGEBRAIC_CALL', algebraic_calls=1):
                pass
            try:
                roots = enumerate_quotient(view, target) if len(active)==4 else enumerate_general_quotient(view, target)
                filtered = filter_bilinear_dual_roots(roots, box=(float(np.min(self.lower)),float(np.max(self.upper))), validate=False) if len(active)==4 else filter_general_roots(roots, box=(float(np.min(self.lower)),float(np.max(self.upper))))
                event.update(status='ENUMERATED', resolved_paths=roots['resolved_paths'],
                             unresolved_paths=roots['unresolved_paths'],
                             distinct_finite_endpoints=roots['distinct_finite_endpoints'],
                             all_paths_numerically_resolved=roots['all_paths_numerically_resolved'],
                             certificate_flag='NOT_RIGOROUSLY_CERTIFIED', counters=roots['counters'])
                self.model.counts['algebraic_resolved_roots'] += roots['resolved_paths']
                self.model.counts['algebraic_unresolved_roots'] += roots['unresolved_paths']
                event['physical_filter'] = []
                for row in filtered['roots']:
                    summary = {key:value for key,value in row.items() if key not in ('x','m')}
                    if row['status'] != 'CANDIDATE':
                        self.rejections[row.get('reason','UNKNOWN_ROOT_REJECTION')] += 1
                    else:
                        origin = 'q'+str(event['expected_degree'])+'_r'+str(round_index)+'_b'+str(block_index)+'_tau'+str(tau)+'_p'+str(row['path'])
                        delta = (row['m']-block.m0)[active]
                        exact_state = block.state(delta)
                        archive_size_before = len(self.archive)
                        candidate = self.add(row['m'], origin, state=exact_state)
                        summary['candidate'] = candidate
                        if candidate is not None:
                            entry = self.archive[candidate]
                            candidate_mask = entry['_y'][np.ix_([receiver], sources)]
                            mask_old = float(np.linalg.norm(current-observed)**2)
                            mask_new = float(np.linalg.norm(candidate_mask-observed)**2)
                            full_change = entry['residual_norm']**2-incumbent['residual_norm']**2
                            expected_mask_gain = -(2*float(tau)-float(tau)**2)*mask_old
                            summary.update(full_loss_change=full_change,
                                           actual_mask_loss_change=mask_new-mask_old,
                                           exact_target_mask_loss_change=expected_mask_gain,
                                           target_equality_relative=float(np.linalg.norm(candidate_mask-target)/max(np.linalg.norm(target),1e-300)),
                                           withheld_loss_change=full_change-(mask_new-mask_old),
                                           material_rms_from_incumbent=float(np.sqrt(np.mean((entry['m']-incumbent['m'])**2))),
                                           box_normalized_rms_from_incumbent=float(np.sqrt(np.mean(((entry['m']-incumbent['m'])/self.width)**2))))
                            if candidate < archive_size_before:
                                summary['correction_eligibility'] = 'EXCLUDED_DUPLICATE_OR_UNCHANGED'
                                self.rejections['CORRECTION_DUPLICATE_OR_UNCHANGED'] += 1
                            elif candidate in self.corrected_candidates:
                                summary['correction_eligibility'] = 'EXCLUDED_ALREADY_CORRECTED'
                                self.rejections['CORRECTION_ALREADY_CORRECTED'] += 1
                            else:
                                # A worse all-data root can enter a better local
                                # basin after correction. Incumbent stays archived.
                                summary['correction_eligibility'] = 'NEW_BRANCH_ELIGIBLE'
                                corrections.append(candidate)
                    event['physical_filter'].append(summary)
            except BudgetExhausted:
                event['status'] = 'BUDGET_STOP'
                raise
            except (ValueError, ArithmeticError, np.linalg.LinAlgError, FloatingPointError) as exc:
                self.rejections['ALGEBRAIC_CALL_FAILURE'] += 1
                self.model.counts['algebraic_unresolved_roots'] += event['expected_degree']
                event.update(status='FAILED_PAID', error=repr(exc))
            finally:
                event['wall_seconds'] = self.elapsed()-event['started_seconds']
                self.events.append(event)
        count = int(self.config['correction_per_block'])
        for candidate in sorted(set(corrections), key=lambda i: self.archive[i]['relative_residual'])[:count]:
            if self.config['correction_max_nfev']:
                self.corrected_candidates.add(candidate)
                self.fit(self.archive[candidate]['m'], self.config['correction_max_nfev'], 'fullj_correction_'+str(candidate))

    def coordinate(self, sweeps=None, initialize=True):
        """Exact rational 1D all-observation minimization on every coordinate.

        Along an exact one-coordinate slice, full loss is a real quadratic
        divided by a real quadratic. Its stationary numerator has degree<=2.
        Every finite real stationary point and both box endpoints are scored;
        incumbent is retained. These are global coordinate choices conditional
        on a nonsingular fixed-complement Schur representation, not global
        material optimization or a multidimensional branch certificate.
        """
        if initialize:
            for index, start in enumerate(self.starts(int(self.config['algebraic_initial_starts']))):
                if self.config['initial_max_nfev']:
                    self.fit(start, self.config['initial_max_nfev'], 'coordinate_initial_'+str(index))
                else:
                    self.add(start, 'coordinate_initial_'+str(index))
            if self.best_seen <= self.config['stop_relative_residual']:
                self.stop_reason = 'DATA_RESIDUAL_GOAL_BEFORE_COORDINATE_SETUP'
                return
            self.kernel = ReducedKernel(self.model, self.midpoint, active=np.arange(self.model.n))
        order = np.random.default_rng(int(self.config['seed'])).permutation(self.model.n)
        for sweep in range(int(self.config['coordinate_sweeps'] if sweeps is None else sweeps)):
            if self.best_seen <= self.config['stop_relative_residual']:
                self.stop_reason = 'DATA_RESIDUAL_GOAL'
                break
            for coordinate in order:
                self.check_budget()
                incumbent_index = self.beam()[0]
                incumbent = self.archive[incumbent_index]
                block = self.kernel.restrict([coordinate], incumbent['m'])
                K = block.K[0,0]
                a = block.y0-self.y
                b = block.L@block.q-K*a
                n0, n1, n2 = float(np.vdot(a,a).real), float(2*np.vdot(a,b).real), float(np.vdot(b,b).real)
                d1, d2 = float(-2*K.real), float(abs(K)**2)
                polynomial = np.array([n1-n0*d1, 2*(n2-n0*d2), n2*d1-n1*d2])
                with self.model.operation('RATIONAL_COORDINATE_STATIONARITY', coordinate_polynomial_calls=1,
                                          algebraic_calls=1, polynomial_root_calls=1):
                    self.root_calls += 1
                    nonzero = np.flatnonzero(np.abs(polynomial)>1e-14*max(np.max(np.abs(polynomial)),1e-300))
                    roots = np.roots(polynomial[:nonzero[-1]+1][::-1]) if len(nonzero) and nonzero[-1]>0 else np.array([])
                lo, hi = self.lower[coordinate]-block.m0[coordinate], self.upper[coordinate]-block.m0[coordinate]
                deltas = [lo, hi, incumbent['m'][coordinate]-block.m0[coordinate]]
                rejected = []
                for root in roots:
                    if abs(root.imag) > 1e-8*(1+abs(root.real)):
                        rejected.append('NONREAL_COORDINATE_STATIONARY_POINT')
                    elif root.real < lo-1e-8 or root.real > hi+1e-8:
                        rejected.append('OUTSIDE_COORDINATE_BOX')
                    else:
                        deltas.append(float(root.real))
                candidates = []
                for delta in deltas:
                    if abs(1-K*delta) <= 1e-10:
                        rejected.append('COORDINATE_POLE')
                        continue
                    state = block.state([delta])
                    candidate = self.add(state['m'], 'coordinate_s'+str(sweep)+'_c'+str(coordinate), state=state)
                    if candidate is not None:
                        candidates.append(candidate)
                for reason in rejected:
                    self.rejections[reason] += 1
                self.events.append(dict(method='EXACT_RATIONAL_COORDINATE', sweep=sweep, coordinate=int(coordinate),
                                        incumbent=incumbent_index, stationary_polynomial=polynomial,
                                        stationary_roots=roots, candidates=candidates, rejections=rejected,
                                        seconds=self.elapsed()))
            if self.config['coordinate_correction_max_nfev']:
                self.fit(self.archive[self.beam()[0]]['m'], self.config['coordinate_correction_max_nfev'],
                         'coordinate_fullj_correction_s'+str(sweep))
            order = np.roll(order,1)

    def globaldual(self):
        """All-receiver variable projection, then original-data validation/polish.

        The eliminated dual norm is a material-dependent metric and is only a
        proposal objective. Near-resolvent poles can make it misleading; every
        endpoint is scored by a fresh physical full forward and online ranking
        always uses the original observation residual.
        """
        q, K, L = self.model.E0, self.model.G, self.model.S
        with self.model.operation('GLOBALDUAL_SETUP', algebraic_SVD=3,
                                  globaldual_setups=1):
            singular = np.linalg.svd(q,compute_uv=False)
            if singular[-1] <= 1e-12*singular[0]:
                raise ValueError('GLOBALDUAL_SOURCE_NOT_FULL_COLUMN_RANK')
            Z = null_space(q.T,rcond=1e-12).T
            W0 = self.y@np.linalg.pinv(q,rcond=1e-12)
            V0, ZK = L+W0@K, Z@K
            if np.linalg.norm(Z@q)>1e-10*max(np.linalg.norm(Z)*np.linalg.norm(q),1e-300):
                raise ValueError('GLOBALDUAL_NULLSPACE_RESIDUAL')
        dimension = Z.shape[0]
        if dimension == 0:
            raise ValueError('GLOBALDUAL_REQUIRES_N_GREATER_THAN_SOURCE_COUNT')
        objective_scale = max(float(np.linalg.norm(W0)**2),1e-300)
        self.events.append(dict(method='GLOBALDUAL_SETUP', nullspace_dimension=dimension,
                                source_singular_values=singular, objective_scale=objective_scale,
                                nullspace_relative_residual=float(np.linalg.norm(Z@q)/max(np.linalg.norm(Z)*np.linalg.norm(q),1e-300)),
                                objective='dual constitutive least-squares; proposal metric only'))
        for index, start in enumerate(self.starts(int(self.config['globaldual_nstarts']))):
            self.check_budget()
            begun = self.elapsed()
            record = dict(method='GLOBALDUAL_VARPRO', start=index, started_seconds=begun)
            last = {}

            def objective(m):
                self.check_budget()
                A = Z-ZK*m[None,:]
                B0 = W0-V0*m[None,:]
                with self.model.operation('GLOBALDUAL_INNER_QR', algebraic_QR=1,
                                          globaldual_inner_factorizations=1,
                                          algebraic_complex_rhs=self.model.nr,
                                          globaldual_inner_complex_rhs=self.model.nr,
                                          globaldual_objective_gradient_calls=1):
                    CT, _, rank, _ = lstsq(A.T,-B0.T,cond=1e-12,lapack_driver='gelsy',check_finite=False)
                if rank != dimension:
                    self.model.counts['globaldual_rank_failures'] += 1
                    raise ValueError('GLOBALDUAL_INNER_RANK_FAILURE')
                C = CT.T
                W, V = W0+C@Z, V0+C@ZK
                error = W-V*m[None,:]
                value = float(np.vdot(error,error).real)/objective_scale
                gradient = -2*np.real(np.sum(error.conj()*V,axis=0))/objective_scale
                if not np.isfinite(value) or not np.all(np.isfinite(gradient)):
                    raise FloatingPointError('NONFINITE_GLOBALDUAL_OBJECTIVE')
                last.update(m=m.copy(), objective=value, rank=int(rank),
                            inner_orthogonality_relative=float(np.linalg.norm(error@A.conj().T)/max(np.linalg.norm(error)*np.linalg.norm(A),1e-300)))
                return value, gradient
            try:
                fit = minimize(objective,start,method='L-BFGS-B',jac=True,
                               bounds=list(zip(self.lower,self.upper)),
                               options=dict(maxiter=int(self.config['globaldual_maxiter']),
                                            maxfun=int(self.config['globaldual_maxfun']),
                                            ftol=self.config['globaldual_ftol'],gtol=self.config['globaldual_gtol']))
                record.update(status=int(fit.status), success=bool(fit.success), message=str(fit.message),
                              nfev=int(fit.nfev), njev=int(fit.njev), iterations=int(fit.nit),
                              dual_objective=float(fit.fun), endpoint=fit.x.copy(),
                              last_inner=last)
                state = self.model.forward(fit.x,label='globaldual_endpoint_'+str(index))
                archive_before = len(self.archive)
                candidate = self.add(fit.x,'globaldual_endpoint_'+str(index),state=state,
                                     metadata=dict(dual_objective=float(fit.fun),proposal_only=True))
                record['candidate'] = candidate
                if candidate is not None and candidate >= archive_before and self.config['globaldual_correction_max_nfev']:
                    self.fit(fit.x,self.config['globaldual_correction_max_nfev'],'globaldual_fullj_correction_'+str(index))
                if self.best_seen <= self.config['stop_relative_residual']:
                    self.stop_reason = 'DATA_RESIDUAL_GOAL'
                    break
            except BudgetExhausted:
                record.update(status='BUDGET_STOP',last_inner=last)
                raise
            except (ValueError,np.linalg.LinAlgError,FloatingPointError) as exc:
                self.rejections['GLOBALDUAL_START_FAILURE'] += 1
                record.update(status='FAILED_PAID',error=repr(exc),last_inner=last)
            finally:
                record['wall_seconds'] = self.elapsed()-begun
                self.events.append(record)

    def finish(self):
        # The globally best paid trial is retained even if an optimizer budget
        # interrupt occurred between fun and jac or before endpoint archival.
        if self.best_seen_material is not None and not any(np.array_equal(self.best_seen_material,row['m']) for row in self.archive):
            self.archive.append(dict(m=self.best_seen_material.copy(), relative_residual=self.best_seen,
                                     residual_norm=self.best_seen*self.scale,
                                     discovery_seconds=self.best_trace[-1]['seconds'],
                                     discovery_full_complex_rhs=self.best_trace[-1]['full_complex_rhs'],
                                     discovery_full_lu=self.best_trace[-1]['full_lu'],
                                     origin='best_paid_trial', metadata={}))
        candidates = [{key:value for key,value in row.items() if not key.startswith('_')} for row in self.archive]
        ranking = sorted(range(len(candidates)), key=lambda i:candidates[i]['relative_residual'])
        counts = dict(self.model.counts)
        accounting = dict(full_factorizations=counts.get('LU',0),
                          full_complex_rhs=counts.get('actual_complex_RHS',0),
                          reduced_factorizations=counts.get('reduced_LU',0),
                          reduced_complex_rhs=counts.get('reduced_complex_RHS',0),
                          counts=counts, operation_seconds=dict(self.model.timings),
                          convention='Full factors include all-N fixed-active setup; reduced factors include exact complementary Schur solves and active rational evaluation. RHS counts are actual complex columns, not LU-equivalent estimates.')
        return json_safe(dict(status=self.status, stop_reason=self.stop_reason, method=self.method,
                              solver_revision='R3_PUBLIC_V6_GLOBALDUAL_RECIPROCAL_CONTROL',
                              config=self.config, material_dimension=self.model.n,
                              receiver_count=self.model.nr, source_count=self.model.ns,
                              all_materials_unknown=True, objective='unweighted_all_observation_complex_L2',
                              selection='minimum all-observation relative residual; no GT',
                              candidates=candidates, ranking=ranking, top1_index=ranking[0] if ranking else None,
                              top1=candidates[ranking[0]] if ranking else None,
                              best_so_far=self.best_trace, fits=self.fits, algebraic_events=self.events,
                              rejections=dict(self.rejections), wall_seconds=self.elapsed(), accounting=accounting,
                              blas_thread_pools=threadpool_info(),
                              scope='Finite discrete scalar-TM nonlinear VIE. Conditional exact coordinate slices; numerical roots only. No global coverage or optimization certificate.'))


def solve_case(public_npz_or_dict, method='algebraic', config=None):
    """Return a JSON-safe result with observation-selected candidates and costs."""
    public = load_public(public_npz_or_dict)
    aliases = {'quotient':'algebraic', 'q6':'algebraic', 'baseline':'fullj'}
    method = aliases.get(method, method)
    if method not in ('algebraic','fullj','coordinate','globaldual'):
        raise ValueError('method must be algebraic, fullj, coordinate, or globaldual')
    with threadpool_limits(limits=1):
        run = Run(public, method, config)
        try:
            {'fullj':run.fullj, 'algebraic':run.algebraic, 'coordinate':run.coordinate,
             'globaldual':run.globaldual}[method]()
        except BudgetExhausted as exc:
            run.status, run.stop_reason = 'BUDGET_EXHAUSTED', str(exc)
        except Exception as exc:
            run.status, run.stop_reason = 'FAILED_PAID', repr(exc)
        return run.finish()
