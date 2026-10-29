"""Deterministic 24-path numerical multiaffine homotopy for four variables.

Numerical path resolution is recorded explicitly.  This module does not turn
floating point path tracking into a rigorous completeness certificate.
"""
import itertools
import time
import numpy as np
from scipy.linalg import null_space


class MultiaffineSystem:
    def __init__(self, coefficients, normalize=True):
        c = np.asarray(coefficients)
        self.r = int(round(np.log2(c.shape[0])))
        if c.ndim != 2 or c.shape[0] != 1 << self.r:
            raise ValueError('Need coefficients[2**r,equations]')
        self.raw_coefficients = c.copy()
        self.scale = np.max(np.abs(c), axis=0) if normalize else np.ones(c.shape[1])
        if np.any(self.scale == 0):
            raise ValueError('An identically zero equation cannot define a square root problem')
        self.c = c/self.scale
        self.supports = [tuple(i for i in range(self.r) if (mask >> i) & 1)
                         for mask in range(1 << self.r)]
        self.evaluations = self.jacobian_evaluations = 0

    def evaluate(self, x, jac=False):
        x = np.asarray(x, dtype=complex)
        if x.shape != (self.r,):
            raise ValueError('Incorrect variable dimension')
        self.evaluations += 1
        f = np.zeros(self.c.shape[1], dtype=complex)
        J = np.zeros((self.c.shape[1], self.r), dtype=complex) if jac else None
        if jac:
            self.jacobian_evaluations += 1
        for mask, support in enumerate(self.supports):
            monomial = np.prod(x[list(support)]) if support else 1.
            f += self.c[mask]*monomial
            if jac:
                for variable in support:
                    rest = [j for j in support if j != variable]
                    J[:, variable] += self.c[mask]*(np.prod(x[rest]) if rest else 1.)
        return (f, J) if jac else f


def _start_evaluate(x, a):
    r = len(x)
    differences = x[None, :]-a
    values = np.prod(differences, axis=1)
    J = np.empty((r, r), dtype=complex)
    for j in range(r):
        J[:, j] = np.prod(differences[:, [k for k in range(r) if k != j]], axis=1)
    return values, J


def enumerate_multiaffine(coefficients, seed=2026102903, max_steps=4000,
                          tolerance=2e-11, initial_step=.01, min_step=1e-9,
                          max_step=.04, max_newton=12, divergence_threshold=1e8,
                          model=None, wall_seconds=None):
    """Track the 4! isolated start roots through one frozen complex homotopy.

    G_i=product_j(x_j-a_ij), H=(1-t)*gamma*G+t*P.  Failed, singular,
    divergent-numerical, and real/complex terminal paths remain separate.
    """
    system = MultiaffineSystem(coefficients)
    r = system.r
    if r != 4 or system.c.shape[1] != r:
        raise ValueError('This frozen enumerator requires four variables and four real equations')
    if np.max(np.abs(np.imag(system.c))) > 1e-14:
        raise ValueError('P must have real coefficients before complex continuation')
    rng = np.random.default_rng(seed)
    # Distinct entries in each column ensure all permutation start roots are simple.
    a = rng.uniform(-.8, 1.8, (r, r))+1j*rng.uniform(.3, 1.5, (r, r))
    gamma = np.exp(1j*rng.uniform(.2, 2*np.pi-.2))
    paths, roots = [], []
    counters = dict(algebraic_linear_solves=0, homotopy_evaluations=0,
                    newton_iterations=0, rejected_steps=0)
    begin = time.perf_counter()

    def homotopy(x, t):
        p, Jp = system.evaluate(x, jac=True)
        g, Jg = _start_evaluate(x, a)
        counters['homotopy_evaluations'] += 1
        return (1-t)*gamma*g+t*p, (1-t)*gamma*Jg+t*Jp, p-gamma*g

    def solve(A, b):
        counters['algebraic_linear_solves'] += 1
        return np.linalg.solve(A, b)

    for path_index, permutation in enumerate(itertools.permutations(range(r))):
        x = np.empty(r, dtype=complex)
        for equation, variable in enumerate(permutation):
            x[variable] = a[equation, variable]
        t, step, accepted, rejected = 0., initial_step, 0, 0
        status, reason, last_condition = 'UNRESOLVED', None, None
        maximum_norm = float(np.linalg.norm(x))
        path_start = time.perf_counter()
        for attempt in range(max_steps):
            if wall_seconds is not None and time.perf_counter()-begin >= wall_seconds:
                status, reason = 'UNRESOLVED_WALL_LIMIT', 'Frozen wall budget exhausted'
                break
            dt = min(step, 1.-t)
            try:
                H, J, dHdt = homotopy(x, t)
                predictor = x-dt*solve(J, dHdt)
                z = predictor.copy()
                converged, corrections = False, 0
                next_t = t+dt
                for corrections in range(1, max_newton+1):
                    Hn, Jn, _ = homotopy(z, next_t)
                    if np.linalg.norm(Hn, ord=np.inf) <= tolerance:
                        converged = True
                        break
                    update = solve(Jn, -Hn)
                    counters['newton_iterations'] += 1
                    z += update
                    if not np.all(np.isfinite(z)):
                        break
                    if np.linalg.norm(update, ord=np.inf) <= tolerance*(1+np.linalg.norm(z, ord=np.inf)):
                        check = homotopy(z, next_t)[0]
                        converged = np.linalg.norm(check, ord=np.inf) <= tolerance*10
                        if converged:
                            break
                predictor_error = np.linalg.norm(z-predictor)/(1+np.linalg.norm(z))
                # Reject a large Newton jump: otherwise tracker can silently hop paths.
                if not converged or predictor_error > .15:
                    raise ArithmeticError('corrector failed or predictor correction exceeded .15')
                x, t = z, next_t
                accepted += 1
                maximum_norm = max(maximum_norm, float(np.linalg.norm(x)))
                if maximum_norm > divergence_threshold:
                    status, reason = 'UNRESOLVED_DIVERGENT_NUMERICAL', 'Affine norm exceeded frozen threshold; infinity not certified'
                    break
                if t >= 1.-1e-15:
                    residual, endpoint_J = system.evaluate(x, jac=True)
                    last_condition = float(np.linalg.cond(endpoint_J))
                    if np.linalg.norm(residual, ord=np.inf) > tolerance*20:
                        status, reason = 'UNRESOLVED_ENDPOINT_RESIDUAL', 'P residual exceeds frozen numerical goal'
                    elif not np.isfinite(last_condition) or last_condition > 1e12:
                        status, reason = 'UNRESOLVED_SINGULAR_ENDPOINT', 'Endpoint condition exceeds 1e12'
                    else:
                        status = 'RESOLVED_FINITE_NUMERICAL'
                        roots.append(x.copy())
                    break
                if corrections <= 4 and predictor_error < .005:
                    step = min(max_step, step*1.35)
                elif corrections >= 8 or predictor_error > .05:
                    step = max(min_step, step*.65)
            except (np.linalg.LinAlgError, ArithmeticError, FloatingPointError) as exc:
                rejected += 1
                counters['rejected_steps'] += 1
                step *= .5
                if step < min_step:
                    status, reason = 'UNRESOLVED_MIN_STEP', str(exc)
                    break
        else:
            status, reason = 'UNRESOLVED_STEP_LIMIT', 'Frozen maximum attempted steps exhausted'
        endpoint_residual = float(np.linalg.norm(system.evaluate(x), ord=np.inf))
        paths.append(dict(path=path_index, permutation=list(permutation), status=status,
                          reason=reason, t=float(t), x=x.copy(), accepted_steps=accepted,
                          rejected_steps=rejected, maximum_affine_norm=maximum_norm,
                          endpoint_residual=endpoint_residual,
                          endpoint_condition=last_condition,
                          wall_seconds=time.perf_counter()-path_start))

    counters.update(polynomial_evaluations=system.evaluations,
                    polynomial_J_evaluations=system.jacobian_evaluations,
                    tracked_paths=len(paths))
    if model is not None:
        elapsed = time.perf_counter()-begin
        model.timings['ALGEBRAIC_HOMOTOPY'] += elapsed
        with model.operation('ALGEBRAIC_HOMOTOPY_SUMMARY', **counters,
                             externally_measured_wall_seconds=elapsed):
            pass
    resolved = sum(p['status'] == 'RESOLVED_FINITE_NUMERICAL' for p in paths)
    return dict(method='MULTIAFFINE_24_PATH_HOMOTOPY', seed=seed, gamma=gamma,
                start_constants=a, paths=paths, roots=np.asarray(roots).reshape(-1, r),
                resolved_paths=resolved, unresolved_paths=len(paths)-resolved,
                all_paths_numerically_resolved=resolved == 24,
                completeness='NOT_RIGOROUSLY_CERTIFIED',
                counters=counters, wall_seconds=time.perf_counter()-begin,
                equation_scales=system.scale,
                limits=dict(max_steps=max_steps,tolerance=tolerance,min_step=min_step,
                            divergence_threshold=divergence_threshold))


def filter_real_box_roots(result, polynomial, model=None, target=None,
                          imaginary_tolerance=1e-7, merge_rms=1e-7,
                          denominator_floor=1e-10, validate=True):
    """Filter numerically real physical roots; preserve every rejection reason."""
    kernel = polynomial['kernel']
    h = polynomial['halfwidth']
    rows, candidates = [], []
    for path in result['paths']:
        if path['status'] != 'RESOLVED_FINITE_NUMERICAL':
            continue
        x = path['x']
        row = dict(path=path['path'], x=x.copy(), status='REJECTED')
        if np.max(np.abs(x.imag)) > imaginary_tolerance:
            row['reason'] = 'COMPLEX'
        elif np.any(x.real < -imaginary_tolerance) or np.any(x.real > 1+imaginary_tolerance):
            row['reason'] = 'OUTSIDE_REAL_BOX'
        else:
            delta = h*x.real
            with kernel.model.operation('ALGEBRAIC_POLE_FILTER', algebraic_determinants=1):
                determinant = np.linalg.det(kernel.coefficients(delta)[0])
            row['denominator_abs'] = float(abs(determinant))
            if abs(determinant) <= denominator_floor:
                row['reason'] = 'DENOMINATOR_NEAR_ZERO'
            else:
                m = kernel.m0+kernel.U@delta
                row['m'] = m
                if any(np.linalg.norm(m-old)/np.sqrt(len(m)) <= merge_rms for old in candidates):
                    row['reason'] = 'DUPLICATE_MATERIAL_ROOT'
                else:
                    candidates.append(m)
                    row.update(status='CANDIDATE', reason=None)
                    if validate:
                        if model is None or target is None:
                            raise ValueError('Independent validation requires model and target')
                        try:
                            from .algebraic_solver import validate_candidate
                        except ImportError:
                            from algebraic_solver import validate_candidate
                        row['validation'] = validate_candidate(model, m, target,
                                                               box=polynomial['box'],
                                                               label='homotopy_path_'+str(path['path']))
                        row['status'] = 'QUALIFIED' if row['validation']['qualified'] else 'REJECTED_INDEPENDENT_FORWARD'
        rows.append(row)
    return dict(roots=rows, qualified=[row for row in rows if row['status']=='QUALIFIED'],
                candidate_count=len(candidates), completeness=result['completeness'])


class BilinearDualSystem:
    """Exact four-equation bidegree(1,1) dual of one-receiver/two-source data.

    Variables are c_1,c_2,d_1,d_2.  d is independent during continuation;
    physical roots must satisfy d=conj(c).  Full conjugacy and node checks are
    performed only by the physical filter.
    """
    def __init__(self, kernel, target, rank_rtol=1e-12):
        if kernel.r != 4 or kernel.model.nr != 1 or kernel.model.ns != 2:
            raise ValueError('Bilinear dual requires active4, receiver1, sources2')
        self.kernel = kernel
        self.target = np.asarray(target, dtype=complex)
        q = kernel.q
        with kernel.model.operation('BILINEAR_DUAL_SETUP',
                                    algebraic_SVD=3, bilinear_coefficient_arrays=1):
            singular = np.linalg.svd(q, compute_uv=False)
            if singular[-1] <= rank_rtol*singular[0]:
                raise ValueError('q does not have complex column rank2')
            self.w0 = ((self.target-kernel.y0)@np.linalg.pinv(q, rcond=rank_rtol)).reshape(4)
            # No conjugation in q.T: w q = target is complex linear in w.
            self.N = null_space(q.T, rcond=rank_rtol).T
            self.W = np.vstack([self.w0, self.N])
            self.V = np.vstack([kernel.L.reshape(4)+self.w0@kernel.K,
                                self.N@kernel.K])
            self.raw_coefficients = np.empty((4, 3, 3), dtype=complex)
            for i in range(4):
                self.raw_coefficients[i] = (np.outer(self.W[:, i], self.V[:, i].conj())
                                             - np.outer(self.V[:, i], self.W[:, i].conj()))
            self.scale = np.max(np.abs(self.raw_coefficients), axis=(1, 2))
            if np.any(self.scale == 0):
                raise ValueError('Degenerate zero dual equation')
            self.coefficients = self.raw_coefficients/self.scale[:, None, None]
        self.evaluations = self.jacobian_evaluations = 0

    def evaluate(self, x, jac=False):
        x = np.asarray(x, dtype=complex)
        c, d = np.r_[1., x[:2]], np.r_[1., x[2:]]
        self.evaluations += 1
        f = np.einsum('a,b,iab->i', c, d, self.coefficients)
        if not jac:
            return f
        self.jacobian_evaluations += 1
        J = np.column_stack([np.einsum('b,ib->i', d, self.coefficients[:, j, :]) for j in (1, 2)]
                            + [np.einsum('a,ia->i', c, self.coefficients[:, :, j]) for j in (1, 2)])
        return f, J


def enumerate_bilinear_dual(kernel, target, seed=2026102904, max_steps=4000,
                            tolerance=2e-11, initial_step=.01, min_step=1e-9,
                            max_step=.04, max_newton=12, divergence_threshold=1e8,
                            wall_seconds=None):
    """Track six generic product-affine starts in P²×P², in an affine chart.

    Affine chart escape is unresolved, never silently discarded.  Complete
    resolution of six distinct nonsingular affine endpoints is numerical
    evidence; rigorous projective completeness is not claimed here.
    """
    system = BilinearDualSystem(kernel, target)
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(4, 3))+1j*rng.normal(size=(4, 3))
    B = rng.normal(size=(4, 3))+1j*rng.normal(size=(4, 3))
    gamma = np.exp(1j*rng.uniform(.2, 2*np.pi-.2))
    counters = dict(algebraic_linear_solves=0, homotopy_evaluations=0,
                    newton_iterations=0, rejected_steps=0)
    begin, paths, roots = time.perf_counter(), [], []

    def start_evaluate(x):
        av, bv = A@np.r_[1., x[:2]], B@np.r_[1., x[2:]]
        return av*bv, np.column_stack([A[:, 1:]*bv[:, None], B[:, 1:]*av[:, None]])

    def homotopy(x, t):
        p, Jp = system.evaluate(x, jac=True)
        g, Jg = start_evaluate(x)
        counters['homotopy_evaluations'] += 1
        return (1-t)*gamma*g+t*p, (1-t)*gamma*Jg+t*Jp, p-gamma*g

    def solve(matrix, vector):
        counters['algebraic_linear_solves'] += 1
        return np.linalg.solve(matrix, vector)

    for path_index, chosen in enumerate(itertools.combinations(range(4), 2)):
        complement = tuple(i for i in range(4) if i not in chosen)
        x = np.r_[solve(A[list(chosen), 1:], -A[list(chosen), 0]),
                  solve(B[list(complement), 1:], -B[list(complement), 0])]
        t, step, accepted, rejected = 0., initial_step, 0, 0
        status, reason, condition = 'UNRESOLVED', None, None
        maximum_norm = float(np.linalg.norm(x))
        path_begin = time.perf_counter()
        for attempt in range(max_steps):
            if wall_seconds is not None and time.perf_counter()-begin >= wall_seconds:
                status, reason = 'UNRESOLVED_WALL_LIMIT', 'Frozen wall budget exhausted'
                break
            dt = min(step, 1.-t)
            try:
                _, J, dHdt = homotopy(x, t)
                prediction = x-dt*solve(J, dHdt)
                z = prediction.copy()
                converged = False
                for corrections in range(1, max_newton+1):
                    residual, Jn, _ = homotopy(z, t+dt)
                    if np.linalg.norm(residual, ord=np.inf) <= tolerance:
                        converged = True
                        break
                    update = solve(Jn, -residual)
                    counters['newton_iterations'] += 1
                    z += update
                    if not np.all(np.isfinite(z)):
                        break
                    if np.linalg.norm(update, ord=np.inf) <= tolerance*(1+np.linalg.norm(z, ord=np.inf)):
                        converged = np.linalg.norm(homotopy(z, t+dt)[0], ord=np.inf) <= tolerance*10
                        if converged:
                            break
                error = np.linalg.norm(z-prediction)/(1+np.linalg.norm(z))
                if not converged or error > .15:
                    raise ArithmeticError('corrector failed or predictor correction exceeded .15')
                x, t = z, t+dt
                accepted += 1
                maximum_norm = max(maximum_norm, float(np.linalg.norm(x)))
                if maximum_norm > divergence_threshold:
                    status, reason = 'UNRESOLVED_DIVERGENT_NUMERICAL', 'Affine chart escape; projective infinity not certified'
                    break
                if t >= 1.-1e-15:
                    residual, endpoint_J = system.evaluate(x, jac=True)
                    condition = float(np.linalg.cond(endpoint_J))
                    if np.linalg.norm(residual, ord=np.inf) > tolerance*20:
                        status, reason = 'UNRESOLVED_ENDPOINT_RESIDUAL', 'Endpoint residual above tolerance'
                    elif not np.isfinite(condition) or condition > 1e12:
                        status, reason = 'UNRESOLVED_SINGULAR_ENDPOINT', 'Endpoint condition exceeds 1e12'
                    else:
                        status = 'RESOLVED_FINITE_NUMERICAL'
                        roots.append(x.copy())
                    break
                if corrections <= 4 and error < .005:
                    step = min(max_step, step*1.35)
                elif corrections >= 8 or error > .05:
                    step = max(min_step, step*.65)
            except (np.linalg.LinAlgError, ArithmeticError, FloatingPointError) as exc:
                rejected += 1
                counters['rejected_steps'] += 1
                step *= .5
                if step < min_step:
                    status, reason = 'UNRESOLVED_MIN_STEP', str(exc)
                    break
        else:
            status, reason = 'UNRESOLVED_STEP_LIMIT', 'Maximum attempted steps exhausted'
        paths.append(dict(path=path_index, start_A_equations=list(chosen),
                          status=status, reason=reason, t=float(t), x=x.copy(),
                          accepted_steps=accepted, rejected_steps=rejected,
                          maximum_affine_norm=maximum_norm,
                          endpoint_residual=float(np.linalg.norm(system.evaluate(x), ord=np.inf)),
                          endpoint_condition=condition,
                          wall_seconds=time.perf_counter()-path_begin))
    distinct = []
    for root in roots:
        if all(np.linalg.norm(root-other) > 1e-6*(1+np.linalg.norm(root)) for other in distinct):
            distinct.append(root)
    counters.update(polynomial_evaluations=system.evaluations,
                    polynomial_J_evaluations=system.jacobian_evaluations,
                    tracked_paths=6)
    elapsed = time.perf_counter()-begin
    kernel.model.timings['BILINEAR_SIX_PATH_HOMOTOPY'] += elapsed
    # Put a summary into op as CountingBook intentionally ignores append.
    with kernel.model.operation('BILINEAR_HOMOTOPY_SUMMARY', **counters,
                                externally_measured_wall_seconds=elapsed):
        pass
    return dict(method='BILINEAR_DUAL_SIX_PATH_HOMOTOPY', seed=seed,
                gamma=gamma, start_A=A, start_B=B, paths=paths,
                roots=np.asarray(roots).reshape(-1, 4),
                resolved_paths=len(roots), unresolved_paths=6-len(roots),
                distinct_finite_endpoints=len(distinct),
                all_paths_numerically_resolved=len(roots)==6 and len(distinct)==6,
                completeness='NOT_RIGOROUSLY_CERTIFIED', counters=counters,
                wall_seconds=elapsed, system=system,
                limits=dict(max_steps=max_steps, tolerance=tolerance,
                            min_step=min_step, divergence_threshold=divergence_threshold))


def filter_bilinear_dual_roots(result, box=(.5, 4.), conjugacy_tolerance=1e-6,
                               material_imaginary_tolerance=1e-7,
                               node_floor=1e-10, denominator_floor=1e-10,
                               constitutive_tolerance=1e-8,
                               merge_rms=1e-7, validate=True):
    """Recover real materials from the dual; independently validate candidates."""
    system = result['system']
    kernel, model = system.kernel, system.kernel.model
    rows, candidates = [], []
    for path in result['paths']:
        if path['status'] != 'RESOLVED_FINITE_NUMERICAL':
            continue
        x = path['x']
        row = dict(path=path['path'], x=x.copy(), status='REJECTED')
        conjugacy = float(np.linalg.norm(x[2:]-x[:2].conj())/(1+np.linalg.norm(x[:2])))
        row['conjugacy_relative'] = conjugacy
        if conjugacy > conjugacy_tolerance:
            row['reason'] = 'NONPHYSICAL_CONJUGACY'
        else:
            w = system.w0+x[:2]@system.N
            v = kernel.L.reshape(4)+w@kernel.K
            row['minimum_node_abs'] = float(np.min(np.abs(v)))
            if np.min(np.abs(v)) <= node_floor:
                row['reason'] = 'UNRESOLVED_DUAL_NODE'
            else:
                complex_delta = w/v
                delta = np.real(w*v.conj())/np.abs(v)**2
                constitutive = np.linalg.norm(w-delta*v)/max(np.linalg.norm(w), 1e-300)
                row['maximum_material_imaginary'] = float(np.max(np.abs(complex_delta.imag)))
                row['dual_constitutive_relative'] = float(constitutive)
                if np.max(np.abs(complex_delta.imag)) > material_imaginary_tolerance or constitutive > constitutive_tolerance:
                    row['reason'] = 'NONREAL_MATERIAL'
                else:
                    m = kernel.m0+kernel.U@delta
                    row['m'] = m
                    with model.operation('ALGEBRAIC_POLE_FILTER', algebraic_determinants=1):
                        determinant = np.linalg.det(kernel.coefficients(delta)[0])
                    row['denominator_abs'] = float(abs(determinant))
                    if np.any(m < box[0]-1e-8) or np.any(m > box[1]+1e-8):
                        row['reason'] = 'OUTSIDE_REAL_BOX'
                    elif abs(determinant) <= denominator_floor:
                        row['reason'] = 'DENOMINATOR_NEAR_ZERO'
                    elif any(np.linalg.norm(m-old)/np.sqrt(len(m)) <= merge_rms for old in candidates):
                        row['reason'] = 'DUPLICATE_MATERIAL_ROOT'
                    else:
                        candidates.append(m)
                        row.update(status='CANDIDATE', reason=None)
                        if validate:
                            try:
                                from .algebraic_solver import validate_candidate
                            except ImportError:
                                from algebraic_solver import validate_candidate
                            row['validation'] = validate_candidate(model, m, system.target,
                                                                   box=box,
                                                                   label='dual_path_'+str(path['path']))
                            row['status'] = 'QUALIFIED' if row['validation']['qualified'] else 'REJECTED_INDEPENDENT_FORWARD'
        rows.append(row)
    return dict(roots=rows, qualified=[r for r in rows if r['status']=='QUALIFIED'],
                candidate_count=len(candidates), completeness=result['completeness'])


def classify_determinantal_endpoints(result, polynomial, relative_floor=1e-7):
    """Audit saved complex24 endpoints against both cleared denominators.

    Complexification retains Re and Im coefficient polynomials separately.
    Thus the second denominator conjugates coefficients but does not conjugate
    x; taking the actual imaginary part of a complex-evaluated response would
    test a different system.  This diagnostic performs no full physics solve.
    """
    kernel, h = polynomial['kernel'], polynomial['halfwidth']
    coefficient_system = MultiaffineSystem(polynomial['coefficients'])
    rows = []
    for path in result['paths']:
        if path['status'] != 'RESOLVED_FINITE_NUMERICAL':
            continue
        x = np.asarray(path['x'], dtype=complex)
        delta = h*x
        A = np.eye(kernel.r)-delta[:, None]*kernel.K
        Abar = np.eye(kernel.r)-delta[:, None]*kernel.K.conj()
        with kernel.model.operation('DETERMINANT_BASE_LOCUS_AUDIT', algebraic_determinants=2,
                                    algebraic_SVD=2):
            det, detbar = np.linalg.det(A), np.linalg.det(Abar)
            singular, singularbar = np.linalg.svd(A, compute_uv=False), np.linalg.svd(Abar, compute_uv=False)
        # Smallest/largest singular value is invariant to a scalar matrix scale.
        gap, gapbar = singular[-1]/singular[0], singularbar[-1]/singularbar[0]
        small, smallbar = bool(gap <= relative_floor), bool(gapbar <= relative_floor)
        rows.append(dict(path=path['path'], determinant=det, conjugate_coefficient_determinant=detbar,
                         resolvent_relative_singular_gap=float(gap),
                         conjugate_coefficient_relative_singular_gap=float(gapbar),
                         denominator_small=small, conjugate_denominator_small=smallbar,
                         cleared_polynomial_residual=float(np.linalg.norm(coefficient_system.evaluate(x), ord=np.inf)),
                         classification='REGULAR_BOTH' if not small and not smallbar
                         else 'BOTH_BASE_LOCUS' if small and smallbar else 'ONE_BASE_LOCUS'))
    return dict(endpoints=rows, counts={name:sum(row['classification']==name for row in rows)
                                      for name in ('REGULAR_BOTH','BOTH_BASE_LOCUS','ONE_BASE_LOCUS')},
                singular_gap_relative_floor=relative_floor, new_physics_RHS=0)
