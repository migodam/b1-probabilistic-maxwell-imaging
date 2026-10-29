"""Balanced n=2s, one-receiver, s-source bilinear quotient eigenroots.

This separate implementation leaves all frozen four-material sources intact.
It constructs maximal minors in sparse polynomials and retains every numerical
gate.  Invertible leading forms are a conditional algebraic construction, not
a floating-point completeness certificate.
"""
import itertools
import math
import time
import numpy as np
from scipy.linalg import null_space


class GeneralQuotientFailure(ValueError):
    def __init__(self, message, diagnostic=None):
        super().__init__(message)
        self.diagnostic = diagnostic or {}


def homogeneous_monomials(variables, degree):
    if variables == 1:
        return [(degree,)]
    return [(first,)+tail for first in range(degree, -1, -1)
            for tail in homogeneous_monomials(variables-1, degree-first)]


class GeneralBilinearDualSystem:
    def __init__(self, kernel, target, rank_rtol=1e-12):
        self.kernel, self.target = kernel, np.asarray(target, dtype=complex)
        self.s, self.n = kernel.model.ns, kernel.r
        if kernel.model.nr != 1 or self.n != 2*self.s:
            raise GeneralQuotientFailure('REQUIRES_ONE_RECEIVER_AND_ACTIVE_N_EQUAL_2S')
        with kernel.model.operation('GENERAL_BILINEAR_SETUP', algebraic_SVD=3):
            singular = np.linalg.svd(kernel.q, compute_uv=False)
            if singular[-1] <= rank_rtol*singular[0]:
                raise GeneralQuotientFailure('SOURCE_MATRIX_NOT_FULL_COLUMN_RANK',
                                             dict(q_singular_values=singular))
            self.w0 = ((self.target-kernel.y0)@np.linalg.pinv(kernel.q, rcond=rank_rtol)).reshape(self.n)
            self.N = null_space(kernel.q.T, rcond=rank_rtol).T
            if self.N.shape != (self.s, self.n):
                raise GeneralQuotientFailure('UNEXPECTED_COMPLEX_NULLSPACE_DIMENSION')
            self.W = np.vstack([self.w0, self.N])
            self.V = np.vstack([kernel.L.reshape(self.n)+self.w0@kernel.K,
                                self.N@kernel.K])
            self.raw_coefficients = np.asarray([
                np.outer(self.W[:, i], self.V[:, i].conj())
                - np.outer(self.V[:, i], self.W[:, i].conj()) for i in range(self.n)])
            self.scale = np.max(np.abs(self.raw_coefficients), axis=(1, 2))
            if np.any(self.scale == 0):
                raise GeneralQuotientFailure('IDENTICALLY_ZERO_BILINEAR_EQUATION')
            self.coefficients = self.raw_coefficients/self.scale[:, None, None]
        self.evaluations = self.jacobian_evaluations = 0

    def evaluate(self, x, jac=False):
        x = np.asarray(x, dtype=complex)
        if x.shape != (self.n,):
            raise ValueError('Wrong complex variable dimension')
        c, d = np.r_[1., x[:self.s]], np.r_[1., x[self.s:]]
        self.evaluations += 1
        f = np.einsum('a,b,iab->i', c, d, self.coefficients)
        if not jac:
            return f
        self.jacobian_evaluations += 1
        J = np.column_stack([
            np.einsum('b,ib->i', d, self.coefficients[:, j, :]) for j in range(1, self.s+1)]
            + [np.einsum('a,ia->i', c, self.coefficients[:, :, j]) for j in range(1, self.s+1)])
        return f, J


def _multiply_linear(polynomial, linear):
    result = {}
    for exponent, value in polynomial.items():
        for variable, coefficient in enumerate(linear):
            if variable == 0:
                key = exponent
            else:
                changed = list(exponent)
                changed[variable-1] += 1
                key = tuple(changed)
            result[key] = result.get(key, 0)+value*coefficient
    return result


def maximal_minors(coefficients):
    n, width, _ = coefficients.shape
    s = width-1
    basis = [exponent for degree in range(s+1) for exponent in homogeneous_monomials(s, degree)]
    leading = homogeneous_monomials(s, s+1)
    monomials = basis+leading
    output = []
    permutations = [(permutation, (-1)**sum(permutation[i] > permutation[j]
                     for i in range(width) for j in range(i+1, width)))
                    for permutation in itertools.permutations(range(width))]
    for selected in itertools.combinations(range(n), width):
        polynomial = {}
        for permutation, sign in permutations:
            term = {(0,)*s: complex(sign)}
            for row, column in zip(selected, permutation):
                term = _multiply_linear(term, coefficients[row, :, column])
            for exponent, value in term.items():
                polynomial[exponent] = polynomial.get(exponent, 0)+value
        output.append([polynomial.get(exponent, 0) for exponent in monomials])
    return np.asarray(output, dtype=complex), basis, leading


def multiplication_matrices(minors, basis, leading, condition_limit=1e12):
    dimension = len(basis)
    scale = np.max(np.abs(minors), axis=1)
    if np.any(scale == 0):
        raise GeneralQuotientFailure('ZERO_MAXIMAL_MINOR')
    normalized = minors/scale[:, None]
    T, R = normalized[:, dimension:], normalized[:, :dimension]
    if T.shape[0] != T.shape[1]:
        raise GeneralQuotientFailure('NON_SQUARE_LEADING_BLOCK', dict(shape=T.shape))
    condition = float(np.linalg.cond(T))
    diagnostic = dict(leading_block_condition=condition, leading_shape=T.shape,
                      quotient_dimension=dimension, maximal_minors=minors)
    if not np.isfinite(condition) or condition > condition_limit:
        raise GeneralQuotientFailure('UNRESOLVED_LEADING_BLOCK', diagnostic)
    reduction = np.linalg.solve(T, -R)
    basis_index, leading_index = {e:i for i,e in enumerate(basis)}, {e:i for i,e in enumerate(leading)}
    matrices = []
    for variable in range(len(basis[0])):
        matrix = np.zeros((dimension, dimension), dtype=complex)
        for column, exponent in enumerate(basis):
            changed = list(exponent)
            changed[variable] += 1
            changed = tuple(changed)
            if changed in basis_index:
                matrix[basis_index[changed], column] = 1.
            else:
                matrix[:, column] = reduction[leading_index[changed]]
        matrices.append(matrix)
    commutators = [float(np.linalg.norm(a@b-b@a)/max(np.linalg.norm(a)*np.linalg.norm(b), 1e-300))
                   for a,b in itertools.combinations(matrices, 2)]
    diagnostic.update(reduction=reduction, multiplication_matrices=matrices,
                      pairwise_relative_commutators=commutators,
                      maximum_relative_commutator=max(commutators, default=0.))
    return matrices, diagnostic


def enumerate_general_quotient(kernel, target, tolerance=2e-10, max_polish=5,
                               separating=(1., .371+.193j, -.227+.419j)):
    begin = time.perf_counter()
    system = GeneralBilinearDualSystem(kernel, target)
    s, n = system.s, system.n
    if s != 3:
        raise GeneralQuotientFailure('THIS_FROZEN_CAMPAIGN_REQUIRES_S3_N6')
    expected = math.comb(n, s)
    paths, roots = [], []
    counts = dict(maximal_minors=math.comb(n, s+1), determinant_permutation_terms=math.comb(n, s+1)*math.factorial(s+1),
                  polynomial_linear_multiplications=math.comb(n, s+1)*math.factorial(s+1)*(s+1),
                  quotient_linear_solves=0, quotient_eigenproblems=0, quotient_recovery_least_squares=0,
                  quotient_polish_linear_solves=0, quotient_condition_SVDs=0)
    with kernel.model.operation('GENERAL_QUOTIENT_ALGEBRA', maximal_minors=counts['maximal_minors'],
                                determinant_permutation_terms=counts['determinant_permutation_terms'],
                                polynomial_linear_multiplications=counts['polynomial_linear_multiplications']):
        minors, basis, leading = maximal_minors(system.coefficients)
        counts['quotient_condition_SVDs'] += 1
        matrices, diagnostic = multiplication_matrices(minors, basis, leading)
        counts['quotient_linear_solves'] += 1
        if diagnostic['maximum_relative_commutator'] > 1e-7:
            raise GeneralQuotientFailure('NONCOMMUTING_MULTIPLICATION_MATRICES', diagnostic)
        linear = sum(coefficient*matrix for coefficient,matrix in zip(separating, matrices))
        eigenvalues, eigenvectors = np.linalg.eig(linear.T)
        counts['quotient_eigenproblems'] += 1
        counts['quotient_condition_SVDs'] += 1
        eigen_condition = float(np.linalg.cond(eigenvectors))
        eigen_gap = min(abs(a-b)/(1+abs(a)+abs(b)) for a,b in itertools.combinations(eigenvalues, 2))
        diagnostic.update(basis=basis, leading_monomials=leading, separating=separating,
                          eigenvector_condition=eigen_condition,
                          minimum_relative_separating_eigengap=float(eigen_gap))
        if eigen_condition > 1e12 or eigen_gap <= 1e-10:
            raise GeneralQuotientFailure('NONSEPARATING_OR_ILL_CONDITIONED_EIGENPROBLEM', diagnostic)
        for index, vector in enumerate(eigenvectors.T):
            c = np.asarray([np.vdot(vector, matrix.T@vector)/np.vdot(vector, vector) for matrix in matrices])
            A = np.einsum('a,iab->ib', np.r_[1., c], system.coefficients)
            d, _, rank, singular = np.linalg.lstsq(A[:, 1:], -A[:, 0], rcond=1e-12)
            counts['quotient_recovery_least_squares'] += 1
            recovery_condition = float(singular[0]/max(singular[-1], 1e-300))
            x, polished = np.r_[c, d], 0
            for iteration in range(max_polish):
                residual, J = system.evaluate(x, jac=True)
                if np.linalg.norm(residual, ord=np.inf) <= tolerance:
                    break
                try:
                    counts['quotient_polish_linear_solves'] += 1
                    x += np.linalg.solve(J, -residual)
                    polished += 1
                except np.linalg.LinAlgError:
                    break
            residual, J = system.evaluate(x, jac=True)
            error, condition = float(np.linalg.norm(residual, ord=np.inf)), float(np.linalg.cond(J))
            counts['quotient_condition_SVDs'] += 1
            status = 'RESOLVED_FINITE_NUMERICAL' if rank == s and recovery_condition < 1e12 and error <= tolerance and condition < 1e12 else 'UNRESOLVED_GENERAL_QUOTIENT_ENDPOINT'
            paths.append(dict(path=index, x=x.copy(), status=status, endpoint_residual=error,
                              endpoint_condition=condition, recovery_rank=int(rank), recovery_condition=recovery_condition,
                              polish_steps=polished, separating_eigenvalue=eigenvalues[index]))
            if status == 'RESOLVED_FINITE_NUMERICAL':
                roots.append(x.copy())
    distinct = []
    for root in roots:
        if all(np.linalg.norm(root-old) > 1e-6*(1+np.linalg.norm(root)) for old in distinct):
            distinct.append(root)
    counts.update(polynomial_evaluations=system.evaluations, polynomial_J_evaluations=system.jacobian_evaluations)
    with kernel.model.operation('GENERAL_QUOTIENT_ENDPOINT_SUMMARY', **{key:value for key,value in counts.items()
                                if key not in ('maximal_minors','determinant_permutation_terms','polynomial_linear_multiplications')}):
        pass
    return dict(method='GENERAL_BALANCED_QUOTIENT_EIGENROOTS', s=s, active_materials=n,
                expected_degree=expected, paths=paths, roots=np.asarray(roots).reshape(-1, n),
                resolved_paths=len(roots), unresolved_paths=expected-len(roots),
                distinct_finite_endpoints=len(distinct),
                all_paths_numerically_resolved=len(roots)==len(distinct)==expected,
                diagnostic=diagnostic, system=system, counters=counts,
                wall_seconds=time.perf_counter()-begin, certificate_flag='NOT_RIGOROUSLY_CERTIFIED')


def filter_general_roots(result, box=(.5,30.), conjugacy_tolerance=1e-6,
                         material_imaginary_tolerance=1e-7, node_floor=1e-10,
                         denominator_floor=1e-10, constitutive_tolerance=1e-8,
                         merge_rms=1e-7):
    system = result['system']
    kernel, model, s = system.kernel, system.kernel.model, system.s
    rows, candidates = [], []
    for path in result['paths']:
        if path['status'] != 'RESOLVED_FINITE_NUMERICAL':
            continue
        x = path['x']
        row = dict(path=path['path'], x=x.copy(), status='REJECTED')
        conjugacy = float(np.linalg.norm(x[s:]-x[:s].conj())/(1+np.linalg.norm(x[:s])))
        row['conjugacy_relative'] = conjugacy
        if conjugacy > conjugacy_tolerance:
            row['reason'] = 'NONPHYSICAL_CONJUGACY'
        else:
            w = system.w0+x[:s]@system.N
            v = kernel.L.reshape(kernel.r)+w@kernel.K
            row['minimum_node_abs'] = float(np.min(np.abs(v)))
            if np.min(np.abs(v)) <= node_floor:
                row['reason'] = 'UNRESOLVED_DUAL_NODE'
            else:
                ratio = w/v
                delta = np.real(w*v.conj())/np.abs(v)**2
                error = np.linalg.norm(w-delta*v)/max(np.linalg.norm(w), 1e-300)
                row.update(maximum_material_imaginary=float(np.max(np.abs(ratio.imag))), dual_constitutive_relative=float(error))
                if np.max(np.abs(ratio.imag)) > material_imaginary_tolerance or error > constitutive_tolerance:
                    row['reason'] = 'NONREAL_MATERIAL'
                else:
                    m = kernel.m0+kernel.U@delta
                    row['m'] = m
                    with model.operation('GENERAL_POLE_FILTER', algebraic_determinants=1):
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
        rows.append(row)
    return dict(roots=rows, candidate_count=len(candidates), certificate_flag='NOT_RIGOROUSLY_CERTIFIED')
