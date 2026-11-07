"""R11 fixed-reference Green material metric, over the unchanged scalar TM CIE.

The reference factor is frozen only within one material readout and refreshed
at the next alternating step. This is not MM for the material-dependent loss,
nor a convexity claim about Maxwell inversion. All real material DOF remain.
"""
from pathlib import Path
import json
import time
import numpy as np
from scipy.optimize import minimize
from research.research_max_struct_r2.src.cie_global import _baseline_readout
from research.research_max_struct_r10.src.lowrank_denominator import (
    LowRankDenominator, matrix_material_readout, prepare_green_strong)


def gradient_parts(operator, J, Q, m, beta, Os, rho, reg):
    """Counterfactual diagnostic at one common state; no physical forward."""
    factor = operator.factor(m, beta)
    w = operator.solve(factor, J-m[:, None]*Q)
    a = operator.solve(factor, w@Os, adjoint=True)
    t = beta*w-operator.action(w)
    constant = -rho*beta**2*np.real(np.sum(a.conj()*Q, axis=1))
    metric = -rho*beta**2*np.real(np.sum(a.conj()*t, axis=1))
    prior, gp = reg(m)
    return dict(constitutive_gradient=constant, metric_gradient=metric,
                prior_gradient=gp, full_gradient=constant+metric+gp,
                lagged_gradient=constant+gp, transported_defect=w,
                raw_defect=J-m[:, None]*Q,
                value=.5*rho*beta**2*float(np.vdot(w, w@Os).real)+float(prior))


def frozen_value_gradient(operator, factor, J, Q, m, beta, Os, rho, reg):
    w = operator.solve(factor, J-m[:, None]*Q)
    a = operator.solve(factor, w@Os, adjoint=True)
    prior, gp = reg(m)
    value = .5*rho*beta**2*float(np.vdot(w, w@Os).real)+float(prior)
    gradient = -rho*beta**2*np.real(np.sum(a.conj()*Q, axis=1))+gp
    operator.cost['material_inner_evaluations'] += 1
    operator.cost['material_prior_evaluations'] += 1
    operator.cost['source_metric_products'] += 1
    return value, gradient


def lagged_material_readout(blocks, J, Q, m, beta, Os, Ls, rho,
                            bounds, reg, config, cost):
    started = time.process_time()
    before = cost['material_inner_evaluations']
    operator = config['matrix_denominator']
    operator.cost = cost
    reference = m.copy()
    factor = operator.factor(reference, beta)
    seed, _ = _baseline_readout(J@Ls, Q@Ls, m, beta, np.ones(J.shape[1]),
                                rho, bounds, reg, False, config, cost)

    def fun(x):
        return frozen_value_gradient(operator, factor, J, Q, x, beta, Os, rho, reg)

    old, _ = fun(m)
    seed_value, _ = fun(seed)
    start = seed if seed_value < old else m
    opt = minimize(fun, start, jac=True, method='L-BFGS-B',
                   bounds=list(zip(*bounds)),
                   options=dict(maxiter=int(config.get('matrix_inner_maxiter', 8)),
                                maxls=20, ftol=1e-12,
                                gtol=float(config.get('material_inner_gtol', 1e-8))))
    best_value = min(old, seed_value)
    accepted = bool(np.isfinite(opt.fun) and opt.fun <= best_value)
    result = opt.x.copy() if accepted else start.copy()
    final = float(opt.fun) if accepted else best_value
    cost['matrix_material_readouts'] += 1
    cost['lagged_material_readouts'] += 1
    cost['matrix_material_LBFGS_iterations'] += int(opt.nit)
    cost['matrix_material_readout_cpu_microseconds'] += int(1e6*(time.process_time()-started))
    return result, dict(material_inner_nit=int(opt.nit), material_inner_nfev=int(opt.nfev),
        material_inner_total_evaluations=cost['material_inner_evaluations']-before,
        material_inner_success=bool(opt.success), material_inner_message=str(opt.message),
        material_old_block_loss=old, material_scalar_seed_block_loss=seed_value,
        material_final_block_loss=final, material_fallback=not accepted,
        matrix_rank=operator.rank, matrix_projection=operator.projection,
        denominator_reference='input material; factor frozen during this readout only',
        reported_loss='fixed-reference loss, not moving-denominator R10 objective',
        monotonicity_scope='fixed-reference material block only; not outer alternating or physical Phi')


def prepare_lagged(model, config):
    supplied = prepare_green_strong(model, config)
    path = Path(config['receipt_path']).parent/'DENOMINATOR_SETUP.json'
    info = json.loads(path.read_text())
    info.update(study='R11', derivative_policy='lagged at input material within each readout',
                inner_objective='convex quadratic plus inherited convex prior on shared real m',
                reference_refresh='every alternating material block; no outer MM guarantee')
    path.write_text(json.dumps(info, indent=2)+'\n')
    return supplied


def prepare_full(model, config):
    return prepare_green_strong(model, config)


full_material_readout = matrix_material_readout


def contracted_lagged_material_readout(blocks, J, Q, m, beta, Os, Ls, rho,
                                       bounds, reg, config, cost):
    """Freeze only T=I-D_z K, preserving ordinary CIE contraction z=m/(1+βm).

The data term is quadratic in z; the original R(m(z)) is retained and need
not be convex in z. At K=0 this is the ordinary CIE material objective.
"""
    started = time.process_time()
    before = cost['material_inner_evaluations']
    operator = config['matrix_denominator']
    operator.cost = cost
    factor = operator.factor(m, beta)
    reference_a = factor[0]
    V = Q+beta*J
    zold = m/(1+beta*m)
    seed, _ = _baseline_readout(J@Ls, Q@Ls, m, beta, np.ones(J.shape[1]),
                                rho, bounds, reg, False, config, cost)
    zseed = seed/(1+beta*seed)

    def fun(z):
        material = z/(1-beta*z)
        # A_ref = D_a_ref T_ref, so T_ref^-1 b=A_ref^-1 D_a_ref b.
        w = operator.solve(factor, reference_a[:, None]*(J-z[:, None]*V))
        a = reference_a[:, None]*operator.solve(factor, w@Os, adjoint=True)
        prior, gp = reg(material)
        value = .5*rho*beta**2*float(np.vdot(w, w@Os).real)+float(prior)
        gradient = -rho*beta**2*np.real(np.sum(a.conj()*V, axis=1))
        gradient += gp/(1-beta*z)**2
        cost['material_inner_evaluations'] += 1
        cost['material_prior_evaluations'] += 1
        cost['source_metric_products'] += 1
        return value, gradient

    old, _ = fun(zold)
    seed_value, _ = fun(zseed)
    start = zseed if seed_value < old else zold
    lo, hi = bounds
    zlo, zhi = lo/(1+beta*lo), hi/(1+beta*hi)
    opt = minimize(fun, start, jac=True, method='L-BFGS-B', bounds=list(zip(zlo, zhi)),
                   options=dict(maxiter=int(config.get('matrix_inner_maxiter', 8)),
                                maxls=20, ftol=1e-12,
                                gtol=float(config.get('material_inner_gtol', 1e-8))))
    best_value = min(old, seed_value)
    accepted = bool(np.isfinite(opt.fun) and opt.fun <= best_value)
    zresult = opt.x.copy() if accepted else start.copy()
    result = zresult/(1-beta*zresult)
    cost['matrix_material_readouts'] += 1
    cost['contracted_lagged_material_readouts'] += 1
    cost['matrix_material_LBFGS_iterations'] += int(opt.nit)
    cost['matrix_material_readout_cpu_microseconds'] += int(1e6*(time.process_time()-started))
    return result, dict(material_inner_nit=int(opt.nit), material_inner_nfev=int(opt.nfev),
        material_inner_total_evaluations=cost['material_inner_evaluations']-before,
        material_inner_success=bool(opt.success), material_inner_message=str(opt.message),
        material_old_block_loss=old, material_scalar_seed_block_loss=seed_value,
        material_final_block_loss=float(opt.fun) if accepted else best_value,
        material_fallback=not accepted, matrix_rank=operator.rank,
        matrix_projection=operator.projection,
        denominator_reference='T=I-D_zref K fixed; ordinary CIE scalar contraction remains active',
        reported_loss='fixed-T loss in z with exact inherited R(m(z)); not moving R10 objective',
        monotonicity_scope='fixed-T block only; prior may be nonconvex in z; no outer MM claim')


def prepare_contracted_lagged(model, config):
    supplied = prepare_green_strong(model, config)
    path = Path(config['receipt_path']).parent/'DENOMINATOR_SETUP.json'
    info = json.loads(path.read_text())
    info.update(study='R11', derivative_policy='freeze only nonlocal T; retain CIE z contraction',
                inner_objective='quadratic data term in z plus unchanged R(m(z))',
                reference_refresh='every material readout; no outer MM guarantee')
    path.write_text(json.dumps(info, indent=2)+'\n')
    return supplied
