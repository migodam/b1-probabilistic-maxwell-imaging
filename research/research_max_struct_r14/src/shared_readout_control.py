"""Same finite chi-coordinate material solver for the four frozen metric pairs.

Only the two material call sites in the inherited step receive private bindings.
Current/PR/line-search/stage/cost orchestration remains the original core code.
This exposed four-case control changes no main or independent-cohort algorithm.
"""
import time
import numpy as np
from scipy.optimize import minimize
from research.research_max_struct_r13.src.full_current import _local_function
from research.research_max_struct_r4.src.cie_onehop import onehop_material_value_gradient
from dualmetric_core import StatefulDualMetric


SOLVER_TAG='shared_LBFGS60'


def shared_scalar_value_gradient(Jw,Qw,m,beta,rho,reg,cost):
    """Full real derivative of P0(m)(J-mQ), with source-whitened Jw/Qw."""
    d=1+float(beta)*m
    residual=float(beta)*(Jw-m[:,None]*Qw)/d[:,None]
    derivative=-float(beta)*(Qw+float(beta)*Jw)/d[:,None]**2
    prior,gprior=reg(m)
    cost['material_inner_evaluations']+=1
    cost['material_prior_evaluations']+=1
    cost['shared_scalar_value_gradient_evaluations']+=1
    cost['scalar_complex_divisions']+=residual.size+derivative.size
    value=.5*rho*float(np.vdot(residual,residual).real)+float(prior)
    gradient=rho*np.real(np.sum(residual.conj()*derivative,axis=1))+np.asarray(gprior,float)
    return value,gradient


def _common_minimize(fun,m,bounds,config,cost):
    """Same m start, chi coordinates, bounds, tolerances and old-m safeguard."""
    start_cpu=time.process_time();start_evaluations=cost['material_inner_evaluations']
    old=np.asarray(m,float).copy();old_value=float(fun(old)[0])
    if not np.isfinite(old_value):
        raise FloatingPointError('Nonfinite material objective at the inherited legal starting material')
    lower,upper=bounds
    result=minimize(fun,old,jac=True,method='L-BFGS-B',bounds=list(zip(lower,upper)),
        options=dict(maxiter=60,ftol=1e-12,gtol=1e-8,maxls=25))
    candidate=np.asarray(result.x,float)
    finite=bool(np.all(np.isfinite(candidate)))
    final_value=float(fun(candidate)[0]) if finite else np.inf
    fallback=not np.isfinite(final_value) or final_value>old_value
    chosen=old.copy() if fallback else candidate.copy()
    accepted_value=old_value if fallback else final_value
    cost['material_readouts']+=1
    cost['shared_material_LBFGS_readouts']+=1
    cost['material_inner_LBFGS_iterations']+=int(result.nit)
    cost['material_inner_LBFGS_evaluations']+=int(result.nfev)
    cost['material_inner_fallbacks']+=int(fallback)
    return chosen,dict(material_solver=SOLVER_TAG,material_solver_coordinates='original real chi',
        material_inner_start='current material m, no scalar seed',material_seed_used=False,
        material_inner_nit=int(result.nit),material_inner_nfev=int(result.nfev),
        material_inner_total_evaluations=cost['material_inner_evaluations']-start_evaluations,
        material_inner_success=bool(result.success),material_inner_message=str(result.message),
        material_old_block_loss=old_value,material_final_block_loss=accepted_value,
        material_candidate_block_loss=final_value if np.isfinite(final_value) else None,
        material_fallback=bool(fallback),material_readout_CPU_seconds=time.process_time()-start_cpu)


def common_onehop_readout(blocks,J,Q,m,beta,Os,Ls,rho,bounds,reg,config,cost):
    """M1 full existing derivative; no seed, Arc, MM or frozen metric."""
    def objective(z):
        return onehop_material_value_gradient(blocks,J,Q,z,beta,Os,rho,reg,cost)
    return _common_minimize(objective,m,bounds,config,cost)


def common_scalar_readout(J,Q,m,beta,weights,rho,bounds,reg,has_reg,config,cost):
    """Core global branch already supplies J@Ls/Q@Ls and unit source weights."""
    if not np.array_equal(np.asarray(weights),np.ones(J.shape[1])):
        raise ValueError('Shared scalar control requires the inherited source-whitened global branch')
    def objective(z):
        return shared_scalar_value_gradient(J,Q,z,beta,rho,reg,cost)
    return _common_minimize(objective,m,bounds,config,cost)


class SharedMaterialLBFGS(StatefulDualMetric):
    """Same inherited complete state/codec, private material-function bindings."""
    def _bind(self,model,Y,config):
        cfg=dict(config)
        if cfg.get('sourcewise_cg',False):
            raise ValueError('The frozen shared-readout mechanism control is global PR only')
        cfg.update(material_solver=SOLVER_TAG,material_solver_coordinates='real_chi',
            material_inner_maxiter=60,material_inner_ftol=1e-12,material_inner_gtol=1e-8,
            material_inner_maxls=25,material_seed_used=False,
            shared_material_old_m_nonincrease_safeguard=True)
        super()._bind(model,Y,cfg)

    step=_local_function(StatefulDualMetric.step,
        mobius_material_readout=common_onehop_readout,
        _baseline_readout=common_scalar_readout)
