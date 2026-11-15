"""R10 projected-Green matrix denominator; fixed-current, real-material only.

K=B U*, A=I+beta D_m-D_m K.  No reduction of material unknowns.
This changes the auxiliary material objective, never the Maxwell forward.
"""
from collections import Counter
import time
import numpy as np
from scipy.linalg import svd, lu_factor, lu_solve
from scipy.optimize import minimize
from pathlib import Path
import json
from research.research_max_struct_r2.src.cie_global import _baseline_readout


def prepare_bases(G, S, threshold=.1):
    """Background-only Two-fold range; thresholds declared before DEV runs."""
    started=time.process_time()
    _,sg,vhg=svd(G,full_matrices=False,check_finite=False)
    _,so,vho=svd(S,full_matrices=False,check_finite=False)
    Vg=vhg[sg>=threshold*sg[0]].conj().T
    Vo=vho[so>=threshold*so[0]].conj().T
    residual=Vg-Vo@(Vo.conj().T@Vg)
    u,s,_=svd(residual,full_matrices=False,check_finite=False)
    extra=u[:,s>1e-10]
    twofold=np.column_stack((Vo,extra))
    r=twofold.shape[1]
    # Equal-rank spectral controls; both fixed without material truth.
    bases=dict(twofold=twofold,green_equal_rank=vhg[:r].conj().T,green_strong=Vg)
    return bases,dict(threshold=threshold,green_strong_rank=Vg.shape[1],
        receiver_strong_rank=Vo.shape[1],twofold_rank=r,
        basis_source='fixed public G and S; right singular vectors; no GT',
        material_dof=G.shape[0],setup_cpu_seconds=time.process_time()-started,
        setup_G_SVD=1,setup_S_SVD=1,setup_projected_basis_SVD=1)


class LowRankDenominator:
    def __init__(self,G,U,projection='galerkin',cost=None):
        self.U=np.asarray(U,complex)
        self.B=G@self.U
        if projection=='galerkin':self.B=self.U@(self.U.conj().T@self.B)
        elif projection!='right':raise ValueError(projection)
        self.projection=projection
        self.n,self.rank=self.U.shape
        self.cost=Counter() if cost is None else cost

    def action(self,X):
        self.cost['lowrank_Green_actions']+=1
        self.cost['lowrank_Green_complex_RHS']+=X.shape[1]
        return self.B@(self.U.conj().T@X)

    def factor(self,m,beta):
        a=1+beta*m
        z=m/a
        F=z[:,None]*self.B
        T=np.eye(self.rank)-self.U.conj().T@F
        lu=lu_factor(T,check_finite=False)
        self.cost['denominator_small_LU']+=1
        self.cost['denominator_small_LU_dimension_sum']+=self.rank
        return a,F,lu

    def solve(self,factor,rhs,adjoint=False):
        a,F,lu=factor
        self.cost['denominator_small_triangular_complex_RHS']+=rhs.shape[1]
        if adjoint:
            v=lu_solve(lu,F.conj().T@rhs,trans=2,check_finite=False)
            return (rhs+self.U@v)/a[:,None]
        b=rhs/a[:,None]
        return b+F@lu_solve(lu,self.U.conj().T@b,check_finite=False)

    def value_gradient(self,J,Q,m,beta,Os,rho,reg):
        factor=self.factor(m,beta)
        w=self.solve(factor,J-m[:,None]*Q)
        weighted=w@Os
        a=self.solve(factor,weighted,adjoint=True)
        h=Q+beta*w-self.action(w)
        prior,gprior=reg(m)
        value=.5*rho*beta**2*float(np.vdot(w,weighted).real)+float(prior)
        gradient=-rho*beta**2*np.real(np.sum(a.conj()*h,axis=1))+gprior
        self.cost['material_inner_evaluations']+=1
        self.cost['material_prior_evaluations']+=1
        self.cost['source_metric_products']+=1
        return value,gradient


def matrix_material_readout(blocks,J,Q,m,beta,Os,Ls,rho,bounds,reg,config,cost):
    """Exact derivative, bounded LBFGS; no four-color independence assumption."""
    started=time.process_time();before=cost['material_inner_evaluations']
    operator=config['matrix_denominator'];operator.cost=cost
    seed,_=_baseline_readout(J@Ls,Q@Ls,m,beta,np.ones(J.shape[1]),rho,
                             bounds,reg,False,config,cost)
    def fun(x):return operator.value_gradient(J,Q,x,beta,Os,rho,reg)
    old,_=fun(m);seed_value,_=fun(seed)
    start=seed if seed_value<old else m
    opt=minimize(fun,start,jac=True,method='L-BFGS-B',bounds=list(zip(*bounds)),
        options=dict(maxiter=int(config.get('matrix_inner_maxiter',8)),
                     maxls=20,ftol=1e-12,gtol=float(config.get('material_inner_gtol',1e-8))))
    best_value=min(old,seed_value)
    accepted=bool(np.isfinite(opt.fun) and opt.fun<=best_value)
    result=opt.x.copy() if accepted else start.copy()
    final=float(opt.fun) if accepted else best_value
    cost['matrix_material_readouts']+=1
    cost['matrix_material_LBFGS_iterations']+=int(opt.nit)
    cost['matrix_material_readout_cpu_microseconds']+=int(1e6*(time.process_time()-started))
    return result,dict(material_inner_nit=int(opt.nit),material_inner_nfev=int(opt.nfev),
        material_inner_total_evaluations=cost['material_inner_evaluations']-before,
        material_inner_success=bool(opt.success),material_inner_message=str(opt.message),
        material_old_block_loss=old,material_scalar_seed_block_loss=seed_value,
        material_final_block_loss=final,material_fallback=not accepted,
        matrix_rank=operator.rank,matrix_projection=operator.projection)


def _prepare(model,config,basis):
    started=time.process_time()
    bases,info=prepare_bases(model.G,model.S,.1)
    operator=LowRankDenominator(model.G,bases[basis],'galerkin')
    supplied=dict(config,matrix_denominator=operator,matrix_inner_maxiter=8)
    out=Path(config['receipt_path']).parent
    info.update(basis=basis,projection='galerkin',matrix_inner_maxiter=8,
                setup_complete_cpu_seconds=time.process_time()-started,
                status='FROZEN_BEFORE_FRESH_GT_SCORE',truth_access='NONE')
    (out/'DENOMINATOR_SETUP.json').write_text(json.dumps(info,indent=2)+'\n')
    np.savez_compressed(out/'DENOMINATOR_BASIS.npz',U=operator.U,B=operator.B)
    return supplied


def prepare_twofold(model,config):return _prepare(model,config,'twofold')
def prepare_green_strong(model,config):return _prepare(model,config,'green_strong')
