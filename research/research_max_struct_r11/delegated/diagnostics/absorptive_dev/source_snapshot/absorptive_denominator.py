"""Root-derived fixed absorptive auxiliary denominator; original Maxwell physics unchanged."""
from pathlib import Path
from collections import Counter
import json,time
import numpy as np
from scipy.linalg import lu_factor,lu_solve
from research.research_max_struct_r10.src.lowrank_denominator import LowRankDenominator,prepare_bases,matrix_material_readout

class AbsorptiveDenominator(LowRankDenominator):
    """A=I+c D_m-D_m K, c=-i beta; residual prefactor remains real beta."""
    def factor(self,m,beta):
        c=-1j*float(beta);a=1+c*m
        self.cost['absorptive_diagonal_builds']+=1
        if self.rank==0:
            self.cost['denominator_scalar_factorizations']+=1
            return a,None,None
        F=(m/a)[:,None]*self.B
        T=np.eye(self.rank,dtype=complex)-self.U.conj().T@F
        lu=lu_factor(T,check_finite=False)
        self.cost['denominator_small_LU']+=1
        self.cost['denominator_small_LU_dimension_sum']+=self.rank
        return a,F,lu
    def solve(self,factor,rhs,adjoint=False):
        a,F,lu=factor
        if self.rank==0:
            self.cost['denominator_diagonal_complex_RHS']+=rhs.shape[1]
            return rhs/(a.conj() if adjoint else a)[:,None]
        self.cost['denominator_small_triangular_complex_RHS']+=rhs.shape[1]
        if adjoint:
            v=lu_solve(lu,F.conj().T@rhs,trans=2,check_finite=False)
            return (rhs+self.U@v)/a.conj()[:,None]
        b=rhs/a[:,None]
        return b+F@lu_solve(lu,self.U.conj().T@b,check_finite=False)
    def action(self,X):
        if self.rank==0:
            self.cost['zero_K_actions']+=1
            self.cost['zero_K_complex_RHS']+=X.shape[1]
            return np.zeros_like(X)
        return super().action(X)
    def value_gradient(self,J,Q,m,beta,Os,rho,reg):
        beta=float(beta);c=-1j*beta
        fac=self.factor(m,beta);w=self.solve(fac,J-m[:,None]*Q);weighted=w@Os
        a=self.solve(fac,weighted,adjoint=True);h=Q+c*w-self.action(w)
        prior,gprior=reg(m)
        value=.5*rho*beta**2*float(np.vdot(w,weighted).real)+float(prior)
        grad=-rho*beta**2*np.real(np.sum(a.conj()*h,axis=1))+gprior
        self.cost['material_inner_evaluations']+=1;self.cost['material_prior_evaluations']+=1;self.cost['source_metric_products']+=1
        return value,grad

def _prepare(model,config,green):
    started=time.process_time()
    if green:
        bases,info=prepare_bases(model.G,model.S,.1);U=bases['green_strong']
    else:U=np.empty((model.n,0),complex);info=dict(setup_G_SVD=0,setup_S_SVD=0,setup_projected_basis_SVD=0,basis_source='K=0 scalar absorptive control; no nonlocal basis')
    op=AbsorptiveDenominator(model.G,U,'galerkin');out=Path(config['receipt_path']).parent
    info.update(study='R11_ABSORPTIVE_EXPOSED_DEV',method='absorptiveGreenBothBlocks' if green else 'absorptiveScalarBothBlocks',rank=op.rank,projection='galerkin' if green else 'zero_K',shift_c='-1j*beta',outer_prefactor='positive real beta',matrix_inner_maxiter=8,
       seed='Unchanged real-CIE scalar seed proposal, accepted only if it lowers the new absorptive material objective',material_rule='Full moving absorptive gradient plus original real-material prior; not mature scalar CIE',truth_access='NONE',setup_complete_cpu_seconds=time.process_time()-started)
    (out/'DENOMINATOR_SETUP.json').write_text(json.dumps(info,indent=2)+'\n');np.savez_compressed(out/'DENOMINATOR_BASIS.npz',U=op.U,B=op.B)
    return dict(config,matrix_denominator=op,matrix_inner_maxiter=8)
def prepare_absorptive_green(model,config):return _prepare(model,config,True)
def prepare_absorptive_scalar(model,config):return _prepare(model,config,False)
absorptive_material_readout=matrix_material_readout
