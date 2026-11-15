"""Parent-derived fixed twofold denominator and optional residual block hop.

Supplied U is fixed (normally joint_spectral U11); no basis/rank selection here.
No all-box contraction guarantee is claimed for the truncated R correction.
The unchanged R10 matrix_material_readout accepts this operator through config.
"""
from collections import Counter
import time
import numpy as np
from research.research_max_struct_r10.src.lowrank_denominator import (
    LowRankDenominator, matrix_material_readout)
from research.research_max_struct_r4.src.cie_onehop import OneHopAbsorption


class TwofoldDenominator(LowRankDenominator):
    """K=PGP, R=N-PNP, A=I+beta Dm-DmK; inherited low-rank solves."""
    @property
    def cost(self):
        return self._cost

    @cost.setter
    def cost(self, value):
        self._cost=value
        if getattr(self,'blocks',None) is not None:self.blocks.cost=value

    def __init__(self,G,U,shape,*,one_hop=False,cost=None):
        start=time.process_time();wall=time.perf_counter()
        counter=Counter() if cost is None else cost
        U=np.asarray(U,complex)
        if U.ndim!=2 or U.shape[0]!=np.asarray(G).shape[0]:
            raise ValueError('Supplied orthonormal U must have N rows (rank0 permitted)')
        super().__init__(G,U,projection='galerkin',cost=counter)
        self.one_hop=bool(one_hop)
        counter['twofold_denominator_setups']+=1
        counter['denominator_setup_G_times_U_complex_columns']+=self.rank
        counter['denominator_setup_basis_products']+=2
        self.blocks=None
        if self.one_hop:
            self.blocks=OneHopAbsorption(G,shape,counter)
            if self.rank:
                NU=self.blocks.off_action(self.U)
                self.retained_neighbor=self.U.conj().T@NU
                counter['residual_neighbor_setup_basis_products']+=1
            else:self.retained_neighbor=np.empty((0,0),complex)
        self.setup_receipt=dict(rank=self.rank,projection=self.projection,one_hop=self.one_hop,
            CPU_seconds=time.process_time()-start,wall_seconds=time.perf_counter()-wall,
            basis_selection='external supplied fixed U; no SVD or padding',
            neighbor='Original disjoint4cell off-diagonal N; residual R=N-PNP',
            qualification='No all-box contraction guarantee for the truncated residual hop')

    def factor(self,m,beta):
        self.cost['denominator_factor_requests']+=1
        if self.rank:return super().factor(m,beta)
        self.cost['denominator_diagonal_factors']+=1
        return 1+beta*np.asarray(m),None,None

    def solve(self,factor,rhs,adjoint=False):
        self.cost['denominator_solve_actions']+=1
        self.cost['denominator_adjoint_solve_actions' if adjoint else 'denominator_forward_solve_actions']+=1
        self.cost['denominator_solve_input_complex_RHS']+=rhs.shape[1]
        if self.rank:return super().solve(factor,rhs,adjoint=adjoint)
        self.cost['denominator_diagonal_complex_divisions']+=rhs.size
        a=factor[0]
        return rhs/(a.conj() if adjoint else a)[:,None]

    def residual_action(self,rhs,adjoint=False):
        if not self.one_hop:raise ValueError('Residual-neighbor action requires one_hop=True')
        self.cost['residual_neighbor_adjoint_actions' if adjoint else 'residual_neighbor_forward_actions']+=1
        self.cost['residual_neighbor_input_complex_RHS']+=rhs.shape[1]
        result=self.blocks.off_action(rhs,adjoint=adjoint)
        if self.rank:
            C=self.retained_neighbor.conj().T if adjoint else self.retained_neighbor
            result-=self.U@(C@(self.U.conj().T@rhs))
            self.cost['residual_neighbor_projection_basis_products']+=2
            self.cost['residual_neighbor_small_matrix_products']+=1
        return result

    def transform(self,m,beta,rhs,adjoint=False):
        """beta [A^-1 + A^-1 Dm R A^-1], or its exact adjoint."""
        rhs=np.asarray(rhs);vector=rhs.ndim==1
        if vector:rhs=rhs[:,None]
        self.cost['denominator_transform_adjoint_calls' if adjoint else 'denominator_transform_forward_calls']+=1
        factor=self.factor(m,beta)
        v=self.solve(factor,rhs,adjoint=adjoint)
        if self.one_hop:
            if adjoint:
                extra=self.solve(factor,self.residual_action(np.asarray(m)[:,None]*v,adjoint=True),adjoint=True)
            else:
                extra=self.solve(factor,np.asarray(m)[:,None]*self.residual_action(v))
            result=beta*(v+extra)
        else:result=beta*v
        return result[:,0] if vector else result

    def adjoint(self,m,beta,rhs):
        return self.transform(m,beta,rhs,adjoint=True)

    transform_adjoint=adjoint

    def value_gradient(self,J,Q,m,beta,Os,rho,reg):
        if not self.one_hop:
            # Preserve the original helper's arithmetic and acceptance inputs.
            return super().value_gradient(J,Q,m,beta,Os,rho,reg)
        m=np.asarray(m);factor=self.factor(m,beta)
        v=self.solve(factor,J-m[:,None]*Q)
        Rv=self.residual_action(v)
        w=self.solve(factor,m[:,None]*Rv)
        z=v+w;b=z@Os
        p=self.solve(factor,b,adjoint=True)
        q=self.solve(factor,b+self.residual_action(m[:,None]*p,adjoint=True),adjoint=True)
        H=Q+beta*v-self.action(v)
        T=Rv-beta*w+self.action(w)
        prior,gprior=reg(m)
        value=.5*rho*beta**2*float(np.vdot(z,b).real)+float(prior)
        gradient=rho*beta**2*np.real(np.sum(-q.conj()*H+p.conj()*T,axis=1))+gprior
        self.cost['material_inner_evaluations']+=1
        self.cost['material_prior_evaluations']+=1
        self.cost['source_metric_products']+=1
        return value,gradient


# Explicit alias only: seed, old-value comparison, LBFGS8, bounds and prior unchanged.
twofold_material_readout=matrix_material_readout
