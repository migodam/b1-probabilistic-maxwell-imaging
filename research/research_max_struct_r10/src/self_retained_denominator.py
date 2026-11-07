"""DEV control: exact original Green diagonal + projected nonlocal remainder.

Separate file preserves the already frozen primary candidate implementation.
"""
import numpy as np
from scipy.linalg import lu_factor,lu_solve
from research.research_max_struct_r10.src.lowrank_denominator import LowRankDenominator


class SelfRetainedDenominator(LowRankDenominator):
    def __init__(self,G,U,projection='galerkin',cost=None):
        super().__init__(G,U,projection,cost)
        self.diagonal=np.diag(G)-np.sum(self.B*self.U.conj(),axis=1)
        self.projection+=' plus exact original self diagonal'

    def action(self,X):
        return super().action(X)+self.diagonal[:,None]*X

    def factor(self,m,beta):
        a=1+(beta-self.diagonal)*m
        F=(m/a)[:,None]*self.B
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
            return (rhs+self.U@v)/a.conj()[:,None]
        b=rhs/a[:,None]
        return b+F@lu_solve(lu,self.U.conj().T@b,check_finite=False)
