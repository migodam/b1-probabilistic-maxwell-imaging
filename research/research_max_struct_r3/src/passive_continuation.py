"""Fixed passive auxiliary medium m+i*gamma over the disclosed ROI.

Reuses frozen disk-pulse G/E/S and its counted retained-LU primitives.
This is not a data transformation: the target remains the measured original Y.
"""
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from research.research_fbg_r3_experiments.src.physics import MaxwellModel

class PassiveModel(MaxwellModel):
    def __init__(self, public, gamma):
        super().__init__(public)
        self.gamma=float(gamma)
    def state(self,m):
        m=np.asarray(m,float).reshape(self.n)
        if self._current is not None and np.array_equal(m,self._current['m']):
            return self._current
        c=m+1j*self.gamma
        with self._physics.operation('FULL_FORWARD',rhs=self.ns,LU=1,source_RHS=self.ns,full_factorization_dimension=self.n):
            L=self._physics.identity-c[:,None]*self.G
            rhs=c[:,None]*self.E
            lu=lu_factor(L,check_finite=False)
            j=lu_solve(lu,rhs,check_finite=False)
            Q=self.E+self.G@j; Y=self.S@j
            st=dict(m=m.copy(),contrast=c,j=j,Q=Q,E=Q,Y=Y,y=Y,L=L,lu=lu,
                state_relative_residual=float(np.linalg.norm(L@j-rhs)/np.linalg.norm(rhs)))
        self._current=st
        self._physics.state=st;self._physics.factor=lu
        return st
