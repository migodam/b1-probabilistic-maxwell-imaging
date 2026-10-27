"""Synthetic matrix checks only; these do not qualify a Maxwell solver."""
import unittest
import numpy as np
from scipy.sparse import csc_matrix,csr_matrix
from bgga.operators import FixedSpaceCache,goal_galerkin

class SyntheticCounter:
    def __init__(self):self.events=[]
    def reserve(self,*args,**kw):self.events.append(('synthetic_not_PDE',args,kw))
    def dense(self,*args,**kw):self.events.append(('synthetic_dense',args,kw))

class Operators(unittest.TestCase):
    def test_synthetic_fixed_source_and_material_rhs(self):
        rng=np.random.default_rng(271027)
        n=23;U=np.arange(6,13)
        A=7*np.eye(n)+.1*(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
        bw=rng.normal(size=(n,2))+1j*rng.normal(size=(n,2))
        R=csr_matrix(rng.normal(size=(5,n))+1j*rng.normal(size=(5,n)))
        direct=rng.normal(size=(2,5))+1j*rng.normal(size=(2,5))
        V=np.zeros((n,n),complex);V[U,U]=np.linspace(.1,1.5,len(U))+.03j
        bx=np.zeros_like(bw);bx[U]=rng.normal(size=(len(U),2))+.7j
        V=csc_matrix(V);book=SyntheticCounter()
        cache=FixedSpaceCache(csc_matrix(A),bw,R,direct,U,book);cache.build_ports()
        u,y,d=cache.full(V,bx)
        for fn in (cache.woodbury,cache.schur,cache.schur_from_paid_ports,cache.krylov):
            v,z,dd=fn(V,bx)
            self.assertLess(np.linalg.norm(v-u)/np.linalg.norm(u),1e-9)
            self.assertLess(np.linalg.norm(z-y)/np.linalg.norm(y),1e-9)
            self.assertLess(dd['backward'],1e-10)
        # Dropping target-dependent scattered-field forcing must fail.
        wrong=np.linalg.solve(A-V.toarray(),bw)
        self.assertGreater(np.linalg.norm(wrong-u),.01)
        q=np.eye(n,dtype=complex)
        v,z,dd=goal_galerkin(csc_matrix(A)-V,bw+bx,R,direct,q,book)
        self.assertLess(np.linalg.norm(v-u),1e-10)
    def test_public_ports_refuse_large_allocation(self):
        n=8;c=FixedSpaceCache(csc_matrix(np.eye(n)),np.ones((n,1)),csr_matrix(np.eye(n)),np.zeros((1,n)),np.arange(n),SyntheticCounter())
        with self.assertRaises(MemoryError):c.build_ports(max_bytes=1)

if __name__=='__main__':unittest.main()
