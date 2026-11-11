"""Public-only real-isotropic adapter of the existing full-vector CM+RR DDA.

Material has N real entries; J/Q have 3*N cell-major x/y/z rows. This module
contains no target descriptions, truth paths or truth loading. It retains the
R8 source-normalized physical objective and quadratic prior exactly.
"""
from __future__ import annotations
from types import SimpleNamespace
import time
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from research.research_max_struct_r8.src.vector_pr3d import VectorPR3D


class DenseVectorModel:
    components = 3

    def __init__(self, side=6, width=.24, *, lower=0., upper=8.):
        self.legacy = VectorPR3D(ns=int(side), width=float(width), lossy=False)
        b = self.legacy
        self.G, self.E, self.S = b.G, b.inc, b.S
        self.n, self.ns, self.nr = b.N, b.P, b.S.shape[0]
        self.N, self.M = b.N, b.M
        self.side, self.width, self.h = int(side), float(width), b.h
        self.shape = (self.side,) * 3
        self.points = self.centers = b.points
        self.lower, self.upper = float(lower), float(upper)
        self.bounds = [(self.lower, self.upper)] * self.n
        self._current = None
        self.reg = SimpleNamespace(evaluate=lambda m: (*b.prior(m), None),smooth=None,tv=None,epsilon=None,
            name='R8_vector_quadratic_edge_mean_lam2e-4_plus_l2mean1e-6')
        self.metadata = dict(components=3, material_unknowns=self.n,
            complex_current_components=self.M, real_isotropic_shared_material=True,
            row_order='cell-major x/y/z, all sources share one real chi per cell',
            material_box=[self.lower,self.upper], width_lambda=self.width,
            h_lambda=self.h, cell_phase_at_upper=float(2*np.pi*self.h*np.sqrt(1+self.upper)),
            backend='R1 TensorDDA reused by R8 VectorPR3D; CM+RR dyadic point DDA',
            objective='mean_s ||S J_s-Y_s||²/||Y_s||² + R8 quadratic prior',
            prior='2e-4 mean nearest-neighbour |chi_i-chi_j|² + 1e-6 mean chi²',
            continuum_accuracy_certificate=False, unknown_support='entire inverse cube',
            legacy_setup=b.setup)

    @property
    def cost(self):
        return self.legacy.counts

    @property
    def counts(self):
        return self.cost

    def state(self, m):
        m = np.asarray(m, float)
        if m.shape != (self.n,) or not np.all(np.isfinite(m)):
            raise ValueError('material must contain N finite real isotropic contrasts')
        if self._current is not None and np.array_equal(m, self._current['m']):
            self.cost['cached_state_hits'] = self.cost.get('cached_state_hits',0)+1
            return self._current
        t,c = time.perf_counter(),time.process_time()
        d = np.repeat(m, 3)
        L = np.eye(self.M)-d[:,None]*self.G
        fac = lu_factor(L, check_finite=False)
        J = lu_solve(fac, d[:,None]*self.E, check_finite=False)
        Q = self.E+self.G@J
        Y = self.S@J
        residual = float(np.linalg.norm(J-d[:,None]*Q)/max(np.linalg.norm(J),np.finfo(float).tiny))
        if not np.all(np.isfinite(J)) or residual > 1e-9:
            raise RuntimeError(f'full-vector physical residual {residual}')
        self.cost['physical_LU'] += 1
        self.cost['source_RHS'] += self.ns
        self.cost['physical_G_columns'] = self.cost.get('physical_G_columns',0)+self.ns
        self.cost['physical_S_columns'] = self.cost.get('physical_S_columns',0)+self.ns
        self._current = dict(m=m.copy(),j=J,J=J,Q=Q,Y=Y,lu=fac,L=L,
            state_relative_residual=residual,cpu_seconds=time.process_time()-c,
            wall_seconds=time.perf_counter()-t)
        return self._current

    def forward(self,m):
        return self.state(m)['Y']

    def _data(self,m,Y):
        st=self.state(m);Y=np.asarray(Y,complex)
        if Y.shape != (self.nr,self.ns):
            raise ValueError('Y must have shape (receivers, sources)')
        norms=np.sum(abs(Y)**2,axis=0)
        if np.any(norms<=0):raise ValueError('source normalization requires nonzero observation columns')
        R=st['Y']-Y
        return st,R,norms,float(np.mean(np.sum(abs(R)**2,axis=0)/norms))

    def physical_score(self,m,Y):
        st,R,norms,data=self._data(m,Y)
        return data+self.reg.evaluate(m)[0]

    def data_value_gradient(self,m,Y):
        st,R,norms,data=self._data(m,Y)
        adj=lu_solve(st['lu'],self.S.conj().T@(R/norms[None,:]/self.ns),
                     trans=2,check_finite=False)
        self.cost['adjoint_RHS'] += self.ns
        self.cost['physical_S_adjoint_columns']=self.cost.get('physical_S_adjoint_columns',0)+self.ns
        gradient=2*np.real(np.sum(np.conj(adj)*st['Q'],axis=1).reshape(self.n,3).sum(1))
        prior,gprior,_=self.reg.evaluate(m)
        return data+prior,gradient+gprior

    def physical(self,m,Y=None,gradient=True):
        """R8-compatible physical return convention, backed by one retained LU."""
        st=self.state(m)
        if Y is None:return st['Y'],st['j'],st['Q']
        if gradient:return self.data_value_gradient(m,Y)
        _,_,_,data=self._data(m,Y)
        return data+self.reg.evaluate(m)[0],data,st['j'],st['Q']

    def data_objective(self,residual,Y):
        """Match the source-normalized R8 vector objective in external ledgers."""
        return float(np.mean(np.sum(abs(residual)**2,axis=0)/np.sum(abs(Y)**2,axis=0)))

    def data_terms(self,m,Y):
        return self.data_objective(self.state(m)['Y']-Y,Y),float(self.reg.evaluate(m)[0])

    def public_operators(self):
        """Only acquisition/operator metadata, with no target or observation."""
        return dict(G=self.G,E0=self.E,S=self.S,grid_shape=np.array(self.shape),
            components=np.array(3),width=np.array(self.width),centers=self.points,
            lower=np.array(self.lower),upper=np.array(self.upper),
            m_init=np.full(self.n,.25),material_objective=np.array('R8_source_normalized'))


class VectorPublicModel(DenseVectorModel):
    """Constructor used by the common public NPZ runner.

    The given matrices must equal the established R8 external geometry. This
    preserves mature vector CIE compatibility and prevents silent stale H/E2.
    """
    def __init__(self,public):
        shape=tuple(map(int,np.asarray(public['grid_shape'])))
        if len(shape)!=3 or len(set(shape))!=1:raise ValueError('cubic 3D grid required')
        super().__init__(shape[0],float(public.get('width',.24)),
            lower=float(public.get('lower',0.)),upper=float(public.get('upper',8.)))
        for label,expected,provided in [('G',self.G,public['G']),
                ('E0',self.E,public['E0']),('S',self.S,public['S'])]:
            if expected.shape!=provided.shape or not np.allclose(expected,provided,rtol=2e-13,atol=1e-14):
                raise ValueError(f'Public {label} is not the declared inherited vector acquisition')
        self.reg.smooth=self.reg.tv=self.reg.epsilon=None
        self.reg.name='R8_vector_quadratic_edge_mean_lam2e-4_plus_l2mean1e-6'


def run_vector_auxiliary(model,Y,m0,method,check,audit_callback,audit_every=50,
                         *,outer=120,inner=12):
    """Existing R8 vector PR-CG/material primitives with caller-owned budgets.

    No corrector lives here. CIE/Dual/Global preserve native 3D projected-IFFT
    mf1->2, material beta6->1, Dual current beta3 late, SVD range coordinates,
    shared real material and full quadratic prior. All vector components are
    in every material gradient. Audit scheduling is explicitly caller-owned.
    Callback signature: (m, paid_state, original_Phi, metadata).
    """
    from scipy.optimize import minimize
    names={'cie':'CIE_sourcewise','dual_cie':'Dual_sourcewise',
           'global_dual':'Dual_global','CIE_sourcewise':'CIE_sourcewise',
           'Dual_sourcewise':'Dual_sourcewise','CIE_global':'CIE_global'}
    if method not in names:
        raise ValueError('Mature vector adapter supports CIE/Dual/Global only; scalar G19/A+ are not vector substitutes')
    method=names[method];b=model.legacy;check();z=np.asarray(m0,float).copy()
    before=dict(model.cost);started=time.process_time();b.prstate=None
    prep=b.prepare(Y);check();steps=[]
    for iteration in range(outer):
        check();mf=1 if iteration<outer//2 else 2
        beta=6 if iteration<outer//2 else 1
        current_beta=3 if method.startswith('Dual_') and iteration>=outer//2 else beta
        J=b.current_pr(z,Y,'CIE_TSOM',current_beta,mf,method.endswith('global'))
        cols=b.bases[mf].shape[1]
        for key,value in dict(current_full_matrix_builds=1,current_KB_columns=cols,
                current_SB_columns=cols,current_KJd_columns=b.P,current_SJd_columns=b.P).items():
            b.counts[key]=b.counts.get(key,0)+value
        check()
        def fun(m):
            check();r=b.material(m,J,'CIE_TSOM',beta)
            b.counts['material_G_columns']=b.counts.get('material_G_columns',0)+b.P
            check();return r
        result=minimize(fun,z,jac=True,method='L-BFGS-B',bounds=model.bounds,
            options=dict(maxiter=inner,ftol=1e-12,gtol=1e-8,maxls=25))
        z=result.x
        steps.append(dict(iteration=iteration+1,mf=mf,beta_material=beta,beta_current=current_beta,
            nfev=int(result.nfev),nit=int(result.nit),message=str(result.message)))
        if (iteration+1)%audit_every==0 or iteration+1 in (outer//2,outer):
            check();st=model.state(z);value=model.physical_score(z,Y)
            audit_callback(z.copy(),st,value,dict(vector_method=method,iteration=iteration+1,
                mf=mf,beta_material=beta,beta_current=current_beta,paid_physical_audit=True))
            check()
    return dict(status='FIXED_VECTOR_AUXILIARY_CYCLES_COMPLETE',m=z,
        method=method,iterations=outer,optimizer=steps,subspace=prep,
        CPU_seconds=time.process_time()-started,
        counts={k:v-before.get(k,0) for k,v in model.cost.items()},
        qualification='R8 mature reduced-vector adaptation; native3D IFFT range mf1/2, not scalar FFT transplant')

# Stable single import point for the common R13 runner; implementations remain
# in their public-only scalar source-specific module.
try:
    from .physics_inputs import SourceSpecificPublicModel,run_fresnel_auxiliary
except ImportError: # Standalone runtime/source folder in the common runner.
    from physics_inputs import SourceSpecificPublicModel,run_fresnel_auxiliary
