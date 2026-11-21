"""Development arm: defer exact Schur setup to first active ball and then cache it."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import engine
from engine import np, la

# Use a private copy of the same generator, retaining the published arm and
# avoiding process-wide changes to the frozen engine's projector factory.
spec = importlib.util.spec_from_file_location('srap_development_engine',HERE.parent/'engine.py')
development = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = development
spec.loader.exec_module(development)


class LazyCachedBall(engine.BallProjector):
    def __init__(self,model,Y,delta,recycled=None):
        start=model.clock.wall()
        self.model,self.E,self.sources=model,model.E,model.ns
        self.cells,self.components=model.n,1
        self.lower,self.upper=np.zeros(model.n),np.full(model.n,20.)
        self.Y=Y if model.ports.valid_mask is None else np.where(model.ports.valid_mask,Y,0.)
        self.delta,self.costs=float(delta),model.clock.counts
        self.cg_tolerance=2e-10
        self.warm=self._coarse=self.cached=None
        self.recycled=recycled
        self.last_A_info={}
        r15=engine.legacy.R15_ROOT/'runtime/field_current.py'
        mod=engine.legacy.module_at('srap_field_current',r15 if r15.is_file() else engine.legacy.R12_ROOT/'src/field_current.py')
        self.material=mod.DenseFieldCurrentProjectors.material.__get__(self)
        self._fields=mod.DenseFieldCurrentProjectors._fields.__get__(self)
        if recycled is not None and recycled.shape[1]:
            gb=self.B(recycled)
            vals,rot=la.eigh(gb.conj().T@gb,check_finite=False)
            basis=recycled@rot
            self._coarse=basis,(vals/(1+vals))[:,None]*basis.conj().T
        self.lazy_setup=model.clock.wall()-start
        self.costs['DR_setup_wall']=self.lazy_setup
        self.costs['lazy_coarse_setup_s']=self.lazy_setup

    def affine_ball(self,x,q):
        if self.cached is not None:
            self.costs['lazy_cached_calls']+=1
            xf,qf=self.cached.affine_ball(x,q)
            self.last_A_info=dict(self.cached.last_A_info)
            return xf,qf
        self.costs['H_iteration_rhs']+=self.sources
        xu=self.solve(x+self.BH(q-self.E),start=self.warm)
        self.warm=xu.copy()
        r=self.scale*self.model.ports.receive(xu)-self.Y
        raw=float(np.linalg.norm(r))
        if raw<=self.delta:
            self.costs['A_projections']+=1
            self.costs['lazy_inactive_no_schur']+=1
            self.last_A_info=dict(multiplier=0.,affine_data_residual=raw,
                                  radius=self.delta,unconstrained_residual=raw,relative_excess=0.)
            return xu,self.E+self.B(xu)
        self.costs['lazy_active_setup_calls']+=1
        self.cached=engine.BallProjector(self.model,self.Y,self.delta,self.recycled)
        self.costs['DR_setup_wall']+=self.lazy_setup
        # The first unconstrained solution already passed the full H residual check.
        self.cached.warm=xu.copy()
        xf,qf=self.cached.affine_ball(x,q)
        self.last_A_info=dict(self.cached.last_A_info)
        return xf,qf


development.BallProjector=LazyCachedBall


def run_cie_arm(model,Y,case,method,m,delta):
    if method not in ('srap_lazy','lazy_srap','srap'):
        raise ValueError(f'Unknown lazy development arm {method}')
    yield from development.run_cie_arm(model,Y,case,'srap',m,delta)
