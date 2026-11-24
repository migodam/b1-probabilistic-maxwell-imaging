"""Frozen development-informed API for matched CIE branch pruning."""

from datetime import datetime, timezone

import numpy as np
from threadpoolctl import threadpool_limits

from branches import expand
from fast_support import TileProjection
from native_box import OPTICS, NativeBox, clock, elapsed
from restricted import score_hard
from uniform import UniformLow

SUPPORT_T=(8.,)
SUPPORT_ALPHA=(2.**-14,)


class Checks:
    def __init__(self,spectrum,S,B,*,radius_scope='native'):
        start=clock()
        self.geometry=NativeBox(spectrum,S,B)
        self.radius_scope=radius_scope
        self.aperture=None
        with np.load(OPTICS/'synthetic/public/fine_geometry.npz',allow_pickle=False) as f:
            if np.array_equal(spectrum,f['spectrum']):
                for arch in ('PW180_10x20','PW360_10x20'):
                    if np.array_equal(S,f[arch+'__S']) and np.array_equal(B,f[arch+'__B']):
                        self.aperture=arch
        self.uniform=UniformLow(self.geometry,self.aperture) if self.aperture else None
        self.tile_projection=TileProjection(self.geometry) if self.aperture else None
        self.setup_cost=elapsed(start)

    def score(self,Y,r,branch):
        """Use PRUNE only under the archived geometry and supplied valid native ball."""
        start=clock()
        y=np.asarray(Y,complex)
        g=self.geometry
        if y.shape!=(g.nr,g.ns) or not np.all(np.isfinite(y)) or not np.isfinite(r) or r<0:
            raise ValueError('Expected complete common-aperture observations and a valid radius value.')
        material_upper=expand(branch,g.side)
        proof_allowed=self.radius_scope=='native' and self.aperture is not None
        if branch.name=='wildcard':
            result=dict(branch=branch.name,decision='RETAIN',residual_lower=0.,radius=r,
                        rigor='declared_wildcard_no_screen')
        elif self.uniform and np.all(material_upper==.5):
            result=self.uniform.score(y,r)
        elif np.any(material_upper==0) and 0<int(np.count_nonzero(material_upper))<=576:
            fast=None
            tile=np.asarray(branch.caps,float)
            if (self.tile_projection is not None and tile.shape==(4,4) and
                    int(np.count_nonzero(tile))==1 and float(tile.max())==4.):
                try:
                    fast=self.tile_projection.score(y,r,branch)
                except (ArithmeticError,ValueError,np.linalg.LinAlgError) as error:
                    fast=dict(decision='RETAIN',error=str(error),**elapsed(start))
            if fast and fast['decision']=='PRUNE':
                result=fast
            else:
                result=score_hard(g,y,r,branch,t_grid=SUPPORT_T,alpha_grid=SUPPORT_ALPHA,max_active=576,max_proofs=1)
                result['translated_tile_check']=fast
        else:
            result=g.score(y,r,branch,max_proofs=1,max_rank=850)
        if not proof_allowed and result['decision']=='PRUNE':
            result['possible_conditional_exclusion']=True
            result['decision']='RETAIN'
        return {**result,'timestamp_utc':datetime.now(timezone.utc).isoformat(),
                'aperture':self.aperture,'model_scope':self.radius_scope,
                'native_physical_pruning_eligible':proof_allowed,'check_time':elapsed(start)}


def limited_checks(*args,**kwargs):
    """Set one BLAS thread around frozen replay if the caller lacks a thread budget."""
    with threadpool_limits(1):
        return Checks(*args,**kwargs)
