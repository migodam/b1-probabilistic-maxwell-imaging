"""Free independent real/loss material FEM model; no shape labels or truth I/O.

Pinned R2 P3/PML mesh and material quadrature match the original digit
experiment. Every LU and solved RHS passes through the supplied shared book.
A weighted adjoint uses four receiver-weighted RHS per frequency, never full J.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from contextlib import contextmanager
import fcntl
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
R2 = ROOT.parent / 'research_w_gq_r2'
_vendor=ROOT/'vendor'
_backend_paths=[_vendor] if _vendor.exists() else [R2/'src',R2/'reuse/previous_gq/src',R2/'reuse/fixed_vendor']
for _path in reversed(_backend_paths+[ROOT.parents[2]/'research-w-continuation/.deps']):
    if str(_path) not in sys.path:
        sys.path.insert(0,str(_path))
from skfem import asm
from scipy.sparse import diags
from r2.physics import Eval, WallGrid, S1_TX, RX16, line_green, _axis, _real_array, _readonly, WALL_BOX
from wslam.fem_backend import _contrast_mass, _forcing


def default_config():
    return dict(model='R4_FREE_MATERIAL_P3_WEIGHTED_ADJOINT_V2_CACHED',
                backend_source='R4_VENDOR' if _vendor.exists() else 'READONLY_R2_FALLBACK',
                assembly_strategy='cached_combined',
                wall_x=[-1.5,-.5,.5,1.5], wall_y=[.5,.875,1.25],
                target_x_edges=np.linspace(-.7,.7,8).tolist(),
                target_y_edges=np.linspace(1.4,2.8,8).tolist(),
                frequencies=[.85,1.,1.15], sigma_real_imag=.01,
                h=.17, targetstep=.10, wall_refine=1,
                max_factor_bytes=1500000000)


@contextmanager
def physical_lock(book):
    """One campaign-wide lock, with the root supplied by the shared book."""
    root=Path(book.root);root.mkdir(parents=True, exist_ok=True)
    with (root/'.fem.lock').open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            book.check()
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


@dataclass(frozen=True, eq=False)
class PublicWall(WallGrid):
    def __post_init__(self):
        x=_axis(self.x,WALL_BOX[:2],'wall x');y=_axis(self.y,WALL_BOX[2:],'wall y')
        shape=(len(y)-1,len(x)-1)
        eps=_real_array(self.epsinf,shape,'epsinf');ell=_real_array(self.ell,shape,'ell')
        if np.any((eps<1)|(eps>6)) or np.any((ell<0)|(ell>.25)):
            raise ValueError('wall requires independent eps [1,6], loss [0,.25]')
        for name,value in [('x',x),('y',y),('epsinf',eps),('ell',ell)]:
            object.__setattr__(self,name,_readonly(value))
    @property
    def support(self):
        return np.ones(self.shape,bool)


@dataclass
class Prediction:
    data: np.ndarray
    diagnostics: list


class FEMModel:
    def __init__(self, config, book, phase='R4_FWI'):
        self.config={**default_config(), **config};self.book=book;self.phase=phase
        c=self.config;self.wx=np.asarray(c['wall_x'],float);self.wy=np.asarray(c['wall_y'],float)
        self.tx_edges=np.asarray(c['target_x_edges'],float);self.ty_edges=np.asarray(c['target_y_edges'],float)
        for edge in (self.tx_edges,self.ty_edges):
            if edge.ndim!=1 or len(edge)<2 or not np.isfinite(edge).all() or np.any(np.diff(edge)<=0):
                raise ValueError('target edges must increase')
        self.wall_shape=(len(self.wy)-1,len(self.wx)-1);self.target_shape=(len(self.ty_edges)-1,len(self.tx_edges)-1)
        self.nwall=int(np.prod(self.wall_shape));self.ntarget=int(np.prod(self.target_shape));self.nparam=2*(self.nwall+self.ntarget)
        self.sigma=float(config.get('sigma_real_imag',config.get('sigma',.01)))
        if self.sigma!=.01:raise ValueError('R4 preserves sigma=.01 independently per real/imag component')
        self.evaluators=[];self.cell_data=[];self.incident=[];self.quadrature_cache=[]
        book.check(new_job=True)
        with physical_lock(book):
            for freq in c['frequencies']:
                ev=Eval(book,phase,float(freq),S1_TX,self.wx,self.wy,
                        wall_refine=c['wall_refine'],h=c['h'],targetstep=c['targetstep'],
                        max_factor_bytes=c['max_factor_bytes'])
                p=ev.p;points=p.material_points
                inside=((points[:,0]>=self.tx_edges[0])&(points[:,0]<=self.tx_edges[-1])&
                        (points[:,1]>=self.ty_edges[0])&(points[:,1]<=self.ty_edges[-1]))
                ix=np.clip(np.searchsorted(self.tx_edges,points[:,0],side='right')-1,0,self.target_shape[1]-1)
                iy=np.clip(np.searchsorted(self.ty_edges,points[:,1],side='right')-1,0,self.target_shape[0]-1)
                self.evaluators.append(ev);self.cell_data.append((inside,iy*self.target_shape[1]+ix))
                self.incident.append(line_green(points,p.tx,p.k))
                from .currents import quadrature_interpolation
                cache_cpu=time.process_time();cache_wall=time.perf_counter()
                active=np.asarray(ev._inside_wall|inside);P=quadrature_interpolation(ev,active)
                self.quadrature_cache.append(dict(P=P,weights=p.material_weights[active].copy(),
                    incident=self.incident[-1][active].copy(),wall=ev._inside_wall[active],
                    target=inside[active],wall_cells=ev._cells[active],target_cells=(iy*self.target_shape[1]+ix)[active]))
                book.event(phase,'r4_exact_quadrature_cache',frequency=float(freq),
                    quadrature_points=int(active.sum()),sparse_interpolation_nnz=int(P.nnz),
                    cpu_seconds=time.process_time()-cache_cpu,wall_seconds=time.perf_counter()-cache_wall,
                    factorizations=0,rhs=0)

    def bounds_and_scale(self):
        nw,nt=self.nwall,self.ntarget
        return (np.r_[np.ones(nw),np.zeros(nw+2*nt)],
                np.r_[np.full(nw,6.),np.full(nw,.25),np.full(nt,8.),np.full(nt,.25)],
                np.r_[np.ones(nw),np.full(nw,.10),np.full(nt,2.),np.full(nt,.10)])

    def initial(self,start=0):
        """Two public uniform starts, independent of observations and world ID."""
        if start not in (0,1):raise ValueError('only frozen starts 0,1')
        return np.r_[np.full(self.nwall,3. if start==0 else 3.8),
                     np.full(self.nwall,.06 if start==0 else .08),
                     np.full(self.ntarget,.10 if start==0 else .40),
                     np.full(self.ntarget,.03 if start==0 else .06)]

    def validate(self,z):
        z=np.asarray(z,float);lo,hi,_=self.bounds_and_scale()
        if z.shape!=(self.nparam,) or not np.isfinite(z).all() or np.any(z<lo-1e-12) or np.any(z>hi+1e-12):
            raise ValueError('illegal free-material coordinates')
        return np.clip(z,lo,hi)

    def _run(self,z,y=None,gradient=False,frequency_indices=None):
        z=self.validate(z);nw,nt=self.nwall,self.ntarget
        wall=PublicWall(self.wx,self.wy,z[:nw].reshape(self.wall_shape),z[nw:2*nw].reshape(self.wall_shape))
        indices=list(range(len(self.evaluators))) if frequency_indices is None else list(frequency_indices)
        if not indices or any(i<0 or i>=len(self.evaluators) for i in indices):raise ValueError('frequency indices')
        if gradient and y is None:raise ValueError('weighted gradient requires observations')
        if y is not None:
            y=np.asarray(y)
            if y.shape==(len(self.evaluators),4,16):y=y[indices]
            if y.shape!=(len(indices),4,16) or not np.isfinite(y).all():raise ValueError('observations must be finite complex [frequency,source,receiver]')
        data=[];diagnostics=[];g=np.zeros(self.nparam)
        with physical_lock(self.book):
            for jj,ii in enumerate(indices):
                self.book.check();ev=self.evaluators[ii];p=ev.p;inside,cells=self.cell_data[ii];inc=self.incident[ii]
                cpu=time.process_time();walltime=time.perf_counter()
                ct=np.zeros(len(cells),complex)
                ct[inside]=z[2*nw:2*nw+nt][cells[inside]]+1j*z[2*nw+nt:][cells[inside]]/p.freq
                strategy=self.config.get('assembly_strategy','cached_combined')
                if strategy=='legacy_two_pass':
                    A,b=ev.background(wall)
                    mass=asm(_contrast_mass,p.material_basis,contrast=ct.reshape(p._material_shape))
                    A=(A-p.k**2*mass[p._free_dofs][:,p._free_dofs]).tocsc()
                    for src in range(4):
                        b[:,src]+=asm(_forcing,p.material_basis,load=(p.k**2*ct*inc[:,src]).reshape(p._material_shape))[p._free_dofs]
                elif strategy=='cached_combined':
                    # Same weak operator/forcing: combine both disjoint public
                    # material regions before one assembly; reuse analytic inc.
                    combined=wall.values(p.material_points,p.freq)-1.+ct
                    mass=asm(_contrast_mass,p.material_basis,contrast=combined.reshape(p._material_shape))
                    A=(p._background_matrix-p.k**2*mass[p._free_dofs][:,p._free_dofs]).tocsc()
                    b=np.column_stack([asm(_forcing,p.material_basis,
                        load=(p.k**2*combined*inc[:,src]).reshape(p._material_shape))[p._free_dofs] for src in range(4)])
                elif strategy=='sparse_quadrature':
                    cache=self.quadrature_cache[ii];P=cache['P'];weights=cache['weights']
                    wm=cache['wall'];tm=cache['target'];wc=cache['wall_cells'];tc=cache['target_cells']
                    chi=np.zeros(len(weights),complex)
                    chi[wm]=z[:nw][wc[wm]]-1.+1j*z[nw:2*nw][wc[wm]]/p.freq
                    chi[tm]+=z[2*nw:2*nw+nt][tc[tm]]+1j*z[2*nw+nt:][tc[tm]]/p.freq
                    A=(p._background_matrix-p.k**2*(P.T@diags(weights*chi)@P)).tocsc()
                    b=p.k**2*np.asarray(P.T@(weights[:,None]*chi[:,None]*cache['incident']))
                else:raise ValueError('unknown assembly strategy')
                factor,A,storage=ev.factor(A,label=self.phase+' free material LU')
                u=factor.solve(b,label=self.phase+' primal4')
                primal_defect=factor.last_defect
                d=ev.direct+np.asarray(ev.R@u).T;data.append(d.copy())
                adjoint_defect=None
                if gradient:
                    adj=factor.solve(np.asarray(ev.R.conj().T@((d-y[jj])/self.sigma**2).T),trans='H',label=self.phase+' weighted adjoint4')
                    adjoint_defect=factor.last_defect
                    gw=np.zeros(nw,complex);gt=np.zeros(nt,complex)
                    if strategy=='sparse_quadrature':
                        total=cache['incident']+np.asarray(P@u);av=np.asarray(P@adj)
                        prod=p.k**2*weights*np.einsum('qs,qs->q',np.conj(av),total)
                        np.add.at(gw,wc[wm],prod[wm]);np.add.at(gt,tc[tm],prod[tm])
                    else:
                        for sl,(scattered,av) in ev._material_chunks(u,adj):
                            total=scattered+inc[sl]
                            prod=p.k**2*p.material_weights[sl]*np.einsum('qs,qs->q',np.conj(av),total)
                            wm=ev._inside_wall[sl];tm=inside[sl]
                            np.add.at(gw,ev._cells[sl][wm],prod[wm]);np.add.at(gt,cells[sl][tm],prod[tm])
                    g[:nw]+=gw.real;g[nw:2*nw]+=(1j/p.freq*gw).real
                    g[2*nw:2*nw+nt]+=gt.real;g[2*nw+nt:]+=(1j/p.freq*gt).real
                diagnostics.append(dict(frequency=p.freq,ndofs=p.ndofs,factor_bytes=storage,assembly_strategy=strategy,
                    primal_backward=primal_defect,adjoint_backward=adjoint_defect,
                    cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-walltime))
        prediction=Prediction(np.asarray(data),diagnostics)
        if y is None:return prediction
        value=.5*float(np.sum(np.abs(prediction.data-y)**2))/self.sigma**2
        return value,g,prediction

    def predict(self,z,frequency_indices=None):
        return self._run(z,frequency_indices=frequency_indices)

    def gradient(self,z,y,frequency_indices=None):
        return self._run(z,y,True,frequency_indices)

    def run(self,z,y=None,gradient=False,frequency_indices=None):
        result=self._run(z,y,gradient,frequency_indices)
        if y is None:return result.data,result.diagnostics
        value,g,pred=result;return value,g,pred.data,pred.diagnostics

    def maps(self,z):
        z=self.validate(z);nw,nt=self.nwall,self.ntarget
        return dict(wall_eps=z[:nw].reshape(self.wall_shape),wall_loss=z[nw:2*nw].reshape(self.wall_shape),
                    target_contrast=z[2*nw:2*nw+nt].reshape(self.target_shape),
                    target_loss=z[2*nw+nt:].reshape(self.target_shape),
                    target_x_edges=self.tx_edges,target_y_edges=self.ty_edges)
