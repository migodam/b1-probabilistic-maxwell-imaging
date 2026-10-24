"""Pure geometry algebra for the R4 paid driver; no physics backend import.

Frames are generalized material-RMS directions of J'J+I. Field definitions
use a Euclidean cluster projector and fixed seed vectors, then normalize in
the actual decoder metric. This removes basis-sign/within-cluster gauge.
No finite Euler circuit is called a reversible MH proposal or exact flow.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from .core import to_probit, from_probit, probit_derivative


class GeometryInvalid(ValueError):pass


@dataclass
class Chart:
    q0: np.ndarray
    dimension: int

    def __post_init__(self):
        self.q0=np.array(self.q0,float,copy=True)
        if self.dimension==64:self.x0=self.q0[:64].copy()
        elif self.dimension==66:
            if np.any((self.q0[64:]<=0)|(self.q0[64:]>=1)):
                raise GeometryInvalid('BOUNDARY_NATIVE_ANCHOR_NO_FINITE_EXACT_PROBIT_CHART')
            self.x0=to_probit(self.q0)
        else:raise ValueError('64 field or 66 probit coordinates only')

    def native(self,x):
        x=np.asarray(x,float)
        if x.shape!=(self.dimension,) or not np.isfinite(x).all():raise GeometryInvalid('nonfinite coordinate')
        if self.dimension==66:return from_probit(x)
        q=self.q0.copy();q[:64]=x;return q

    def derivative(self,prior,x):
        D=prior.derivative(self.native(x))
        return D[:,:64] if self.dimension==64 else D*probit_derivative(x)[None,:]

    def jacobian(self,Jq,x):
        return Jq[:,:64] if self.dimension==64 else Jq*probit_derivative(x)[None,:]

    def prior_energy(self,x):return .5*float(np.asarray(x)@x)


def material_real(D):
    D=np.asarray(D)
    return np.r_[D.real,D.imag]/np.sqrt(len(D))


def material_frame(J,D,rank,*,relative_cutoff=1e-8,absolute_cutoff=1e-14):
    J=np.asarray(J,float);D=np.asarray(D);M=material_real(D)
    if J.ndim!=2 or D.ndim!=2 or J.shape[1]!=D.shape[1] or not np.isfinite(J).all() or not np.isfinite(D).all():
        raise GeometryInvalid('invalid frame arrays')
    _,s,Vh=np.linalg.svd(M,full_matrices=False)
    keep=s>max(float(s[0])*relative_cutoff if len(s) else 0.,absolute_cutoff)
    material_rank=int(keep.sum())
    if material_rank<=rank:raise GeometryInvalid(f'insufficient material rank or missing cluster complement: {material_rank} <= {rank}')
    T=Vh[keep].T/s[keep]
    H=(J@T).T@(J@T)+T.T@T
    lam,U=eigh((H+H.T)/2)
    selected=T@U[:,:rank]
    Q=np.linalg.qr(selected,mode='complete')[0]
    P=Q[:,:rank]@Q[:,:rank].T
    gap=float(lam[rank]-lam[rank-1])
    relative_gap=gap/max(abs(float(lam[rank])),abs(float(lam[rank-1])),np.finfo(float).tiny)
    return dict(projector=P,Vw=Q[:,:rank],Vs=Q[:,rank:],directions=selected,eigenvalues=lam,
        material_singular_values=s,material_rank=material_rank,rank=int(rank),
        eigen_gap=gap,relative_gap=float(relative_gap),near_crossing=bool(relative_gap<=1e-3),
        metric_condition=float(s[keep].max()/s[keep].min()),
        generalized_condition=float(max(abs(lam[-1]),1e-30)/max(abs(lam[0]),1e-30)))


def projector_field(frame,D,seed):
    raw=frame['projector']@np.asarray(seed,float)
    scale=float(np.linalg.norm(material_real(D)@raw))
    if not np.isfinite(scale) or scale<=1e-10:raise GeometryInvalid('NO_OP_PROJECTED_MATERIAL_FIELD')
    return raw/scale


def central_bracket(XplusY,XminusY,YplusX,YminusX,h):
    """[X,Y]=DY X-DX Y, including derivative of material normalization."""
    return (np.asarray(YplusX)-YminusX-np.asarray(XplusY)+XminusY)/(2*h)


def rotation(frame,base):
    singular=np.linalg.svd(base['Vw'].T@frame['Vw'],compute_uv=False)
    return dict(projector_operator=float(np.linalg.norm(frame['projector']-base['projector'],2)),
        largest_principal_angle=float(np.max(np.arccos(np.clip(singular,0,1)))),
        relative_gap=frame['relative_gap'],material_rank=frame['material_rank'],near_crossing=frame['near_crossing'])


def material_normal(vector,D,frame):
    M=material_real(D);mapped=M@vector;weak=M@frame['Vw']
    U,s,_=np.linalg.svd(weak,full_matrices=False)
    keep=s>max(float(s[0])*1e-10 if len(s) else 0.,1e-14)
    Q=U[:,keep];normal=mapped-Q@(Q.T@mapped)
    amplitude=float(np.linalg.norm(mapped))
    return dict(latent_norm=float(np.linalg.norm(vector)),material_amplitude=amplitude,
        material_normal_amplitude=float(np.linalg.norm(normal)),
        material_normal_fraction=float(np.linalg.norm(normal)/amplitude) if amplitude else 0.,
        latent_normal_norm=float(np.linalg.norm(vector-frame['projector']@vector)),
        weak_material_rank=int(keep.sum()))


def distance_match(prior,chart,origin_chi,start,direction,radius,*,latent_cap=8.):
    """First dyadic decoder-only distance root, never a paid-response search."""
    start=np.asarray(start,float);direction=np.asarray(direction,float)
    if direction.shape!=start.shape or np.linalg.norm(direction)<=1e-14:raise GeometryInvalid('NO_OP_DIRECT_DIRECTION')
    origin=chart.x0
    def fun(t):
        x=start+t*direction
        if np.linalg.norm(x-origin)>latent_cap*(1+1e-12):raise GeometryInvalid('DIRECT_LATENT_CAP')
        return prior.material_rms(prior.decode(chart.native(x))[0],origin_chi)-radius
    if fun(0)>1e-10:raise GeometryInvalid('distance root begins outside requested radius')
    # Conservative cap respects the complete displacement, also for 2nd segment.
    cap=(latent_cap-np.linalg.norm(start-origin))/np.linalg.norm(direction)
    if cap<=0:raise GeometryInvalid('DIRECT_LATENT_CAP')
    hi=min(cap,max(radius/8,1e-6))
    for _ in range(24):
        if fun(hi)>=0:
            root=brentq(fun,0.,hi,xtol=1e-13)
            x=start+root*direction
            actual=prior.material_rms(prior.decode(chart.native(x))[0],origin_chi)
            return x,dict(scale=float(root),actual_RMS=actual,radius_matched=bool(abs(actual-radius)<=1e-8))
        if hi>=cap:break
        hi=min(cap,2*hi)
    raise GeometryInvalid('PREDECLARED_DIRECTION_CANNOT_REACH_RADIUS_WITHIN_CAP_NO_REDRAW')


def phase_diagnostic(prior,q,base_q):
    _,a=prior.decode(q);_,b=prior.decode(base_q)
    phase=lambda f:np.where(f<=-prior.temperature,-1,np.where(f>=prior.temperature,1,0))
    s=a['field']/prior.temperature
    return dict(seam_crossed_cells=int(np.count_nonzero((phase(a['field'])!=phase(b['field']))&prior.mask)),
        transition_cells=int(np.count_nonzero((abs(s)<1)&prior.mask)),
        minimum_abs_standardized_seam_distance=float(np.min(abs(abs(s[prior.mask])-1))),
        continuous_path_verified=False)


def finite_material_normal(delta_chi,D,frame):
    """Actual finite material displacement normal to the BASE material tangent."""
    delta=np.asarray(delta_chi);mapped=np.r_[delta.real,delta.imag]/np.sqrt(len(delta))
    weak=material_real(D)@frame['Vw'];U,s,_=np.linalg.svd(weak,full_matrices=False)
    keep=s>max(float(s[0])*1e-10 if len(s) else 0.,1e-14)
    Q=U[:,keep];normal=mapped-Q@(Q.T@mapped);total=float(np.linalg.norm(mapped))
    return dict(actual_material_RMS=total,actual_normal_RMS=float(np.linalg.norm(normal)),
        actual_normal_fraction=float(np.linalg.norm(normal)/total) if total else 0.,
        reference='base material tangent; finite displacement, not decoder-linearized')


def heun_corrector(start,start_field,predictor_field,signed_step):
    """Second-order Taylor-consistent segment; C1 gives o(h²), not global O(h³)."""
    return np.asarray(start)+.5*signed_step*(np.asarray(start_field)+np.asarray(predictor_field))
