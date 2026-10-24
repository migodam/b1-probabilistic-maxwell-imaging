"""Frozen SPD matrix pCN kernels; pure algebra, never a physics backend."""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq

KERNELS=('isotropic','weak','weak_bracket','weak_random')
EPSILON=.05
RMS_TARGET=.002
Q_CAP=.05


def real_decoder_metric(Dchi):
    D=np.asarray(Dchi)
    if D.shape!=(576,66) or not np.iscomplexobj(D) or not np.isfinite(D).all():raise ValueError('original Gaussianized Dchi[576,66] required')
    return np.concatenate((D.real,D.imag),axis=0)/24


@dataclass
class FrozenMatrixPCN:
    Q:np.ndarray

    def __post_init__(self):
        self.Q=np.asarray(self.Q,float)
        if self.Q.shape!=(66,66) or not np.isfinite(self.Q).all() or not np.allclose(self.Q,self.Q.T,rtol=0,atol=1e-13):raise ValueError('symmetric finite full66 Q')
        self.values,self.vectors=np.linalg.eigh(self.Q)
        if self.values.min()<=0 or self.values.max()>=1:raise ValueError('all66 eigenvalues strictly between0and1 required')
        self.sqrtQ=(self.vectors*np.sqrt(self.values))@self.vectors.T
        self.rho=(self.vectors*np.sqrt(1-self.values))@self.vectors.T
        self.logdet=float(np.log(self.values).sum())

    def propose(self,x,innovation):
        x,innovation=np.asarray(x,float),np.asarray(innovation,float)
        if x.shape!=(66,) or innovation.shape!=(66,) or not np.isfinite(x).all() or not np.isfinite(innovation).all():raise ValueError('full66 finite current and fixed innovation')
        return self.rho@x+self.sqrtQ@innovation

    def log_proposal(self,destination,origin):
        r=self.vectors.T@(np.asarray(destination)-self.rho@np.asarray(origin))
        return float(-.5*(66*np.log(2*np.pi)+self.logdet+np.sum(r*r/self.values)))

    def material_linear_expected_rms(self,D,x):
        drift=D@((self.rho-np.eye(66))@x)
        variance=float(np.sum((D@self.sqrtQ)**2))
        return float(np.sqrt(drift@drift+variance))


def normalized_kernel(C,D,x,rms_target=RMS_TARGET,q_cap=Q_CAP):
    C=np.asarray(C,float);D=np.asarray(D,float);x=np.asarray(x,float)
    values,vectors=np.linalg.eigh((C+C.T)/2)
    if C.shape!=(66,66) or D.shape!=(1152,66) or x.shape!=(66,) or values.min()<=0:raise ValueError('SPD C and source full66 material metric required')
    maximum=q_cap/values.max()
    DC=D@vectors;projected_x=vectors.T@x
    variance_slope=float(np.sum(DC*DC*values[None]))
    def squared(alpha):
        drift=DC@((np.sqrt(1-alpha*values)-1)*projected_x)
        return float(alpha*variance_slope+drift@drift)
    cap_rms=float(np.sqrt(squared(maximum)));binds=cap_rms<rms_target
    scale=maximum if binds else brentq(lambda a:squared(a)-rms_target*rms_target,0.,maximum,xtol=1e-16,rtol=1e-13)
    Q=(vectors*(scale*values))@vectors.T;Q=(Q+Q.T)/2;kernel=FrozenMatrixPCN(Q)
    actual=kernel.material_linear_expected_rms(D,x)
    if kernel.values.max()>q_cap*(1+1e-10):raise AssertionError('Q eigenvalue cap')
    return kernel,dict(scale=float(scale),lambda_min=float(kernel.values.min()),lambda_max=float(kernel.values.max()),
        material_linear_expected_RMS=actual,requested_RMS=rms_target,cap_binds=bool(binds),matched_RMS=not binds and abs(actual-rms_target)<1e-10,
        drift_squared=float(np.linalg.norm(D@((kernel.rho-np.eye(66))@x))**2),noise_trace=float(np.sum((D@kernel.sqrtQ)**2)),
        Q_cap=q_cap,full_rank=True,metric='source Gaussianized realified Dg/24; drift squared plus noise trace at fixed anchor')


def freeze_kernel_family(Vw,bracket,Dchi,x0,seed):
    V=np.asarray(Vw,float);b=np.asarray(bracket,float);D=real_decoder_metric(Dchi)
    if V.shape!=(66,4) or b.shape!=(66,) or not np.isfinite(V).all() or not np.isfinite(b).all() or not np.allclose(V.T@V,np.eye(4),atol=1e-10):raise ValueError('savedrank4weak frame and full66 bracket required')
    norm=float(np.linalg.norm(b))
    if norm<=1e-14:raise ValueError('MISSING_OR_ZERO_PAID_BRACKET; no replacement/recomputation')
    b=b/norm;random=np.random.default_rng(seed).normal(size=66);random/=np.linalg.norm(random)
    weak=EPSILON*np.eye(66)+V@V.T
    shapes={'isotropic':np.eye(66),'weak':weak,'weak_bracket':weak+np.outer(b,b),'weak_random':weak+np.outer(random,random)}
    result={}
    for name in KERNELS:
        kernel,metadata=normalized_kernel(shapes[name],D,np.asarray(x0,float))
        result[name]=dict(Q=kernel.Q,rho=kernel.rho,sqrtQ=kernel.sqrtQ,metadata=metadata)
    return result,dict(random_direction=random,bracket_direction=b,epsilon=EPSILON,enrichment_relative_weight=1.,direction_normalization='Euclidean unit; no state-dependent adaptation')


def finite_run_diagnostics(chains):
    """Descriptive initial-positive ACF ESS/MCSE; split Rhat across parents invalid.

    Each parent is a different target. No cross-parent Rhat is manufactured.
    One trajectory per kernel/parent cannot establish convergence or Rhat.
    """
    a=np.asarray(chains,float)
    if a.ndim==1:a=a[:,None]
    n=len(a);means=a.mean(0);variance=a.var(0,ddof=1) if n>1 else np.full(a.shape[1],np.nan);ess=[]
    for col in range(a.shape[1]):
        if not np.isfinite(a[:,col]).all() or variance[col]<=0:ess.append(np.nan);continue
        z=a[:,col]-means[col];cov=np.correlate(z,z,mode='full')[n-1:]/np.arange(n,0,-1);acf=cov/cov[0];total=0.
        for lag in range(1,min(n-1,512),2):
            pair=acf[lag]+acf[lag+1]
            if pair<=0:break
            total+=pair
        ess.append(float(min(n,n/max(1.,1+2*total))))
    ess=np.array(ess);mcse=np.sqrt(variance/ess)
    half=n//2
    if half>=2:
        split=np.stack((a[:half],a[-half:]));within=split.var(1,ddof=1).mean(0);between=half*split.mean(1).var(0,ddof=1)
        split_rhat=np.sqrt(np.divide((half-1)/half*within+between/half,within,out=np.full_like(within,np.nan),where=within>0))
    else:split_rhat=np.full(a.shape[1],np.nan)
    return dict(status='FINITE_RUN_UNQUALIFIED',n_steps=n,ESS_initial_positive_pair=ess.tolist(),MCSE_mean=mcse.tolist(),
                constant_function='ESS/MCSE NOT_EVALUATED representedNaN; never MCSE0 certification',Rhat_split_one_chain_descriptive=split_rhat.tolist(),Rhat_interpretation='two chronological halves of one target trajectory; not independent converged chains; never pool different parents',
                common_reference='NONE',qualification='NOT_ESTABLISHED')
