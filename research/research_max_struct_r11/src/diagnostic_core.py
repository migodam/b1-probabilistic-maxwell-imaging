"""Offline fixed-current mechanism measurements, with explicit split and costs."""
from collections import Counter
import numpy as np
from scipy.linalg import svdvals

def cosine(x,y):
    norm=float(np.linalg.norm(x)*np.linalg.norm(y))
    return None if norm==0 else float(np.dot(x,y)/norm)

def projected_gradient(m,g,lo,hi):
    """Box projected-gradient mapping at unit step; log scale separately."""
    return m-np.clip(m-g,lo,hi)

def defect_stats(e,w,Os):
    return dict(raw_defect_fro=float(np.linalg.norm(e)),
        source_weighted_defect_squared=float(np.vdot(e,e@Os).real),
        transformed_defect_fro=float(np.linalg.norm(w)),
        transformed_source_weighted_squared=float(np.vdot(w,w@Os).real))

def measure_start(operator,J,Q,m,beta,Os,rho,reg,physical_gradient,bounds):
    factor=operator.factor(m,beta)
    e=J-m[:,None]*Q;w=operator.solve(factor,e)
    a=operator.solve(factor,w@Os,adjoint=True)
    kw=operator.action(w)
    gconst=-rho*beta**2*np.real(np.sum(a.conj()*Q,axis=1))
    gmetric=-rho*beta**2*np.real(np.sum(a.conj()*(beta*w-kw),axis=1))
    prior,gprior=reg(m);gfull=gconst+gmetric+gprior;glag=gconst+gprior
    K=operator.B@operator.U.conj().T
    A=np.diag(1+beta*m)-m[:,None]*K
    singular=svdvals(A,check_finite=False)
    lo,hi=bounds
    R=beta*e/(1+beta*m[:,None])
    dR=-beta*(Q+beta*J)/(1+beta*m[:,None])**2
    gordinary=rho*np.real(np.sum(dR.conj()*(R@Os),axis=1))+gprior
    gradients=dict(g_const=gconst,g_metric=gmetric,g_prior=gprior,g_full=gfull,g_lagged=glag,g_ordinary=gordinary)
    meta=dict(prior_value=float(prior),beta=beta,cos_const_metric=cosine(gconst,gmetric),
        metric_const_norm_ratio=float(np.linalg.norm(gmetric)/max(np.linalg.norm(gconst),1e-300)),
        denominator_smax=float(singular[0]),denominator_smin=float(singular[-1]),denominator_cond2=float(singular[0]/singular[-1]),
        box_lower_active=int(np.count_nonzero(m<=lo+1e-12)),box_upper_active=int(np.count_nonzero(m>=hi-1e-12)),
        defect=defect_stats(e,w,Os),gradients={})
    for key,g in gradients.items():
        pg=projected_gradient(m,g,lo,hi)
        meta['gradients'][key]=dict(norm=float(np.linalg.norm(g)),inf=float(np.linalg.norm(g,np.inf)),
            descent_cos_negative_physical_gradient=cosine(-g,-physical_gradient),
            box_projected_gradient_norm=float(np.linalg.norm(pg)),box_projected_gradient_inf=float(np.linalg.norm(pg,np.inf)),
            projected_descent_cos_negative_physical_gradient=cosine(-pg,-physical_gradient))
    return meta,dict(**gradients,e=e,w=w,a=a,physical_gradient=physical_gradient)
