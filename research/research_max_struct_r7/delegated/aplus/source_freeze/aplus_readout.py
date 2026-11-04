"""R7 low-cost readouts of the unchanged fixed-current onehop material loss.

Möbius four-color conditional MM uses the actual normalized supplied edge prior.
The adjoint projected control uses the inherited exact derivative. Neither makes
an original Maxwell descent or basin claim. Historical modules stay untouched.
"""
import time
import numpy as np
from research.research_max_struct_r2.src.cie_global import _baseline_readout
from research.research_max_struct_r4.src.cie_onehop import onehop_material_value_gradient


def _seed(blocks,J,Q,m,beta,Os,Ls,rho,bounds,reg,config,cost):
    seed,_=_baseline_readout(J@Ls,Q@Ls,m,beta,np.ones(J.shape[1]),rho,
                             bounds,reg,False,config,cost)
    def value(z):
        R=blocks.apply(z,beta,J-z[:,None]*Q)
        cost['material_inner_evaluations']+=1
        cost['material_prior_evaluations']+=1
        cost['source_metric_products']+=1
        return .5*rho*float(np.vdot(R,R@Os).real)+float(reg(z)[0])
    old,seed_value=value(m),value(seed)
    return (seed.copy() if seed_value<old else m.copy()),old,seed_value,value


def _prior_context(config,blocks):
    prior=config.get('aplus_regularizer')
    if prior is None:
        raise ValueError('Möbius MM requires explicit original aplus_regularizer object')
    cached=config.get('_aplus_prior_context')
    if cached is not None and cached['prior'] is prior:return cached
    D=prior.D.tocsr()
    if np.any(np.diff(D.indptr)!=2):raise ValueError('Prior must have two-cell edges')
    endpoints=D.indices.reshape(-1,2);coeff=D.data.reshape(-1,2)
    if not np.allclose(coeff[:,0],-coeff[:,1],rtol=0,atol=1e-15):
        raise ValueError('MM prior must be weighted material differences')
    n=blocks.n;neighbors=np.full((n,4),n,int);squared=np.zeros((n,4));degree=np.zeros(n,int)
    for pair,values in zip(endpoints,coeff):
        for i,j,v in ((pair[0],pair[1],values[0]),(pair[1],pair[0],values[1])):
            if degree[i]>=4:raise ValueError('Only nearest-neighbor grid prior is supported')
            neighbors[i,degree[i]]=j;squared[i,degree[i]]=v*v;degree[i]+=1
    colors=np.empty(n,int);colors[blocks.indices]=np.arange(4)[None,:]
    if np.any(colors[endpoints[:,0]]==colors[endpoints[:,1]]):
        raise ValueError('Prior edges break the independent four-color assumption')
    cached=dict(prior=prior,neighbors=neighbors,squared=squared,edges=int(prior.edges))
    config['_aplus_prior_context']=cached
    return cached


def _conditional_candidates(q,c,b,h,s,lo,hi,old,cost):
    """Real stationary quintic candidates plus endpoints/current, batch roots.

    f(m)= (n1*m+n2*m²)/|1+b*m|² + h*m²/2-s*m+const.
    Compare every numerically real box root, never choose by truth.
    """
    d1=2*b.real;d2=abs(b)**2
    n1=c.real;n2=.5*q+(c*b).real
    d4=np.stack((np.ones_like(d1),2*d1,d1*d1+2*d2,2*d1*d2,d2*d2),axis=1)
    p=np.zeros((len(q),6));p[:,:5]-=s[:,None]*d4;p[:,1:]+=h[:,None]*d4
    p[:,0]+=n1;p[:,1]+=2*n2;p[:,2]+=n2*d1-n1*d2
    scale=np.maximum(abs(lo),abs(hi));scale=np.maximum(scale,1.)
    p*=scale[:,None]**np.arange(6)[None,:]
    # Nonzero prior curvature and b imply positive fifth coefficient.
    if np.any(p[:,5]<=0):raise ValueError('Degenerate MM polynomial outside enabled prior scope')
    companion=np.zeros((len(q),5,5));companion[:,0,:]=-p[:,4::-1]/p[:,5,None]
    companion[:,np.arange(1,5),np.arange(4)]=1.
    roots=np.linalg.eigvals(companion)
    legal=abs(roots.imag)<=1e-7*(1+abs(roots.real))
    roots=roots.real*scale[:,None]
    legal&=(roots>=lo[:,None])&(roots<=hi[:,None])
    candidates=np.concatenate((old[:,None],lo[:,None],hi[:,None],np.where(legal,roots,old[:,None])),axis=1)
    d=abs(1+b[:,None]*candidates)**2
    objective=(n1[:,None]*candidates+n2[:,None]*candidates**2)/d+.5*h[:,None]*candidates**2-s[:,None]*candidates
    chosen=candidates[np.arange(len(q)),np.argmin(objective,axis=1)]
    cost['aplus_conditional_polynomial_roots']+=roots.size
    cost['aplus_conditional_real_box_roots']+=int(legal.sum())
    cost['aplus_conditional_candidates']+=candidates.size
    return chosen


def mobius_material_readout(blocks,J,Q,m,beta,Os,Ls,rho,bounds,reg,config,cost):
    """Drop-in exact-loss four-color Möbius/MM material readout.

    Default two sweeps, fixed before any full development result. Root candidates
    minimize the conditional MM surrogate up to polynomial floating-point error;
    old material and exact final-loss fallback enforce observed nonincrease.
    """
    started=time.process_time();start_evals=cost['material_inner_evaluations']
    current,old,seed_value,value=_seed(blocks,J,Q,m,beta,Os,Ls,rho,bounds,reg,config,cost)
    context=_prior_context(config,blocks);prior=context['prior'];idx=blocks.indices
    lo,hi=bounds;b=float(beta)-blocks.diag_g;bb=b[idx]
    # Original source metric is represented by its exact supplied Cholesky Ls.
    Jw,Qw=J@Ls,Q@Ls
    cost['source_metric_products']+=2
    H=(Qw+b[:,None]*Jw)[idx];Jb=Jw[idx];G=blocks.goff
    z=current[idx]/(1+bb*current[idx]);V=Jb-z[:,:,None]*H
    R=float(beta)*(V+z[:,:,None]*(G@V))
    sweeps=int(config.get('aplus_sweeps',2));updates=0
    for sweep in range(sweeps):
        before=current.copy()
        for color in range(4):
            cells=idx[:,color]
            # Residual affine coefficient includes every affected row/source.
            B=-float(beta)*z[:,:,None]*G[:,:,color,None]*H[:,color,None,:]
            B[:,color,:]=float(beta)*(-H[:,color,:]+(G@V)[:,color,:])
            A=R-z[:,color,None,None]*B
            q=rho*np.sum(abs(B)**2,axis=(1,2))
            c=rho*np.sum(B.conj()*A,axis=(1,2))
            padded=np.r_[current,0.];neighbors=context['neighbors'][cells]
            other=padded[neighbors];csq=context['squared'][cells]
            weights=(2*prior.smooth+prior.tv/np.sqrt(csq*(current[cells,None]-other)**2+prior.epsilon**2))*csq/context['edges']
            h=weights.sum(axis=1);s=(weights*other).sum(axis=1)
            new=_conditional_candidates(q,c,b[cells],h,s,lo[cells],hi[cells],current[cells],cost)
            newz=new/(1+b[cells]*new)
            R+=(newz-z[:,color])[:,None,None]*B
            current[cells]=new;z[:,color]=newz
            V[:,color,:]=Jb[:,color,:]-newz[:,None]*H[:,color,:]
            updates+=len(cells)
            cost['aplus_conditional_color_batches']+=1
        if np.max(abs(current-before))<=float(config.get('aplus_step_tol',1e-8)):break
    final=value(current);fallback=not np.isfinite(final) or final>min(old,seed_value)
    if fallback:
        current=(m.copy() if old<=seed_value else _baseline_readout(Jw,Qw,m,beta,np.ones(J.shape[1]),rho,bounds,reg,False,config,cost)[0])
        final=min(old,seed_value)
    cost['aplus_readouts']+=1;cost['aplus_sweeps']+=sweep+1;cost['aplus_coordinate_updates']+=updates
    cost['material_inner_fallbacks']+=int(fallback)
    cost['aplus_readout_cpu_microseconds']+=int(1e6*(time.process_time()-started))
    return current,dict(material_inner_nit=sweep+1,material_inner_nfev=3,
        material_inner_total_evaluations=cost['material_inner_evaluations']-start_evals,
        material_inner_success=True,material_inner_message='four-color Möbius conditional MM; numerical stationary roots/endpoints/current',
        material_old_block_loss=old,material_scalar_seed_block_loss=seed_value,
        material_final_block_loss=final,material_fallback=bool(fallback),aplus_kind='mobius_mm',
        aplus_sweeps=sweep+1,aplus_readout_cpu_seconds=time.process_time()-started)


def projected_material_readout(blocks,J,Q,m,beta,Os,Ls,rho,bounds,reg,config,cost):
    """Exact-adjoint projected descent control with positive BB diagonal scale."""
    started=time.process_time();start_evals=cost['material_inner_evaluations']
    current,old,seed_value,value=_seed(blocks,J,Q,m,beta,Os,Ls,rho,bounds,reg,config,cost)
    lo,hi=bounds;previous=None;step=1.;accepted=0;backtracks=0
    for iteration in range(int(config.get('aplus_projected_steps',8))):
        f,g=onehop_material_value_gradient(blocks,J,Q,current,beta,Os,rho,reg,cost)
        if previous is not None:
            dx=current-previous[0];dg=g-previous[1];den=np.dot(dx,dg)
            if den>0:step=float(np.clip(np.dot(dx,dx)/den,1e-5,1e8))
        else:
            step=1/max(float(np.max(abs(g))),1e-10)
        previous=(current.copy(),g.copy())
        for trial in range(16):
            candidate=np.clip(current-step*g,lo,hi);delta=candidate-current
            if np.max(abs(delta))<1e-10:break
            candidate_value=value(candidate)
            if candidate_value<=f+1e-4*np.dot(g,delta):
                current=candidate;accepted+=1;break
            step*=.5;backtracks+=1
        else:break
        if np.max(abs(delta))<1e-10:break
    final=value(current);fallback=not np.isfinite(final) or final>min(old,seed_value)
    if fallback:current=m.copy();final=old
    cost['aplus_readouts']+=1;cost['aplus_projected_accepted_steps']+=accepted
    cost['aplus_projected_backtracks']+=backtracks
    cost['aplus_readout_cpu_microseconds']+=int(1e6*(time.process_time()-started))
    return current,dict(material_inner_nit=accepted,material_inner_nfev=cost['material_inner_evaluations']-start_evals,
        material_inner_total_evaluations=cost['material_inner_evaluations']-start_evals,
        material_inner_success=True,material_inner_message='exact-adjoint projected BB Armijo',
        material_old_block_loss=old,material_scalar_seed_block_loss=seed_value,
        material_final_block_loss=final,material_fallback=bool(fallback),aplus_kind='projected_adjoint',
        aplus_readout_cpu_seconds=time.process_time()-started)
