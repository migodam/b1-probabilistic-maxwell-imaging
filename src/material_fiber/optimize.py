"""Full-physics, GT-blind inversion and directly constrained GN subproblems."""
import numpy as np
from scipy.optimize import minimize,lsq_linear
from .core import OUT,write,BlockStop

PDIAG=np.r_[np.ones(64),np.zeros(2)]

def projected_gradient(q,g):
    pg=g.copy()
    pg[64:][((q[64:]<=1e-10)&(g[64:]>=0))|((q[64:]>=1-1e-10)&(g[64:]<=0))]=0
    return pg

def step(target,p,y,method,radius,iteration=0):
    A=target.jacobian(p);q=p.q;r=target.white(p.state.data-y)
    g=A.T@r+PDIAG*q;H=A.T@A+np.diag(PDIAG)
    damping=max(1e-8,.001*float(np.median(np.diag(H))))
    lo=np.r_[np.full(64,-np.inf),np.maximum(-q[64:],-.2)];hi=np.r_[np.full(64,np.inf),np.minimum(1-q[64:],.2)]
    if method=='fixed_LM':
        MM=np.r_[A,np.diag(PDIAG),np.sqrt(damping)*np.eye(66)];b=np.r_[-r,-PDIAG*q,np.zeros(66)]
        lo0=np.r_[np.full(64,-np.inf),-q[64:]];hi0=np.r_[np.full(64,np.inf),1-q[64:]]
        d=lsq_linear(MM,b,bounds=(lo0,hi0),tol=1e-10,max_iter=300).x
        d*=min(1.,.5/max(1e-30,np.linalg.norm(d[:64])),.2/max(1e-30,np.max(abs(d[64:]))))
        status='bounded_LM_postscale'
    else:
        # Coarse-to-fine only restricts optimizer iterates; the target/prior is unchanged.
        free=np.arange(66)
        if method=='coarse_TR' and iteration<4:
            n=4 if iteration<2 else 6
            modes=[(k,l) for k in range(8) for l in range(8) if k or l]
            free=np.array([0]+[i+1 for i,(k,l) in enumerate(modes) if k<n and l<n]+[64,65])
        hs=H[np.ix_(free,free)]+damping*np.eye(len(free));gg=g[free]
        # Diagonal scaling is numerical only: constraints are transformed exactly.
        scale=1/np.sqrt(np.maximum(np.diag(hs),1e-8));B=np.diag(scale)
        hs=B@hs@B;gg=B@gg;bounds=list(zip(lo[free]/scale,hi[free]/scale))
        def fun(t):return float(gg@t+.5*t@hs@t)
        def jac(t):return gg+hs@t
        if method=='material_TR':
            D=target.prior.derivative(q);M=np.r_[D.real,D.imag]/24;metric=(M[:,free]*scale).T@(M[:,free]*scale)
            cap=radius*.08
        else:
            metric=np.diag((free<64)*scale*scale);cap=radius
        def sphere(t):return float(cap*cap-t@metric@t)
        def sj(t):return -2*metric@t
        fit=minimize(fun,np.zeros(len(free)),jac=jac,bounds=bounds,constraints=[dict(type='ineq',fun=sphere,jac=sj)],method='SLSQP',options=dict(maxiter=120,ftol=1e-9))
        d=np.zeros(66);d[free]=scale*fit.x;status=str(fit.message)
        if sphere(fit.x)<-max(1e-9,1e-5*cap*cap):return None,dict(subproblem=status,invalid_constraint=True)
    pred=float(-g@d-.5*d@H@d)
    pg=projected_gradient(q,g);kkt=float(np.max(abs(pg)/np.sqrt(np.maximum(np.diag(H),1))))
    return d,dict(predicted_reduction=pred,raw_gradient_inf=float(max(abs(pg))),scaled_KKT=kkt,subproblem=status,H_condition=float(np.linalg.cond(H+1e-10*np.eye(66))),Gaussian_step=float(np.linalg.norm(d[:64])),property_step=d[64:].tolist(),radius=radius)

def reconstruct(target,parent,method,max_states=16,initial=None,run=None,iteration_offset=0):
    book=target.book;y=target.observation(parent);u0=book.used[book.stage].copy();run=run or f'{parent}_{method}'
    history=[];radius=.5
    with book.scope(max_states,40*max_states+32):
        p=initial if initial is not None else target.evaluate(np.r_[np.zeros(64),.5,.5],run+'_independent_mean')
        best=p;best_obj=target.objective(p,y);reason='STATE_LIMIT'
        target.persist_cache(p)
        for k in range(max_states+1):
            d,info=step(target,p,y,method,radius,k+iteration_offset)
            row=dict(iteration=k,key=p.key,Phi=target.phi(p,y),objective=target.objective(p,y),qualified=target.qualified(p,y),**info)
            history.append(row);write(OUT/f'{book.stage}/histories/{run}.json',history)
            if info.get('scaled_KKT',1)<1e-3:reason='NUMERICAL_STATIONARY';break
            if d is None:reason='SUBPROBLEM_FAILURE';break
            if book.used[book.stage]['STATE']-u0['STATE']>=max_states:break
            if np.linalg.norm(d)<1e-10 or info.get('predicted_reduction',0)<=1e-12:reason='NO_PREDICTED_DESCENT';break
            new=None;oldobj=target.objective(p,y)
            for alpha in ([1.,.5,.25] if method=='fixed_LM' else [1.]):
                if book.used[book.stage]['STATE']-u0['STATE']>=max_states:break
                q=p.q+alpha*d
                violation=max(0.,float(-min(q[64:])),float(max(q[64:])-1))
                if violation and violation<=64*np.finfo(float).eps:
                    # Only numerical active-bound restoration in latent coordinates.
                    # The original prior and physical material bounds are unchanged.
                    old=q[64:].copy();q[64:]=np.maximum(0.,np.minimum(1.,q[64:]))
                    book.append('ROUND_OFF_BOUND_RESTORATION',max_latent_change=float(max(abs(q[64:]-old))),threshold=float(64*np.finfo(float).eps))
                if not target.prior.legal(q):raise RuntimeError('optimizer outside original P1 support')
                try:new=target.evaluate(q,run+f'_trial{k}_alpha{alpha}')
                except BlockStop:break
                if target.objective(new,y)<oldobj:break
            if new is None:break
            oldobj=target.objective(p,y);newobj=target.objective(new,y);rho=(oldobj-newobj)/max(info['predicted_reduction'],1e-30)
            row.update(trial_key=new.key,trial_objective=newobj,actual_predicted_ratio=float(rho),accepted=bool(newobj<oldobj))
            if method!='fixed_LM':
                if rho<.25:radius=max(.03125,radius/2)
                elif rho>.75:radius=min(4.,radius*2)
            if newobj<oldobj:
                p=new
                if target.qualified(p,y) and not target.qualified(best,y) or (target.qualified(p,y)==target.qualified(best,y) and newobj<best_obj):best=p;best_obj=newobj
            else:
                if method=='fixed_LM':radius/=2  # fixed LM remains a fixed-cap baseline; rejection is retained.
                if radius<=.03125:reason='TR_RADIUS_LIMIT';break
        target.jacobian(best);target.persist_cache(best)
        final=target.record(best,parent,method,phase='reconstruction',run=run,stop_reason=reason)
    final.update(history=str((OUT/f'{book.stage}/histories/{run}.json').relative_to(OUT)),incremental_states=book.used[book.stage]['STATE']-u0['STATE'],incremental_RHS=book.used[book.stage]['RHS']-u0['RHS'])
    write(OUT/f'{book.stage}/runs/{run}.json',final)
    return best,final
