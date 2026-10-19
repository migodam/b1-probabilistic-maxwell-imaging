"""Material-task exploration. These constrained maps are candidate-only."""
import json,time
import numpy as np
from scipy.optimize import minimize
from p1_material.proposals import frozen_latent_step
from .core import OUT,write,save,mpack,BlockStop

METHODS=['random_matched','native','hard_prior','soft_Maxwell','soft_frozen','soft_dynamic','soft_adaptive','full_J']

def geometry(t,p,s,beta=1):
    t.jacobian(p);D=t.prior.derivative(p.q);M=np.r_[D.real,D.imag]/24
    U,d,Vh=np.linalg.svd(M,full_matrices=False);keep=d>d[0]*1e-8
    T=Vh[keep].T/d[keep];Q=U[:,keep];K=p.A@T;R=T[:64].T@T[:64]
    lam,V=np.linalg.eigh(K.T@K+beta*R)
    if lam[0]<-1e-7*max(1,lam[-1]):raise RuntimeError('material observability not PSD')
    lam=np.maximum(lam,0);weak=1/(1+lam*s*s)
    return dict(T=T,Q=Q,K=K,lam=lam,V=V,weak=weak,material_singular=d,material_rank=int(sum(keep)),cutoff=1/(s*s))

def feasible_direction(t,p,v):
    v=np.array(v,copy=True)
    for k in [64,65]:
        if (p.q[k]>=1-1e-10 and v[k]>0) or (p.q[k]<=1e-10 and v[k]<0):v[k]=0
    w=t.prior.derivative(p.q)@v;n=np.linalg.norm(w)/24
    if n<1e-14:return None,None
    return v/n,w/n

def direction(t,p,g,method,rng):
    if method=='random_matched':
        v=t.prior.local_prior_reversible(p.q,rng,beta=.15,property_sd=.04)-p.q
    elif method in ['native','full_J']:
        v=rng.normal(size=66);v[64:]*=.04
    else:
        z=rng.normal(size=len(g['lam']));weight=(g['lam']<=g['cutoff']).astype(float) if method=='hard_prior' else g['weak']
        v=g['T']@(g['V']@(np.sqrt(weight)*z))
    return feasible_direction(t,p,v)

def constrained_fit(t,model,origin,velocity,s,seed,objective='data',g=None,corrector=False,progress_required=True):
    """Same exact decoder/progress constraints; no arbitrary material clipping."""
    P=t.prior;x0=origin.state.chi;goal=x0+s*velocity;q0=origin.q;local=max(1.,s/.02)
    y=t.observation(t.parent);r=t.white(model.state.data-y);scale=max(1.,.5*r@r)
    cache={}
    def values(q):
        if 'q' not in cache or not np.array_equal(q,cache['q']):
            x,_=P.decode(q);D=P.derivative(q);delta=x-x0
            cache.update(q=q.copy(),x=x,D=D,delta=delta)
        return cache
    def f(q):
        z=values(q)
        if objective=='native':return float(np.vdot(z['x']-goal,z['x']-goal).real/(576*s*s)+1e-8*np.sum((q-q0)**2))
        rr=r+t.material_action(model,z['x']-model.state.chi)
        value=.5*rr@rr+.5*q[:64]@q[:64]
        if corrector and g is not None:
            c=g['V'].T@g['Q'].T@mpack(z['x']-model.state.chi)
            value+=2*np.dot(g['weak']*c,c)/(s*s)
        return float(value/scale)
    def jf(q):
        z=values(q);D=z['D']
        if objective=='native':return 2*np.real(D.conj().T@(z['x']-goal))/(576*s*s)+2e-8*(q-q0)
        rr=r+t.material_action(model,z['x']-model.state.chi);Z=model.C@D
        A=np.column_stack([t.white(Z[:,k].reshape(32,8)) for k in range(66)])
        grad=A.T@rr+np.r_[q[:64],0.,0.]
        if corrector and g is not None:
            m=np.r_[D.real,D.imag]/24;c=g['V'].T@g['Q'].T@mpack(z['x']-model.state.chi)
            grad+=4*m.T@g['Q']@g['V']@(g['weak']*c)/(s*s)
        return grad/scale
    def sphere(q):d=values(q)['delta'];return float(np.vdot(d,d).real/(576*s*s)-1)
    def js(q):z=values(q);return 2*np.real(z['D'].conj().T@z['delta'])/(576*s*s)
    def prog(q):return float(np.vdot(velocity,values(q)['delta']).real/(576*s)-.5)
    def jp(q):return np.real(values(q)['D'].conj().T@velocity)/(576*s)
    def ball(q):return float(1-np.sum((q[:64]-q0[:64])**2)/(local*local))
    def jb(q):return np.r_[-2*(q[:64]-q0[:64])/(local*local),0.,0.]
    constraints=[dict(type='eq',fun=sphere,jac=js),dict(type='ineq',fun=ball,jac=jb)]
    if progress_required:constraints.append(dict(type='ineq',fun=prog,jac=jp))
    fit=minimize(f,seed,jac=jf,bounds=[(None,None)]*64+[(0,1)]*2,constraints=constraints,method='SLSQP',options=dict(ftol=1e-8,maxiter=70))
    valid=lambda q:P.legal(q) and abs(sphere(q))<=.0404 and (not progress_required or prog(q)>=-1e-7) and ball(q)>=-1e-7
    q=fit.x
    fallback=False
    if not valid(q) and valid(seed):q=seed.copy();fallback=True
    return q,dict(optimizer_success=bool(fit.success),message=str(fit.message),iterations=int(fit.nit),radius_matched=bool(valid(q)),actual_RMS=P.material_rms(P.decode(q)[0],x0),signed_progress=float((prog(q)+.5)*s),decoder_fallback_to_legal_seed=fallback)

def block(t,parent,anchor,method,seed,radii=(.02,.04),version='v1',maxstates=4,maxrhs=128):
    book=t.book;t.parent=parent;y=t.observation(parent);start=time.perf_counter();before=book.used[book.stage].copy();rng=np.random.default_rng(seed)
    t.jacobian(anchor);model=anchor;current=anchor;rows=[];lastg=None;accepted_since=0;previous_v=None;previous_w=None
    with book.scope(maxstates,maxrhs):
        for s in radii:
            refresh=method=='soft_dynamic' or (method=='soft_adaptive' and (lastg is None or accepted_since>=3))
            if refresh:
                t.jacobian(current);model=current
            basis=current if method in ['soft_dynamic','soft_adaptive'] else anchor
            if method=='soft_adaptive' and lastg is not None and not refresh:model=lastg['model']
            # Current decoder geometry can be rebuilt with cached full-material J without a new solve.
            if basis is not model:
                from .core import Point
                proxy=Point(basis.q,basis.state,basis.key,T=model.T)
                proxy.C=model.C;Z=model.C@t.prior.derivative(basis.q);proxy.A=np.column_stack([t.white(Z[:,k].reshape(32,8)) for k in range(66)])
                g=geometry(t,proxy,s,beta=0 if method=='soft_Maxwell' else 1)
            else:g=geometry(t,basis,s,beta=0 if method=='soft_Maxwell' else 1)
            if previous_w is not None and method in ['soft_dynamic','soft_adaptive']:
                c=g['Q'].T@mpack(previous_w)
                weight=g['weak'];v=g['T']@(g['V']@(np.sqrt(weight)*(g['V'].T@c)))
                v,w=feasible_direction(t,basis,v)
            elif previous_v is not None:v,w=previous_v.copy(),previous_w.copy()
            else:v,w=direction(t,basis,g,method,rng)
            if v is None:
                rows.append(dict(target_RMS=s,status='NO_FEASIBLE_MATERIAL_VELOCITY'));continue
            # Destination remains measured from the original block anchor, not accumulated latent length.
            v,w=feasible_direction(t,anchor,v)
            if v is None:continue
            tangent_alignment=None if previous_w is None else float(np.vdot(previous_w,w).real/(24*24))
            previous_v=v.copy();previous_w=w.copy()
            qseed,matched=frozen_latent_step(t.prior,anchor.q,v,s)
            fit={};q=qseed
            if method in ['native','full_J','soft_dynamic','soft_adaptive']:
                q,fit=constrained_fit(t,model,anchor,w,s,qseed,objective='native' if method=='native' else 'data',g=g)
                matched=fit['radius_matched']
            if not matched:
                rows.append(dict(target_RMS=s,status='DECODER_DISTANCE_OR_PROGRESS_FAILURE',fit=fit));continue
            try:
                pred=t.evaluate(q,f'{parent}_{method}_{version}_s{s}_seed{seed}_predictor')
                actual=t.white(pred.state.data-model.state.data);mlin=t.material_action(model,pred.state.chi-model.state.chi);latent=model.A@(pred.q-model.q)
                eta=float(np.linalg.norm(actual-mlin)/max(1.,np.linalg.norm(actual)));coordinate=float(np.linalg.norm(mlin-latent)/max(1.,np.linalg.norm(actual)))
                rr=t.record(pred,parent,method,anchor,seed=int(seed),version=version,target_RMS=s,role='predictor',fit=fit,Maxwell_secant_relative_error=eta,decoder_secant_relative_error=coordinate,tangent_alignment=tangent_alignment)
                rows.append(rr);best=pred
                if book.used[book.stage]['STATE']-before['STATE']<maxstates:
                    # Dynamic corrector uses its own paid receiver J; frozen controls retain the anchor model.
                    cm=pred if method in ['soft_dynamic','soft_adaptive'] else model
                    t.jacobian(cm)
                    cq,ci=constrained_fit(t,cm,anchor,w,s,pred.q,g=g if method.startswith('soft_') or method=='hard_prior' else None,corrector=True)
                    if ci['radius_matched'] and np.linalg.norm(cq-pred.q)>1e-10:
                        corrected=t.evaluate(cq,f'{parent}_{method}_{version}_s{s}_seed{seed}_corrector')
                        rows.append(t.record(corrected,parent,method,anchor,seed=int(seed),version=version,target_RMS=s,role='corrector',fit=ci))
                        if t.objective(corrected,y)<t.objective(best,y):best=corrected
                current=best;t.persist_cache(current);accepted_since+=int(t.qualified(current,y))
                if method=='soft_adaptive':
                    if version=='v2':trigger=eta>.25
                    else:trigger=float(np.linalg.norm(actual-latent)/max(1.,np.linalg.norm(actual)))>.25
                    actual_drop=t.objective(model,y)-t.objective(best,y)
                    residual=t.white(model.state.data-y);ml=t.material_action(model,best.state.chi-model.state.chi)
                    predicted_drop=t.objective(model,y)-(.5*np.linalg.norm(residual+ml)**2+.5*best.q[:64]@best.q[:64])
                    ratio=actual_drop/predicted_drop if abs(predicted_drop)>1e-12 else 1.
                    if trigger or not .5<=ratio<=1.5:
                        t.jacobian(current);model=current;accepted_since=0
                lastg=dict(model=model)
            except BlockStop:
                rows.append(dict(target_RMS=s,status='BLOCK_PHYSICS_LIMIT'));break
    summary=dict(parent=parent,method=method,seed=int(seed),version=version,rows=rows,states=book.used[book.stage]['STATE']-before['STATE'],RHS=book.used[book.stage]['RHS']-before['RHS'],wall_seconds=time.perf_counter()-start,GT_used=False)
    write(OUT/f'{book.stage}/blocks/{parent}_{method}_{version}_{seed}.json',summary)
    return summary

def stage_B(t):
    parents=json.loads((OUT/'configs/PARENTS.json').read_text());anchors=json.loads((OUT/'A/ANCHORS.json').read_text());summary=[]
    dev=[r for r in parents if r['split']!='heldout'];test=[r for r in parents if r['split']=='heldout']
    # One counted finite-difference/JVP audit on the first qualified training anchor.
    audit=None
    for r in dev:
        p=t.load(anchors[r['parent']]);y=t.observation(r['parent'])
        if not t.qualified(p,y):continue
        t.jacobian(p);v=np.random.default_rng(2026101902).normal(size=66);v[64:]=0;v/=np.linalg.norm(v);h=1e-5
        plus=t.evaluate(p.q+h*v,'B_counted_derivative_plus');minus=t.evaluate(p.q-h*v,'B_counted_derivative_minus')
        with t.book.physical('counted_derivative_JVP'):
            native=t.backend().anchor_jvp(t.prior.derivative(p.q)@v,p.state)
        fd=t.white((plus.state.data-minus.state.data)/(2*h));a=p.A@v;load=np.random.default_rng(19).normal(size=(32,8))+1j*np.random.default_rng(20).normal(size=(32,8))
        dot=t.backend().cached_material_dot_identity(t.prior.derivative(p.q)@v,load,p.state)
        audit=dict(parent=r['parent'],FD_relative_error=float(np.linalg.norm(fd-a)/max(1,np.linalg.norm(a))),native_JVP_relative_error=float(np.linalg.norm(t.white(native)-a)/max(1,np.linalg.norm(a))),adjoint=dot,new_states=2,new_RHS=24)
        write(OUT/'B/DERIVATIVE_AUDIT.json',audit)
        if audit['FD_relative_error']>1e-5 or audit['native_JVP_relative_error']>1e-10:raise RuntimeError('new backend derivative numerical mismatch')
        break
    for rr in dev:
        parent=rr['parent'];p=t.load(anchors[parent]);y=t.observation(parent)
        if not t.qualified(p,y):summary.append(dict(parent=parent,status='INVALID_BASE',methods_not_run=METHODS));continue
        for m in METHODS:summary.append(block(t,parent,p,m,2026101910))
        print(json.dumps(dict(stage='B',parent=parent,consumed=dict(t.book.used['B']))),flush=True)
    # Predeclared single mechanism revision: decoder/Maxwell-split refresh trigger.
    revise=[r for r in dev if (r['split']=='train' and r['index']<4) or r['split']=='validation']
    for rr in revise:
        p=t.load(anchors[rr['parent']]);y=t.observation(rr['parent'])
        if not t.qualified(p,y):continue
        for m in ['soft_dynamic','soft_adaptive']:summary.append(block(t,rr['parent'],p,m,2026101910,version='v2'))
    write(OUT/'B/METHODS_FROZEN_BEFORE_HELDOUT.json',dict(revision='adaptive trigger isolates Maxwell from decoder secant error',versions={m:'v2' if m in ['soft_dynamic','soft_adaptive'] else 'v1' for m in METHODS},heldout_GT_read=False))
    for rr in dev:
        if rr['split']!='validation':continue
        p=t.load(anchors[rr['parent']]);y=t.observation(rr['parent'])
        if not t.qualified(p,y):continue
        for m in ['native','soft_frozen','soft_dynamic','full_J']:summary.append(block(t,rr['parent'],p,m,2026101912,radii=(.08,),version='stress',maxstates=2,maxrhs=64))
    for rr in test:
        parent=rr['parent'];p=t.load(anchors[parent]);y=t.observation(parent)
        if not t.qualified(p,y):summary.append(dict(parent=parent,status='INVALID_BASE',methods_not_run=METHODS));continue
        for seed in ([2026101910,2026101911] if rr['index']<4 else [2026101910]):
            for m in METHODS:summary.append(block(t,parent,p,m,seed,version='v2' if m in ['soft_dynamic','soft_adaptive'] else 'v1'))
        print(json.dumps(dict(stage='B',parent=parent,consumed=dict(t.book.used['B']))),flush=True)
    write(OUT/'B/SUMMARY.json',dict(blocks=summary,derivative_audit=audit,GT_online=False,posterior_samples=False))
