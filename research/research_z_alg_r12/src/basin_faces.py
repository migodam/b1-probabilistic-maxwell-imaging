"""Validate the first box-face event of a measured-source continuation."""
import argparse
import json
from time import perf_counter

from basin_sources import OUT, Wave, bounds, hermitian, load_endpoint, np, source_directions
from scipy.linalg import eigh


def face_state(w,x,W,free,full=False):
    r,E,P,K=w.state(x,with_kernel=True)
    rw=r@W
    _,gx,_=w.reg.evaluate(x)
    grad=((P.T@rw.conj())*E).sum(axis=1).real/w.den+gx
    jj=np.flatnonzero(free)
    D=w.reg.D; de=D@x
    mult=(2*w.reg.smooth+w.reg.tv*w.reg.epsilon**2/(de*de+w.reg.epsilon**2)**1.5)/w.reg.edges
    target=np.arange(w.n) if full else jj
    gn=(P[:,target].conj().T@P[:,jj])*((E@W)[jj]@E[target].conj().T).T
    rr=P.T@rw.conj()
    v=rr[target]@E[jj].T*K[np.ix_(target,jj)]
    u=rr[jj]@E[target].T*K[np.ix_(jj,target)]
    h=(gn+v+u.T).real/w.den+(D[:,target].T@D[:,jj].multiply(mult[:,None])).toarray()
    return grad,h


def normal_point(w,seed,W,free):
    x=seed.copy()
    for it in range(14):
        g,h=face_state(w,x,W,free)
        gn=np.linalg.norm(g[free],np.inf)
        if gn<1.e-13:
            break
        step=np.linalg.solve(h,-g[free])
        # Nearby nondegenerate stationary branch; not a free optimization restart.
        if np.max(abs(step))>2.:
            step*=2./np.max(abs(step))
        x[free]+=step
    low=np.flatnonzero(~free & (seed < (w.lo+w.hi)/2))
    high=np.flatnonzero(~free & (seed >= (w.lo+w.hi)/2))
    ii=np.flatnonzero(free)
    val=np.r_[g[low],-g[high],x[ii]-w.lo,w.hi-x[ii]]
    idx=np.r_[low,high,ii,ii]
    kinds=['release_lower']*len(low)+['release_upper']*len(high)+['hit_lower']*len(ii)+['hit_upper']*len(ii)
    return x,dict(values=val,index=idx,kinds=kinds,free_gradient=float(gn),
                  free_min_eigenvalue=float(eigh(h,eigvals_only=True,subset_by_index=[0,0])[0]))


def trace_case(label):
    t0=perf_counter()
    _,pub,x=load_endpoint(label); w=Wave(pub); I=np.eye(w.s)
    free=(x>w.lo+1e-8)&(x<w.hi-1e-8)
    x0=x.copy()
    changes=[]
    for _ in range(16):
        x0,z0=normal_point(w,x0,I,free)
        bad=np.flatnonzero(z0['values'] < 0)
        if not len(bad):break
        for j in bad:
            index=int(z0['index'][j]); kind=z0['kinds'][j]
            changes.append(dict(index=index,kind=kind,initial_margin=float(z0['values'][j])))
            if kind.startswith('hit'):
                x0[index]=w.lo if kind.endswith('lower') else w.hi
                free[index]=False
            else:
                free[index]=True
    if min(z0['values'])<=0 or z0['free_min_eigenvalue']<=0:
        raise ValueError("The polished state does not satisfy strict-face hypotheses")
    geometry,weights,intern=source_directions(w,x0)
    g,h=face_state(w,x0,I,free,True)
    low=np.flatnonzero(~free & (x0 < (w.lo+w.hi)/2))
    high=np.flatnonzero(~free & (x0 >= (w.lo+w.hi)/2))
    rows={}
    for name,a in weights.items():
        coherent=np.outer(a,a.conj()); delta=coherent-I
        _,ga=w.objective(x0,coherent)
        v=np.linalg.solve(h[free],-(ga-g)[free])
        dg=ga-g+h@v
        rate=np.r_[dg[low],-dg[high],v,-v]
        eta=np.where(rate<0,z0['values']/np.maximum(-rate,1e-300),np.inf)
        chosen=int(np.argmin(eta)); tlin=float(eta[chosen])
        tiny=min(1.e-7,tlin/20)
        xt,zt=normal_point(w,x0,I+tiny*delta,free)
        slope=(zt['values'][chosen]-z0['values'][chosen])/tiny
        lo=0.; hi=max(tlin,1e-8)
        while hi<.25:
            xh,zh=normal_point(w,x0,I+hi*delta,free)
            if np.min(zh['values'])<=0:break
            hi*=2
        for _ in range(20):
            mid=(lo+hi)/2
            xm,zm=normal_point(w,x0,I+mid*delta,free)
            if np.min(zm['values'])<0: hi=mid
            else:lo=mid
        xh,zh=normal_point(w,x0,I+hi*delta,free)
        j=int(np.argmin(zh['values']))
        rows[name]=dict(pred_t=tlin,pred_index=int(z0['index'][chosen]),pred_kind=z0['kinds'][chosen],
            pred_rate=float(rate[chosen]),fd_rate=float(slope),rate_relative_error=float(abs(slope/rate[chosen]-1)),
            observed_t=[lo,hi],observed_index=int(zh['index'][j]),observed_kind=zh['kinds'][j],
            all_margins_valid_below=bool(np.min(z0['values'])>0),free_min_eigenvalue_at_event=zh['free_min_eigenvalue'],
            free_gradient=zh['free_gradient'])
    record={'label':label,'needed_to_polish_active_set':changes,'polish_max_material':float(np.max(abs(x0-x))),
            'starting_free_gradient':z0['free_gradient'],'starting_hessian_min':z0['free_min_eigenvalue'],
            'rows':rows,'cost':{'seconds':perf_counter()-t0,'LU':w.lu,'RHS':w.rhs}}
    (OUT/f'face_{label}.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print(label,[(k,round(v['pred_t'],8),round(v['observed_t'][1],8),v['observed_kind']) for k,v in rows.items()],flush=True)


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--labels',nargs='+',default=['ring','shell'])
    for name in cli.parse_args().labels:trace_case(name)
