"""Re-fit actual data using a two-ball trust-region step in 130 real measurement coordinates."""
from run import (HERE, SIDE, Fit, controller, elapsed, np, proof, read_pair, real, timer)
from pathlib import Path
import json
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits


def compensation_step(j,resid,lo,hi,trust,price=.002):
    """Dual of boxed quadratic minimum displacement with two distance-to-ball losses."""
    lower=np.maximum(lo,-trust); upper=np.minimum(hi,trust)
    design=np.ascontiguousarray(j)
    def obj(y):
        z=-(y@design)
        e=np.clip(z,lower,upper)
        val=float(z@e-.5*(e@e)-resid@y+.5*price*(y@y))
        grad=-resid-design@e+price*y
        for sl in (slice(0,128),slice(128,130)):
            size=np.sqrt(y[sl]@y[sl]+1e-14)
            val+=.85*size
            grad[sl]+=.85*y[sl]/size
        return val,grad
    ans=minimize(obj,np.zeros(130),method='L-BFGS-B',jac=True,options={'maxiter':400,'ftol':1e-15,'gtol':2e-7,'maxcor':24})
    step=np.clip(-(ans.x@design),lower,upper)
    return step,dict(nfev=int(ans.nfev),dual_dimension=130,gradient_norm=float(np.linalg.norm(ans.jac)),stop=str(ans.message))


def run(tag, resume=False):
    clock=timer()
    prop,dec,observed,tx,rx,ap,_=controller()
    _,_,roi_mask,r,_=read_pair()
    roi=np.flatnonzero(roi_mask)
    y=observed[tag]
    with np.load(HERE/f'{tag}_candidate.npz',allow_pickle=False) as p:
        x=p['material_opposite_class']
        old=p['searched_territory']
    if resume:
        with np.load(HERE/f'dual_{tag}.npz',allow_pickle=False) as p:
            x=p['material']
        old=np.unique(np.r_[old,np.flatnonzero(x)])
    fit=Fit(y,r,tx,rx,ap,prop['amplitude'])
    history=[];cost=[];territory=old;mask=np.ones(x.size,bool);mask[roi]=False
    free=np.flatnonzero(mask)
    trust=.70;last=np.inf;stop='iteration budget'
    for i in range(18 if resume else 10):
        rs,blocks=fit.values(x)
        rr=np.concatenate(rs);jac=np.vstack(blocks)
        start=[float(np.linalg.norm(z)) for z in rs]
        val=max(start)
        if val<.995:stop='accepted two-ball point';break
        direction,all_rec=compensation_step(jac[:,free],rr,-x[free],4.5-x[free],trust)
        external=~np.isin(free,territory)
        # New sites are chosen by their bounded, simultaneous contribution in all measured coordinates.
        impact=np.abs(direction)*np.linalg.norm(jac[:,free],axis=0)
        ranked=np.argsort(np.where(external,impact,-1))[::-1]
        count=144 if i<4 else 64
        expanding=free[ranked[:count]]
        territory=np.unique(np.r_[territory,expanding])
        idx=territory[mask[territory]]
        delta,local_rec=compensation_step(jac[:,idx],rr,-x[idx],4.5-x[idx],trust)
        solvers=fit.trial_wave(territory)
        trials=[]
        for alpha in (1.,.5,.25,.125):
            xn=x.copy()
            xn[idx]=np.clip(xn[idx]+alpha*delta,0.,4.5)
            wave=[w.fields(xn[territory],jac=False)[0] for w in solvers]
            gaps=fit.distance(wave)
            trials.append((max(gaps),alpha,gaps))
            if max(gaps)<.995:break
        score,alpha,gaps=min(trials,key=lambda z:z[0])
        history.append(dict(iteration=i,initial=start,target_radius=.85,trust=trust,added_cells=int(len(expanding)),
                active_domain=int(len(idx)),whole_problem=all_rec,constrained_problem=local_rec,candidates=trials))
        print(tag,i,start,gaps,flush=True)
        if score<val:
            x[idx]=np.clip(x[idx]+alpha*delta,0.,4.5)
            if score<.995:stop='accepted two-ball point';break
            trust=min(1.,trust*(1.15 if alpha>.5 else .7))
        else:trust*=.45
        if trust<.06 or (i>=3 and abs(last-min(score,val))<2e-4):
            stop='local step ceased improving';break
        last=min(score,val)
    strict,fields,replay=proof(x,y,r,tx,rx,ap,prop['amplitude'])
    good=all(v['residual_upper_over_r']<=1. for v in strict.values())
    notes=dict(sample=tag,method='two-ball response-space dual for locally boxed compensation, priced linear subproblem then independent nonlinear replay',
           scope='full original material box, contiguous fixed-amplitude core; inner witness search only',
           positive=good,core=float(x[roi[0]]),core_cells=len(roi),same_uniform_core=bool(np.all(x[roi]==x[roi[0]])),
           status=stop,history=history,strict=strict,cost_full_jac=[q.cost for q in fit.jac],cost_line_search=fit.trial_cost,**elapsed(clock))
    dest=HERE/(('continued_' if resume else 'dual_')+tag)
    np.savez_compressed(dest.with_suffix('.npz'),material=x,active=np.flatnonzero(x),roi_indices=roi,
                      observed_400=y[0],observed_460=y[1],**fields)
    dest.with_suffix('.json').write_text(json.dumps(notes,indent=2)+'\n')
    print('one-sided replay',tag,good,[v['residual_upper_over_r'] for v in strict.values()],flush=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('samples',nargs='+',choices=('unit_415','unit_862'))
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    with threadpool_limits(1):
        for id_ in args.samples:run(id_,args.resume)
