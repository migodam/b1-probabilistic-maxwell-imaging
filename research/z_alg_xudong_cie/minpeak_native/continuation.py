"""Reduce the peak of a same-ball material with native-differential box steps."""

import argparse
import json

import numpy as np
from threadpoolctl import threadpool_limits

from search import HERE, Model, full_ring_inputs, bounded_support, clocks, elapsed
from gn_recovery import gauss_newton, damped_box_step


def move(model, m, j, y, upper, active):
    start=clocks()
    incident=model.b+model.green(j)
    adj=model.linear_solve(m,model.s.conj().T,True,None)
    jac=np.einsum('ir,is->rsi',adj[active].conj(),incident[active]).reshape(-1,len(active))
    model.count['analytic_jacobians']+=1
    residual=(model.y-y).ravel()
    best=None
    attempts=[]
    for damping in [.0001,.001,.01]:
        direction=damped_box_step(m[active],jac,residual,upper,damping)
        model.count['projected_linear_steps']+=144
        tester=m.copy()
        tester[active]+=direction
        tester=np.clip(tester,0,upper)
        r,p,current=model.full(tester,gradient=False)
        attempts.append(dict(damping=damping,physical_residual=float(r)))
        if best is None or r<best[0]:
            best = float(r),tester,p,current
    return best,dict(new_cap=upper,trials=attempts,**elapsed(start))


def one(s,p,b,obs,args):
    start=clocks()
    model=Model(s,p,b,obs.y)
    path=args.input/(obs.id+'.npz')
    with np.load(path,allow_pickle=False) as f:
        m,j,y=f['m_selected'],f['j_selected_numeric'],f['Y_selected_numeric']
    norm=float(np.linalg.norm(y-obs.y))
    if norm >= obs.padded_radius:
        print(obs.id,'outside; no certified branch to continue',flush=True)
        return
    eta=obs.padded_radius*.997
    positive=np.flatnonzero(m)
    if len(positive)>args.max_cells:
        raise ValueError('Use a support budget covering the feasible starting material.')
    nearby=bounded_support(m,96,args.max_cells)
    extra=nearby[~np.isin(nearby,positive)][:args.max_cells-len(positive)]
    active=np.sort(np.concatenate([positive,extra]))
    first,interior=gauss_newton(model,m,96,active,float(m.max()),args.steps,.65*obs.padded_radius)
    if first['residual']<eta:
        m,j,y,norm=first['m'],first['j'],first['ports'],first['residual']
    stage=[]
    increment=args.increment
    for k in range(args.peaks):
        cap=max(.5,float(m.max())-increment)
        (r,try_m,try_p,try_j),linear=move(model,m,j,y,cap,active)
        refinement=None
        if r>=eta:
            check,refinement=gauss_newton(model,try_m,96,active,cap,args.steps,obs.padded_radius*.90)
            if check['residual']<r:
                r,try_m,try_p,try_j=check['residual'],check['m'],check['ports'],check['j']
        accepted=r<eta
        if accepted:
            m,y,j,norm=try_m,try_p,try_j,r
        else:
            increment*=.5
        stage.append(dict(cap=cap,final_residual=r,accepted=accepted,linear=linear,refinement=refinement))
        print(obs.id,cap,r,accepted,flush=True)
        payload=dict(m_selected=m,j_selected_numeric=j,Y_selected_numeric=y,Y=obs.y)
        np.savez_compressed(args.output/(obs.id+'.npz'),**payload)
        receipt=dict(observation=obs.id,radius=obs.padded_radius,selected_numeric_residual=norm,
                     selected_peak=float(m.max()),selected_positive_cells=int(np.count_nonzero(m)),
                     eligible_for_independent_validation=bool(norm<eta),interior=interior,
                     stages=stage,counts=dict(model.count),**elapsed(start))
        (args.output/(obs.id+'.json')).write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=type(HERE),default=HERE/'gn')
    parser.add_argument('--output',type=type(HERE),default=HERE/'continuation')
    parser.add_argument('--object',action='append')
    parser.add_argument('--noise',action='append',type=float)
    parser.add_argument('--increment',type=float,default=.15)
    parser.add_argument('--peaks',type=int,default=6)
    parser.add_argument('--steps',type=int,default=8)
    parser.add_argument('--max-cells',type=int,default=2800)
    args=parser.parse_args()
    args.output.mkdir(exist_ok=True,parents=True)
    with threadpool_limits(1):
        s,p,b,data=full_ring_inputs()
        for obs in data:
            if obs.noise_pct in (args.noise or [1.,2.]) and (not args.object or obs.object_id in args.object):
                one(s,p,b,obs,args)
