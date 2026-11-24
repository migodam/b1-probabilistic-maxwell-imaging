"""Multistart projected Gauss–Newton and bounded minimum-peak continuation."""

import argparse
import json

import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from search import (HERE, Model, SIDE, N, full_ring_inputs, groups, coarsen, bounded_support,
                    box_fit, clocks, elapsed)


def damped_box_step(m, jac, data_rhs, cap, ridge_relative):
    j = np.vstack([jac.real,jac.imag])
    rhs = np.concatenate([data_rhs.real,data_rhs.imag])
    s0 = j @ j.T
    lips = eigh(s0,subset_by_index=(s0.shape[0]-1,s0.shape[0]-1),eigvals_only=True)[0]
    ridge = ridge_relative*lips
    x = z = m.copy()
    alpha = 1.
    for _ in range(144):
        dz = z-m
        g = j.T @ (j@dz-rhs)+ridge*dz
        x_next = np.clip(z-g/(lips+ridge),0,cap)
        alpha_next = .5*(1+np.sqrt(1+4*alpha*alpha))
        z = x_next+((alpha-1)/alpha_next)*(x_next-x)
        x,alpha = x_next,alpha_next
    return x-m


def gauss_newton(model, material, side, active, cap, steps, radius):
    start=clocks()
    index=groups(side)
    dims=side*side
    mat=np.clip(coarsen(material,index,dims),0,cap)
    restrict=np.ones(dims,bool)
    restrict[active]=False
    mat[restrict]=0.
    model.count['native_full_wave_calls']+=1
    state= model.linear_solve(mat[index],mat[index,None]*model.b,False,None)
    ports=model.s @ state
    fr=float(np.linalg.norm(model.y-ports))
    hist=[dict(iteration=0,residual=fr,peak=float(mat.max()))]
    adj_0=None
    stop='Gauss–Newton budget'
    for step in range(1,steps+1):
        if fr < .998*radius:
            stop='reached numerical interior of noise ball'
            break
        m = mat[index]
        total=model.b+model.green(state)
        adj_0=model.linear_solve(m,model.s.conj().T,True,adj_0)
        jac_full=np.einsum('ir,is->rsi',adj_0.conj(),total).reshape(-1,N)
        if side < SIDE:
            jac=np.stack([np.bincount(index,weights=z.real,minlength=dims)+
                          1j*np.bincount(index,weights=z.imag,minlength=dims) for z in jac_full])[:,active]
        else:
            jac=jac_full[:,active]
        model.count['analytic_jacobians']+=1
        before=fr
        rho= max(.00025, .013*max(.025, fr/np.linalg.norm(model.y))**1.5)
        for damping in [rho,12*rho,120*rho]:
            direction=damped_box_step(mat[active],jac,(model.y-ports).ravel(),cap,damping)
            model.count['projected_linear_steps']+=144
            moved=False
            for damping_step in [1.,.5,.25,.125]:
                trial=mat.copy()
                trial[active]+=damping_step*direction
                trial=np.clip(trial,0,cap)
                tm=trial[index]
                model.count['native_full_wave_calls']+=1
                try:
                    jv=model.linear_solve(tm,tm[:,None]*model.b,False,state)
                except ArithmeticError:
                    model.count['rejected_line_solves']+=1
                    continue
                y_trial=model.s@jv
                new_residual=float(np.linalg.norm(model.y-y_trial))
                hist.append(dict(iteration=step,residual=new_residual,damping=damping,
                                 fraction=damping_step,accepted=new_residual<fr,peak=float(trial.max())))
                if new_residual<fr:
                    mat,state,ports,fr=trial,jv,y_trial,new_residual
                    moved=True
                    break
                model.count['rejected_line_solves']+=1
            if moved:
                break
        if before-fr <= 2e-8*np.linalg.norm(model.y):
            stop='local step stalled'
            break
    return {'m':mat[index],'ports':ports,'j':state,'residual':fr}, {
        'stop':stop,'level':side,'cap':cap,'best_peak':float(mat.max()),
        'best_residual':fr,'positive_cells':int(np.count_nonzero(mat[index])),
        'iterations':hist,**elapsed(start)}


def candidate_starts(obs,args):
    if args.cie is not None:
        path=args.cie/(obs.id+'.npz')
        if path.exists():
            with np.load(path,allow_pickle=False) as saved:
                yield 'historical_dual_cie', saved['m_seed']
    y,x=np.indices((96,96))/48.-1.
    shape=np.exp(-.5*((x/.38)**2+(y/.38)**2)).ravel()
    for amplitude in args.smooth:
        yield f'center_gaussian_{amplitude:g}',amplitude*shape


def run_one(s,p,b,obs,args):
    start=clocks()
    model=Model(s,p,b,obs.y)
    histories=[]
    best=None
    staged=[]
    for name,seed in candidate_starts(obs,args):
        result,hist=gauss_newton(model,seed,24,np.arange(24**2),args.cap,args.starts_iterations,obs.padded_radius)
        hist['name']=name
        histories.append(hist)
        staged.append((name,result))
        print(obs.id,name,hist['best_residual'],hist['best_peak'],hist['wall_seconds'],flush=True)
        if best is None or result['residual']<best['residual']:
            best=result
    target=obs.padded_radius*.998
    selected=best if best['residual']<target else None
    for name,st in sorted(staged,key=lambda p:p[1]['residual'])[:args.refine_seeds]:
        state = st
        for side in (48,96):
            support=bounded_support(state['m'],side,args.max_cells)
            state,row = gauss_newton(model,state['m'],side,support,args.cap,args.refinement_iterations,obs.padded_radius)
            row['name']=name
            histories.append(row)
            print(obs.id,name,side,row['best_residual'],row['best_peak'],row['wall_seconds'],flush=True)
            if state['residual']<best['residual']:
                best=state
            if state['residual']<target:
                direction_rank = (state['m'].max(),np.count_nonzero(state['m']))
                prior_rank = (selected['m'].max(),np.count_nonzero(selected['m'])) if selected else None
                if prior_rank is None or direction_rank[0]<prior_rank[0] or (
                        args.cap_path=='near' and direction_rank<prior_rank):
                    selected=state
    continuation=[]
    if selected is not None:
        high = float(selected['m'].max())
        decrement = min(.2, .06*high)
        low = .5
        active=bounded_support(selected['m'],96,args.max_cells)
        for _ in range(args.peaks):
            mid = high-.35*(high-low) if args.cap_path=='wide' else max(.5,high-decrement)
            trial,row=box_fit(model,selected['m'],active,96,mid,target,
                             args.peak_calls,args.peak_cpu,True)
            if trial['residual']<best['residual']:
                best=trial
            if trial['residual']<target:
                selected=trial
                high=float(selected['m'].max())
            else:
                if args.cap_path=='wide':
                    low = mid
                else:
                    decrement *= .5
            continuation.append(row)
            print(obs.id,'minpeak',row,flush=True)
    eligible=selected is not None
    if not eligible:
        selected=best
    np.savez_compressed(args.out/(obs.id+'.npz'),m_selected=selected['m'],
                        Y_selected_numeric=selected['ports'],j_selected_numeric=selected['j'],
                        m_best_residual=best['m'],Y_best_residual_numeric=best['ports'],Y=obs.y)
    receipt=dict(observation=obs.id,radius=obs.padded_radius,selected_numeric_residual=selected['residual'],
                 selected_peak=float(selected['m'].max()),selected_positive_cells=int(np.count_nonzero(selected['m'])),
                 eligible_for_independent_validation=eligible,stages=histories,minpeak_continuation=continuation,
                 counts=dict(model.count),**elapsed(start))
    (args.out/(obs.id+'.json')).write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print('DONE',obs.id,receipt['selected_peak'],receipt['selected_numeric_residual'],receipt['wall_seconds'],flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=type(HERE),default=HERE/'gn')
    parser.add_argument('--cie',type=type(HERE))
    parser.add_argument('--object',action='append')
    parser.add_argument('--noise',type=float,action='append')
    parser.add_argument('--smooth',type=float,default=[.9,2.4],nargs='*')
    parser.add_argument('--cap',type=float,default=4.)
    parser.add_argument('--refine-seeds',type=int,default=1)
    parser.add_argument('--starts-iterations',type=int,default=16)
    parser.add_argument('--refinement-iterations',type=int,default=16)
    parser.add_argument('--max-cells',type=int,default=2700)
    parser.add_argument('--peaks',type=int,default=5)
    parser.add_argument('--peak-calls',type=int,default=48)
    parser.add_argument('--peak-cpu',type=float,default=24.)
    parser.add_argument('--cap-path',choices=('near','wide'),default='near')
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    s,p,b,data=full_ring_inputs()
    for obs in data:
        if obs.noise_pct in (args.noise or [1.,2.]) and (not args.object or obs.object_id in args.object):
            run_one(s,p,b,obs,args)


if __name__=='__main__':
    with threadpool_limits(1):
        main()
