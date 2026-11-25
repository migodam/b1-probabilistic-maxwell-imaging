"""Two-weight reduced SDP proposal followed by original Fourier interval verification."""
import cvxpy as cp
from flint import ctx
import numpy as np
from scipy.linalg import eigh,solve

from search import (HERE, MAPS, Geometry, SectorGeometry, box, clock, elapsed, geom_path,
                    hermitian, read_spec, save_json, threadpool_limits)
from exact_replay import ExactPadded,cutoff


def pencil(sector,node,y):
    cells,profile,lo,hi=box(node)
    d=len(cells)
    gamma=sector.covariance(cells)
    cross,v=sector.selected(cells[None,:]); cross,v=cross[0],v[0]
    cov=np.block([[gamma,cross.conj().T],[cross,sector.p0]])
    values,vectors=eigh(hermitian(cov),check_finite=False)
    kept=values>np.max(values)*2.e-13
    root=vectors[:,kept]*np.sqrt(values[kept])[None,:]
    ns=len(root[0]); cols=y.shape[1]
    tv=np.column_stack((root[:2*d,:],v))
    sensor=np.column_stack((root[2*d:,:],-(y-sector.y0)))
    base=np.block([[np.eye(ns),np.zeros((ns,cols))],[np.zeros((cols,ns)),-sector.t0]])
    ident,zero=np.eye(d),np.zeros((d,d))
    dp=np.block([[-ident,.5*sector.b*ident],[.5*sector.b*ident,zero]])
    mid=np.diag((lo+hi)*profile*.5)
    prod=np.diag(lo*hi*profile**2)
    dt=np.block([[ident,-mid],[-mid,prod]])
    return (hermitian(sensor.conj().T@sensor),hermitian(base+tv.conj().T@dp@tv),
            hermitian(tv.conj().T@dt@tv),ns,int(np.count_nonzero(~kept)))


def run():
    spec=read_spec()
    node=next(x for x in spec['nodes'] if x['side']==32 and x['index']==510)
    start=clock()
    g=Geometry(geom_path(32)); sector=SectorGeometry(g,spec['b'])
    with np.load(MAPS/'public'/f'{node["oid"]}.npz',allow_pickle=False) as data:
        y,r=data['Y'],float(data['r'])
    C,Q,V,n,dropped=pencil(sector,node,y)
    prep=elapsed(start)
    # Freeze the source projection at the independently verified positive-state seed.
    mu0=.013619899148313314
    ratio0=1/128
    old=C+mu0*Q+mu0*ratio0*V
    schur=old[n:,n:]-old[n:,:n]@solve(old[:n,:n],old[:n,n:],assume_a='her',check_finite=False)
    eig,w=eigh(hermitian(schur),check_finite=False)
    F=w[:,eig>0]
    rank=F.shape[1]
    W=np.block([[np.eye(n),np.zeros((n,rank))],[np.zeros((len(eig),n)),F]])
    z=cp.Variable((rank,rank),hermitian=True)
    a=cp.Variable(nonneg=True); b=cp.Variable(nonneg=True)
    # The two weights remain affine and independent. The background is scaled
    # by 2^-5 for the solver's variable units, and the footprint by 2^-10.
    mu=a/32; beta=b/1024
    def restrict(mat):
        return hermitian(W.conj().T@mat@W)
    form=restrict(C)+mu*restrict(Q)+beta*restrict(V)
    zborder=cp.bmat([[np.zeros((n,n)),np.zeros((n,rank))],[np.zeros((rank,n)),z]])
    state=C[:n,:n]+mu*Q[:n,:n]+beta*V[:n,:n]
    constraints=[form-zborder>>0,state-1.e-6*mu*np.eye(n)>>0,
                 a>=2.**-15,a<=4,beta>=2.**-15*mu,beta<=.125*mu]
    prob=cp.Problem(cp.Maximize(cp.real(cp.trace(z))),constraints)
    started=clock()
    try:
        prob.solve(solver='SCS',eps=2.e-6,max_iters=2400,verbose=False)
    except cp.SolverError as exc:
        save_json(HERE/'convex_probe.json',{'node':node,'error':str(exc),'setup':prep,'solver_cost':elapsed(started)})
        return
    cost=elapsed(started)
    stats=prob.solver_stats.extra_stats['info']
    stats={k:(None if isinstance(v,float) and not np.isfinite(v) else v) for k,v in stats.items()}
    result=dict(node=node,method='fixed source projection SDP with both sector weights affine',
                full_source_columns=y.shape[1],projection_columns=rank,compressed_current_variables=n,
                truncated_small_covariance_modes=dropped,numeric_solver=prob.status,
                numeric_objective=float(prob.value) if np.isfinite(prob.value) else None,
                setup_cost=prep,solver_cost=cost,solver_info=stats,
                exact_acceptance_requires_fresh_original_operator_replay=True)
    if a.value is not None and b.value is not None and float(a.value)>0:
        started=clock()
        m=float(a.value)/32; ratio=float(b.value)/(1024*m)
        result.update(total_scale=m,ratio=ratio)
        try:
            exact=ExactPadded(g,spec['b'])
            check=cutoff(exact.local(node,ratio),y,r,m)
        except (ArithmeticError,ValueError,ZeroDivisionError) as exc:
            check={'status':'UNRESOLVED','reason':str(exc)}
        result.update(exact_replay=check,exact_replay_cost=elapsed(started))
    save_json(HERE/'convex_probe.json',result)
    print(prob.status,result.get('ratio'),result.get('total_scale'),result.get('exact_replay'),flush=True)


if __name__=='__main__':
    ctx.prec,ctx.threads=192,1
    with threadpool_limits(limits=1):run()
