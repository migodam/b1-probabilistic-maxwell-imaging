"""Compute a source-only obstruction and compare to the data-dependent conditional operator."""
from pathlib import Path
import json,time,sys
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits
from flint import arb,acb,acb_mat,ctx
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]; STUDY=ROOT.parent
sys.path[:0]=[str(ROOT/'src'),str(STUDY/'ip_guided_v4/boxes')]
from fixed_trace import block,diag,mat,adj,positive_dimension
from run import CONTROL_ID,load,herm
from verify import ProofGrid,ball,diagonal_times,gersh
def exact_sign(a):
    p,n=positive_dimension(a)
    return [p,n] if p+n==a.nrows() else None
def large_sign(a,approx):
    ev,v=eigh(herm(approx),driver='evd',check_finite=False)
    n=len(v);nu=int(sum(ev<0))
    for interval,polarity in ((range(0,nu),-1),(range(nu,n),1)):
        if not len(interval):continue
        cols=ball(v[:,interval])
        projected=polarity*(adj(cols)*(a*cols))
        # Gershgorin in bounded memory; definite restricted forms imply independent columns.
        if not gersh(projected,upper=False)>0:raise ArithmeticError('full sign unresolved')
    return nu
def hcat(a,b):return acb_mat(a.nrows(),a.ncols()+b.ncols(),[z for x,y in zip(a.tolist(),b.tolist()) for z in x+y])
def persist(result): (HERE/'source_index.json').write_text(json.dumps(result,indent=2))
def run():
    spec=json.loads((HERE/'contract.json').read_text()); base=STUDY/'ip_guided_v4/boxes/results'
    loaded=load(base,[CONTROL_ID])[CONTROL_ID]; r=np.load(STUDY/'ip_final_v5/strong_control/strong_ring_arrays.npz')
    readings={'disk':(loaded.Y,loaded.radius),'hollow ring':(r['Y'],float(r['radius'][0]))}
    ids=r['ids'];n=len(ids); order=sorted(range(32),key=lambda i:int(f'{i:05b}'[::-1],2))
    numeric,(g,s,b)=ProofGrid(loaded.geometry).restricted(ids);gn,sn,bn=numeric
    ss=adj(s)*s; nn=sn.conj().T@sn
    out={'method':'source-only congruence vs directly conditioned receiver; both strict sign counts','n':n,'rows':[]}
    if '--resume' in sys.argv:
        out=json.loads((HERE/'source_index.json').read_text())
    for name,width,mu in spec['cases']:
        if any(q['class']==name for q in out.get('cost',[])):continue
        t0=time.perf_counter(); t1=time.process_time()
        if name.startswith('solid'):
            l=np.maximum(0,r['c_solid']-width);u=r['c_solid']+width
        else:l=np.zeros(n);u=np.full(n,float(width))
        # Expand the sector directly, minimizing the cost of dense one-sided Green products.
        mid=[(arb(x)+arb(y))/2 for x,y in zip(l,u)]
        prod=[arb(x)*arb(y) for x,y in zip(l,u)]
        mg=diagonal_times(mid,g)
        physical=acb(mu)*(diag([1]*n)-mg-adj(mg)+adj(g)*diagonal_times(prod,g))
        nf=acb(mu)*(diagonal_times(mid,b)-adj(g)*diagonal_times(prod,b))
        cb=acb(mu)*adj(b)*diagonal_times(prod,b)
        gl=((l+u)/2)[:,None]*gn
        physnp=mu*(np.eye(n)-gl-gl.conj().T+gn.conj().T@((l*u)[:,None]*gn))
        k=physical+ss;np0=large_sign(physical,physnp);nu=large_sign(k,physnp+nn)
        qphys=cb-adj(nf)*physical.solve(nf)
        mix=hcat(nf,adj(s)); d=acb_mat(64,64)
        for i in range(32):
            d[i+32,i+32]=1
            for j in range(32):d[i,j]=cb[i,j]
        q=d-adj(mix)*k.solve(mix)
        for obj,(y,rad) in readings.items():
            for count in spec['old_source_counts']:
                pick=order[:count]; v=np.eye(32)[:,pick];z=ball(np.vstack([v,y[:,pick]]))
                omega=diag([arb(rad)**2]*count)
                sp=block(qphys,pick,pick)-omega
                sc=exact_sign(sp); po=adj(z)*q*z-omega; oc=exact_sign(po)
                row=dict(box=name,object=obj,s=count,nu=nu,physical_negative=np0,source_inertia=sc,old_inertia=oc,receiver_directions=32)
                if sc and oc:
                    h=np0-sc[0]; expected=32-h+nu-oc[0]
                    ez=block(q,range(32,64),range(64))*z
                    ry=block(q,range(32,64),range(32,64))-ez*po.solve(adj(ez))
                    rc=exact_sign(ry)
                    row.update(deficit_source_only=h,old_saturated=oc[0]==nu,predicted_receiver_positive=expected,direct_receiver_inertia=rc)
                    if rc:assert rc[0]==expected
                out['rows'].append(row)
        out.setdefault('cost',[]).append({'class':name,'wall':time.perf_counter()-t0,'cpu':time.process_time()-t1})
        persist(out);print(name,np0,nu,out['cost'][-1],flush=True)
if __name__=='__main__':
    ctx.prec=128;ctx.threads=3
    with threadpool_limits(limits=2):run()
