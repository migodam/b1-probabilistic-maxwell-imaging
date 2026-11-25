"""Restrict the constitutive quadratic to actual affine Fourier current spaces."""
from pathlib import Path
import json,sys
from time import perf_counter,process_time
import numpy as np
from flint import acb_mat,arb,acb,ctx
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits
HERE=Path(__file__).resolve().parent
STUDY=HERE.parents[1]
sys.path.insert(0,str(STUDY/'ip_unified/numerics/native'))
from replay import configure,ball,star,low,high,mid,gershgorin,herm
def clock(): return perf_counter(),process_time()
def elapsed(t): return dict(zip(('wall','cpu'),[a-b for a,b in zip(clock(),t)]))
def restrict_array(geo,Y,U,Jd,mu):
    n,l=U.shape;s=Y.shape[1];span=np.column_stack((U,Jd))
    field=geo.numerical_green(span);field[:,l:]+=geo.incident
    resid=geo.sensor@span;resid[:,l:]-=Y
    cross=span.conj().T@field
    return herm(resid.conj().T@resid+mu*(span.conj().T@span-herm(cross)+(cross-cross.conj().T)/1j)),resid
def propose(form,l,s,r):
    K=form[:l,:l];k,vec=eigh(K)
    neg=int(np.count_nonzero(k<0))
    if np.min(abs(k))<1e-22: return None
    L=form[:l,l:];H=herm(form[l:,l:]-L.conj().T@(vec@((vec.conj().T@L)/k[:,None])))
    e=eigh(H,eigvals_only=True)[::-1]
    return dict(inertia=neg,state_float_min_abs=float(min(abs(k))),rank=float(np.linalg.matrix_rank(K)),point=float(np.sqrt(np.maximum(0,e[neg:]).sum())/r))
def enclosure(geo,Y,U,Jd,mu):
    x=ball(np.column_stack((U,Jd)));l=U.shape[1];s=Y.shape[1]
    ef=geo.exact_green(x);rx=geo.sensor_ball*x
    y=ball(Y)
    for i in range(geo.n):
        for j in range(s):ef[i,l+j]+=geo.incident_ball[i,j]
    for i in range(y.nrows()):
        for j in range(s):rx[i,l+j]-=y[i,j]
    cross=star(x)*ef
    return star(rx)*rx+arb(mu)*(star(x)*x-(cross+star(cross))/2+acb(0,-1)*(cross-star(cross)))
def positive_count_lower(mat):
    a=mid(mat);v,w=eigh(herm(a))
    # All eigenvectors form an independent trial basis. Gershgorin on the
    # trailing subspace alone gives a lower bound on its positive dimension.
    q=ball(w);t=star(q)*mat*q
    for i in range(len(v)-1,-1,-1):
        a0=t[i,i].real
        if not (a0>0 and a0>sum((abs(t[i,j]) for j in range(len(v)) if j!=i),arb(0))): return len(v)-1-i
    return len(v)
def main():
    ctx.prec=112
    g,native,data,_=configure('u100_best')
    out={'scope':'exact represented saved CIE affine maps; not a zero-current condition per mode','cases':{}}
    for name in ('optfresh_01','optfresh_05'):
        meta,y=next(x for x in data if x[0]['id'].startswith(name));r=meta['padded_radius']
        tdir=HERE.parent/'routing/cie_diagnostic'
        path=tdir/(meta['id']+'__standard.npz')
        with np.load(path) as f: jd=f['Jd'];std=f['standard_U']
        with np.load(tdir/(meta['id']+'__equal_width_negative_swap.npz')) as f: sw=f['swapped_U']
        rows=[]
        for tag,u in [('fourier100',std),('swapped100',sw)]:
            t0=clock();mat,rx=restrict_array(g,y,u,jd,native.w.scale)
            l=u.shape[1];s=y.shape[1];numeric=propose(mat,l,s,r)
            er=rx[:,l:]-rx[:,:l]@np.linalg.lstsq(rx[:,:l],rx[:,l:],rcond=None)[0]
            rec=dict(space=tag,numeric_full_source=numeric,unconstrained_fit_over_radius=float(np.linalg.norm(er)/r),
                     unconstrained_fit='floating least squares; near-null direction conditioning matters',
                     numerical_cost=elapsed(t0))
            t0=clock(); exact=enclosure(g,y,u,jd,native.w.scale)
            tests=[]
            # One threshold in the true original physical units for each count.
            for j in (1,2,4,8,16,24,32):
                threshold=arb(float(r))**2*(1+arb(1)/65536)/j
                shifted=exact.copy() if hasattr(exact,'copy') else acb_mat(exact)
                for i in range(l,l+s):shifted[i,i]-=threshold
                pos=positive_count_lower(shifted)
                tests.append(dict(j=j,positive_count_lower=pos,required=l+j,excluded=bool(pos>=l+j)))
            rec.update(verified_tests=tests,strict_excludes=any(z['excluded'] for z in tests),outward_cost=elapsed(t0))
            rows.append(rec);print(meta['id'],tag,numeric,rec['strict_excludes'],rec['outward_cost'],flush=True)
            out['cases'][meta['id']]=rows
            (HERE/'restricted.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':
    with threadpool_limits(1): main()
