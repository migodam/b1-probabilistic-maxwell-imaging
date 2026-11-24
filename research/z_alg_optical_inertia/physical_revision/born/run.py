"""Distance from measured complex ports to the first-Born data of any nonnegative image."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name] = '1'
import csv
import json
from pathlib import Path
from time import perf_counter, process_time
import numpy as np
from scipy.linalg import solve
from flint import acb, acb_mat, arb, ctx

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
ctx.prec = 80
ctx.threads = 1

def clocks():
    return np.array([perf_counter(), process_time()])
def mat(x):
    return acb_mat([[acb(complex(v)) for v in row] for row in x])
def realvector(x):
    return np.r_[x.real.ravel(),x.imag.ravel()]
def complexmatrix(x, shape):
    nr = int(np.prod(shape))
    return (x[:nr]+1j*x[nr:]).reshape(shape)
def polar_fit(matrix, y, interior, margins, cap, iterations):
    # Projected accelerated gradient on the primal yields useful dual residuals.
    # Rescaling the negative interior direction makes all cell inner products nonpositive.
    lip = np.linalg.norm(matrix @ matrix.T, ord=2)
    z = np.zeros(matrix.shape[1]); x=z.copy(); t=1.
    best = (-np.inf,None)
    for i in range(iterations+1):
        if i in (0,10,40,iterations):
            v=y-matrix@x
            g=matrix.T@v
            if cap is None:
                if np.max(margins)>=0: raise RuntimeError('No strictly negative cone direction')
                shift=max(0.,float(np.max(g/(-margins))))
                v=v+(shift+max(1.e-9,abs(shift)*1.e-8))*interior
                obj=float(y@v) / np.linalg.norm(v)
            else:
                obj=(float(y@v)-cap*np.maximum(g,0.).sum())/np.linalg.norm(v)
            if obj>best[0]:
                best=(obj,v.copy())
        if i==iterations: break
        grad=matrix.T@(matrix@z-y)
        next_x=np.maximum(0., z-grad/lip)
        if cap is not None: np.minimum(cap,next_x,out=next_x)
        nt=(1+np.sqrt(1+4*t*t))/2
        z=next_x+((t-1)/nt)*(next_x-x);x=next_x;t=nt
    return best
def interval_bounds(direction,y,bs,bb,cap):
    zw=mat(direction); yw=mat(y)
    dot=(zw.conjugate().transpose()*yw)
    pair=sum((dot[i,i].real for i in range(dot.nrows())),arb(0))
    norm=sum((abs(zw[i,j])**2 for i in range(zw.nrows()) for j in range(zw.ncols())),arb(0)).sqrt()
    pos=arb(0); maximum=None
    for s,b in zip(bs,bb):
        v=s*zw
        for i in range(v.nrows()):
            h=sum(((v[i,j]*b[i,j]).real for j in range(v.ncols())),arb(0))
            maximum=h.upper() if maximum is None else maximum.max(h.upper())
            if cap is not None: pos+=h.max(0)
    if cap is None and not maximum<0: return 0.,float(maximum.upper())
    value=(pair-(cap or 0)*pos)/norm
    return max(0.,np.nextafter(float(value.lower()),-np.inf)),np.nextafter(float(maximum.upper()),np.inf)
def main():
    start=clocks()
    with np.load(STUDY/'synthetic/public/fine_geometry.npz') as f:
        S,B=(f['PW360_10x20__'+key] for key in ('S','B'))
    # Real coordinates are fixed as receiver-major/source-minor.
    complex_a=np.einsum('pi,ij->pji',S,B).reshape(S.shape[0]*B.shape[1],S.shape[1])
    ar=np.ascontiguousarray(np.r_[complex_a.real,complex_a.imag]); del complex_a
    gram=ar@ar.T
    interior=solve(gram+1.e-9*np.linalg.norm(gram,2)*np.eye(len(gram)), -ar.sum(axis=1),assume_a='pos')
    margins=ar.T@interior
    print('interior span',margins.min(),margins.max(),flush=True)
    bs=[mat(S[:,i:i+384].conj().T) for i in range(0,len(B),384)]
    bb=[mat(B[i:i+384].conj()) for i in range(0,len(B),384)]
    setup=clocks()-start
    curves={}
    with (STUDY/'submission/curve/noise_curve.csv').open() as f:
        for r in csv.DictReader(f): curves[r['id'] if 'id' in r else next(iter(r.values()))]=r
    out=[];directions={}
    paths=sorted((STUDY/'submission/curve/observations').glob('*pair[12]pct.npz'))
    for path in paths:
        y=np.load(path)['Y']; yv=realvector(y)
        entry=curves[path.stem]
        rad=float(next(v for k,v in entry.items() if k in ('radius','r','joint_radius','noise_radius')))
        for cap in (4.,None):
            key=path.stem+('__cone' if cap is None else '__box4')
            t=clocks()
            estimate,v=polar_fit(ar,yv,interior,margins,cap,64)
            search=clocks()-t;z=complexmatrix(v,y.shape)
            t=clocks();d,maxg=interval_bounds(z,y,bs,bb,cap);audit=clocks()-t
            out.append(dict(observation=path.stem,cap='nonnegative, unbounded' if cap is None else cap,radius=rad,
                      rigorous_lower=d,lower_over_radius=d/rad,max_cell_pair_upper=maxg,excludes=d>rad,
                      search_wall=search[0],search_cpu=search[1],interval_wall=audit[0],interval_cpu=audit[1]))
            directions[key]=z
        print(path.stem,[(a['cap'],round(a['lower_over_radius'],3)) for a in out[-2:]],flush=True)
        with (HERE/'born_bounds.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
    np.savez_compressed(HERE/'directions.npz',**directions)
    (HERE/'cost.json').write_text(json.dumps(dict(precision_bits=ctx.prec,thread_count=1,setup_wall_cpu=setup.tolist(),total_wall_cpu=(clocks()-start).tolist(),n_observations=len(paths),interior_range=[float(margins.min()),float(margins.max())],using_truth=False),indent=2)+'\n')
if __name__=='__main__': main()
