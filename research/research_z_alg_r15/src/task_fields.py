"""Full-wave derivatives for a blank plate and translated localized material features."""
from itertools import product
from pathlib import Path
from time import perf_counter, process_time
import json

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from threadpoolctl import threadpool_limits
from ndt import Air, BACKGROUND, FREQUENCIES_GHZ, grid

OUT = Path(__file__).resolve().parents[1] / 'results'
CENTERS = np.array(list(product((-6., 0., 6.), repeat=2)))  # mm, public before inversion
AMPLITUDE = -1.
RADIUS_MM = 1.55
TRANSITION_MM = .9


def feature(xy, center):
    """A compact, differentiable virtual flaw; center and derivatives use mm."""
    off = xy * 1000. - center
    radius = np.linalg.norm(off, axis=1)
    t = np.clip((RADIUS_MM + .5 * TRANSITION_MM - radius) / TRANSITION_MM, 0., 1.)
    weight = t*t*(3. - 2.*t)
    dw = (6.*t*(1.-t)/TRANSITION_MM)[:, None] * off / np.maximum(radius[:, None], 1.e-30)
    return weight, dw


def local_field(q, derivatives, g, e, s):
    """Exact source and target derivatives near a defect of compact support."""
    a = np.eye(len(q)) - q[:,None]*g
    fact = lu_factor(a, check_finite=False)
    current = lu_solve(fact, q[:,None]*e, check_finite=False)
    total = e + g@current
    receiver = lu_solve(fact, s.T, trans=1, check_finite=False).T
    jac = np.einsum('ri,it,ik->rtk', receiver,total,derivatives,optimize=True)
    return s@current, jac


def collect(n, f, virtual, check):
    air = Air(n, f)
    green = air.dense()
    fac = lu_factor(np.eye(n*n)-BACKGROUND*green, check_finite=False)
    eh = lu_solve(fac, air.E0, check_finite=False)
    sh = lu_solve(fac, air.S.T, trans=1, check_finite=False).T
    # No n^2-by-n^2 inverted host Green is needed for these features.
    index = np.flatnonzero(np.any([v[0]>0 for v in virtual], axis=0))
    gh = lu_solve(fac,green[:,index],check_finite=False)[index]
    jac_pixel = np.einsum('ri,it->rti',sh,eh)
    if n!=24:
        ratio=n//24
        jac_pixel=jac_pixel.reshape(16,8,24,ratio,24,ratio).sum(axis=(3,5)).reshape(16,8,576)
    fields=[]; derivatives=[]
    for weight, dw in virtual:
        a=np.flatnonzero(weight[index]>0)
        ids=index[a]
        dirs=np.column_stack((AMPLITUDE*dw[ids],weight[ids]))
        y,j=local_field(AMPLITUDE*weight[ids],dirs,gh[np.ix_(a,a)],eh[ids],sh[:,ids]);fields.append(y);derivatives.append(j)
    error=None
    if check:
        # The translated test stays in the existing union, checked explicitly.
        h=1.e-5; sec=[]
        for axis in range(2):
            diff=[]
            for sign in (-1,1):
                v=np.zeros(2);v[axis]=h*sign
                wt,_=feature(air.xy,CENTERS[0]+v); ids=np.flatnonzero(wt>0);a=np.searchsorted(index,ids)
                assert np.array_equal(index[a],ids)
                y,_=local_field(AMPLITUDE*wt[ids],np.ones((len(ids),1)),gh[np.ix_(a,a)],eh[ids],sh[:,ids]);diff.append(y)
            d=(diff[1]-diff[0])/(2*h)
            sec.append(float(np.linalg.norm(d-derivatives[0][...,axis])/np.linalg.norm(d)))
        error=sec
    return jac_pixel,np.stack(fields),np.stack(derivatives),error


def main():
    timing=[];checks={};arrays={'centers_mm':CENTERS,'frequencies_ghz':FREQUENCIES_GHZ}
    for n in (24,48):
        virtual=[feature(grid(n),pos) for pos in CENTERS]
        pixels=[];ys=[];js=[]
        for f in FREQUENCIES_GHZ:
            t,c=perf_counter(),process_time()
            p,y,j,err=collect(n,f,virtual,check=f==16.)
            pixels.append(p);ys.append(y);js.append(j)
            timing.append(dict(grid=n,f=float(f),wall_s=perf_counter()-t,cpu_s=process_time()-c))
            if err: checks[str(n)]=err
        arrays.update({f'pixels{n}':np.array(pixels), f'y{n}':np.array(ys), f'jac{n}':np.array(js)})
    np.savez_compressed(OUT/'task_fields.npz',**arrays)
    (OUT/'task_fields_cost.json').write_text(json.dumps(dict(timing=timing,translation_derivative_relative_error=checks,amplitude=AMPLITUDE,radius_mm=RADIUS_MM,transition_mm=TRANSITION_MM),indent=2)+'\n')
    print('Derivative checks',checks,'CPU seconds',sum(x['cpu_s'] for x in timing),flush=True)

if __name__=='__main__':
    with threadpool_limits(1): main()
