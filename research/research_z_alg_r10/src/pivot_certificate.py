"""Rigorous low-dimensional elimination using Bernstein coefficient lower bounds."""
import argparse
from itertools import permutations
import json
import math
import time

import numpy as np
from scipy.linalg import qr
from flint import arb, acb, acb_mat
from bounded_certificate import E, GAMMA, ROOT, ball, matrix, row, square, dot, feasible_witness


def determinant(cols):
    d=len(cols)
    total=acb(0)
    for perm in permutations(range(d)):
        inversions=sum(perm[i]>perm[j] for i in range(d) for j in range(i+1,d))
        term=acb((-1)**inversions)
        for j,i in enumerate(perm):
            term*=cols[j][i]
        total+=term
    return total


def cell_tables(H,R,B,host,y,box):
    k,m=B.shape
    idx=qr(B.T,pivoting=True,mode="economic")[2][:m].tolist()
    free=[i for i in range(k) if i not in idx]
    inv=matrix(B[idx]).inv()
    N=[]
    for g in free:
        nx=row(B[g])*inv
        N.append([acb(int(i==g)) if i in free else -nx[0,idx.index(i)] for i in range(k)])
    H=matrix(H);R=list(row(R).entries())
    hx,yy=list(row(host).entries()),list(row(y).entries())
    lo,hi=map(arb,box)
    table=[]
    for i in range(k):
        line=[]
        for q in (lo,(lo+hi)/2,hi):
            u=[acb(int(j==i))-q*H[j,i] for j in range(k)]
            xx=list((inv*acb_mat([[u[j]] for j in idx])).entries())
            c=dot(hx,xx)+q*R[i]
            line.append(([dot(n,u) for n in N], xx, c, dot(yy,xx)-c))
        table.append(line)
    return table


def tensor(table, ids):
    dim=len(ids)
    grid=np.empty((3,)*dim, dtype=object)
    for index in np.ndindex(grid.shape):
        terms=[table[i][level] for i,level in zip(ids,index)]
        col=[t[0] for t in terms]
        weights=[(-1)**a*determinant(col[:a]+col[a+1:]) for a in range(dim)]
        res=dot(weights,[t[3] for t in terms])
        pred=dot(weights,[t[2] for t in terms])
        power=sum((square(dot(weights,items)) for items in zip(*(t[1] for t in terms))),arb(0))
        grid[index]=square(res)-E**2*power-GAMMA**2*square(pred)
    # Degree-two Bernstein controls from values at endpoints and midpoint, in each coordinate.
    for axis in range(dim):
        g=np.moveaxis(grid,axis,0)
        g[1]=2*g[1]-(g[0]+g[2])/2
    largest=max(float(abs(v).upper()) for v in grid.flat)
    if largest==0:
        return np.zeros(grid.shape,dtype=np.int64)
    factor=arb(2)**(45-math.frexp(largest)[1])
    vals=[int((v*factor).lower().floor().unique_fmpz()) for v in grid.flat]
    if max(abs(a) for a in vals)>=1<<48:
        raise ValueError("Coefficient exceeds signed subdivision scale")
    return np.asarray(vals,dtype=np.int64).reshape(grid.shape)


def halves(array,axis):
    a,b,c=np.moveaxis(array,axis,0)
    ab=(a+b)//2
    bc=(b+c)//2
    mid=(ab+bc)//2
    return np.stack((a,ab,mid),axis=axis),np.stack((mid,bc,c),axis=axis)


def uniformly_positive(controls):
    # Positive quadratic on [0,1] iff a,c>0 and b>−sqrt(ac).
    # Integer rescaling below rounds favorably to retaining, never to excluding.
    a,b,c=np.moveaxis(controls,-1,0)
    bound=np.maximum(np.maximum(abs(a),abs(b)),abs(c))
    shift=np.maximum(np.frexp(bound.astype(float))[1].astype(np.int64)-29,0)
    up=(-b+(1<<shift)-1)>>shift
    ac=(a>>shift)*(c>>shift)
    positive=(a>0)&(c>0)&((b>=0)|(ac>up*up))
    return np.any(np.all(positive,axis=tuple(range(2,positive.ndim))),axis=1)


def cover(table, pivots, depth, coefficients=None):
    dim=len(pivots)
    ids=[[*pivots,j] for j in range(len(table)) if j not in pivots]
    if coefficients is None:
        base=np.stack([tensor(table,ix) for ix in ids])
    else:
        # A permutation of the cofactor columns only changes a common sign.
        base=np.stack([coefficients[tuple(sorted(ix))].transpose([sorted(ix).index(j) for j in ix]) for ix in ids])
    current=base[None,...]
    coords=np.zeros((1,dim),dtype=int)
    total=0
    for _ in range(depth):
        collected=[];addresses=[]
        for ib in range(0,len(coords),128):
            pieces=[(current[ib:ib+128],2*coords[ib:ib+128])]
            for axis in range(dim):
                nxt=[]
                direction=np.eye(dim,dtype=int)[axis]
                for vals,addr in pieces:
                    l,r=halves(vals,2+axis)
                    nxt.extend(((l,addr),(r,addr+direction)))
                pieces=nxt
            for val,addr in pieces:
                total+=len(addr)
                keep=~uniformly_positive(val)
                if np.any(keep):
                    collected.append(val[keep]);addresses.append(addr[keep])
        if not addresses:
            return np.empty((0,dim),dtype=int),total
        current=np.concatenate(collected);coords=np.concatenate(addresses)
    return coords,total


def intervals(indices,box,depth):
    mask=np.zeros(2**depth,dtype=bool)
    mask[np.asarray(indices,dtype=int)]=True
    segments=np.flatnonzero(np.diff(np.r_[False,mask,False].astype(int))).reshape(-1,2)
    return (box[0]+(box[1]-box[0])*segments/2**depth).tolist()


def record(H,R,B,host,y,box,pivots,depth,roots):
    tic=time.process_time()
    tab=cell_tables(H,R,B,host,y,box)
    pixels,n=cover(tab,pivots,depth)
    marginals=[intervals(pixels[:,i],box,depth) for i in range(len(pivots))]
    witnesses=[feasible_witness(z,H,R,B,host,y) for z in roots]
    return dict(depth=depth,pivots=[int(i) for i in pivots],box=box,outer=marginals,
                tested_boxes=n,retained_boxes=len(pixels),witnesses=witnesses,
                witness_coordinates=[[float(z[i]) for i in pivots] for z in roots],
                cpu_seconds=time.process_time()-tic),pixels


def short_experiment():
    results={}
    for k in (12,24,18):
        file=ROOT/"results"/f"short_instance_{k}_noisy.npz"
        if not file.exists():continue
        d=np.load(file)
        # Selected by the noisily inferred alternatives, independent of the true endpoints.
        roots=d['old_candidates']
        pivot=int(np.argmax(np.ptp(roots,axis=0))) if len(roots)>=2 else 0
        for index,label in enumerate(d['protocols']):
            if label=='repeat_old': continue
            B=d['B'][index]; y=d['observed'][index]; host=d['host'][index]
            result,_=record(d['H'],d['R'],B,host,y,(.5,12.),[pivot],12,d['z'])
            result['old_witnesses']=[feasible_witness(z,d['H'],d['R'],B[:,:-1],host[:-1],y[:-1]) for z in d['z']]
            results[f'{k}/{label}']=result
            print(k,label,pivot,result['outer'][0],flush=True)
    (ROOT/'results'/'short_certificate.json').write_text(json.dumps(results,indent=2)+'\n')


def natural_experiment(depth,only,include_old):
    f=np.load(ROOT/'results'/'natural_inputs.npz')
    results={};maps={}
    valid=f['case_ids'].tolist()
    names=valid if not only else only
    for name in names:
        H=f[name+'/H'];R=f[name+'/R'];pv=f[name+'/pivots'].astype(int).tolist()
        truth=f[name+'/truth']
        base=f[name+'/clean_witnesses']
        # First witness is the actual sample, followed by distinct old-data witnesses.
        zs=[truth,*[z for z in base if np.linalg.norm(z-truth)>1e-4]]
        for arm in ('active','random','fisher'):
            key=name+'/'+arm;B=f[key+'/B'];y=f[key+'/y'];hy=f[key+'/y_host']
            if 'fisher'==arm and np.array_equal(B,f[name+'/random/B']) and name+'/random' in results:
                results[key]=dict(results[name+'/random'],same_design_as='random')
                maps[key]=maps[name+'/random']
                continue
            rec,pix=record(H,R,B,hy,y,(.5,30.),pv,depth,zs)
            rec['old_witnesses']=[feasible_witness(z,H,R,B[:,:3],hy[:3],y[:3]) for z in zs]
            results[key]=rec;maps[key]=pix
            print(key,'pivots',pv,'boxes',rec['retained_boxes'],'outer',rec['outer'][0],flush=True)
        if include_old:
            oldpv=[*pv,next(i for i in range(len(R)) if i not in pv)]
            B=f[name+'/B0']; oldy=f[name+'/y_old']; h0=np.zeros(3,dtype=complex)
            rec,pix=record(H,R,B,h0,oldy,(.5,30.),oldpv,depth,zs)
            results[name+'/old']=rec;maps[name+'/old']=pix
    tag='natural_old_certificate' if include_old else 'natural_certificate'
    (ROOT/'results'/(tag+'.json')).write_text(json.dumps(results,indent=2)+'\n')
    np.savez_compressed(ROOT/'results'/(tag+'_cells.npz'),**maps)


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--short',action='store_true')
    cli.add_argument('--case',action='append')
    cli.add_argument('--depth',type=int,default=8)
    cli.add_argument('--include-old',action='store_true')
    a=cli.parse_args()
    if a.short:
        short_experiment()
    else:
        natural_experiment(a.depth,a.case,a.include_old)

