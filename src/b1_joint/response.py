"""Explicit tensor-Chebyshev response approximation of paid full-Mie arrays.

Its density target is separate from the numerical full-Mie target. Nothing in
this module imports a physical solver. Cartesian contractions avoid forming
quadrature-node-by-polynomial-feature design matrices.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
from numpy.polynomial.chebyshev import chebvander
from .material import RMAX,geometry

PHYSICAL_BOUNDS=np.array([[0.,RMAX[0]],[0.,RMAX[1]],[2.2,2.8],[.04,.16]])

def nodes(n=5):
    t=(1-np.cos(np.pi*np.arange(n)/(n-1)))/2
    return [lo+t*(hi-lo) for lo,hi in PHYSICAL_BOUNDS]

def fit(axes,values):
    out=np.asarray(values).copy()
    for j,a in enumerate(axes):
        lo,hi=PHYSICAL_BOUNDS[j];v=chebvander(2*(a-lo)/(hi-lo)-1,len(a)-1)
        out=np.moveaxis(np.tensordot(np.linalg.inv(v),out,axes=(1,j)),0,j)
    return out

class Response:
    def __init__(self,path):
        data=np.load(path);self.coeff={float(k):data[f'coeff_{k:g}'] for k in data['frequencies']}
        self.bounds=PHYSICAL_BOUNDS;self.evaluations=0
    def mesh(self,axes,k=1.6):
        out=self.coeff[float(k)]
        for j,a in enumerate(axes):
            lo,hi=self.bounds[j];v=chebvander(2*(np.asarray(a)-lo)/(hi-lo)-1,out.shape[j]-1)
            out=np.moveaxis(np.tensordot(v,out,axes=(1,j)),0,j)
        self.evaluations+=int(np.prod([len(a) for a in axes]))
        return out.reshape(-1,out.shape[-1])
    def predict(self,x,frequencies=(1.6,)):
        x=np.atleast_2d(x);p=np.column_stack((geometry(x),x[:,2:]));out=[]
        for k in frequencies:
            c=self.coeff[float(k)];n=c.shape[0];flat=c.reshape(n*n,n*n,-1);rows=[]
            for start in range(0,len(p),4096):
                z=p[start:start+4096];v=[chebvander(2*(z[:,j]-a)/(b-a)-1,n-1) for j,(a,b) in enumerate(self.bounds)]
                g=(v[0][:,:,None]*v[1][:,None,:]).reshape(len(z),n*n)
                e=(v[2][:,:,None]*v[3][:,None,:]).reshape(len(z),n*n)
                q=(g@flat.reshape(n*n,-1)).reshape(len(z),n*n,-1)
                rows.append(np.einsum('ni,nio->no',e,q))
            out.append(np.concatenate(rows) if rows else np.empty((0,c.shape[-1])))
        self.evaluations+=len(x)*len(frequencies)
        return np.concatenate(out,axis=1)

def save_bank(path,axes,fields,frequencies=(1.6,2.4)):
    coeff={f'coeff_{k:g}':fit(axes,fields[i]) for i,k in enumerate(frequencies)}
    # Embedded three-node independent response order, not an integration order.
    ids=[0,2,4];a3=[a[ids] for a in axes]
    for i,k in enumerate(frequencies):
        coarse=fields[i][np.ix_(ids,ids,ids,ids)]
        coeff[f'coarse_coeff_{k:g}']=fit(a3,coarse)
    np.savez_compressed(path,frequencies=np.array(frequencies),**{f'axis_{i}':a for i,a in enumerate(axes)},**{f'fields_{k:g}':f for k,f in zip(frequencies,fields)},**coeff)
