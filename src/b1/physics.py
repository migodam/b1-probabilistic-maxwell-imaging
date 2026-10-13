"""3D vector multilayer Mie by tangential-field boundary matching.

Nonmagnetic isotropic layers. Global matching is the main implementation;
scalar admittance propagation is a separately assembled audit implementation.
All material-frequency public calls must reserve the cumulative ledger first.
"""
from __future__ import annotations
import json, importlib.util
from pathlib import Path
import numpy as np
from scipy.special import spherical_jn, spherical_yn
from .budget import Ledger

BASE=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('a13_mie_reused',BASE/'vendor/a13_mie_functions.py')
a13=importlib.util.module_from_spec(spec); spec.loader.exec_module(a13)

def _psi(n,z):
    j=spherical_jn(n,z); y=spherical_yn(n,z)
    return np.array([z*j,z*y]),np.array([j+z*spherical_jn(n,z,True),y+z*spherical_yn(n,z,True)])

def topology_layers(H,eps,config):
    voids=config['topologies'][int(H)]['void_intervals']; bounds=sorted({b for interval in voids for b in interval if b>0}|{config['outer_radius']})
    radii=[]; perm=[]; lo=0.
    for hi in bounds:
        mid=(lo+hi)/2; vacuum=any(a<=mid<=b for a,b in voids)
        radii.append(hi); perm.append(1.+0j if vacuum else complex(eps)); lo=hi
    return np.asarray(radii,float),np.asarray(perm,complex)

def _validate(radii,eps,k,nmax):
    if len(radii)!=len(eps) or not len(radii): raise ValueError('layer dimensions')
    if not np.all(np.diff(np.r_[0.,radii])>0): raise ValueError('radii increasing')
    if k<=0 or nmax<1 or np.any(np.real(eps)<=0) or np.any(np.imag(eps)<0): raise ValueError('passive positive dielectric only')

def _matching(radii,eps,k,nmax):
    _validate(radii,eps,k,nmax); m=np.sqrt(eps); L=len(radii)
    a=np.zeros(nmax,complex); b=a.copy(); maxres=0.
    for n in range(1,nmax+1):
        for polarization,out in [('TM',a),('TE',b)]:
            w=1/m if polarization=='TM' else m
            A=np.zeros((2*L,2*L),complex); rhs=np.zeros(2*L,complex)
            for j,r in enumerate(radii):
                vals,ders=_psi(n,m[j]*k*r)
                columns=[0] if j==0 else [2*j-1,2*j]
                for q,col in enumerate(columns):
                    A[2*j,col]=vals[q]; A[2*j+1,col]=w[j]*ders[q]
                if j+1<L:
                    rv,rd=_psi(n,m[j+1]*k*r)
                    for q,col in enumerate([2*(j+1)-1,2*(j+1)]):
                        A[2*j,col]=-rv[q]; A[2*j+1,col]=-w[j+1]*rd[q]
                else:
                    v,d=_psi(n,k*r); xi=v[0]+1j*v[1]; dxi=d[0]+1j*d[1]
                    A[2*j,-1]=xi; A[2*j+1,-1]=dxi
                    rhs[2*j]=v[0]; rhs[2*j+1]=d[0]
            # Exact invertible equilibration. No regularization or jitter.
            col=np.maximum(np.max(np.abs(A),axis=0),1e-300)
            Ae=A/col; row=np.maximum(np.max(np.abs(Ae),axis=1),1e-300)
            y=np.linalg.solve(Ae/row[:,None],rhs/row); sol=y/col
            backward=np.linalg.norm(A@sol-rhs)/(np.linalg.norm(A)*np.linalg.norm(sol)+np.linalg.norm(rhs)+1e-300)
            maxres=max(maxres,float(backward)); out[n-1]=sol[-1]
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)): raise FloatingPointError('Mie coefficients nonfinite')
    return a,b,dict(boundary_backward_residual=maxres,dense_boundary_solves=2*nmax,layer_count=L)

def _recursive(radii,eps,k,nmax):
    _validate(radii,eps,k,nmax); m=np.sqrt(eps); arrays=[]
    for polarization in ['TM','TE']:
        w=1/m if polarization=='TM' else m; coeffs=[]
        for n in range(1,nmax+1):
            c=np.array([1.+0j,0.+0j])
            for j in range(len(radii)-1):
                v,d=_psi(n,m[j]*k*radii[j]); val=c@v; der=w[j]*(c@d)
                v2,d2=_psi(n,m[j+1]*k*radii[j])
                matrix=np.vstack((v2,w[j+1]*d2)); scale=np.maximum(np.max(abs(matrix),axis=0),1e-300)
                c=np.linalg.solve(matrix/scale,np.array([val,der]))/scale
                c/=np.max(abs(c))
            v,d=_psi(n,m[-1]*k*radii[-1]); Y=w[-1]*(c@d)/(c@v)
            ext,dp=_psi(n,k*radii[-1]); xi=ext[0]+1j*ext[1]; dx=dp[0]+1j*dp[1]
            coeffs.append((dp[0]-Y*ext[0])/(dx-Y*xi))
        arrays.append(np.asarray(coeffs))
    return *arrays,dict(implementation='independent_admittance_recursion')

def angular_data(a,b,angles,k):
    # Known range phase removed: vector far-field AMPLITUDE A = i(S1,S2)/k.
    # Two orthogonal unit incident polarizations / scattering-plane projections.
    z=np.asarray([a13.mie_s1_s2(np.cos(np.deg2rad(th)),a,b) for th in angles],complex).ravel()*1j/k
    return np.r_[z.real,z.imag]

class Forward:
    def __init__(self,config,ledger,cache_dir):
        self.config=config; self.ledger=ledger; self.path=Path(cache_dir); self.path.mkdir(parents=True,exist_ok=True)
        self.index=self.path/'index.json'; self.entries=json.loads(self.index.read_text()) if self.index.exists() else {}
    def run(self,H,er,ei,*,kind='3D_VECTOR_MIE',nmax=None,implementation='matching',cache=True,radii=None,eps_layers=None,k=None):
        c=self.config; n=c['nmax'] if nmax is None else nmax; kw=c['wavenumber'] if k is None else k
        rr,ee=topology_layers(H,complex(er,ei),c) if radii is None else (np.asarray(radii,float),np.asarray(eps_layers,complex))
        params=dict(H=int(H),er=float(er),ei=float(ei),radii=rr.tolist(),eps_layers=[[z.real,z.imag] for z in ee],k=kw,nmax=n,implementation=implementation,angles=c['angles_degrees'])
        key=json.dumps(params,sort_keys=True,separators=(',',':'))
        if cache and key in self.entries:
            self.ledger.cache_hit(key); q=np.load(self.path/self.entries[key]); return q['data'],q['a'],q['b'],json.loads(str(q['meta']))
        def calc():
            if implementation=='matching': a,b,info=_matching(rr,ee,kw,n)
            elif implementation=='recursive': a,b,info=_recursive(rr,ee,kw,n)
            elif implementation=='a13_homogeneous':
                if len(rr)!=1: raise ValueError('homogeneous oracle needs one layer')
                a,b=a13.mie_coefficients(n,np.sqrt(ee[0]),kw*rr[-1]); info=dict(implementation='unmodified_A13_homogeneous')
            else: raise ValueError('unknown solver')
            f=angular_data(a,b,c['angles_degrees'],kw)
            qe,qs,qa=a13.mie_efficiencies(a,b,kw*rr[-1]); info.update(Qext=qe,Qsca=qs,Qabs=qa,physical_dimension='3D_VECTOR')
            return f,a,b,info
        out=self.ledger.evaluate(kind,params,calc)
        name=f'forward_{len(self.entries):05d}.npz'
        # Noncached audits also retained but never overwrite an existing cache.
        if not cache: name=f'audit_{self.ledger.count():05d}.npz'
        np.savez_compressed(self.path/name,data=out[0],a=out[1],b=out[2],meta=json.dumps(out[3]))
        if cache:
            self.entries[key]=name; self.index.write_text(json.dumps(self.entries,indent=2)+'\n')
        return out
