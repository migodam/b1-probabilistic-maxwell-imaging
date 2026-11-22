"""Outward-rounded verification using the actual fine-grid Fourier multiplier.

Keeps the omitted spatial Fourier coefficients as a lower operator envelope, then
uses exact finite aperture Fourier/receiver/source Grams. The material has 9216
free coordinates in every test. Only the proof matrix is small.
"""
import argparse
import json
import time
from collections import defaultdict
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import cholesky, eigh
from flint import acb, acb_mat, arb, ctx
from threadpoolctl import threadpool_limits
from fourier_bound import STUDY, PUBLIC, DATA, grid, gram, read, h
sys.path.insert(0,str(STUDY/'core'))
from certificate import ball, star, frobenius_squared, gershgorin, safe_float

def hi(x): return safe_float(x,upper=True)
def lo(x): return safe_float(x,upper=False)

def sub(a,rows,cols):
    rows=list(rows);cols=list(cols)
    return acb_mat(len(rows),len(cols),[a[i,j] for i in rows for j in cols])
def stack(a,b):
    return acb_mat(a.nrows()+b.nrows(),a.ncols(),a.entries()+b.entries())

class ExactGeometry:
    def __init__(self,S,B,Q,p,qmax):
        n=int(np.sqrt(S.shape[1]));T=np.column_stack((S.conj().T,B))
        qs=np.arange(-qmax,qmax+1)
        root=[acb(arb(2)*k/p).exp_pi_i() for k in range(p)]
        dirichlet=[sum((root[(k*j)%p] for j in range(n)),acb(0)) for k in range(p)]
        f=acb_mat(len(qs),n,[root[(-int(k)*j)%p] for k in qs for j in range(n)])
        # Column-major amplitudes indexed by (y mode, x mode), with a single 1/p for 2D.
        columns=[]
        for z in T.T:
            out=f*ball(z.reshape(n,n))*f.transpose()
            columns.append([v/p for v in out.entries()])
        up=acb_mat(T.shape[1],len(Q),[v for column in columns for v in column]).transpose()
        g0=star(ball(T))*ball(T)
        dd=(Q[None,:,:]-Q[:,None,:])%p
        unit=acb_mat(len(Q),len(Q),[dirichlet[int(x)]*dirichlet[int(y)]/(p*p) for x,y in dd.reshape(-1,2)])
        self.unit,self.up,self.g0=unit,up,g0
        self.double_gram=gram(S,B,Q,p)

    def projected(self,idx,rec):
        r=len(idx);n=r+rec; ns=self.g0.nrows()-rec
        unit,up,g0=self.unit,self.up,self.g0
        FF=acb_mat(n,n)
        for i,ii in enumerate(idx):
            for j,jj in enumerate(idx): FF[i,j]=unit[ii,jj]
            for j in range(rec):
                FF[i,r+j]=up[ii,j];FF[r+j,i]=up[ii,j].conjugate()
        for i in range(rec):
            for j in range(rec):FF[r+i,r+j]=g0[i,j]
        FB=stack(sub(up,idx,range(rec+ns)),sub(g0,range(rec),range(rec+ns)))
        return FF,FB,g0

def profile_proof(g,spec,Q,rec,ns,u,t,mu,c):
    start=time.perf_counter()
    u0,t0,mu0,c0=map(arb,(u,t,mu,c))
    alpha=mu0*c0
    masked=np.zeros(spec.shape,dtype=bool); masked[Q[:,0],Q[:,1]]=True
    minom=arb(1000000)
    for z in spec[~masked]:
        d=1-u0*arb(float(z.real))+t0*arb(float(z.imag))-c0
        if d<minom:minom=d
        if not d>=0: raise ArithmeticError('Omitted Fourier entry invalidates lower envelope')
    ix=[];coeff=[]
    for i,(x,y) in enumerate(Q):
        z=spec[x,y];v=mu0*(1-u0*arb(float(z.real))+t0*arb(float(z.imag))-c0)
        if v==0:continue
        if not (v>0 or v<0):raise ArithmeticError('unresolved retained coefficient')
        ix.append(i);coeff.append(v)
    ff,fb,gg=g.projected(ix,rec)
    dd=coeff+[arb(1)]*rec;n=len(dd)
    M=ff/alpha
    for i,di in enumerate(dd):M[i,i]+=1/di
    idx=np.r_[ix,len(Q)+np.arange(rec+ns)];gd=g.double_gram[np.ix_(idx,idx)]
    md=gd[:n,:n]/(mu*c);md.flat[::n+1]+=1/np.array([float(v) for v in dd])
    bd=gd[:n,len(ix):];de,V=eigh(h(md),driver='evd',check_finite=False)
    neg=int(np.count_nonzero(de<0));nu=int(n-neg-sum(d>0 for d in dd))
    if not 0<=nu<ns:raise ArithmeticError('No surviving illuminations.')
    gap=float(min(abs(de))*.5)
    negbasis=V[:,:neg];W=negbasis*np.sqrt(-2*de[:neg])
    lifted=md+W@W.conj().T;lifted.flat[::n+1]-=gap
    F=ball(cholesky(h(lifted),lower=True)); Wb=ball(W); Vb=ball(negbasis)
    rr=M+Wb*star(Wb)-F*star(F)
    for i in range(n):rr[i,i]-=arb(gap)
    eps=frobenius_squared(rr).sqrt().upper()
    gaplow=(arb(gap)-eps).lower()
    if not gaplow>0:raise ArithmeticError('Positive shift is unresolved.')
    if neg and not gershgorin(star(Vb)*(M*Vb+Vb*arb(gap)),upper=True)<0:
        raise ArithmeticError('Negative subspace unresolved.')
    X=ball(V@((V.conj().T@bd)/de[:,None]))
    defect=M*X-fb
    loss=frobenius_squared(defect)/(gaplow*alpha**2)
    # Upper bound on U, in Loewner order, sufficient for any actual observation.
    Ub=gg/alpha-(star(fb)*X+star(X)*fb-star(X)*(M*X))/alpha**2
    for i in range(gg.nrows()):Ub[i,i]+=loss.upper()
    return Ub,nu,dict(selected_fourier_coeffs=len(ix),fine_cells=9216,small_matrix_dimension=n,
                     certified_surrogate_nu=nu,positive_gap_of_small_system_lower=lo(gaplow),
                     minimum_omitted_symbol_slack=lo(minom),source_elimination_uncertainty_upper=hi(loss),
                     wall_seconds=time.perf_counter()-start)

def one_y(U,nu,Y,u,t,mu):
    ns=Y.shape[1]; y=ball(Y); factors=acb_mat(ns,ns)
    for i in range(ns):factors[i,i]=acb(arb(mu)*arb(u)/2,arb(mu)*arb(t)/2)
    A=stack(y,factors); Hb=star(y)*y-star(A)*U*A
    Hc=np.array([[complex(z) for z in row] for row in Hb.tolist()]);v,V=eigh(h(Hc))
    V=V[:,::-1];lower_list=[]
    for k in range(nu+1,ns+1):
        W=ball(V[:,:k]);wx=gershgorin(star(W)*Hb*W,upper=False)
        gr=gershgorin(star(W)*W,upper=True)
        if wx<=0:break
        lower_list.append((wx/gr).lower())
    total=sum(lower_list,arb(0));return (max(0.,lo(total.sqrt())) if total>0 else 0.,
                                             [lo(x) for x in lower_list])

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry',choices=('both','180','360'),default='both')
    arg=parser.parse_args()
    ctx.prec,ctx.threads=96,2
    meta={d['id']:d for d in json.loads((PUBLIC/'manifest.json').read_text())['observations']}
    scan=json.loads((DATA/'fine_scan.json').read_text())['rows']; todo={}
    for r in scan:
        u=r['u']
        if u==.25 or r['d']<=r['radius'] or (arg.geometry!='both' and arg.geometry not in r['id']):continue
        # Verify the best crossed cap per type/source setting, including the actual real-only controls.
        k=(r['id'],u,'real' if r['t']==0 else 'optical')
        if k not in todo or r['d']>todo[k]['d']: todo[k]=r
    profiles=defaultdict(list)
    for r in todo.values():
        geom=meta[r['id']]['geometry'];key=(geom,r['grid'],r['u'],r['t'],r['mu'],r['c'])
        profiles[key].append(r)
    geometry_key=None; report={'profiles':[],'cases':[],'failures':[]};start=time.perf_counter()
    for key,rows in sorted(profiles.items()):
        name,q,u,t,mu,c=key
        if geometry_key!=(name,q):
            spec,S,B,_=read(rows[0]['id'],name);Q=grid(q,len(spec))
            g=ExactGeometry(S,B,Q,len(spec),q);geometry_key=(name,q)
        try:
            U,nu,proof=profile_proof(g,spec,Q,len(S),B.shape[1],u,t,mu,c)
            report['profiles'].append(dict(key=list(key),**proof))
            for r in rows:
                with np.load(PUBLIC/meta[r['id']]['file']) as f:y=f['Y']
                bound,terms=one_y(U,nu,y,u,t,mu)
                report['cases'].append(dict(id=r['id'],profile=list(key),cap=u,nu=nu,
                    distance_lower=bound,descending_positive_tail_lower=terms,noise_radius=r['radius'],
                    margin_lower=lo(arb(bound)-arb(r['radius'])),status='excluded' if bound>r['radius'] else 'unresolved'))
        except (ArithmeticError,ZeroDivisionError,ValueError) as e:
            report['failures'].append(dict(key=list(key),reason=str(e),ids=[r['id'] for r in rows]))
        report['wall_seconds']=time.perf_counter()-start
        suffix='' if arg.geometry=='both' else '_'+arg.geometry
        (DATA/('fine_verified'+suffix+'.json')).write_text(json.dumps(report,indent=2)+'\n')
        print(key,'passed',sum(x['status']=='excluded' for x in report['cases']),'failures',len(report['failures']),flush=True)
if __name__=='__main__':
    with threadpool_limits(1):main()
