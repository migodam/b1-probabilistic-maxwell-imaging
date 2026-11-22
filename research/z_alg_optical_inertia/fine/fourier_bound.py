"""Direct finite-grid material test. Fourier truncation gives a lower matrix bound;
all 96x96 independent material coefficients remain free."""
from pathlib import Path
import json, sys, time
import numpy as np
from scipy.linalg import eigh, solve
from threadpoolctl import threadpool_limits
from scipy.fft import fft2
STUDY=Path(__file__).resolve().parents[1]
ROOT=STUDY.parents[1]
PUBLIC=STUDY/'synthetic/public'
DATA=STUDY/'results/theory'

def h(a): return (a+a.conj().T)/2

def read(case,geometry):
    with np.load(PUBLIC/'fine_geometry.npz') as f:
        spec=f['spectrum'];S=f[geometry+'__S'];B=f[geometry+'__B']
    with np.load(PUBLIC/(case+'.npz')) as f:Y=f['Y']
    return spec,S,B,Y

def export_geometry():
    sys.path.insert(0,str(ROOT/'research/srap_cie_array_study_20261121'))
    from data2d import SyntheticOperators
    names=['PW180_10x20','PW360_10x20']; ops=SyntheticOperators(96,names)
    fields={'spectrum':ops.green.spectrum}
    for name in names:
        fields[name+'__S']=ops.fine_s[ops.rx_ids[name]]
        fields[name+'__B']=ops.fine_e[:,ops.tx_ids[name]]
    np.savez_compressed(PUBLIC/'fine_geometry.npz',**fields)

def grid(qmax,pad):
    a,b=np.meshgrid(np.arange(-qmax,qmax+1),np.arange(-qmax,qmax+1),indexing='ij')
    return np.stack((a.ravel()%pad,b.ravel()%pad),axis=1)

def gram(S,B,Q,size):
    n=int(np.sqrt(S.shape[1]));p=size
    T=np.column_stack((S.conj().T,B))
    z=fft2(T.reshape(n,n,-1),s=(p,p),axes=(0,1))/p
    U_T=z[Q[:,0],Q[:,1]]
    g=np.fft.ifft(np.r_[np.ones(n),np.zeros(p-n)])*p
    delta=(Q[None,:,:]-Q[:,None,:])%p
    UU=(g[delta[:,:,0]]*g[delta[:,:,1]])/p**2
    return h(np.block([[UU,U_T],[U_T.conj().T,T.conj().T@T]]))

def profile(spec,u,t,Q):
    z=1-u*spec.real+t*spec.imag
    outside=np.ones(z.shape,dtype=bool);outside[Q[:,0],Q[:,1]]=False
    # Downward dyadic rounding provides an independently checkable Loewner envelope.
    c=np.floor(32*z[outside].min())/32.0
    if c<=0: raise ValueError('Retained Fourier square must be enlarged or cap reduced.')
    return c,z[Q[:,0],Q[:,1]]-c

def solve_small(gram,u,t,Y,nrec,mu,c,coeff):
    ns=Y.shape[1];n=len(coeff)+nrec
    gramF=gram[:n,:n]
    Lc=np.zeros((n+ns,ns),complex)
    Lc[len(coeff):n]=Y
    Lc[n:]=np.eye(ns)*mu*(u+1j*t)/2
    FL=gram[:n,:]@Lc
    LL=h(Lc.conj().T@gram@Lc)
    d=np.r_[mu*coeff,np.ones(nrec)]
    M=gramF/(mu*c);M.flat[::n+1]+=1./d
    de,V=eigh(M,driver='evd',check_finite=False)
    nu=int((de>0).sum()-(d>0).sum())
    sol=V@((V.conj().T@FL)/de[:,None])
    H=h(Y.conj().T@Y-LL/(mu*c)+FL.conj().T@sol/(mu*c)**2)
    hd=eigh(H,eigvals_only=True)[::-1]
    return dict(nu=nu,h=hd.tolist(),distance=float(np.sqrt(np.maximum(0.,hd[max(0,nu):]).sum())),
                cap=u,t=t,mu=mu,c=c,fourier_terms=len(coeff))

def main():
    start=time.time()
    if not (PUBLIC/'fine_geometry.npz').exists():export_geometry()
    meta=json.loads((PUBLIC/'manifest.json').read_text())['observations'];rows=[]
    grids=(7,12) # nested lower models; all hypotheses retained
    for geom in ('PW180_10x20','PW360_10x20'):
        views=[x for x in meta if x['geometry']==geom]
        spec,S,B,_=read(views[0]['id'],geom)
        yy={v['id']:np.load(PUBLIC/v['file'])['Y'] for v in views}
        for qmax in grids:
            q=grid(qmax,len(spec));G=gram(S,B,q,len(spec));nr=len(S);ns=B.shape[1]
            for u in (.25,.5,.75):
                for t in ((0.,.3,.4,.55,.75,1.) if u==.5 else (0.,.3,1.)):
                    try:c,co=profile(spec,u,t,q)
                    except ValueError:continue
                    # Factor once for every object and every noise level.
                    idx=np.flatnonzero(abs(co)>1e-12)
                    valid=np.r_[idx,len(co)+np.arange(nr+ns)]
                    GG=G[np.ix_(valid,valid)]
                    co=co[idx]; r=len(co);n=r+nr
                    for mu in (.00001,.0001,.00035,.001,.0025,.005):
                        d=np.r_[mu*co,np.ones(nr)]
                        M=GG[:n,:n]/(mu*c);M.flat[::n+1]+=1/d
                        de,V=eigh(M,driver='evd',check_finite=False)
                        nu=int((de>0).sum()-(d>0).sum())
                        if nu>=ns:continue
                        FL=GG[:n,r:];P=V.conj().T@FL
                        U=h(GG[r:,r:]/(mu*c)-(P.conj().T/de)@P/(mu*c)**2)
                        for v in views:
                            y=yy[v['id']]
                            compact=np.vstack((y,np.eye(ns)*mu*(u+1j*t)/2))
                            H=h(y.conj().T@y-compact.conj().T@U@compact)
                            hd=eigh(H,eigvals_only=True)[::-1]
                            rows.append(dict(id=v['id'],radius=v['radius_noise_fro'],grid=qmax,
                                             u=u,t=t,mu=mu,c=c,nu=nu,h=hd.tolist(),
                                             d=float(np.sqrt(np.maximum(hd[nu:],0).sum()))))
            print(geom,qmax,'pairs above noise',len({(r['id'],r['u']) for r in rows if r['d']>r['radius']}),flush=True)
            (DATA/'fine_scan.json').write_text(json.dumps(dict(rows=rows,wall_seconds=time.time()-start),separators=(',',':'))+'\n')
if __name__=='__main__':
    with threadpool_limits(1):main()
