"""Nondispersive two-band Gram lift: verify a rank obstruction and test partial coherence."""
from pathlib import Path
import argparse,csv,json,time,os
for k in ('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.linalg import eigh,lu_factor,lu_solve,block_diag
from threadpoolctl import threadpool_limits
from verified_wave import P,VerifiedWave,up,fro,rd,opnorm
DATA=P/'ports/public';OUT=P/'results'
def opfile(side,f,q):
    root=DATA/'operators'
    path=root/f'n{side}_f{f}_a{q:02}_c048.npz'
    return path if path.exists() else root/f'n{side}_f{f}_a{q:02}.npz'
def run_obstruction():
    rows=[]
    for side in [32,96]:
        # Geometry-only, deterministic widely spread cells; no material or measurement is read.
        idx=np.random.default_rng(271135).choice(side**2,80,replace=False)
        for q in [6,8,12,16]:
            tic=time.perf_counter(); cpu=time.process_time()
            mats=[];errors=[]
            for f in [400,460]:
                with np.load(opfile(side,f,q)) as d:
                    g=VerifiedWave(d)
                a=g.matrix(idx)
                errors.append(up(len(idx)*g.lagerr+rd(g.R.shape[0])*fro(g.L[idx])*fro(g.R[:,idx])+rd(2)*fro(a)))
                mats.append(a)
            D=mats[0]-mats[1];dE=up(sum(errors)+rd(1)*(fro(mats[0])+fro(mats[1])))
            W=lu_solve(lu_factor(D),np.eye(len(idx)),check_finite=False)
            wn=opnorm(W)
            defect=up(fro(np.eye(len(idx))-D@W)+rd(len(idx))*(np.sqrt(len(idx))+fro(D)*fro(W))+dE*wn)
            row=dict(side=side,q=q,verified_rank_lower=len(idx) if defect<1 else 0,inverse_defect_upper=defect,
                     singular_value_lower=float(max(0.,(1-defect)/wn)*(1-1e-12)),negative_inertia_lower=max(0,len(idx)-2*q) if defect<1 else 0,
                     independent_source_columns=2*q,tail_forced_zero=bool(defect<1 and len(idx)-2*q>=2*q),
                     seconds_wall=time.perf_counter()-tic,seconds_cpu=time.process_time()-cpu)
            rows.append(row);print('rank',row,flush=True)
    (OUT/'frequency_coherence_rank.json').write_text(json.dumps(dict(cells=80,full_field_independent_4p5=True,
       freqs_MHz=[400,460],uniform_common_sector_weight=True,rows=rows),indent=2))
def h(A):return (A+A.conj().T)/2
def run_scan(q):
    side=32;n=side**2; disks=[]
    for f in [400,460]:
        p=opfile(side,f,q)
        with np.load(p) as a:
            g=VerifiedWave(a);G=g.matrix(np.arange(n));S=g.S;B=g.B/np.sqrt(2)
        disks.append((G,S,B))
    # Reflectors/loads may vary; plain row count of 2q across bands is recorded.
    I=np.eye(n);Q=q*2;G1,S1,B1=disks[0];G2,S2,B2=disks[1]
    st=[S1.conj().T@S1,S2.conj().T@S2]
    ya=[];radii=[]
    for j in range(1,9):
        z=np.load(DATA/'cases'/f'n32_c{j:02}_mimo.npz')
        keys=[f'f{f}_a{q:02}_c048_e000' for f in (400,460)]
        ya.append([z[k+'_Y']/np.sqrt(2) for k in keys]);radii.append(np.linalg.norm([float(z[k+'_radius']) for k in keys])/np.sqrt(2))
    rows=[]
    for cap in [.25,.5,1.]:
        # The inputs and eigendecompositions are shared across all eight observations.
        for gamma in [0.,.60,.92]:
            mu=.02
            A1=st[0]+mu*(I-cap*h(G1));A2=st[1]+mu*(I-cap*h(G2))
            C=mu*gamma*(I-cap*(G2+G1.conj().T)/2)
            K=np.block([[A1,C],[C.conj().T,A2]])
            t=time.perf_counter();cpu=time.process_time()
            vals,Qv=eigh(K,driver='evd',check_finite=False);nu=int(np.sum(vals<0))
            split=[eigh(a,eigvals_only=True,check_finite=False) for a in (A1,A2)] if not gamma else []
            for j,(Y1,Y2) in enumerate(ya):
                H0=block_diag(Y1.conj().T@Y1,Y2.conj().T@Y2)
                L=-np.block([[S1.conj().T@Y1+mu*cap/2*B1,mu*cap/2*gamma*B2],
                            [mu*cap/2*gamma*B1,S2.conj().T@Y2+mu*cap/2*B2]])
                T=Qv.conj().T@L
                H=h(H0-T.conj().T@(T/vals[:,None]));spec=eigh(H,eigvals_only=True)[::-1]
                val=np.sqrt(np.maximum(spec[nu:],0).sum())/radii[j]
                separate=None
                if not gamma:
                    tot=0.
                    for k in range(2):
                        block=H[k*q:(k+1)*q,k*q:(k+1)*q];eig=eigh(block,eigvals_only=True)[::-1]
                        tot+=np.maximum(eig[int(np.sum(split[k]<0)):],0).sum()
                    separate=float(np.sqrt(tot)/radii[j])
                rows.append(dict(case=f'n32_c{j+1:02}',q=q,cap=cap,coherence=gamma,mu=mu,nu=nu,s=Q,
                             ratio_numeric=float(val),separate_ratio_numeric=separate,minimum_abs_state_eigenvalue=float(np.min(abs(vals))),
                             geometry_wall_s=time.perf_counter()-t,geometry_cpu_s=time.process_time()-cpu))
            print('q',q,'u',cap,'gamma',gamma,'nu',nu,'first ratio',rows[-8]['ratio_numeric'],flush=True)
    with (OUT/f'frequency_coherence_scan_{q}.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--rank',action='store_true');ap.add_argument('--scan',nargs='*',type=int,default=[])
    opt=ap.parse_args()
    with threadpool_limits(1):
        if opt.rank:run_obstruction()
        for q in opt.scan:run_scan(q)
