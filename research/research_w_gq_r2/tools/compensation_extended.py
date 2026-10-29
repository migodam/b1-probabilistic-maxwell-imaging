"""Reoptimize the old compensation against added channels, then compare null secants."""
from pathlib import Path
from dataclasses import replace
import json,numpy as np
from scipy.optimize import lsq_linear
from r2.accounting import BudgetBook,plain
from r2.physics import Eval,WallGrid,S1_TX
from r2.worlds import white,whiteJ
from compensation import base
ROOT=Path(__file__).resolve().parents[1]
def main():
    out=ROOT/'results/compensation_extended';out.mkdir(exist_ok=False);old=ROOT/'results/compensation';b=BudgetBook();w0,t0=base();target=replace(t0,pose=t0.pose+[.05,0]);d=np.load(old/'BASE.npz');ref=d['data'];n=w0.epsinf.size;scale=np.r_[np.full(n,.5),np.full(n,.05)]
    evs=[Eval(b,'E2',f,S1_TX,w0.x,w0.y,h=.15 if f>1 else .17) for f in [.85,1.,1.15]];curves=[]
    for contract,indices in [('A1',[1]),('A2',[0,1,2])]:
        e=np.load(old/'endpoint_b0_dx0.05.npz');wall=WallGrid(w0.x,w0.y,e['eps'],e['ell']);hist=[]
        for it in range(5):
            ss=[evs[i].evaluate(wall,target,label=f'new-channel corrector {contract} it{it}') for i in indices]
            data=np.asarray([s.data for s in ss]);r=white(data-ref[indices]);norm=float(np.linalg.norm(r));hist.append(dict(iteration=it,d=norm))
            if norm<1e-6:break
            As=[evs[i].wall_jacobian(s) for i,s in zip(indices,ss)];A=whiteJ(np.concatenate(As))*scale[None,:]
            z=np.r_[wall.epsinf.ravel(),wall.ell.ravel()];lo=(np.r_[np.full(n,1.5),np.zeros(n)]-z)/scale;hi=(np.r_[np.full(n,6.),np.full(n,.25)]-z)/scale
            sol=lsq_linear(np.vstack((A,1e-4*np.eye(2*n))),np.r_[-r,np.zeros(2*n)],bounds=(lo,hi),tol=1e-10,max_iter=100)
            ok=False
            for alpha in [1.,.5]:
                zz=z+alpha*scale*sol.x;w=WallGrid(w0.x,w0.y,zz[:n].reshape(w0.shape),zz[n:].reshape(w0.shape))
                xx=np.asarray([evs[i].evaluate(w,target,label=f'{contract} nonlinear safeguard').data.copy() for i in indices]);nr=float(np.linalg.norm(xx-ref[indices])/.01)
                hist[-1].setdefault('trials',[]).append(dict(alpha=alpha,d=nr))
                if nr<norm:wall=w;ok=True;break
            if not ok:break
        nested=np.asarray([ev.evaluate(wall,target,label=f'{contract} reoptimized endpoint all data').data.copy() for ev in evs])
        np.savez_compressed(out/f'{contract}_ENDPOINT.npz',eps=wall.epsinf,ell=wall.ell,data=nested)
        curves.append(dict(contract=contract,dx=.05,history=hist,d_A0=float(np.linalg.norm((nested-ref)[1,:2])/.01),d_A1=float(np.linalg.norm((nested-ref)[1])/.01),d_A2=float(np.linalg.norm(nested-ref)/.01),material_RMS=float(np.sqrt(np.mean((wall.epsinf-w0.epsinf)**2))),status='five-step bounded local optimization, not global profile infimum'))
        (out/'RESULTS.json').write_text(json.dumps(dict(curves=curves),default=plain,indent=2));print(json.dumps(curves[-1]),flush=True)
    # Two distinct first-order compensators; no claim that two samples optimize ker(A).
    A=whiteJ(d['A'][1].reshape(4,16,-1)[:2].reshape(32,-1))*scale[None,:];B=whiteJ(d['B'][1].reshape(4,16,-1)[:2].reshape(32,-1));u,s,vh=np.linalg.svd(A,full_matrices=True)
    rank=int(np.sum(s>np.finfo(float).eps*max(A.shape)*s[0]));a=np.linalg.lstsq(A,-B[:,0],rcond=None)[0];null=vh[rank:];sec=[]
    U=u[:,:rank]
    for k in range(2):
        direction=a+(null[0]*5 if k and len(null) else 0);physical=direction*scale
        curvature=[]
        for h in [.0001,.0002]:
            dat=[]
            for sign in [-1,1]:
                q=sign*h*physical;w=WallGrid(w0.x,w0.y,w0.epsinf+q[:n].reshape(w0.shape),w0.ell+q[n:].reshape(w0.shape));t=replace(t0,pose=t0.pose+[sign*h,0]);dat.append(evs[1].evaluate(w,t,label='second-order compensator secant').data.copy())
            second=white((dat[0]+dat[1]-2*ref[1])[:2])/h**2;res=second-U@(U.T@second)
            curvature.append(dict(h=h,second_norm=float(np.linalg.norm(second)),projected_second_norm=float(np.linalg.norm(res)),scope='float-range projection; no uniform derivative error or ker-optimum claim'))
        sec.append(dict(alternative=k,first_order_residual=float(np.linalg.norm(A@direction+B[:,0])),physical_direction_max=float(np.max(abs(physical))),curvature=curvature))
    (out/'SECOND_ORDER.json').write_text(json.dumps(sec,indent=2))
if __name__=='__main__':main()
