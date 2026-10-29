"""Bounded continuation: two additional A2 correctors, then independent endpoint."""
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
    b=BudgetBook();out=ROOT/'results/compensation_final';out.mkdir(exist_ok=False)
    (out/'PLAN.json').write_text(json.dumps(dict(reason='five-step A2 local fit still decreasing at1.224; diagnose optimization ceiling with exactly two additional steps, no threshold change',max_outer=2,max_trial_per_outer=2,max_FE_extra=18,max_RHS_extra=180,remaining_phase_borrowing=False)))
    w0,t0=base();t=replace(t0,pose=t0.pose+[.05,0]);n=w0.epsinf.size;scale=np.r_[np.full(n,.5),np.full(n,.05)];old=ROOT/'results/compensation_extended';d=np.load(old/'A2_ENDPOINT.npz');w=WallGrid(w0.x,w0.y,d['eps'],d['ell']);ref=np.load(ROOT/'results/compensation/BASE.npz')['data'];fs=[.85,1.,1.15]
    evs=[Eval(b,'E2',f,S1_TX,w.x,w.y,h=.15 if f>1 else .17) for f in fs];history=[]
    for it in range(2):
        ss=[e.evaluate(w,t,label='last bounded A2 state') for e in evs];r=white(np.asarray([s.data for s in ss])-ref);norm=float(np.linalg.norm(r));A=whiteJ(np.concatenate([e.wall_jacobian(s) for e,s in zip(evs,ss)]))*scale[None,:]
        cur=np.r_[w.epsinf.ravel(),w.ell.ravel()];lo=(np.r_[np.full(n,1.5),np.zeros(n)]-cur)/scale;hi=(np.r_[np.full(n,6.),np.full(n,.25)]-cur)/scale
        sol=lsq_linear(np.vstack((A,1e-4*np.eye(2*n))),np.r_[-r,np.zeros(2*n)],bounds=(lo,hi),tol=1e-10,max_iter=100)
        row=dict(iteration=it,start=norm,trials=[],linear_optimality=float(sol.optimality));accepted=False
        for alpha in [1.,.5]:
            z=cur+alpha*scale*sol.x;ww=WallGrid(w.x,w.y,z[:n].reshape(w.shape),z[n:].reshape(w.shape));new=np.asarray([e.evaluate(ww,t,label='last bounded A2 trial').data.copy() for e in evs]);nd=float(np.linalg.norm(new-ref)/.01);row['trials'].append(dict(alpha=alpha,d=nd))
            if nd<norm:w=ww;accepted=True;break
        history.append(row)
        if not accepted:break
    result=np.asarray([e.evaluate(w,t,label='sealed final A2 endpoint').data.copy() for e in evs]);np.savez_compressed(out/'ENDPOINT.npz',eps=w.epsinf,ell=w.ell,data=result);rows=[]
    for i,f in enumerate(fs):
        e=Eval(b,'E2',f,S1_TX,w.x,w.y,wall_refine=2,h=.145,targetstep=.1,max_dofs=120000);s=e.evaluate(w,t,label='final A2 independent endpoint')
        baseref=np.load(old/f'INDEPENDENT_A2_f{f}.npz')['base'];rows.append(dict(f=f,error=float(np.linalg.norm(s.data-result[i])/.01),distance=float(np.linalg.norm(s.data-baseref)/.01)))
        np.savez_compressed(out/f'INDEPENDENT_f{f}.npz',base=baseref,endpoint=s.data)
    record=dict(history=history,coarse_distance=float(np.linalg.norm(result-ref)/.01),refined_distance=float(np.sqrt(sum(r['distance']**2 for r in rows))),rows=rows,dx=.05,material_RMS=float(np.sqrt(np.mean((w.epsinf-w0.epsinf)**2))),bounds=dict(eps=[float(w.epsinf.min()),float(w.epsinf.max())],loss=[float(w.ell.min()),float(w.ell.max())]),status='bounded seven-step development profile; no global minimum claim')
    (out/'RESULTS.json').write_text(json.dumps(record,default=plain,indent=2));print(json.dumps(record,default=plain))
if __name__=='__main__':main()
