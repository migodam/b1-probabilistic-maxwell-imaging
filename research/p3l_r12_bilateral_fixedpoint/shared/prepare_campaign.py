"""Read-only historical observation/anchor adapter; three paid RML branch entries.

All strong anchors use original P1 and exact Maxwell, randomized observations
and prior targets. Endpoint frequencies are never called posterior masses.
"""
from .physics import *
from .budget import BudgetStop
from scipy.optimize import least_squares
import argparse

def run(parent):
    out=ROOT/'results/P1'/parent/'COMMON'
    if (out/'SEALED.json').exists():return
    out.mkdir(parents=True,exist_ok=True)
    online=np.load(ROOT/'data/online'/f'{parent}.npz');y=online['observation']
    pr=Problem('P1');ph=Physics(pr,y,'P1/'+parent+'/COMMON',out/'states',local_rhs_cap=1200)
    historical=R9/'results/possibility_volume/parents'/parent/'anchor/ANCHOR.npz'
    x0=np.load(historical)['x'];anchor=ph.evaluate(x0,True);states=[anchor];metadata=[]
    rng=np.random.default_rng(np.random.SeedSequence([20261030,120012,101,int(parent[-3:])]))
    for j in range(3):
        eta=rng.normal(size=512);xi=rng.normal(size=66)
        last=[None]
        def fun(x):
            last[0]=ph.evaluate(x,True)
            return np.r_[last[0]['residual']-eta,x-xi]
        def jac(x):
            s=ph.evaluate(x,True)
            return np.vstack((s['jac'],np.eye(66)))
        before=ph.rhs
        try:
            sol=least_squares(fun,x0,jac=jac,max_nfev=8,ftol=1e-7,xtol=1e-7,gtol=1e-7,method='trf',tr_solver='exact')
            s=ph.evaluate(sol.x,True);states.append(s)
            metadata.append(dict(branch=j,status=int(sol.status),nfev=sol.nfev,Phi=s['phi'],RHS=ph.rhs-before,prior=float(.5*sol.x@sol.x),material_RMS_from_original=material_rms(s['material'],anchor['material'])))
        except BudgetStop as e:
            metadata.append(dict(branch=j,status='BUDGET_STOP',error=str(e),RHS=ph.rhs-before));break
    np.savez_compressed(out/'ANCHORS.npz',xs=np.stack([s['x'] for s in states]),materials=np.stack([s['material'] for s in states]),Fs=np.stack([s['F'] for s in states]),phi=np.array([s['phi'] for s in states]),jac=np.stack([s['jac'] for s in states]),dmaterial=np.stack([s['dmaterial'] for s in states]),L=np.stack([s['L'] for s in states]),R=np.stack([s['R'] for s in states]))
    atomic_json(out/'SEALED.json',dict(parent=parent,scope='exposedR9DEV; multipleRMLbranchentries; noGT',new_RHS=ph.rhs,new_LU=ph.lu,branch_search=metadata,historical_original_anchor_cost='read-onlyR9; standalonesharedinitialcostdisclosedseparately',randomized_target_is_posterior=False,utc=time.time()))
    print(parent,'ANCHORS',len(states),'Phi',[round(s['phi'],2) for s in states],'RHS',ph.rhs,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--parents',default='0,1,2,3');a=ap.parse_args()
    for i in map(int,a.parents.split(',')):run(f'R9V_FRESH_{i:03d}')
