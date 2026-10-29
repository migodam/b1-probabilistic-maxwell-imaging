"""Small continuous off-grid imaging control; no global-optimum or broad-support claim."""
from pathlib import Path
from dataclasses import replace
import argparse,json,time,gc,numpy as np
from scipy.optimize import minimize
from r2.accounting import BudgetBook,plain
from r2.physics import Eval,WallGrid,Target,S1_TX
from r2.worlds import micro_wall,public_wall,target_defs,make_target
ROOT=Path(__file__).resolve().parents[1]
class CallLimit(Exception):pass
def main():
    out=ROOT/'results/continuous_fwi_third';out.mkdir(exist_ok=False);book=BudgetBook()
    truewall=micro_wall(1,1/16,0);truth=Target([.28,2.13],.13,4.2,.10,'ellipse',orientation=0.)
    y=np.load(ROOT/'results/continuous_fwi/OBSERVATION.npz')['y']
    (out/'PLAN.json').write_text(json.dumps(dict(status='DEVELOPMENT_ADAPTATION',reason='original public starts omitted allowed ellipse orientation0; add same neutral ROI-center start to all three methods',no_truth_in_initialization=True,max_actual_evaluations_each=10,extra_FE_hard_cap=30)))
    # Initialization is selected only from the already paid public prior responses.
    library=target_defs()[1];scores=[];bo=ROOT/'results/imaging'
    for g in range(6):
        for j,t in enumerate(library):
            data=np.load(bo/f'macro_g{g}_f1_q{j}.npz')['data'];scores.append((float(np.sum(abs(data-y)**2)),g,j))
    scores.sort();initial=[(scores[0][1],0)]
    records=[]
    for method in ['KNOWN_WALL_ORACLE','JOINT_MACRO_CELL_FWI','EFFECTIVE_MACRO_FWI']:
        for start,(g,j) in enumerate(initial):
            wall= truewall if method=='KNOWN_WALL_ORACLE' else public_wall(g)
            mask=wall.support;nc=int(mask.sum());npar=0 if method=='KNOWN_WALL_ORACLE' else (2 if method=='EFFECTIVE_MACRO_FWI' else 2*nc)
            tar=Target([0.,2.1],.12,3.,.065,'ellipse',orientation=0.);ev=Eval(book,'E1',1.,S1_TX,wall.x,wall.y)
            if npar==2:wv=np.array([float(wall.epsinf[mask].mean()),float(wall.ell[mask].mean())])
            elif npar:wv=np.r_[wall.epsinf[mask],wall.ell[mask]]
            else:wv=np.zeros(0)
            z0=np.r_[wv,tar.pose,tar.contrast,tar.loss,tar.radius]
            bounds=([(1.5,6),(0,.25)] if npar==2 else [(1.5,6)]*nc+[(0,.25)]*nc if npar else [])+[(-1,1),(1.6,2.6),(0,8),(0,.25),(.08,.16)]
            calls=[];best=[float('inf'),None,None,None];last=[None,None]
            def fun(z):
                if last[0] is not None and np.array_equal(z,last[0]):return last[1]
                if len(calls)>=10:raise CallLimit('10 actual full objective evaluations per public shape start')
                e=wall.epsinf.copy();l=wall.ell.copy()
                if npar==2:e[mask]=z[0];l[mask]=z[1]
                elif npar:e[mask]=z[:nc];l[mask]=z[nc:2*nc]
                w=WallGrid(wall.x,wall.y,e,l);v=z[npar:];t=replace(tar,pose=v[:2],contrast=v[2],loss=v[3],radius=v[4])
                state=ev.evaluate(w,t,label=f'offgrid {method} start{start} evaluation{len(calls)}');res=state.data-y
                grad=ev.objective_gradient(state,res/.01**2)
                value=.5*float(np.sum(abs(res)**2))/.01**2
                if npar:
                    # Same full-grid quadratic material penalty evaluated on each class.
                    pe=(e[mask]-3.35)/.7;pl=(l[mask]-.065)/.08
                    value+=.5*float(np.mean(pe**2+pl**2))
                    ge=grad['epsinf'][mask]+pe/(.7*nc);gl=grad['ell'][mask]+pl/(.08*nc)
                    gw=np.array([ge.sum(),gl.sum()]) if npar==2 else np.r_[ge,gl]
                else:gw=np.zeros(0)
                gg=np.r_[gw,grad['target']]
                row=dict(call=len(calls)+1,objective=value,gradient_norm=float(np.linalg.norm(gg)),target=t.metadata(),wall_id_initial=g,selected_by='public finite likelihood only',primal_backward=state.residual)
                calls.append(row);np.savez_compressed(out/f'{method}_s{start}_call{len(calls)}.npz',z=z,data=state.data,eps=e,ell=l,gradient=gg)
                if value<best[0]:best[:]=[value,z.copy(),w,t]
                last[:]=[z.copy(),(value,gg)];return value,gg
            try:
                result=minimize(fun,z0,jac=True,bounds=bounds,method='L-BFGS-B',options=dict(maxiter=100,maxls=8,gtol=1e-6,ftol=1e-12));stop=result.message
            except CallLimit as error:stop='NONCONVERGED_CALL_CAP: '+str(error)
            rec=dict(method=method,start=start,shape=tar.shape,best_objective=best[0],best_target=best[3].metadata(),history=calls,stop=stop,calls=len(calls),same_actual_call_budget=True)
            records.append(rec);(out/'ONLINE_OUTPUTS.json').write_text(json.dumps(records,default=plain,indent=2))
            del ev;gc.collect();print(method,start,stop,flush=True)
    (out/'ONLINE_SEALED.json').write_text(json.dumps(dict(time_unix=time.time(),oracle_wall_only_for_labeled_oracle=True,online_initialization='public finite macro/target likelihood, same shape starts',scope='one exposed-development offgrid target; fixed selected support per start, no claim general unknown support')))
    scores=[]
    for method in sorted(set(r['method'] for r in records)):
        r=min((r for r in records if r['method']==method),key=lambda a:a['best_objective']);q=r['best_target']
        scores.append(dict(method=method,chosen_start_by_objective=r['start'],position_error=float(np.linalg.norm(np.asarray(q['pose'])-truth.pose)),radius_error=abs(q['radius']-truth.radius),eps_real_error=abs(q['contrast']-truth.contrast),loss_error=abs(q['loss']-truth.loss),shape_correct=q['shape']==truth.shape,objective=r['best_objective'],stop=r['stop'],predicted=q))
    (out/'RESULTS.json').write_text(json.dumps(scores,indent=2));print(json.dumps(scores))
if __name__=='__main__':main()
