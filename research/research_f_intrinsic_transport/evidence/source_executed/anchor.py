"""F1: prior-only initialization; online loader cannot access truth arrays."""
import json, sys, time, traceback
import numpy as np
from scipy.optimize import least_squares
from .common import ROOT,OUT,Book,Prior,BudgetStop,whiten,response_jac,save_endpoint,dump,TAU
from prior_constrained import prior
sys.path.insert(0,str(prior.A30/'src'))
from a30.physics import FrozenMaxwellAdapter

def run():
    protocol=json.loads((OUT/'configs/EXECUTION_OVERLAY.json').read_text())
    gate=json.loads((OUT/'GATE1.json').read_text())
    if not gate['paid_stage_allowed']:raise RuntimeError('F0 gate blocks paid stage')
    book=Book('F1');summaries=[];error=None
    try:
        for case in protocol['cases']:
            meta=json.loads((OUT/f'inputs/online/{case}.json').read_text());data=np.load(OUT/f'inputs/online/{case}.npz')
            if set(data.files)!={'observation','sigma'}:raise RuntimeError('online fields forbidden')
            y=data['observation'];sigma=float(data['sigma']);P=Prior(meta)
            model=FrozenMaxwellAdapter.from_background_npz(prior.INPUTS/'background'/f"{meta['background_id']}.npz",book,purpose_prefix='RESEARCH_F_F1')
            states={};jacobians={};history=[];st=time.perf_counter()
            def full(p):
                key=np.asarray(p).tobytes()
                if key not in states:
                    if len(states)>=12:raise BudgetStop('per-condition twelve F1 states')
                    chi,row=P.decode(p);cert=prior.verify_exact_endpoint(chi,row,background_id=meta['background_id'])
                    if cert['prior_status']!='P0_EXACT':raise RuntimeError('prior failure')
                    i=len(states);book.context=dict(case=case,operation='F1_likelihood',iteration=i)
                    ss=model.full_state(chi,purpose=f'{case}_{i}')
                    path=OUT/f'results/F1/{case}/state_{i:02d}.npz';save_endpoint(ss,p,path,case=case,iteration=i,GT_used=False)
                    states[key]=(ss,np.array(p),path);r=whiten(ss.data-y,sigma)
                    history.append(dict(id=i,p=p.tolist(),residual=float(np.linalg.norm(r)),likelihood_cost=float(r@r/2),backward_residual=ss.residual,prior='P0_EXACT',array=str(path.relative_to(ROOT))))
                    dump(OUT/f'results/F1/{case}/HISTORY.json',history)
                return states[key][0]
            def fun(p):return whiten(full(p).data-y,sigma)
            def jac(p):
                key=np.asarray(p).tobytes();ss=full(p)
                if key not in jacobians:
                    book.context=dict(case=case,operation='F1_full_material_Jacobian')
                    T=model.receiver_transfer(ss);C=response_jac(T,ss);D=P.derivative(p);V=C@D
                    A=np.column_stack([whiten(V[:,i].reshape(32,8),sigma) for i in range(8)])
                    jacobians[key]=(T,C,A)
                    np.savez_compressed(states[key][2].with_name(states[key][2].stem+'_JAC.npz'),transfer=T,C=C,A=A,D=D)
                return jacobians[key][2]
            fit=least_squares(fun,P.p0,jac=jac,bounds=(np.full(8,1e-10),np.full(8,1-1e-10)),
                method='trf',ftol=1e-8,xtol=1e-8,gtol=1e-6,max_nfev=12,x_scale='jac')
            ss=full(fit.x);A=jac(fit.x);g=A.T@fun(fit.x)
            lo=fit.x<=1e-7;hi=fit.x>=1-1e-7;pg=g.copy();pg[lo & (g>0)]=0;pg[hi & (g<0)]=0
            geo=P.frame(fit.x,jacobians[fit.x.tobytes()][1],sigma)
            normalized_stationarity=float(np.linalg.norm(np.linalg.pinv(A,rcond=1e-7)@fun(fit.x)))
            qualified=bool(np.linalg.norm(fun(fit.x))<=TAU and geo['rank']==7 and ss.residual<1e-10 and fit.optimality<=protocol['anchor_projected_gradient_tolerance'])
            final=states[fit.x.tobytes()][2]
            info=dict(case=case,role='INDEPENDENT_PHYSICS_ONLY_CONDITIONAL_RECONSTRUCTION',GT_used=False,N3_used=False,
                known_conditional_geometry=meta['known_conditioned_parameters'],same_parent_two_backgrounds=True,
                status=int(fit.status),message=fit.message,nfev=fit.nfev,njev=fit.njev,states=len(states),
                initial_residual=history[0]['residual'],final_residual=float(np.linalg.norm(fun(fit.x))),tau=TAU,
                optimality=float(fit.optimality),gradient_inf=float(np.max(abs(g))),projected_gradient_inf=float(np.max(abs(pg))),
                projected_GN_coordinate_step=normalized_stationarity,material_rank=geo['rank'],response_per_material_spectrum=geo['sv'],
                material_singular_values=geo['material_sv'],qualified=qualified,final_array=str(final.relative_to(ROOT)),
                final_p=fit.x,wall_seconds=time.perf_counter()-st,posterior_samples=False,prior_status='P0_EXACT',
                observable_information='local full-rank quotient response, stationarity and original likelihood; not GT-recovery certificate')
            summaries.append(info);dump(OUT/'results/F1/SUMMARY.json',summaries);print(json.dumps(dict(case=case,residual=info['final_residual'],qualified=qualified,optimality=info['optimality'],cost=book.totals())),flush=True)
        dump(OUT/'GATE2.json',dict(qualified_anchors=[x['case'] for x in summaries if x['qualified']],passed=any(x['qualified'] for x in summaries),scope='conditional exposed P0, no GT/N3 initialization'))
    except BaseException as e:
        error=repr(e);dump(OUT/'results/F1/FAILURE.json',dict(error=error,traceback=traceback.format_exc(),cost=book.totals()));raise
    finally:book.end(error)

if __name__=='__main__':run()
