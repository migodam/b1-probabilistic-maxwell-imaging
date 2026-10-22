"""Forward-of-estimator diagnostic, with no truth or posterior claim.

A material risk action can lie outside the single-decoder material image.
Evaluating its Maxwell response does not make it a posterior draw or a PPC.
"""
from pathlib import Path
import json, sys, time
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT,Book,Target,write,save
from infer_p3r2_locked import frozen_contract


def main():
    manifest,parents,_,_=frozen_contract(OUT,OUT/'METHODS_FROZEN_FOR_LOCKED.json')
    for folder in ('estimates','trajectory_estimates'):
        if json.loads((OUT/'L'/folder/'COMPLETE.json').read_text()).get('status')!='COMPLETE':
            raise RuntimeError('All online risk actions must be complete before forward diagnostics')
    methods=manifest['forward_decision_methods']
    if methods!=['fused','residual','projected_direct','classical_postprocess','dynamic_curvature__arc_length_trapezoid']:
        raise ValueError('unregistered forward decision diagnostic')
    b=Book('L','locked_decision_forward',max_wall=3600);status='COMPLETE';error=None
    rows=[];t=Target(b)
    try:
        write(b.jobdir/'PROTOCOL.json',dict(methods=methods,parents=parents,
            target_id=t.cfg['target_id'],GT_read=False,likelihood_target_changed=False,
            interpretation='Forward of deterministic material risk actions, not posterior-predictive checks',
            original_prior_image='Unconstrained material estimates need not be representable by one latent',
            new_RHS_cap=8*len(methods)*len(parents),shared_acquisition=True))
        for index,parent in enumerate(parents):
            y=t.observation(parent)
            with np.load(OUT/'L/estimates'/f'{parent}.npz') as z:
                estimates={name:z[name].copy() for name in methods if '__' not in name}
            with np.load(OUT/'L/trajectory_estimates'/f'{parent}.npz') as z:
                estimates['dynamic_curvature__arc_length_trapezoid']=z['dynamic_curvature__arc_length_trapezoid'].copy()
            for method,chi in estimates.items():
                passive=bool(chi.shape==(576,) and np.isfinite(chi).all() and np.min(chi.imag)>=0 and np.min(1+chi.real)>0)
                row=dict(parent=parent,method=method,passive_positive_epsilon=passive,
                    single_P1_decoder_legality='NOT_CLAIMED',posterior_predictive=False)
                if not passive:
                    rows.append(dict(**row,status='NOT_EVALUATED_INVALID_DECISION',new_RHS=0));continue
                key=parent+'__'+method;b.context=dict(endpoint=key,purpose='RISK_ESTIMATE_FORWARD_NO_GT')
                save(b.jobdir/'attempts'/f'{key}.npz',chi=chi,purpose=np.array('risk action only'))
                try:
                    before=time.perf_counter()
                    with b.physical('decision_full_forward'):
                        state=t.backend().full_state(chi,purpose='RISK_ESTIMATE_'+key)
                        if state.residual>1e-10:raise RuntimeError('Decision backward residual contract')
                    phi=float(np.sum(abs(state.data-y)**2)/t.sigma**2)
                    save(b.jobdir/'endpoints'/f'{key}.npz',chi=chi,prediction=state.data,
                         backward_residual=state.residual,Phi=phi)
                    rows.append(dict(**row,status='COMPLETE',Phi=phi,new_RHS=8,wall_seconds=time.perf_counter()-before))
                    del state
                except BaseException as exception:
                    write(b.jobdir/'failures'/f'{key}.json',dict(error=repr(exception),cost_preserved=True));raise
            write(b.jobdir/'RESULTS.json',dict(status='IN_PROGRESS',rows=rows,GT_read=False))
            if (index+1)%32==0:print(json.dumps(dict(parents=index+1,physics='forward of risk actions only')),flush=True)
        write(b.jobdir/'RESULTS.json',dict(status='COMPLETE',rows=rows,GT_read=False,
            interpretation='Full Maxwell prediction of a risk decision is not a posterior predictive distribution'))
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:b.close(status,error)


if __name__=='__main__':main()
