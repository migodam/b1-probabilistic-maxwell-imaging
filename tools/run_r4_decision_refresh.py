#!/usr/bin/env python3
"""Parent-release-only fixed-star decision-driven finite-solver pilot.

No GT, no training, no chain. Without matching release JSON this script cannot
construct Book/Target or perform paid work. See AUDIT/risk/DECISION_REFRESH.md.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
import argparse,json,sys,time,resource,gc,traceback
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.decision_refresh import DecisionRefresh,charged_inverse

PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(8));N_PROPOSALS=128;BETA=.002;SEED=20261025089
MAX_RHS=250000;MAX_WALL=1200;JOB='decision_refresh_star_v1'


def json_finite(value):
    if isinstance(value,dict):return {k:json_finite(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [json_finite(v) for v in value]
    if isinstance(value,(float,np.floating)) and not np.isfinite(value):return None
    if isinstance(value,np.generic):return value.item()
    return value


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--release',type=Path,required=True);a=p.parse_args()
    release=json.loads(a.release.read_text())
    if release.get('released') is not True or release.get('experiment')!=JOB or release.get('campaign_id')!='P3LR4-20261025':
        raise PermissionError('parent release for this exact prototype is required')
    # Imports/physics construction are deferred until after the explicit gate.
    from p3r4.core import OUT,Book,Target,BudgetStop,save,write,from_probit
    result=OUT/'results/II'/JOB
    if result.exists():raise FileExistsError('no implicit paid rerun; preserve prior evidence')
    for parent in PARENTS:
        for path in (OUT/'results/contexts/dev'/f'{parent}.npz',OUT/'results/contexts/cache'/f'{parent}.npz'):
            if not path.exists():raise FileNotFoundError(path)
    book=Book('II',JOB,max_wall=MAX_WALL);status='COMPLETE';error=None;rows=[];summaries=[];costs=dict(RHS=0,states=0,LU=0)
    operation_wall={'physics':0.}
    original=book._update
    def checked_update(event,**kw):
        if costs['RHS']+kw.get('rhs',0)>MAX_RHS:raise BudgetStop('fixed star pilot250kRHS cap')
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
        if rss>3*1024**3:raise BudgetStop('fixed star pilot3GiB memory cap')
        answer=original(event,**kw)
        operation_wall['physics']+=kw.get('physics',0.)
        for key,arg in [('RHS','rhs'),('states','states'),('LU','lu')]:costs[key]+=kw.get(arg,0)
        return answer
    book._update=checked_update
    start=time.perf_counter()
    try:
        result.mkdir(parents=True);target=Target(book,keep_fields=True)
        write(result/'FROZEN_DESIGN.json',dict(experiment=JOB,parents=PARENTS,proposals_per_anchor=N_PROPOSALS,seed=SEED,beta=BETA,
          proposal='full66 Gaussianized pCN from fixed online x0; exact u=ndtr(eta), unchanged original prior',same_U=True,no_GT=True,
          base='same exact paid native/probit96+64 context for all arms',refinement_max=4,
          inverse='576 actual charged solve_anchor identity-column RHS once per anchor; approximate left inverse defect retained',
          krylov='ordinary anchored-LU preconditioned GMRES, restart8,maxiter1,rtol1e-10; each actual preconditioner LU-RHS charged',
          full_audit='every proposal independently paid8-source Full; this is research cost, not deployment-free',
          max_RHS=MAX_RHS,max_wall=MAX_WALL,chains=False,posterior_MCSE='NOT_RUN',floating_certificate=False,
          all_point_arrays_saved=True,release=str(a.release.resolve())))
        for parent_index,parent in enumerate(PARENTS):
            folder=result/parent;folder.mkdir();parent_start=time.perf_counter();initial_cost=costs.copy()
            anchor=target.load_cache(OUT/'results/contexts/cache'/f'{parent}.npz',owner='R4_newDEV_common96+64_paid_anchor')
            with np.load(OUT/'results/contexts/dev'/f'{parent}.npz',allow_pickle=False) as z:
                x0=z['x'].copy();observation=z['observation'].copy()
            if x0.shape!=(66,) or not np.isfinite(x0).all():raise ValueError('finite online Gaussianized anchor x')
            native=from_probit(x0);anchor_difference=float(np.max(np.abs(native-anchor.q)))
            if anchor_difference>1e-12:raise ValueError('online x differs from exact cached anchor q')
            model=target.backend();base_phi=target.phi(anchor,observation);rng=np.random.default_rng(np.random.SeedSequence([SEED,parent_index]))
            proposals=np.sqrt(1-BETA**2)*x0[None]+BETA*rng.normal(size=(N_PROPOSALS,66));logu=np.log(rng.uniform(size=N_PROPOSALS))
            qs=np.array([from_probit(x) for x in proposals]);chis=np.array([target.prior.decode(q)[0] for q in qs])
            save(folder/'FROZEN_PROPOSALS.npz',x=proposals,q=qs,chi=chis,logu=logu,base_x=x0,base_q=anchor.q,base_Phi=np.array(base_phi))
            setup_start=time.perf_counter()
            def paid_inverse(rhs):
                before=costs['RHS']
                with book.physical('decision_refresh_inverse576'):
                    answer=model.solve_anchor(rhs,anchor.state,role='R4_decision_refresh_inverse_columns')
                return answer,costs['RHS']-before
            R0=charged_inverse(paid_inverse,model.n);inverse_setup_seconds=time.perf_counter()-setup_start
            bound_setup_start=time.perf_counter()
            cache=DecisionRefresh(anchor.state.chi,anchor.state.current,anchor.state.field,model.G,model.S,model.incident,R0,target.sigma,model.n,True)
            bound_setup_seconds=time.perf_counter()-bound_setup_start
            save(folder/'PAID_R0.npz',R0=R0,base_q=anchor.q,base_chi=anchor.state.chi,
                 left_inverse_defect_F=np.array(cache.left_inverse_defect_F),right_inverse_residual_F=np.array(cache.inverse_residual))
            write(folder/'SETUP.json',dict(inverse_RHS=model.n,inverse_setup_seconds=inverse_setup_seconds,bound_setup_seconds=bound_setup_seconds,bound_setup_dense_matrix_products=2,
                R0F=cache.R0F,SF=cache.SF,left_inverse_defect_F=cache.left_inverse_defect_F,right_inverse_residual_F=cache.inverse_residual,
                anchor_probit_difference=anchor_difference,cached_R0_bytes=R0.nbytes,base_Phi=base_phi))
            def paid_preconditioner(rhs):
                before=costs['RHS']
                with book.physical('ordinary_GMRES_actual_LU_preconditioner'):
                    answer=model.solve_anchor(rhs,anchor.state,role='R4_ordinary_GMRES_preconditioner')
                return answer,costs['RHS']-before
            parent_rows=[]
            for candidate,(q,chi,uniform) in enumerate(zip(qs,chis,logu)):
                book.context=dict(parent=parent,candidate=candidate,purpose='fixed_star_decision_audit')
                rr=cache.richardson_decision(chi,observation,base_phi,float(uniform),4)
                # Count the eight iterative material-system right-hand sides
                # in addition to every actual cached-LU preconditioner action.
                book._update('ITERATIVE_MATERIAL_STATE_ATTEMPT',rhs=8,states=1,lu=0,
                             role='ordinary_GMRES_source_systems',parent=parent,candidate=candidate)
                gm=cache.ordinary_krylov(chi,observation,base_phi,float(uniform),paid_preconditioner=paid_preconditioner)
                # Save approximate currents/fields before any Full-audit result.
                persistence_start=time.perf_counter()
                for iteration,point in enumerate(rr['approximate_points']):
                    save(folder/'approximate'/f'c{candidate:03d}_richardson_i{iteration}.npz',q=q,chi=chi,**point)
                save(folder/'approximate'/f'c{candidate:03d}_gmres.npz',q=q,chi=chi,**{k:gm['point'][k] for k in ('current','field','prediction','residual')})
                approximate_persistence_seconds=time.perf_counter()-persistence_start
                physics_before=operation_wall['physics'];tick=time.perf_counter();full=target.evaluate(q,purpose=f'{parent}_decision_refresh_full_audit_c{candidate:03d}');full_seconds=time.perf_counter()-tick
                full_physics_seconds=operation_wall['physics']-physics_before
                phi=target.phi(full,observation);full_accept=bool(uniform<min(0.,base_phi-phi))
                row=dict(parent=parent,candidate=candidate,full_endpoint=full.key,Full_Phi=phi,Full_accept=full_accept,Full_seconds_including_persistence=full_seconds,Full_physics_seconds=full_physics_seconds,approximate_persistence_seconds=approximate_persistence_seconds,
                         RSS_peak_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024),logu=float(uniform),
                         material_RMS=target.prior.material_rms(chi,anchor.state.chi),richardson={},gmres={})
                for name,arm in [('richardson',rr),('gmres',gm)]:
                    point=arm['point'];decision=arm['decision'];fallback=decision=='UNKNOWN';resolved=full_accept if fallback else decision=='ACCEPT'
                    row[name]=dict(decision_before_Full=decision,fallback=fallback,final_accept=resolved,decision_mismatch=resolved!=full_accept,
                        bound_contains_Full=point['phi_lower']<=phi<=point['phi_upper'],
                        phi_lower=point['phi_lower'],phi_upper=point['phi_upper'],approximate_phi=point['approximate_phi'],
                        eta_upper=point['eta_upper'],inverse_upper=point['inverse_upper'],rho_F=point['rho_F'],
                        relative_prediction_error=float(np.linalg.norm(point['prediction']-full.state.data)/max(np.linalg.norm(full.state.data),1e-300)),
                        wall_seconds=arm['wall_seconds'],counts=arm['counts'],
                        equivalent_independent_arm_Full_fallback_RHS=8 if fallback else 0,
                        refinement_history=arm.get('history'),iterations_per_source=arm.get('iterations_per_source'),gmres_info=arm.get('gmres_info'))
                with (folder/'PAIRED_ROWS.jsonl').open('a') as f:f.write(json.dumps(json_finite(row),allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
                for arm in (rr,gm):
                    for key,value in arm['counts'].items():book.counts['decision_refresh_'+key]+=value
                rows.append(row);parent_rows.append(row);del full,rr,gm
            summary=dict(parent=parent,proposals=len(parent_rows),costs={k:costs[k]-initial_cost[k] for k in costs},
                         inverse_setup_RHS=model.n,inverse_setup_seconds=inverse_setup_seconds,bound_setup_seconds=bound_setup_seconds,elapsed_seconds=time.perf_counter()-parent_start)
            for name in ('richardson','gmres'):
                arms=[r[name] for r in parent_rows];fallbacks=sum(r['fallback'] for r in arms)
                rhs=sum(r['counts'].get('krylov_preconditioner_LU_RHS',0) for r in arms)
                summary[name]=dict(unknown_fallbacks=fallbacks,decision_mismatches=sum(r['decision_mismatch'] for r in arms),
                    Full_bound_misses=sum(not r['bound_contains_Full'] for r in arms),evaluation_wall_seconds=sum(r['wall_seconds'] for r in arms),
                    paid_preconditioner_RHS=rhs,equivalent_standalone_RHS_with_bound_setup=model.n+rhs+8*fallbacks,
                    shared_audit_deployment_accounting='standalone cost is accounting counterfactual; audit Full was actually paid once for every proposal')
            summaries.append(summary);write(folder/'SUMMARY.json',summary)
            write(result/'PROGRESS.json',dict(completed=len(summaries),costs=costs,last=summary));print(json.dumps(summary),flush=True)
            del cache,R0,anchor,chis;model._transfers.clear();gc.collect()
        write(result/'RESULTS.json',dict(status='COMPLETE',parents=summaries,n_proposals=len(rows),actual_paid_costs=costs,
             timing_scope='Fullwall includes endpoint persistence; approximate compute and persistence are separately saved; wholejob is actual end-to-end cost',actual_Full_audit_RHS=8*len(rows),inverse_setup_RHS=576*len(PARENTS),wall_seconds=time.perf_counter()-start,
             no_GT=True,chains=False,posterior_MCSE='NOT_RUN',floating_certificate=False,
             limitations=['one-step star decisions are not chain speedup','accepted interval does not create exact cache','norms and arithmetic are not directed interval certification','all Full audits are research costs','ordinary GMRES LU applications actually charged']))
    except BaseException as e:
        status='STOPPED' if isinstance(e,BudgetStop) else 'ERROR';error=repr(e)
        write(result/'ERROR.json',dict(status=status,error=error,traceback=traceback.format_exc(),costs=costs));raise
    finally:
        book.close(status,error);write(result/'FINAL_JOB_COST.json',dict(status=status,error=error,costs=costs,receipt=str(book.jobdir/'CLOSED.json')))

if __name__=='__main__':main()
