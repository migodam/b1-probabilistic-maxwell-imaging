#!/usr/bin/env python3
"""Authorized first R4 paired finite-material evaluator; no test GT access.

One fixed exposed DEV set, fixed decoder-only proposals, two genuine paid Full
passes. W is acquired only through the booked adapter. No automatic new seeds,
no favorable response selection, no unbooked inverse solve or fresh LU.
"""
import os
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
from pathlib import Path
import sys, json, time, resource, traceback, gc
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.optimize import brentq
from p3r4.core import OUT, Book, Target, BudgetStop, save, write
from p3r4.reuse import FiniteMaterialCache, changed_support, broyden_update

JOB='finite_feedback_reuse_first'
RESULT=OUT/'results/II/reuse_first'
OLD=ROOT.parent/'b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration'
PARENTS=tuple(f'DEV_{i:03d}' for i in range(10,18))
SEED=202610250401
MAX_RHS=100000
RSS_CAP=3*1024**3


def rss_check():
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform!='darwin':rss*=1024
    if rss>=RSS_CAP:raise BudgetStop('pilot RSS cap')
    return rss


def observation(target,parent):
    # Exact filename allowlist: never enumerate or read offline truth arrays.
    if parent not in PARENTS:raise ValueError('parent not allowlisted')
    path=OLD/'inputs/online'/f'{parent}.npz'
    with np.load(path,allow_pickle=False) as z:
        if set(z.files)!={'observation','sigma','target_id'} or str(z['target_id'])!=target.cfg['target_id'] or float(z['sigma'])!=target.sigma:
            raise ValueError('original online target metadata')
        return z['observation'].copy()


def freeze_proposals(prior,anchor,parent_index):
    rng=np.random.default_rng(np.random.SeedSequence([SEED,parent_index]))
    q0=anchor.q;chi0=anchor.state.chi
    rows=[]
    for i in range(32):
        direction=np.r_[rng.normal(size=64),0.,0.]
        direction[:64]/=np.linalg.norm(direction[:64])
        radius=(.02,.04)[i%2]
        def distance(t):return prior.material_rms(prior.decode(q0+t*direction)[0],chi0)-radius
        hi=.0625
        while distance(hi)<0 and hi<8:hi*=2
        if distance(hi)<0:raise RuntimeError('frozen direction cannot reach preset radius within 8 latent SD; no redraw')
        scale=brentq(distance,0,hi,xtol=1e-13)
        q=q0+scale*direction
        chi,aux=prior.decode(q)
        rows.append(dict(index=i,q=q,chi=chi,direction=direction,scale=scale,radius=radius,
            actual_RMS=prior.material_rms(chi,chi0),logu=float(np.log(rng.uniform())),
            endpoint_transition_count=int(np.count_nonzero((abs(aux['field']/prior.temperature)<1)&prior.mask)),
            changed=changed_support(chi0,chi)))
    if len({r['q'].tobytes() for r in rows})!=32:raise RuntimeError('proposal duplicate')
    return rows


def append(path,obj):
    with path.open('a') as f:
        f.write(json.dumps(obj,allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())


def run():
    if RESULT.exists():raise FileExistsError('pilot output exists; never rerun paid job implicitly')
    RESULT.mkdir(parents=True)
    write(RESULT/'FROZEN_DESIGN.json',dict(parents=PARENTS,exposure='ALL_HISTORICAL_EXPOSED',seed=SEED,
        distinct_proposals_per_anchor=32,warm_passes=2,material_RMS=[.02,.04],properties='fixed original u',
        proposal_rule='independent Gaussian field64 directions normalized in latent space; decoder-only distance root; no redraw',
        full_repeats_paid=True,no_GT=True,new_RHS_cap=MAX_RHS,max_wall_seconds=1200,RSS_cap_bytes=RSS_CAP,
        Krylov='NOT_RUN; full cost reservation for iterative operator actions not implemented in this pilot',
        primary='matched evaluator cost and numerical decision replay, not posterior/MCSE'))
    book=Book('II',JOB,max_wall=1200)
    jobcounts=dict(RHS=0,states=0,LU=0)
    original_update=book._update
    def bounded_update(event,**kwargs):
        if jobcounts['RHS']+kwargs.get('rhs',0)>MAX_RHS:raise BudgetStop('pilot RHS cap')
        rss_check()
        result=original_update(event,**kwargs)
        for key,argument in [('RHS','rhs'),('states','states'),('LU','lu')]:jobcounts[key]+=kwargs.get(argument,0)
        return result
    book._update=bounded_update
    target=Target(book,keep_fields=True)
    anchor_summaries=[];all_rows=[];error=None;status='COMPLETE'
    started=time.perf_counter()
    try:
      for parent_index,parent in enumerate(PARENTS):
        folder=RESULT/parent;folder.mkdir()
        anchor_start=time.perf_counter();prior_cost=jobcounts.copy()
        load_start=time.perf_counter()
        anchor=target.load_cache(OLD/'A/anchors'/f'{parent}.npz',owner='R2_A_exposed_anchor_common_to_both_arms')
        y=observation(target,parent);model=target.backend()
        load_seconds=time.perf_counter()-load_start
        freeze_start=time.perf_counter();proposals=freeze_proposals(target.prior,anchor,parent_index)
        active=np.unique(np.concatenate([row['changed'] for row in proposals])).astype(int)
        save(folder/'FROZEN_PROPOSALS.npz',q=np.stack([row['q'] for row in proposals]),
            chi=np.stack([row['chi'] for row in proposals]),radius=np.array([row['radius'] for row in proposals]),
            logu=np.array([row['logu'] for row in proposals]),direction=np.stack([row['direction'] for row in proposals]),active=active)
        freeze_seconds=time.perf_counter()-freeze_start
        # Honest setup: inverse columns paid exactly once by adapter callback.
        book.context=dict(parent=parent,purpose='active_inverse_columns',active_size=len(active))
        setup_start=time.perf_counter()
        with book.physical('active_column_setup'):
            W=model.solve_anchor(np.eye(model.n,dtype=complex)[:,active],anchor.state,role='R4_active_inverse_columns')
        column_seconds=time.perf_counter()-setup_start
        setup_start=time.perf_counter()
        cache=FiniteMaterialCache(anchor.state.chi,anchor.state.current,anchor.state.field,anchor.state.data,
            model.G,model.S,model.incident,active,W,parent,len(active),'R4_'+JOB+'_'+parent)
        save(folder/'ACTIVE_COLUMN_CACHE.npz',active=active,W=W,K=cache.K,O=cache.O,B=cache.B,GW=cache.GW,
            base_q=anchor.q,base_chi=anchor.state.chi)
        cache_setup_seconds=time.perf_counter()-setup_start
        phi0=target.phi(anchor,y)
        jac_start=time.perf_counter();J0=target.jacobian(anchor);jac_seconds=time.perf_counter()-jac_start
        f0=target.white(anchor.state.data)
        arm_times={'Full':[],'reuse':[]};parent_rows=[]
        for repeat in range(2):
          Jsec=J0.copy()
          for candidate in proposals:
            i=candidate['index'];q=candidate['q'];chi=candidate['chi'];h=q-anchor.q
            frozen_prediction=f0+J0@h
            secant_prediction=f0+Jsec@h
            # Frozen design order avoids selecting order from measured response.
            order=('Full','reuse') if (i+repeat)%2==0 else ('reuse','Full')
            p=None;reuse=None
            durations={}
            for arm in order:
              tick=time.perf_counter()
              if arm=='Full':
                p=target.evaluate(q,purpose=f'{parent}_paired_p{repeat}_c{i:02d}')
                full_key=p.key
              else:
                reuse=cache.evaluate(chi)
                save(folder/'reuse_endpoints'/f'p{repeat}_c{i:02d}.npz',q=q,chi=chi,material_increment=chi-anchor.state.chi,
                    current=reuse['current'],field=reuse['field'],prediction=reuse['prediction'],
                    backward_residual=reuse['backward_residual'])
              durations[arm]=time.perf_counter()-tick
              arm_times[arm].append(durations[arm])
            fullphi=target.phi(p,y);reusephi=float(np.linalg.norm(reuse['prediction']-y)**2/target.sigma**2)
            prior_ratio=-.5*float(q[:64]@q[:64]-anchor.q[:64]@anchor.q[:64])
            fullratio=phi0-fullphi+prior_ratio;reuseratio=phi0-reusephi+prior_ratio
            full_decision=candidate['logu']<min(0.,fullratio)
            reuse_decision=candidate['logu']<min(0.,reuseratio)
            fullwhite=target.white(p.state.data)
            frozenphi=.5*float(np.linalg.norm(frozen_prediction-target.white(y))**2)
            secantphi=.5*float(np.linalg.norm(secant_prediction-target.white(y))**2)
            row=dict(parent=parent,repeat=repeat,candidate=i,Full_endpoint=full_key,active_size=len(active),
                increment_rank=len(candidate['changed']),endpoint_transition_count=candidate['endpoint_transition_count'],
                radius=candidate['radius'],actual_RMS=candidate['actual_RMS'],logu=candidate['logu'],
                Full_Phi=fullphi,reuse_Phi=reusephi,delta_Phi=reusephi-fullphi,
                relative_prediction_error=float(np.linalg.norm(reuse['prediction']-p.state.data)/max(np.linalg.norm(p.state.data),1e-30)),
                whitened_prediction_error=float(np.linalg.norm(target.white(reuse['prediction']-p.state.data))),
                Full_logratio=fullratio,reuse_logratio=reuseratio,Full_accept=bool(full_decision),reuse_accept=bool(reuse_decision),
                decision_mismatch=bool(full_decision!=reuse_decision),decision_margin=abs(candidate['logu']-min(0.,fullratio)),
                reuse_backward_residual=reuse['backward_residual'],small_condition_monitor=reuse['small_condition_monitor'],
                Full_wall=durations['Full'],reuse_wall=durations['reuse'],reuse_internal_wall=reuse['profile']['wall_seconds'],
                frozenJ_white_error=float(np.linalg.norm(frozen_prediction-fullwhite)),frozenJ_delta_Phi=frozenphi-fullphi,
                secant_white_error=float(np.linalg.norm(secant_prediction-fullwhite)),secant_delta_Phi=secantphi-fullphi,
                frozenJ_decision_mismatch=bool((candidate['logu']<min(0.,phi0-frozenphi+prior_ratio))!=full_decision),
                secant_decision_mismatch=bool((candidate['logu']<min(0.,phi0-secantphi+prior_ratio))!=full_decision),
                reduced_factorizations=reuse['profile']['reduced_factorizations'],reduced_rhs_columns=reuse['profile']['reduced_solve_rhs_columns'],
                cache_bytes=cache.bytes,RSS_peak_bytes=rss_check())
            append(folder/'PAIRED_ROWS.jsonl',row);parent_rows.append(row);all_rows.append(row)
            # Existing paid responses update a diagnostic only, after scoring current candidate.
            Jsec=broyden_update(Jsec,h,fullwhite-f0)
            del p,reuse
          gc.collect()
        fullwall=sum(arm_times['Full']);reusewall=sum(arm_times['reuse'])
        setupwall=column_seconds+cache_setup_seconds
        unique_fullwall=sum(r['Full_wall'] for r in parent_rows if r['repeat']==0)
        unique_reusewall=sum(r['reuse_wall'] for r in parent_rows if r['repeat']==0)
        summary=dict(parent=parent,exposure='EXPOSED',completed_distinct=32,passes=2,active_size=len(active),
            increment_rank_min=min(r['increment_rank'] for r in parent_rows),increment_rank_max=max(r['increment_rank'] for r in parent_rows),
            Full_new_RHS=64*8,reuse_new_setup_RHS=len(active),Full_new_LU=64,reuse_reduced_LU=64,
            column_setup_wall=column_seconds,cache_setup_and_persistence_wall=cache_setup_seconds,Full_endpoint_wall=fullwall,reuse_endpoint_wall=reusewall,
            reuse_with_setup_wall=reusewall+setupwall,Full_over_reuse_warm_ratio=fullwall/reusewall,
            Full_over_reuse_with_setup_ratio=fullwall/(reusewall+setupwall),
            first_pass_Full_wall=unique_fullwall,first_pass_reuse_with_setup_wall=unique_reusewall+setupwall,
            common_cache_load_and_online_read_wall=load_seconds,common_proposal_freezing_and_persistence_wall=freeze_seconds,
            cached_J_composition_wall=jac_seconds,common_historical_anchor_cold_min_RHS=8,
            historical_receiver_J_acquisition_min_RHS=32,actual_original_anchor_optimization_cost='RETAINED_IN_R2_LEDGER_NOT_FREE_OR_RECONSTRUCTED_HERE',
            cache_bytes=cache.bytes,anchor_total_wall=time.perf_counter()-anchor_start,
            job_new_costs={k:jobcounts[k]-prior_cost[k] for k in jobcounts},
            max_relative_prediction_error=max(r['relative_prediction_error'] for r in parent_rows),
            max_absolute_Phi_error=max(abs(r['delta_Phi']) for r in parent_rows),decision_mismatches=sum(r['decision_mismatch'] for r in parent_rows),
            break_even_two_pass_wall=bool(reusewall+setupwall<fullwall),fallback_count=0,rebase_count=0)
        write(folder/'SUMMARY.json',summary);anchor_summaries.append(summary)
        write(RESULT/'PROGRESS.json',dict(completed_parents=len(anchor_summaries),last=summary,job_new_costs=jobcounts))
        print(json.dumps(dict(parent=parent,active=len(active),Full_wall=fullwall,reuse_with_setup_wall=reusewall+setupwall,
            max_rel_error=summary['max_relative_prediction_error'],decision_mismatches=summary['decision_mismatches'],cost=jobcounts)),flush=True)
        del cache,W,J0,anchor,proposals,Jsec
        gc.collect()
      fullwall=sum(s['Full_endpoint_wall'] for s in anchor_summaries)
      warmwall=sum(s['reuse_endpoint_wall'] for s in anchor_summaries)
      setupwall=sum(s['column_setup_wall']+s['cache_setup_and_persistence_wall'] for s in anchor_summaries)
      write(RESULT/'RESULTS.json',dict(status='COMPLETE',parents=anchor_summaries,independent_original_objects=8,exposure='ALL_EXPOSED',
          new_costs=jobcounts,Full_endpoint_wall=fullwall,reuse_warm_endpoint_wall=warmwall,reuse_setup_wall=setupwall,
          reuse_with_setup_wall=warmwall+setupwall,Full_over_reuse_ratio=fullwall/(warmwall+setupwall),
          controlled_elapsed_at_summary=time.perf_counter()-started,RSS_peak_bytes=rss_check(),
          max_relative_prediction_error=max(r['relative_prediction_error'] for r in all_rows),max_absolute_Phi_error=max(abs(r['delta_Phi']) for r in all_rows),
          decision_mismatches=sum(r['decision_mismatch'] for r in all_rows),accepted_Full=sum(r['Full_accept'] for r in all_rows),
          frozenJ_decision_mismatches=sum(r['frozenJ_decision_mismatch'] for r in all_rows),secant_decision_mismatches=sum(r['secant_decision_mismatch'] for r in all_rows),
          reduced_LU_count=sum(r['reduced_factorizations'] for r in all_rows),reduced_rhs_columns=sum(r['reduced_rhs_columns'] for r in all_rows),
          decision='STOP_REUSE_BRANCH_NO_WALL_BREAK_EVEN' if warmwall+setupwall>=fullwall else 'EVALUATOR_SPEEDUP_ONLY_PARENT_REVIEW_REQUIRED',
          posterior_MCSE='NOT_RUN',prior_speedup_claim=False,deterministic_certificate=False,Krylov='NOT_RUN',
          proposals_are_posterior_samples=False,final_scientific_review='PARENT_OWNER'))
    except BaseException as exc:
      status='STOPPED' if isinstance(exc,BudgetStop) else 'ERROR';error=repr(exc)
      write(RESULT/'ERROR.json',dict(status=status,error=error,traceback=traceback.format_exc(),job_new_costs=jobcounts,completed_parents=len(anchor_summaries)))
      raise
    finally:
      book.close(status=status,error=error)
      write(RESULT/'FINAL_JOB_COST.json',dict(job=JOB,status=status,error=error,new_costs=jobcounts,
          total_wall_seconds=time.perf_counter()-started,RSS_peak_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
          receipt=str(book.jobdir/'CLOSED.json')))


if __name__=='__main__':
    run()
