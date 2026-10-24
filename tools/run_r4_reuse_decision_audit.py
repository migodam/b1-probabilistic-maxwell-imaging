#!/usr/bin/env python3
"""Predeclared one-step star pCN decision audit; never a speedup rescue."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
import sys,json,time,traceback,gc
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,Target,BudgetStop,save,write
from p3r4.reuse import FiniteMaterialCache,changed_support
from run_r4_reuse_pilot import observation,OLD,PARENTS,rss_check,append

JOB='finite_feedback_pcn_decision_audit'
RESULT=OUT/'results/II/reuse_pcn_decision_audit'
SEED=20261025078
BETA=.002
MAX_RHS=5000


def run():
    if RESULT.exists():raise FileExistsError('no implicit rerun')
    RESULT.mkdir(parents=True)
    write(RESULT/'FROZEN_DESIGN.json',dict(parents=PARENTS,seed=SEED,beta=BETA,proposals_per_anchor=32,
        proposal='field64 pCN from fixed anchor: sqrt(1-beta^2)*z0+beta*nu; original u fixed',
        prior_reversible=True,log_MH_ratio='Phi0-Phi_candidate',same_uniform=True,
        no_GT=True,all_exposed=True,amplitude_adaptation=False,max_RHS=MAX_RHS,max_wall=300,
        purpose='one-step star acceptance decision audit; not chains or speedup rescue',
        prior_speedup_stop_preserved='STOP_REUSE_BRANCH_NO_WALL_BREAK_EVEN',
        fallback='if support outside paid active columns, separately paid Full fallback'))
    book=Book('II',JOB,max_wall=300);original_update=book._update
    costs=dict(RHS=0,states=0,LU=0)
    def limited(event,**kw):
        rss_check()
        if costs['RHS']+kw.get('rhs',0)>MAX_RHS:raise BudgetStop('decision-audit RHS cap')
        result=original_update(event,**kw)
        for key,arg in [('RHS','rhs'),('states','states'),('LU','lu')]:costs[key]+=kw.get(arg,0)
        return result
    book._update=limited
    target=Target(book,keep_fields=True);status='COMPLETE';error=None;rows=[];summaries=[];started=time.perf_counter()
    try:
      for parent_index,parent in enumerate(PARENTS):
        folder=RESULT/parent;folder.mkdir()
        anchor=target.load_cache(OLD/'A/anchors'/f'{parent}.npz','R2_common_exposed_anchor')
        model=target.backend();y=observation(target,parent);phi0=target.phi(anchor,y)
        path=OUT/'results/II/reuse_first'/parent/'ACTIVE_COLUMN_CACHE.npz'
        with np.load(path,allow_pickle=False) as a:
            if not np.array_equal(a['base_q'],anchor.q) or not np.array_equal(a['base_chi'],anchor.state.chi):raise ValueError('old W base mismatch')
            active=a['active'].copy();W=a['W'].copy()
        cache=FiniteMaterialCache(anchor.state.chi,anchor.state.current,anchor.state.field,anchor.state.data,
            model.G,model.S,model.incident,active,W,parent,len(active),'finite_feedback_reuse_first_'+parent)
        rng=np.random.default_rng(np.random.SeedSequence([SEED,parent_index]));proposals=[]
        for i in range(32):
            q=anchor.q.copy();q[:64]=np.sqrt(1-BETA*BETA)*q[:64]+BETA*rng.normal(size=64)
            chi,_=target.prior.decode(q)
            proposals.append((q,chi,float(np.log(rng.uniform()))))
        save(folder/'FROZEN_PROPOSALS.npz',q=np.stack([z[0] for z in proposals]),chi=np.stack([z[1] for z in proposals]),logu=np.array([z[2] for z in proposals]),active=active,beta=np.array(BETA))
        parent_rows=[]
        for i,(q,chi,logu) in enumerate(proposals):
            support=changed_support(anchor.state.chi,chi);fits=bool(np.isin(support,active).all())
            # Support-based fallback is decided before either candidate response.
            full=target.evaluate(q,purpose=f'{parent}_decision_Full_c{i:02d}')
            fullphi=target.phi(full,y)
            if fits:
                reuse=cache.evaluate(chi)
                raw=dict(current=reuse['current'],field=reuse['field'],prediction=reuse['prediction'],backward_residual=reuse['backward_residual'])
                key=None
            else:
                fallback=target.evaluate(q,purpose=f'{parent}_decision_reuse_Full_fallback_c{i:02d}')
                raw=dict(current=fallback.state.current,field=fallback.state.field,prediction=fallback.state.data,backward_residual=fallback.state.residual)
                key=fallback.key
            save(folder/'reuse_endpoints'/f'c{i:02d}.npz',q=q,chi=chi,material_increment=chi-anchor.state.chi,**raw)
            reusephi=float(np.linalg.norm(raw['prediction']-y)**2/target.sigma**2)
            fullratio=phi0-fullphi;reuseratio=phi0-reusephi
            fullaccept=bool(logu<min(0.,fullratio));reuseaccept=bool(logu<min(0.,reuseratio))
            row=dict(parent=parent,candidate=i,beta=BETA,Full_endpoint=full.key,reuse_fallback_endpoint=key,
                active_size=len(active),increment_rank=len(support),support_fits=bool(fits),Full_fallback=not fits,
                material_RMS=target.prior.material_rms(chi,anchor.state.chi),logu=logu,Full_Phi=fullphi,reuse_Phi=reusephi,
                Full_logratio=fullratio,reuse_logratio=reuseratio,Full_accept=fullaccept,reuse_accept=reuseaccept,
                decision_mismatch=fullaccept!=reuseaccept,decision_margin=abs(logu-min(0.,fullratio)),
                relative_prediction_error=float(np.linalg.norm(raw['prediction']-full.state.data)/max(np.linalg.norm(full.state.data),1e-30)),
                delta_Phi=reusephi-fullphi,backward_residual=float(raw['backward_residual']),RSS_peak_bytes=rss_check())
            append(folder/'PAIRED_ROWS.jsonl',row);rows.append(row);parent_rows.append(row)
            del full,raw
            if not fits:del fallback
        summary=dict(parent=parent,Full_accepts=sum(z['Full_accept'] for z in parent_rows),rejections=sum(not z['Full_accept'] for z in parent_rows),
            decision_mismatches=sum(z['decision_mismatch'] for z in parent_rows),fallbacks=sum(z['Full_fallback'] for z in parent_rows),
            minimum_decision_margin=min(z['decision_margin'] for z in parent_rows))
        write(folder/'SUMMARY.json',summary);summaries.append(summary)
        write(RESULT/'PROGRESS.json',dict(completed=len(summaries),costs=costs,last=summary))
        print(json.dumps(summary),flush=True)
        del cache,W,anchor,proposals
        gc.collect()
      write(RESULT/'RESULTS.json',dict(status='COMPLETE',seed=SEED,beta=BETA,parents=summaries,
          distinct_candidates=len(rows),Full_accepts=sum(z['Full_accept'] for z in rows),Full_rejections=sum(not z['Full_accept'] for z in rows),
          decision_mismatches=sum(z['decision_mismatch'] for z in rows),fallbacks=sum(z['Full_fallback'] for z in rows),
          min_decision_margin=min(z['decision_margin'] for z in rows),max_relative_prediction_error=max(z['relative_prediction_error'] for z in rows),
          max_absolute_Phi_error=max(abs(z['delta_Phi']) for z in rows),max_backward_residual=max(z['backward_residual'] for z in rows),
          material_RMS_min=min(z['material_RMS'] for z in rows),material_RMS_max=max(z['material_RMS'] for z in rows),
          new_costs=costs,earlier_paid_inverse_column_RHS=sum(len(np.load(OUT/'results/II/reuse_first'/p/'ACTIVE_COLUMN_CACHE.npz')['active']) for p in PARENTS),
          new_inverse_column_RHS=0,elapsed_seconds=time.perf_counter()-started,RSS_peak_bytes=rss_check(),
          prior_reversible_likelihood_only=True,speedup_stop_preserved='STOP_REUSE_BRANCH_NO_WALL_BREAK_EVEN',
          chains=False,posterior_MCSE='NOT_RUN',deterministic_certificate=False,scientific_review='PARENT_OWNER'))
    except BaseException as exc:
        status='STOPPED' if isinstance(exc,BudgetStop) else 'ERROR';error=repr(exc)
        write(RESULT/'ERROR.json',dict(status=status,error=error,traceback=traceback.format_exc(),costs=costs))
        raise
    finally:
        book.close(status,error)
        write(RESULT/'FINAL_JOB_COST.json',dict(status=status,error=error,new_costs=costs,elapsed_seconds=time.perf_counter()-started,receipt=str(book.jobdir/'CLOSED.json')))

if __name__=='__main__':run()
