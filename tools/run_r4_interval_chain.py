#!/usr/bin/env python3
"""Independent Full/interval arms, then after-close path audit. Never autoruns."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,json,sys,time,resource,gc,traceback,subprocess
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.interval_chain import IntervalState,approximate_state,resolve_step
from p3r4.decision_refresh import DecisionRefresh,charged_inverse

EXPERIMENT='interval_chain_v1';PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(16,24))
STEPS=1024;BETA=.002;SEED=20261025090
# Conservative disjoint reservations enforce200000 combined RHS without either
# arm reading the other's live outputs or shared counters.
RHS_CAP={'full':65536,'interval':134464};MAX_WALL=1200
OUT=ROOT/'research/p3l_r4_conditional_transport';RESULT=OUT/'results/II'/EXPERIMENT


def finite(value):
    if isinstance(value,dict):return {k:finite(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [finite(v) for v in value]
    if isinstance(value,(float,np.floating)) and not np.isfinite(value):return None
    if isinstance(value,np.generic):return value.item()
    return value


def append(path,value):
    with path.open('a') as f:f.write(json.dumps(finite(value),allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())


def point_from_full(x,q,point,phi):
    return IntervalState(np.array(x,copy=True),np.array(q,copy=True),point.state.chi.copy(),
        dict(current=point.state.current.copy(),field=point.state.field.copy(),prediction=point.state.data.copy(),
             residual=point.state.chi[:,None]*point.state.field-point.state.current,
             phi_lower=phi,phi_upper=phi,approximate_phi=phi,inverse_upper=None),True,0,point.key)


def persist(path,state,save):
    save(path,x=state.x,q=state.q,chi=state.chi,current=state.point['current'],field=state.point['field'],
         prediction=state.point['prediction'],residual=state.point['residual'],
         phi_interval=np.array([state.lower,state.upper]),exact=np.array(state.exact),refinements=np.array(state.refinements),endpoint=np.array(state.endpoint))


def compare():
    from p3r4.core import write
    if not all((RESULT/arm/'ARM_COMPLETE.json').exists() for arm in ('full','interval')):raise ValueError('both arms must close before any other-arm output is read')
    gates={arm:json.loads((RESULT/arm/'ARM_COMPLETE.json').read_text()) for arm in ('full','interval')}
    if any(v['status']!='COMPLETE' for v in gates.values()):raise ValueError('both independent arms must close COMPLETE before comparison')
    if gates['full']['input_identities']!=gates['interval']['input_identities']:raise ValueError('online initial cache changed between independent arms')
    if gates['full']['pid']==gates['interval']['pid']:raise ValueError('independent arm processes required')
    for gate in gates.values():
        receipt=json.loads(Path(gate['receipt']).read_text())
        if receipt['status']!='COMPLETE':raise ValueError('both physics Books must be closed')
    if sum(g['costs']['RHS'] for g in gates.values())>200000:raise ValueError('combinedRHS cap')
    summaries=[];mismatches=0;misses=0
    for parent in PARENTS:
        with np.load(RESULT/'full'/parent/'CHAIN.npz') as z:f={k:z[k] for k in z.files}
        with np.load(RESULT/'interval'/parent/'CHAIN.npz') as z:i={k:z[k] for k in z.files}
        flags=np.flatnonzero(f['accept']!=i['accept']);paths=np.flatnonzero(np.any(f['x_path']!=i['x_path'],axis=1))
        # Candidate and current potentials were recorded before each decision;
        # exact Full current/candidate values are opened only at this audit.
        candidate_miss=(f['candidate_phi']<i['candidate_interval'][:,0])|(f['candidate_phi']>i['candidate_interval'][:,1])
        current_miss=(f['current_phi_before']<i['current_interval'][:,0])|(f['current_phi_before']>i['current_interval'][:,1])
        history_misses=[]
        with (RESULT/'interval'/parent/'STEPS.jsonl').open() as stream:
            for line in stream:
                saved=json.loads(line);step=saved['step']
                for h in saved['history']:
                    for endpoint,value in [('current',f['current_phi_before'][step]),('candidate',f['candidate_phi'][step])]:
                        lo=h[endpoint+'_lower'];hi=h[endpoint+'_upper'];lo=0. if lo is None else lo;hi=float('inf') if hi is None else hi
                        if not lo<=value<=hi:history_misses.append(dict(step=step,attempt=h['attempt'],endpoint=endpoint))
        count=len(flags)+len(paths);miss=int(candidate_miss.sum()+current_miss.sum()+len(history_misses));mismatches+=count;misses+=miss
        summaries.append(dict(parent=parent,accept_flag_mismatches=flags.tolist(),state_path_mismatches=paths.tolist(),
                              candidate_bound_misses=np.flatnonzero(candidate_miss).tolist(),current_bound_misses=np.flatnonzero(current_miss).tolist(),all_refinement_history_bound_misses=history_misses,
                              Full_accepts=int(f['accept'].sum()),interval_accepts=int(i['accept'].sum())))
    verdict='PASS_FIXED_KERNEL_PATH_COST_AUDIT' if mismatches==0 and misses==0 else 'NO_GO_PATH_OR_ENCLOSURE_MISMATCH'
    result=dict(status=verdict,parents=summaries,total_path_or_flag_mismatches=mismatches,total_bound_misses=misses,
        arms=gates,combined_RHS=sum(g['costs']['RHS'] for g in gates.values()),
        claim_scope='fixed kernel decision and wholejob cost only; not posterior mixing/MCSE qualification',floating_certificate=False)
    if (RESULT/'COMPARISON.json').exists():raise FileExistsError('preserve closed audit')
    write(RESULT/'COMPARISON.json',result);print(json.dumps(result))
    if mismatches or misses:raise SystemExit(2)


def run(arm,release_path):
    release=json.loads(release_path.read_text())
    if release.get('released') is not True or release.get('experiment')!=EXPERIMENT or release.get('campaign_id')!='P3LR4-20261025':raise PermissionError('parent release required')
    from p3r4.core import Book,Target,BudgetStop,save,write,from_probit
    folder=RESULT/arm
    if folder.exists():raise FileExistsError('fresh arm only; no implicit replay')
    input_identities=[]
    for parent in PARENTS:
        for p in (OUT/'results/contexts/dev'/f'{parent}.npz',OUT/'results/contexts/cache'/f'{parent}.npz'):
            if not p.exists():raise FileNotFoundError(p)
            st=p.stat();input_identities.append(dict(path=str(p.resolve()),bytes=st.st_size,mtime_ns=st.st_mtime_ns))
    book=Book('II',EXPERIMENT+'_'+arm,max_wall=MAX_WALL);original=book._update;costs=dict(RHS=0,states=0,LU=0);status='COMPLETE';error=None;start=time.perf_counter();parent_summaries=[]
    def limited(event,**kw):
        if costs['RHS']+kw.get('rhs',0)>RHS_CAP[arm]:raise BudgetStop('combined200000RHS disjoint arm reservation')
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)>3*1024**3:raise BudgetStop('3GiB arm cap')
        answer=original(event,**kw)
        for k,arg in [('RHS','rhs'),('states','states'),('LU','lu')]:costs[k]+=kw.get(arg,0)
        return answer
    book._update=limited
    try:
        folder.mkdir(parents=True);target=Target(book,keep_fields=True)
        write(folder/'FROZEN_ARM.json',dict(experiment=EXPERIMENT,arm=arm,parents=PARENTS,steps=STEPS,beta=BETA,seed=SEED,
           proposal='full66Gaussianized pCN, independent process SeedSequence([seed,parent_index]); draw normal66 thenuniform eachstep',
           input_identities=input_identities,other_arm_read=False,Full_RHS_cap=RHS_CAP[arm],combined_cap=200000,max_wall=MAX_WALL,no_GT=True,
           persistence='both arms candidate+poststep full current/field/prediction/residual and x/q/chi, including rejection repeats; interval refinements additionally persisted',
           exact_cache='only charged Full gives exact potential; acceptance itself never changes originalR0 cache',floating_certificate=False))
        for parent_index,parent in enumerate(PARENTS):
            d=folder/parent;d.mkdir();initial=costs.copy();tick=time.perf_counter()
            base=target.load_cache(OUT/'results/contexts/cache'/f'{parent}.npz',owner='shared_newDEV96+64_initial_only')
            with np.load(OUT/'results/contexts/dev'/f'{parent}.npz',allow_pickle=False) as z:x=z['x'].copy();y=z['observation'].copy()
            if x.shape!=(66,) or not np.isfinite(x).all() or np.max(np.abs(from_probit(x)-base.q))>1e-12:raise ValueError('online Gaussianized anchor/cache mismatch')
            current=point_from_full(x,base.q,base,target.phi(base,y));model=target.backend();cache=None
            if arm=='interval':
                def inverse(rhs):
                    before=costs['RHS']
                    with book.physical('interval_original_R0_inverse_columns'):v=model.solve_anchor(rhs,base.state,role='R4_interval_original_inverse')
                    return v,costs['RHS']-before
                R0=charged_inverse(inverse,576);cache=DecisionRefresh(base.state.chi,base.state.current,base.state.field,model.G,model.S,model.incident,R0,target.sigma,576,True)
                save(d/'PAID_ORIGINAL_R0.npz',R0=R0,base_q=base.q,base_chi=base.state.chi,left_inverse_defect_F=np.array(cache.left_inverse_defect_F))
            persist(d/'initial.npz',current,save)
            rng=np.random.default_rng(np.random.SeedSequence([SEED,parent_index]));accept=[];paths=[current.x.copy()];candidate_phi=[];current_phi=[];candidate_intervals=[];current_intervals=[];fallbacks=0
            for step in range(STEPS):
                step_start=time.perf_counter();book.context=dict(parent=parent,step=step,arm=arm,purpose='fixed_kernel_interval_chain')
                proposed=np.sqrt(1-BETA*BETA)*current.x+BETA*rng.normal(size=66);logu=float(np.log(rng.uniform()));q=from_probit(proposed);chi=target.prior.decode(q)[0]
                current_before=current.point['approximate_phi'];record_index=[0]
                def record(state,role,event):
                    persist(d/'refinements'/f's{step:04d}_{role}_{record_index[0]:02d}_{event}.npz',state,save);record_index[0]+=1
                def paid(state,role):
                    before=costs['RHS'];p=target.evaluate(state.q,purpose=f'{parent}_chain_{arm}_{step:04d}_{role}')
                    return point_from_full(state.x,state.q,p,target.phi(p,y)),costs['RHS']-before
                if arm=='full':
                    candidate=IntervalState(proposed,q,chi,{})
                    candidate,charged=paid(candidate,'candidate')
                    if charged!=8:raise AssertionError('Full arm each candidate must charge8RHS')
                    decision=logu<min(0.,current.lower-candidate.lower);history=[];full_calls=1
                else:
                    candidate=approximate_state(cache,proposed,q,chi,y);record(candidate,'candidate','initial_approx')
                    answer=resolve_step(cache,current,candidate,y,logu,paid,4,record)
                    current,candidate,decision=answer['current'],answer['candidate'],answer['accepted'];history=answer['history'];full_calls=answer['Full_calls'];fallbacks+=full_calls
                persist(d/'candidates'/f's{step:04d}.npz',candidate,save)
                current_intervals.append([current.lower,current.upper]);candidate_intervals.append([candidate.lower,candidate.upper]);current_phi.append(current.lower if current.exact else current.point['approximate_phi']);candidate_phi.append(candidate.lower if candidate.exact else candidate.point['approximate_phi'])
                if decision:current=candidate
                # Rejections are explicit repeated states with full arrays on both arms.
                persist(d/'states'/f's{step:04d}.npz',current,save);accept.append(bool(decision));paths.append(current.x.copy())
                append(d/'STEPS.jsonl',dict(step=step,accept=bool(decision),logu=logu,current_exact_after=current.exact,
                    current_phi_interval=[current.lower,current.upper],candidate_phi_interval=[candidate.lower,candidate.upper],
                    Full_calls=full_calls,history=history,wall_seconds=time.perf_counter()-step_start,
                    RSS_peak_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024)))
            save(d/'CHAIN.npz',accept=np.array(accept),x_path=np.array(paths),candidate_phi=np.array(candidate_phi),current_phi_before=np.array(current_phi),
                 candidate_interval=np.array(candidate_intervals),current_interval=np.array(current_intervals))
            summary=dict(parent=parent,steps=STEPS,accepts=sum(accept),Full_fallback_calls=fallbacks,wall_seconds=time.perf_counter()-tick,
                         costs={k:costs[k]-initial[k] for k in costs},matrix_actions=dict(cache.counts) if cache is not None else {})
            parent_summaries.append(summary);write(d/'SUMMARY.json',summary);print(json.dumps(summary),flush=True)
            del current,base,cache;model._transfers.clear();gc.collect()
        write(folder/'RESULTS.json',dict(status='COMPLETE',arm=arm,parents=parent_summaries,costs=costs,wall_seconds=time.perf_counter()-start,no_GT=True,other_arm_read=False,
              cost_scope='actual wholejob including matched endpoint persistence,repeated states,setup,refinements,failures; not compute-only timing',mixing='NOT_QUALIFIED'))
    except BaseException as e:
        status='STOPPED' if isinstance(e,BudgetStop) else 'ERROR';error=repr(e)
        write(folder/'ERROR.json',dict(status=status,error=error,traceback=traceback.format_exc(),costs=costs));raise
    finally:
        book.close(status,error)
        write(folder/'ARM_COMPLETE.json',dict(status=status,error=error,arm=arm,pid=os.getpid(),costs=costs,wall_seconds=time.perf_counter()-start,
              receipt=str(book.jobdir/'CLOSED.json'),input_identities=input_identities,other_arm_read=False,no_GT=True,parents=parent_summaries))


def pair(release_path):
    # Coordinator owns no Book/Target and never opens GT or an active arm's
    # arrays. Each arm and the post-close comparison run in separate processes.
    release=json.loads(release_path.read_text())
    if release.get('released') is not True or release.get('experiment')!=EXPERIMENT or release.get('campaign_id')!='P3LR4-20261025':raise PermissionError('parent release required')
    compatibility=OUT/'results/II/interval_chain_pair'
    if compatibility.exists() or RESULT.exists():raise FileExistsError('pair already has artifacts; parent must reconcile rather than replay paid work')
    compatibility.mkdir(parents=True)
    started=time.monotonic();records=[]
    (compatibility/'PAIR_STARTED.json').write_text(json.dumps(dict(status='STARTED',pid=os.getpid(),experiment=EXPERIMENT,release=str(release_path.resolve()),coordinator_physics=False,GT_read=False),indent=2)+'\n')
    for arm in ('interval','full','compare'):
        command=[sys.executable,'-u',str(Path(__file__).resolve()),arm]
        if arm!='compare':command+=['--release',str(release_path.resolve())]
        child=subprocess.run(command,cwd=ROOT)
        records.append(dict(role=arm,returncode=child.returncode))
        if child.returncode:
            failure=dict(status='NO_GO' if arm=='compare' else 'FAILED',experiment=EXPERIMENT,children=records,coordinator_physics=False,GT_read=False,comparison=str(RESULT/'COMPARISON.json'),scope='no paid retry; preserve independent-arm evidence')
            name='PAIR_COMPLETE.json' if arm=='compare' else 'PAIR_FAILED.json'
            (compatibility/name).write_text(json.dumps(failure,indent=2)+'\n')
            raise SystemExit(child.returncode)
    comparison=json.loads((RESULT/'COMPARISON.json').read_text())
    if comparison['status']!='PASS_FIXED_KERNEL_PATH_COST_AUDIT':raise RuntimeError('unexpected comparison status')
    marker=dict(status='COMPLETE',experiment=EXPERIMENT,children=records,comparison=str(RESULT/'COMPARISON.json'),
                comparison_status=comparison['status'],combined_RHS=comparison['combined_RHS'],
                coordinator_physics=False,GT_read=False,coordinator_wall_seconds=time.monotonic()-started,
                scope='compatibility completion marker for independently closed interval/full arms and third-process comparison; fixed-kernel decision/cost only, no posterior qualification')
    (compatibility/'PAIR_COMPLETE.json').write_text(json.dumps(marker,indent=2)+'\n')
    print(json.dumps(marker),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('arm',nargs='?',default='pair',choices=('interval','full','compare','pair'));p.add_argument('--release',type=Path);a=p.parse_args()
    if a.arm=='compare':compare()
    else:
        if a.release is None:p.error('--release required; no autorun')
        if a.arm=='pair':pair(a.release)
        else:run(a.arm,a.release)

if __name__=='__main__':main()
