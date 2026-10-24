#!/usr/bin/env python3
"""Opt-in common stored-I/O/endpoint-once cost revision; no automatic release.

Independent interval and optimized Full arms; after-close whole-path/proof audit.
Closed original/tight pairs are untouched. Ordinary float64, not certification.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import argparse,json,sys,time,resource,subprocess,gc,traceback
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.interval_chain import IntervalState,approximate_state,resolve_step
from p3r4.decision_refresh import charged_inverse
from p3r4.decision_refresh_fast import DecisionRefreshFast
from p3r4.evaluation import checkpoint_identity
OUT=ROOT/'research/p3l_r4_conditional_transport';EXPERIMENT='interval_chain_fast_v1';RESULT=OUT/'results/II'/EXPERIMENT
PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(40,48));STEPS=1024;BETA=.002;SEED=20261025131
CAP={'interval':134464,'full':65536};MAX_WALL=1200
# Explicit future wrappers may select a new ledger stage/cohort. Defaults and
# closed source snapshots retain the original StageII DEV protocol.
LEDGER_STAGE='II'
CONTEXT_DIR=OUT/'results/contexts/dev'
CACHE_DIR=OUT/'results/contexts/cache'
# Prospective grid wrappers only; ordinary existing defaults stay unchanged.
CASE_DESIGNS=None
REUSE_ORIGINAL_R0_PER_PARENT=False
MATERIAL_MOVE_AUDIT=False
EXACT_SCALAR_COPY_BYTES=0
# Opt-in extension hooks; the closed fast-pair defaults remain unchanged.
CACHE_CLASS_BY_ARM={}
CACHE_CONFIGURE_BY_ARM={}
RESOLVER_BY_ARM={}
FROZEN_DESIGN_EXTENSION={}


def clean(v):
    if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    if isinstance(v,np.ndarray):return clean(v.tolist())
    if isinstance(v,(float,np.floating)) and not np.isfinite(v):return None
    if isinstance(v,np.generic):return v.item()
    return v


def append(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('a') as f:f.write(json.dumps(clean(value),allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())


def exact_state(x,q,p,phi,endpoint_id):
    global EXACT_SCALAR_COPY_BYTES
    EXACT_SCALAR_COPY_BYTES+=np.asarray(x).nbytes+np.asarray(q).nbytes
    # Native Full arrays already have one paid immutable on-disk endpoint. Keep
    # views for scalar decisions; do not duplicate/copy their full fields.
    state=IntervalState(np.array(x,copy=True),np.array(q,copy=True),p.state.chi,
        dict(current=p.state.current,field=p.state.field,prediction=p.state.data,
             residual=None,
             phi_lower=phi,phi_upper=phi,approximate_phi=phi,inverse_upper=None),True,0,p.key)
    state.endpoint_id=endpoint_id
    return state


class EndpointWriter:
    def __init__(self,folder,book,save,incident):
        self.incident=np.asarray(incident)
        self.folder=Path(folder);self.book=book;self.save=save;self.latest={};self.versions={};self.paths=set();self.stats=dict(array_writes=0,array_overwrites=0,array_bytes_written=0,array_write_seconds=0.,Target_endpoint_references=0,attempt_writes=0,proof_events=0,rhs_norm_entries=0,rhs_norm_temporary_bytes=0)
    def attempt(self,step,x,q,chi,logu):
        path=self.folder/'attempts'/f'c{step:04d}.npz';tick=time.perf_counter();self.save(path,x=x,q=q,chi=chi,logu=np.array(logu))
        self.paths.add(path);self.stats['attempt_writes']+=1;self.stats['array_bytes_written']+=path.stat().st_size;self.stats['array_write_seconds']+=time.perf_counter()-tick
    def proof(self,state,role,event,step):
        self.stats['proof_events']+=1
        rho=state.point.get('rho_F')
        if rho is None and not state.exact:rho=float(np.linalg.norm(state.point['residual']))
        append(self.folder/'ENDPOINT_PROOFS.jsonl',dict(step=step,endpoint_id=state.endpoint_id,role=role,event=event,
           accuracy_iteration=state.refinements,exact=state.exact,Phi_interval=[state.lower,state.upper],
           rho_F=rho,residual_white_norm=state.point.get('residual_white_norm'),
           error_white_upper=state.point.get('error_white_upper'),eta_upper=state.point.get('eta_upper'),inverse_upper=state.point.get('inverse_upper'),
           reconstruction='originalR0/base/chi; jinit=j0+R0(dchi*E0), then indicated cumulative Richardson recurrences; Full uses paid endpoint',same_physical_endpoint=True))
    def final(self,state):
        endpoint_id=state.endpoint_id;previous=self.latest.get(endpoint_id)
        signature=(state.exact,state.refinements,state.endpoint)
        if previous is not None and previous['signature']==signature:return
        version=self.versions.get(endpoint_id,0)+1;self.versions[endpoint_id]=version
        if state.exact and state.endpoint!='initial_cache':
            path=self.book.jobdir/'endpoints'/f'{state.endpoint}.npz'
            if not path.exists():raise FileNotFoundError('paidFull endpoint must already be immediately persisted')
            self.stats['Target_endpoint_references']+=1
            self.paths.add(self.book.jobdir/'attempts'/f'{state.endpoint}.npz')
        elif state.endpoint=='initial_cache':
            path=Path(state.initial_cache_path)
        else:
            path=self.folder/'endpoints'/f'{endpoint_id}.npz';tick=time.perf_counter();overwrite=path.exists()
            rhs=state.chi[:,None]*self.incident;rhs_norm=float(np.linalg.norm(rhs));rho=state.point['rho_F'] if 'rho_F' in state.point else float(np.linalg.norm(state.point['residual']))
            backward=rho/(rhs_norm if rhs_norm else 1.)
            self.stats['rhs_norm_entries']+=rhs.size;self.stats['rhs_norm_temporary_bytes']+=rhs.nbytes
            self.save(path,q=state.q,chi=state.chi,current=state.point['current'],field=state.point['field'],prediction=state.point['prediction'],
                      backward_residual=np.array(backward))
            self.stats['array_writes']+=1;self.stats['array_overwrites']+=int(overwrite);self.stats['array_bytes_written']+=path.stat().st_size;self.stats['array_write_seconds']+=time.perf_counter()-tick
        self.paths.add(path)
        entry=dict(signature=signature,path=str(path.resolve()),array_version=version,exact=state.exact,accuracy_iteration=state.refinements)
        self.latest[endpoint_id]=entry
        append(self.folder/'ENDPOINT_VERSIONS.jsonl',dict(endpoint_id=endpoint_id,**entry,Phi_interval=[state.lower,state.upper],
            physical_material_unchanged_on_overwrite=previous is not None,paidFull_never_overwritten=True))
    def seal(self):return [dict(endpoint_id=k,**v,identity=checkpoint_identity(v['path'])) for k,v in self.latest.items()]


def run(arm,release_path):
    release=json.loads(release_path.read_text())
    if release.get('released') is not True or release.get('experiment')!=EXPERIMENT or release.get('campaign_id')!='P3LR4-20261025':raise PermissionError('matching parentrelease required')
    # Only this new revision opts in. Defaults in core.save stay compressed.
    if os.environ.get('R4_ARRAY_CODEC')!='stored':raise ValueError('both fastarms require explicit R4_ARRAY_CODEC=stored')
    from p3r4.core import Book,Target,BudgetStop,save,write,from_probit,array_io_profile
    folder=RESULT/arm
    if folder.exists():raise FileExistsError('no paidreplay')
    identities=[]
    for parent in PARENTS:
        for path in (CACHE_DIR/f'{parent}.npz',CONTEXT_DIR/f'{parent}.npz'):
            identities.append(checkpoint_identity(path))
    book=Book(LEDGER_STAGE,EXPERIMENT+'_'+arm,max_wall=MAX_WALL);original=book._update;costs=dict(RHS=0,states=0,LU=0);status='COMPLETE';error=None;summaries=[];all_seals=[];all_array_seals=[];started=time.perf_counter()
    def limited(event,**kw):
        if costs['RHS']+kw.get('rhs',0)>CAP[arm]:raise BudgetStop('pair200000RHS reservationcap')
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)>3*1024**3:raise BudgetStop('3GiB armcap')
        value=original(event,**kw)
        for k,a in [('RHS','rhs'),('states','states'),('LU','lu')]:costs[k]+=kw.get(a,0)
        return value
    book._update=limited
    try:
        folder.mkdir(parents=True);target=Target(book,keep_fields=True)
        design=dict(experiment=EXPERIMENT,arm=arm,parents=PARENTS,steps=STEPS,beta=BETA,seed=SEED,array_codec='stored',input_identities=identities,
          same_kernel=True,two_endpoint_intervals=True,tight_bound_unchanged=True,algebra='rowG2 exact norm identity+immutable2entry boundcache',
          evidence_contract='one canonical latest-array file per approximate material endpoint; paidFullnative endpoint referenced, never duplicated; same-endpoint accuracy replacement versioned with append-only proof scalars; repeatstates references+exactx',
          common_attempts='x/q/chi immediately stored for every proposal including failures',intermediate_refinements='scalarenclosure/residualhistory only, deterministically reconstructible from originalR0/base+chi',
          caps=CAP,max_wall=MAX_WALL,GT_read=False,floating_certificate=False,posterior_reference='NONE',success='0flag/path/allproofenclosuremiss and >=20%wholejobwall saving including coldR0/setup/allpersistence vscommonstoredoptimizedFull')
        design.update(FROZEN_DESIGN_EXTENSION)
        write(folder/'FROZEN_ARM.json',design);write(book.jobdir/'COMMAND.json',dict(argv=sys.argv,array_codec=os.environ['R4_ARRAY_CODEC'],noGT=True))
        cases=CASE_DESIGNS if CASE_DESIGNS is not None else [dict(parent=p,case_id=p,parent_index=i,seed=SEED,beta=BETA,steps=STEPS) for i,p in enumerate(PARENTS)]
        inverse_pool={}
        for case in cases:
            parent=case['parent'];parent_index=case['parent_index'];case_id=case['case_id'];case_steps=case['steps'];case_beta=case['beta'];case_seed=case['seed']
            candidate=None
            if CASE_DESIGNS is not None:book.context=dict(parent=parent,case_id=case_id,beta=case_beta,seed=case_seed,arm=arm,purpose='original_anchor_case_setup')
            d=folder/case_id;d.mkdir();before=costs.copy();tick=time.perf_counter();cachepath=CACHE_DIR/f'{parent}.npz'
            base=target.load_cache(cachepath,owner='paid online96+64 initialanchor, commonbotharms')
            with np.load(CONTEXT_DIR/f'{parent}.npz') as z:x0=z['x'];y=z['observation']
            if not np.isfinite(x0).all() or np.max(abs(from_probit(x0)-base.q))>1e-12:raise ValueError('onlineexactprobitbase')
            current=exact_state(x0,base.q,base,target.phi(base,y),'initial');current.endpoint='initial_cache';current.initial_cache_path=str(cachepath)
            model=target.backend();writer=EndpointWriter(d,book,save,model.incident);writer.final(current);cache=None;setup={}
            if arm!='full':
                setup_start=time.perf_counter()
                def acquire(rhs):
                    old=costs['RHS']
                    with book.physical('fast_original_R0_576RHS'):R=model.solve_anchor(rhs,base.state,role='R4_fast_R0')
                    return R,costs['RHS']-old
                reused=REUSE_ORIGINAL_R0_PER_PARENT and parent in inverse_pool
                if reused:R0=inverse_pool[parent]['R0']
                else:
                    R0=charged_inverse(acquire,576)
                    if REUSE_ORIGINAL_R0_PER_PARENT:
                        R0.flags.writeable=False;inverse_pool[parent]=dict(R0=R0,first_case_id=case_id)
                cache=CACHE_CLASS_BY_ARM.get(arm,DecisionRefreshFast)(base.state.chi,base.state.current,base.state.field,model.G,model.S,model.incident,R0,target.sigma,576,True)
                if arm in CACHE_CONFIGURE_BY_ARM:CACHE_CONFIGURE_BY_ARM[arm](cache,book,costs)
                reconstruction=(folder/'shared_R0'/f'{parent}.npz') if REUSE_ORIGINAL_R0_PER_PARENT else d/'RECONSTRUCTION_R0.npz'
                if not reused:save(reconstruction,R0=R0,base_q=base.q,base_chi=base.state.chi,base_current=base.state.current,base_field=base.state.field)
                writer.paths.add(reconstruction)
                setup=dict(wall_including_R0_boundcopy_and_persistence=time.perf_counter()-setup_start,RHS=0 if reused else 576,fast_setup=cache.fast_setup_profile,
                           weak_background_identity=checkpoint_identity(OUT/'inputs/background/weak.npz'),base_cache_identity=checkpoint_identity(cachepath))
                if REUSE_ORIGINAL_R0_PER_PARENT:setup.update(shared_original_R0_cache_hit=reused,shared_original_R0_first_case=inverse_pool[parent]['first_case_id'],shared_original_R0_path=str(reconstruction),cold_owner_RHS=576)
            if CASE_DESIGNS is not None:write(d/'SETUP.json',clean(setup))
            rng=np.random.default_rng(np.random.SeedSequence([case_seed,parent_index]));xs=[x0.copy()];flags=[];cand_phi=[];cur_phi=[];cand_intervals=[];cur_intervals=[];fallbacks=0;moves=[];proposals=[];anchor_moves=[]
            for step in range(case_steps):
                candidate=None
                endpoint_id=f'c{step:04d}';proposed=np.sqrt(1-case_beta*case_beta)*current.x+case_beta*rng.normal(size=66);logu=float(np.log(rng.uniform()));q=from_probit(proposed);chi=target.prior.decode(q)[0]
                proposal_move=float(np.linalg.norm(chi-current.chi)/24) if MATERIAL_MOVE_AUDIT else None
                writer.attempt(step,proposed,q,chi,logu);book.context=dict(parent=parent,step=step,arm=arm,purpose='commonstored_fixedkernel_fastpair')
                if CASE_DESIGNS is not None:book.context.update(case_id=case_id,beta=case_beta,seed=case_seed)
                def paid(state,role):
                    old=costs['RHS'];p=target.evaluate(state.q,purpose=f'{parent}_{arm}_s{step:04d}_{role}')
                    return exact_state(state.x,state.q,p,target.phi(p,y),state.endpoint_id),costs['RHS']-old
                def record(state,role,event):writer.proof(state,role,event,step)
                if arm=='full':
                    candidate=IntervalState(proposed,q,chi,{});candidate.endpoint_id=endpoint_id;candidate,charged=paid(candidate,'candidate')
                    if charged!=8:raise AssertionError('Fullsource8RHS')
                    accepted=logu<min(0.,current.lower-candidate.lower);history=[];full_calls=1
                    record(candidate,'candidate','paid_Full')
                else:
                    candidate=approximate_state(cache,proposed,q,chi,y);candidate.endpoint_id=endpoint_id;record(candidate,'candidate','initial_approx')
                    answer=RESOLVER_BY_ARM.get(arm,resolve_step)(cache,current,candidate,y,logu,paid,4,record)
                    current,candidate,accepted=answer['current'],answer['candidate'],answer['accepted'];history=answer['history'];full_calls=answer['Full_calls'];fallbacks+=full_calls
                writer.final(current);writer.final(candidate)
                cur_intervals.append([current.lower,current.upper]);cand_intervals.append([candidate.lower,candidate.upper]);cur_phi.append(current.point['approximate_phi']);cand_phi.append(candidate.point['approximate_phi'])
                before_endpoint=current.endpoint_id
                if accepted:current=candidate
                flags.append(bool(accepted));xs.append(current.x.copy())
                if MATERIAL_MOVE_AUDIT:moves.append(proposal_move if accepted else 0.);proposals.append(proposal_move);anchor_moves.append(float(np.linalg.norm(current.chi-base.state.chi)/24))
                reference=dict(step=step,accepted=bool(accepted),x=current.x,logu=logu,current_before=before_endpoint,candidate=endpoint_id,current_after=current.endpoint_id,
                    current_exact=current.exact,current_array_version=writer.versions[current.endpoint_id],candidate_array_version=writer.versions[candidate.endpoint_id],
                    current_interval=[current.lower,current.upper],history=history,Full_calls=full_calls)
                if MATERIAL_MOVE_AUDIT:reference.update(material_move_RMS=moves[-1],proposal_material_RMS=proposals[-1],state_material_RMS_from_anchor=anchor_moves[-1])
                append(d/'CHAIN_REFERENCES.jsonl',reference)
            move_arrays=dict(material_move_RMS=np.array(moves),proposal_material_RMS=np.array(proposals),state_material_RMS_from_anchor=np.array(anchor_moves)) if MATERIAL_MOVE_AUDIT else {}
            save(d/'CHAIN.npz',x_path=np.array(xs),accept=np.array(flags),candidate_phi=np.array(cand_phi),current_phi_before=np.array(cur_phi),
                candidate_interval=np.array(cand_intervals),current_interval=np.array(cur_intervals),**move_arrays)
            writer.paths.add(d/'CHAIN.npz')
            seals=writer.seal();all_seals.extend(seals);all_array_seals.extend(checkpoint_identity(path) for path in sorted(writer.paths))
            summary=dict(parent=parent,steps=case_steps,accepts=sum(flags),fallbacks=fallbacks,costs={k:costs[k]-before[k] for k in costs},wholeparentwall=time.perf_counter()-tick,
                setup=setup,io=writer.stats,cumulative_core_array_io=array_io_profile(),cumulative_exact_xq_copy_bytes=EXACT_SCALAR_COPY_BYTES,matrix_actions=dict(cache.counts) if cache else {},latest_endpoint_arrays=seals)
            if CASE_DESIGNS is not None:summary.update(case_id=case_id,beta=case_beta,seed=case_seed)
            if MATERIAL_MOVE_AUDIT:summary.update(mean_material_move_RMS=float(np.mean(moves)),mean_proposal_material_RMS=float(np.mean(proposals)),max_proposal_material_RMS=float(np.max(proposals)),terminal_material_RMS_from_anchor=anchor_moves[-1])
            summaries.append(summary);write(d/'SUMMARY.json',clean(summary));print(json.dumps(clean({k:v for k,v in summary.items() if k!='latest_endpoint_arrays'})),flush=True)
            del current,base,cache,writer;model._transfers.clear();gc.collect()
        write(folder/'RESULTS.json',clean(dict(status='COMPLETE',arm=arm,costs=costs,parents=summaries,wall_seconds=time.perf_counter()-started,array_codec='stored',GT_read=False)))
    except BaseException as exc:
        if MATERIAL_MOVE_AUDIT:
            # Keep any already-computed last fields at a hard-cap stop, while
            # retaining UNKNOWN status and never making an unverified decision.
            try:
                if 'writer' in locals():
                    for state in (locals().get('current'),locals().get('candidate')):
                        if state is not None and state.point.get('current') is not None:writer.final(state)
                    all_seals.extend(writer.seal());all_array_seals.extend(checkpoint_identity(path) for path in sorted(writer.paths))
            except BaseException as publication_error:
                write(folder/'PARTIAL_ARRAY_PUBLICATION_ERROR.json',dict(error=repr(publication_error),original_error=repr(exc)))
        status='STOPPED' if isinstance(exc,BudgetStop) else 'ERROR';error=repr(exc);write(folder/'ERROR.json',dict(status=status,error=error,traceback=traceback.format_exc(),costs=costs));raise
    finally:
        book.close(status,error)
        write(folder/'ARM_COMPLETE.json',clean(dict(status=status,pid=os.getpid(),arm=arm,costs=costs,receipt=str(book.jobdir/'CLOSED.json'),
            wall_seconds=time.perf_counter()-started,exact_xq_copy_bytes=EXACT_SCALAR_COPY_BYTES,latest_array_seals=all_seals,all_array_seals=all_array_seals,core_array_io=array_io_profile(),input_identities=identities,array_codec='stored',GT_read=False)))


def compare():
    if not all((RESULT/a/'ARM_COMPLETE.json').exists() for a in ('full','interval')):raise ValueError('botharmsmustclosebeforeother-armreads')
    gates={a:json.loads((RESULT/a/'ARM_COMPLETE.json').read_text()) for a in ('full','interval')}
    if any(v['status']!='COMPLETE' or v['array_codec']!='stored' for v in gates.values()):raise ValueError('bothCOMPLETEcommonstoredarmsrequired')
    if gates['full']['pid']==gates['interval']['pid'] or gates['full']['input_identities']!=gates['interval']['input_identities']:raise ValueError('independentprocess/sameinitialinputsrequired')
    for gate in gates.values():
        if json.loads(Path(gate['receipt']).read_text())['status']!='COMPLETE':raise ValueError('closedBookrequired')
        for identity in gate['all_array_seals']:
            if checkpoint_identity(identity['path'])!=identity:raise ValueError('sealed arrayartifact changed')
        for record in gate['latest_array_seals']:
            if checkpoint_identity(record['path'])!=record['identity']:raise ValueError('latestcanonicalendpointarraychanged')
    comparisons=[];misses=0
    for parent in PARENTS:
        with np.load(RESULT/'full'/parent/'CHAIN.npz') as z:f={k:z[k] for k in z.files}
        with np.load(RESULT/'interval'/parent/'CHAIN.npz') as z:i={k:z[k] for k in z.files}
        flag=np.flatnonzero(f['accept']!=i['accept']);path=np.flatnonzero(np.any(f['x_path']!=i['x_path'],axis=1));proof_misses=[]
        candidate_miss=np.flatnonzero((f['candidate_phi']<i['candidate_interval'][:,0])|(f['candidate_phi']>i['candidate_interval'][:,1]))
        current_miss=np.flatnonzero((f['current_phi_before']<i['current_interval'][:,0])|(f['current_phi_before']>i['current_interval'][:,1]))
        # Every intermediate proof, including initial approx and later CURRENT
        # refinements at the SAME old material endpoint, gets post-close audit.
        exact_phi={'initial':float(f['current_phi_before'][0])};exact_phi.update({f'c{k:04d}':float(v) for k,v in enumerate(f['candidate_phi'])})
        for line in (RESULT/'interval'/parent/'ENDPOINT_PROOFS.jsonl').read_text().splitlines():
            record=json.loads(line);lo,hi=record['Phi_interval'];lo=0. if lo is None else lo;hi=float('inf') if hi is None else hi
            if not lo<=exact_phi[record['endpoint_id']]<=hi:proof_misses.append(dict(endpoint_id=record['endpoint_id'],step=record['step'],iteration=record['accuracy_iteration']))
        count=len(flag)+len(path)+len(candidate_miss)+len(current_miss)+len(proof_misses);misses+=count
        comparisons.append(dict(parent=parent,flag_mismatches=flag.tolist(),path_mismatches=path.tolist(),candidate_enclosure_misses=candidate_miss.tolist(),current_enclosure_misses=current_miss.tolist(),all_intermediate_proof_misses=proof_misses))
    process_timing=json.loads((RESULT/'ARM_PROCESS_WALL.json').read_text())
    if any(process_timing[a]['returncode']!=0 for a in ('full','interval')):raise ValueError('complete externalprocess timingrequired')
    wall={a:process_timing[a]['seconds'] for a in ('full','interval')};saving=1-wall['interval']/wall['full'];combined=sum(g['costs']['RHS'] for g in gates.values())
    if combined>200000:raise ValueError('combinedhardRHS cap')
    status='NO_GO_NUMERICAL_OR_PATH_MISMATCH' if misses else 'PASS_FIXED_KERNEL_GE20_PERCENT_WALL_SAVING' if saving>=.2 else 'NO_GO_WHOLE_WALL_SAVING_BELOW20_PERCENT'
    from p3r4.core import write
    result=dict(status=status,misses=misses,parents=comparisons,whole_arm_wall_seconds=wall,whole_wall_relative_saving=saving,combined_RHS=combined,
        codec='stored',cost_scope='actualwholejob inclcoldR0,immutablecopy/setup,proofs,allarraywrites/rewrites,Fullnativeendpoints; no compute-only claim',
        evidence_contract='canonical latestsame-endpoint approximatearray,versionedscalarproofhistory,allproposalattempts+rejectstatexreferences; closedoriginalpairsunaltered',
        float64_interval_certification=False,posterior_reference='NONE',posterior_qualification='UNQUALIFIED')
    if (RESULT/'COMPARISON.json').exists():raise FileExistsError('no auditoverwrite')
    write(RESULT/'COMPARISON.json',result);print(json.dumps(result),flush=True)
    if misses or saving<.2:raise SystemExit(2)


def pair(release):
    marker=OUT/'results/II/interval_chain_fast_pair'
    if marker.exists() or RESULT.exists():raise FileExistsError('no paidreplay')
    check=json.loads(release.read_text())
    if check.get('released') is not True or check.get('experiment')!=EXPERIMENT or check.get('campaign_id')!='P3LR4-20261025':raise ValueError('parentrelease required')
    marker.mkdir(parents=True);children=[];timings={};env=dict(os.environ,R4_ARRAY_CODEC='stored')
    for arm in ('interval','full','compare'):
        process_start=time.perf_counter()
        process=subprocess.run([sys.executable,'-u',str(Path(__file__).resolve()),arm,'--release',str(release.resolve())],cwd=ROOT,env=env)
        timings[arm]=dict(seconds=time.perf_counter()-process_start,returncode=process.returncode)
        children.append(dict(arm=arm,returncode=process.returncode))
        if arm=='full' and process.returncode==0:(RESULT/'ARM_PROCESS_WALL.json').write_text(json.dumps(timings,indent=2)+'\n')
        if process.returncode:
            (marker/'FAILED.json').write_text(json.dumps(dict(status='NO_GO' if arm=='compare' else 'FAILED',children=children),indent=2)+'\n');raise SystemExit(process.returncode)
    (marker/'COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',children=children,comparison=str(RESULT/'COMPARISON.json'),array_codec='stored',scope='implementationoptimization/fixedkernelcostonly'),indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('arm',nargs='?',default='pair',choices=('pair','interval','full','compare'));p.add_argument('--release',type=Path);a=p.parse_args()
    if a.arm=='compare':compare()
    elif a.release is None:p.error('explicitparentrelease required; no autorun')
    elif a.arm=='pair':pair(a.release)
    else:run(a.arm,a.release)

if __name__=='__main__':main()
