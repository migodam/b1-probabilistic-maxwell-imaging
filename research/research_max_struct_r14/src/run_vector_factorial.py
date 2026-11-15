"""Fixed qualified vector-DDA factorial and unchanged mature sourcewise controls.

Use --task-id for one child, or --run-all --workers 2 for the frozen 24 tasks.
All CPU caps are absolute whole-process CPU. No truth/scoring is accessed.
"""
import os
for _name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
              'VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
import argparse
from copy import deepcopy
import json
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
import traceback


class BudgetEnd(BaseException):
    pass


def json_value(value):
    if isinstance(value, dict): return {str(k):json_value(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)): return [json_value(v) for v in value]
    if isinstance(value, Path): return str(value)
    if hasattr(value, 'tolist'): return json_value(value.tolist())
    if isinstance(value, complex): return dict(real=value.real,imag=value.imag)
    return value


def dump(path, value):
    Path(path).write_text(json.dumps(json_value(value),indent=2,allow_nan=False)+'\n',encoding='utf8')


def rss_bytes():
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return int(usage.ru_maxrss if sys.platform == 'darwin' else usage.ru_maxrss*1024)


class VectorTask:
    def __init__(self,args,manifest,task):
        self.args,self.manifest,self.task=args,manifest,task
        self.root=args.root.resolve();self.protocol=manifest['protocol']
        self.out=self.root/'RESULTS'/task['campaign']/task['task_id']/('attempt_'+args.attempt)
        self.out.mkdir(parents=True,exist_ok=False)
        self.wall=time.perf_counter();self.deadline=float(self.protocol['proposal_CPU_cap'])
        self.cap=float(self.protocol['CPU_cap']);self.next_audit=2.
        self.phase='initial';self.model=None;self.b=None
        self.events=[];self.history=[];self.pool=[];self.native=[]
        self.best=self.secondary=None;self.phase_best={};self.snapshots={}
        self.saved=set();self.completed=None;self.pending_J=None;self.pending_cycle=0
        self.status='RUNNING';self.termination={};self.error=None
        self.last_corrector_m=None

    def check(self):
        if time.process_time() >= self.deadline: raise BudgetEnd()

    def event(self,label,function,metadata=None):
        before=dict(self.model.cost) if self.model is not None else {}
        cpu,wall=time.process_time(),time.perf_counter();entry=dict(label=label,metadata=metadata)
        try:
            result=function();entry['status']='COMPLETE';return result
        except BaseException as exc:
            entry.update(status='INTERRUPTED' if isinstance(exc,BudgetEnd) else 'FAILED',error=repr(exc));raise
        finally:
            after=dict(self.model.cost) if self.model is not None else {}
            delta={k:after.get(k,0)-before.get(k,0) for k in set(after)|set(before)
                   if isinstance(after.get(k,0),(int,float)) and isinstance(before.get(k,0),(int,float))
                   and after.get(k,0)!=before.get(k,0)}
            entry.update(CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
                         process_CPU_seconds=time.process_time(),operation_delta=delta)
            self.events.append(entry)

    def record(self,m,state,extra=None):
        np=self.np;residual=state['Y']-self.Y
        D=self.model.data_objective(residual,self.Y);R=float(self.model.reg.evaluate(m)[0])
        elapsed=time.process_time();legal=bool(np.isfinite(m).all() and np.all(m>=self.lo-1e-12) and np.all(m<=self.hi+1e-12))
        feasible=bool(np.linalg.norm(residual)<=self.delta)
        meta=dict(candidate=len(self.history),phase=self.phase,CPU_seconds=elapsed,
            wall_seconds=time.perf_counter()-self.wall,D=D,R=R,Phi=D+R,
            relative_residual=float(np.linalg.norm(residual)/self.normY),
            residual_over_delta=float(np.linalg.norm(residual)/self.delta),noise_feasible=feasible,
            noise_qualification='inherited raw CN99 radius; crossgrid mismatch not covered' if 'crossgrid' in self.task['case_id'] else 'inherited matched-discrete known-sigma CN99',
            legal=legal,eligible=bool(legal and elapsed<=self.deadline and elapsed<=self.cap),
            deadline=self.deadline,lower_active=int(np.count_nonzero(m<=self.lo+1e-8)),
            upper_active=int(np.count_nonzero(m>=self.hi-1e-8)),
            physical_constitutive_relative_residual=float(state.get('state_relative_residual',np.linalg.norm(state['j']-np.repeat(m,3)[:,None]*state['Q'])/max(np.linalg.norm(state['j']),1e-300))),
            physical_cost=dict(self.model.cost),extra=extra)
        self.history.append(meta);self.pool.append(np.array(m,copy=True))
        if meta['eligible']:
            better=self.best is None or meta['Phi']<self.best['metadata']['Phi']
            phase_better=self.phase not in self.phase_best or meta['Phi']<self.phase_best[self.phase]['metadata']['Phi']
            secondary_better=feasible and (self.secondary is None or (R,D+R)<(self.secondary['metadata']['R'],self.secondary['metadata']['Phi']))
            candidate=None
            if better or phase_better or secondary_better:
                candidate=dict(m=np.array(m,copy=True),j=state['j'].copy(),Q=state['Q'].copy(),Y_prediction=state['Y'].copy(),metadata=meta)
            if better:self.best=candidate
            if phase_better:self.phase_best[self.phase]=candidate
            if secondary_better:self.secondary=candidate
        return meta

    def physical(self,m,extra=None):
        self.check()
        def action():
            state=self.model.state(m);return self.record(m,state,extra)
        result=self.event('paid original vector physical audit',action,extra)
        self.check();return result

    def audit_due(self):return time.process_time()>=self.next_audit

    def advance_audit_clock(self):
        self.next_audit=(int(time.process_time()/2)+1)*2.

    def full_payload(self,complete):
        np=self.np;b=self.b;J=complete['J'];Q=b.inc+b.G@J
        b.counts['snapshot_auxiliary_G_columns']=b.counts.get('snapshot_auxiliary_G_columns',0)+b.P
        pr=complete['prstate'];cycle=complete['cycle'];mf=complete['mf']
        return dict(checkpoint_version=1,state_family='R14_vector4_existing_VectorPR3D',
            m=complete['m'].copy(),J=J.copy(),Q=Q,A=None if pr is None else pr['gam'].copy(),
            old_grad=None if pr is None else pr['w'].copy(),
            old_direction=None if pr is None else pr['direction'].copy(),prstate=deepcopy(pr),
            stage=0 if cycle<=60 else 1,global_iteration=cycle,iteration=cycle if cycle<=60 else cycle-60,
            mf=mf,stage_boundary_next_reset=bool(cycle==60),
            Vd=b.V.copy(),JD=b.Jd.copy(),basis_mf1=b.bases[1].copy(),basis_mf2=b.bases[2].copy(),
            E2=b.E2.copy(),data_weights=1/np.sum(abs(self.Y)**2,axis=0),
            state_weights=1/b.E2,rho=b.rho,material_box=[self.lo,self.hi],
            components=3,material_unknowns=b.N,current_rows=b.M,material_kind='one shared real chi per voxel',
            Hop_data=b.Hop.data.copy(),Hop_indices=b.Hop.indices.copy(),Hop_indptr=b.Hop.indptr.copy(),Hop_shape=b.Hop.shape,
            prior=dict(lam=b.lam,l2=b.l2,edges=b.edges.copy()),
            operation_counts=dict(b.counts),method=self.task['method'],
            frozen_protocol=self.protocol,public_input=str((self.root/self.task['input_file']).resolve()),
            phase=self.phase,complete_material_readout=True,
            source_version='vector4_v1; source snapshot retained before first run',
            wall_seconds=time.perf_counter()-self.wall,process_CPU_seconds=time.process_time())

    def save_completed(self,name,complete=None):
        complete=self.completed if complete is None else complete
        if complete is None:return
        self.event('full vector state export',lambda:self.save_npz(self.out/(name+'.npz'),self.full_payload(complete)),dict(cycle=complete['cycle'],name=name))

    def capture(self,m,J,cycle,mf,prstate=None):
        if self.completed is not None and cycle<=self.completed['cycle']:return
        pr=self.b.prstate if prstate is None else prstate
        self.completed=dict(m=self.np.array(m,copy=True),J=self.np.array(J,copy=True),
            cycle=int(cycle),mf=int(mf),prstate=deepcopy(pr))
        if cycle in self.protocol['fixed_states'] and cycle not in self.saved:
            self.save_completed(f'CHECKPOINT_{cycle:05d}');self.saved.add(cycle)

    def install_observers(self,pair_mode):
        b=self.b;original_prepare=b.prepare;original_current=b.current_pr;original_material=b.material
        def prepare(Y):
            result=self.event('native vector DPCS and projected-IFFT range preparation',lambda:original_prepare(Y))
            b.counts['projected_basis_range_SVD']=b.counts.get('projected_basis_range_SVD',0)+2
            b.counts['basis_spatial_IFFT_calls']=b.counts.get('basis_spatial_IFFT_calls',0)+72
            b.counts['basis_component_columns_before_range_SVD']=b.counts.get('basis_component_columns_before_range_SVD',0)+216
            self.subspace=result
            self.completed=dict(m=self.m0.copy(),J=b.Jd.copy(),cycle=0,mf=1,prstate=None)
            return result
        def current(z,Y,metric,beta,mf,global_pr):
            if not pair_mode and self.pending_J is not None:
                self.capture(z,self.pending_J,self.pending_cycle,1 if self.pending_cycle<=60 else 2)
            self.pending_cycle+=1
            result=self.event('native vector exact-line PR current update',lambda:original_current(z,Y,metric,beta,mf,global_pr),
                dict(cycle=self.pending_cycle,metric=metric,beta=beta,mf=mf,global_PR=bool(global_pr)))
            if pair_mode:
                columns=b.bases[mf].shape[1]
                for key,value in dict(current_full_matrix_builds=1,current_KB_columns=columns,
                    current_SB_columns=columns,current_KJd_columns=b.P,current_SJd_columns=b.P).items():
                    b.counts[key]=b.counts.get(key,0)+value
            if metric=='A_TSOM':
                b.counts['current_P1_Hop_forward_calls']=b.counts.get('current_P1_Hop_forward_calls',0)+2
                b.counts['current_P1_Hop_complex_RHS']=b.counts.get('current_P1_Hop_complex_RHS',0)+b.M+b.P
                b.counts['current_P1_diagonal_scaled_entries']=b.counts.get('current_P1_diagonal_scaled_entries',0)+b.M*(b.M+b.P)
            self.pending_J=result.copy()
            return result
        def material(z,J,kind,beta):
            result=original_material(z,J,kind,beta)
            if pair_mode:b.counts['material_G_columns']=b.counts.get('material_G_columns',0)+b.P
            if kind=='A_TSOM':
                b.counts['material_P1_Hop_forward_calls']=b.counts.get('material_P1_Hop_forward_calls',0)+1
                b.counts['material_P1_Hop_adjoint_calls']=b.counts.get('material_P1_Hop_adjoint_calls',0)+1
                b.counts['material_P1_Hop_complex_RHS']=b.counts.get('material_P1_Hop_complex_RHS',0)+2*b.P
            return result
        b.prepare,b.current_pr,b.material=prepare,current,material

    def run_pairs(self,pair):
        b=self.b;self.check();z=self.m0.copy();b.prstate=None;b.prepare(self.Y);self.check()
        kinds={'0':'CIE_TSOM','1':'A_TSOM'}
        for index in range(120):
            self.check();mf=1 if index<60 else 2;betas={'0':6 if index<60 else 3,'1':6 if index<60 else 1}
            J=b.current_pr(z,self.Y,kinds[pair[0]],betas[pair[0]],mf,True);self.check()
            def fun(m):
                self.check();result=b.material(m,J,kinds[pair[1]],betas[pair[1]]);self.check();return result
            result=self.event('full analytic vector material L-BFGS-B',lambda:self.minimize(fun,z,jac=True,
                method='L-BFGS-B',bounds=self.model.bounds,
                options=dict(maxiter=12,ftol=1e-12,gtol=1e-8,maxls=25)),dict(cycle=index+1,kind=kinds[pair[1]],beta=betas[pair[1]]))
            z=result.x
            self.native.append(dict(cycle=index+1,mf=mf,beta_current=betas[pair[0]],beta_material=betas[pair[1]],
                nit=int(result.nit),nfev=int(result.nfev),success=bool(result.success),message=str(result.message),
                lower_active=int(self.np.count_nonzero(z<=self.lo+1e-8)),upper_active=int(self.np.count_nonzero(z>=self.hi-1e-8)),
                CPU_seconds=time.process_time()))
            self.capture(z,J,index+1,mf)
            if self.audit_due() or index+1 in (60,120):
                self.physical(z,dict(source='pair completed material update',cycle=index+1,mf=mf));self.advance_audit_clock()
            self.check()
        return dict(status='FIXED_VECTOR_AUXILIARY_CYCLES_COMPLETE',iterations=120,m=z)

    def run_control(self,method):
        from research.research_max_struct_r13.src.vector_adapter import run_vector_auxiliary
        def audit(m,state,value,metadata):
            self.record(m,state,metadata)
            self.capture(m,self.pending_J,int(metadata['iteration']),int(metadata['mf']))
            self.advance_audit_clock();self.check()
        result=self.event('unchanged R13 mature sourcewise auxiliary loop',lambda:run_vector_auxiliary(
            self.model,self.Y,self.m0,method,self.check,audit,audit_every=50,audit_due=self.audit_due,outer=120,inner=12))
        self.capture(result['m'],self.pending_J,int(result['iterations']),2)
        self.native=result.get('optimizer',[])
        return result

    def corrector(self):
        self.deadline=min(self.cap,time.process_time()+float(self.protocol['corrector_CPU_allowance']))
        self.phase='corrector';start=self.best['m'].copy() if self.best is not None else self.m0.copy()
        def objective(m):
            self.check();value,gradient=self.model.data_value_gradient(m,self.Y)
            self.record(m,self.model.state(m),dict(source='common original physical corrector'))
            self.last_corrector_m=self.np.array(m,copy=True);self.check()
            return value*1e6,gradient*1e6
        result=self.event('same R13 original vector physical corrector',lambda:self.minimize(objective,start,jac=True,
            method='L-BFGS-B',bounds=self.model.bounds,
            options=dict(maxiter=3000,maxfun=12000,ftol=1e-12,gtol=.01,maxls=30)))
        return dict(message=str(result.message),iterations=int(result.nit),evaluations=int(result.nfev),success=bool(result.success))

    def execute(self):
        try:
            sys.path[:0]=[str(self.root/'src'),str(self.args.legacy_repo.resolve())]
            import numpy as np
            from scipy.optimize import minimize
            from scipy.stats import chi2
            from dualmetric_core import _save_state_npz
            from research.research_max_struct_r13.src.vector_adapter import VectorPublicModel
            self.np=np;self.minimize=minimize;self.save_npz=_save_state_npz
            def load():
                keys={'G','E0','S','grid_shape','components','width','centers','lower','upper','m_init','Y','sigma','delta','material_bounds','material_objective'}
                with np.load(self.root/self.task['input_file'],allow_pickle=False) as archive:
                    self.public_keys_read=[key for key in archive.files if key in keys]
                    return {key:archive[key].copy() for key in self.public_keys_read}
            self.public=self.event('load public-only vector acquisition and observations',load)
            def setup():self.model=VectorPublicModel(self.public);self.b=self.model.legacy
            self.event('existing full-vector model/operator setup and declared matrix agreement',setup)
            self.Y=self.public['Y'];self.m0=self.public['m_init'].copy();self.normY=float(np.linalg.norm(self.Y))
            self.lo,self.hi=float(self.public['lower']),float(self.public['upper'])
            if (self.lo,self.hi)!=(0.,8.) or self.model.n!=216 or self.Y.shape!=(180,60):
                raise ValueError('Frozen vector4 input dimensions/box changed')
            self.sigma=float(self.public['sigma']);self.delta=float(self.public.get('delta',self.sigma*np.sqrt(chi2.ppf(.99,2*self.Y.size)/2)))
            self.physical(self.m0,dict(source='original initial material'))
            self.snapshots['initial']=self.phase_best.get('initial')
            self.phase='proposal';pair_mode=self.task['method'].startswith('pair')
            self.install_observers(pair_mode)
            dump(self.out/'ACTUAL_CONFIG.json',dict(task=self.task,protocol=self.protocol,model=self.model.metadata,
                auxiliary='full analytic material L-BFGS-B for every pair; no Arc',truth_access='NONE',
                native_mf_transition='reset gam and PR at mf change; no 2D zero padding',
                actual_pair=self.task['method'][4:] if pair_mode else None,
                public_keys_read=self.public_keys_read,source_snapshot='audit/vector4_v1/source_snapshot'))
            try:
                result=self.run_pairs(self.task['method'][4:]) if pair_mode else self.run_control(self.task['method'])
                self.termination['proposal']=result['status']
            except BudgetEnd:self.termination['proposal']='WHOLE_PROCESS_CPU_PROPOSAL'
            # Raw is the paid incumbent including the original initialization.
            # Attempt01's preserved source snapshot used proposal-only here;
            # RAW_SELECTION_AUDIT found all24 saved arrays already identical.
            self.snapshots['raw']=self.best
            self.save_completed('CHECKPOINT_LAST_COMPLETE')
            if self.pending_J is not None:
                self.event('in-flight PR state export',lambda:_save_state_npz(self.out/'PROPOSAL_IN_FLIGHT.npz',dict(
                    cycle=self.pending_cycle,prstate=deepcopy(self.b.prstate),J=self.pending_J.copy(),
                    completed_cycle=None if self.completed is None else self.completed['cycle'],
                    note='PR may belong to an in-flight readout if proposal interrupted; not mislabeled complete material state')))
            try:self.termination['corrector']=self.corrector()
            except BudgetEnd:self.termination['corrector']='WHOLE_PROCESS_CPU'
            self.snapshots['corrected']=self.phase_best.get('corrector')
            self.status='COMPLETE'
        except BudgetEnd:
            self.status='INCOMPLETE_SETUP';self.termination['setup']='WHOLE_PROCESS_CPU_PROPOSAL'
        except Exception as exc:
            self.status='FAILED';self.error=repr(exc);traceback.print_exc()
            (self.out/'FAILURE.txt').write_text(traceback.format_exc())
        finally:
            self.finish()
        return 0 if self.status=='COMPLETE' else 1

    def finish(self):
        if self.model is not None and hasattr(self,'np'):
            np=self.np
            self.snapshots['primary']=self.best
            secondary_fallback=self.secondary is None
            self.snapshots['secondary']=self.secondary if self.secondary is not None else self.best
            for name,value in self.snapshots.items():
                if value is None:continue
                self.event('selected physical candidate export',lambda:np.savez_compressed(self.out/(name+'.npz'),
                    **{key:array for key,array in value.items() if key!='metadata'}),dict(view=name))
                dump(self.out/(name+'.json'),dict(value['metadata'],secondary_fallback_to_primary=secondary_fallback if name=='secondary' else False))
            if self.pool:self.event('full paid material candidate pool export',lambda:np.savez_compressed(self.out/'CANDIDATE_POOL.npz',m=np.asarray(self.pool)))
            if self.completed is not None and 'CHECKPOINT_LAST_COMPLETE.npz' not in {p.name for p in self.out.glob('*.npz')}:
                try:self.save_completed('CHECKPOINT_LAST_COMPLETE')
                except Exception as exc:self.error=(self.error or '')+'; final checkpoint: '+repr(exc)
            counts=dict(self.model.cost)
        else:counts={};secondary_fallback=True
        dump(self.out/'HISTORY.json',self.history);dump(self.out/'NATIVE_HISTORY.json',self.native)
        dump(self.out/'COST_LEDGER.json',self.events)
        receipt=dict(status=self.status,error=self.error,task=self.task,protocol=self.protocol,termination=self.termination,
            whole_process_CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-self.wall,
            operation_counts_all=counts,physical_cost={key:value for key,value in counts.items() if key.startswith('physical') or key in ('source_RHS','adjoint_RHS')},
            actual_peak_RSS_bytes=rss_bytes(),peak_memory_scope='OS getrusage ru_maxrss high-water for this complete child process through final exports',
            completed_material_cycles=None if self.completed is None else self.completed['cycle'],
            saved_full_state_cycles=sorted(self.saved),candidates=len(self.history),
            eligible_candidates=sum(int(row['eligible']) for row in self.history),
            out_of_budget_candidates=sum(int(not row['eligible']) for row in self.history),
            secondary_fallback_to_primary=secondary_fallback,
            truth_access='NONE',scoring_performed=False,
            physical_qualification='inherited finite CM+RR vector-DDA; no new Mie/PDE/continuum claim',
            material_unknowns=216,current_components_per_source=648,
            source_provenance=self.manifest['provenance'],
            timer_scope='absolute whole-process CPU includes imports/setup/proposal/corrector/checkpoint/selected exports/failed attempts; controller records OS complete child-lifetime aggregate',
            note='native control material summaries may be unavailable if interrupted inside unchanged R13 loop; full last-completed and in-flight PR states retained')
        dump(self.out/'RECEIPT.json',receipt)
        print(json.dumps(dict(task_id=self.task['task_id'],status=self.status,CPU_seconds=receipt['whole_process_CPU_seconds'],out=str(self.out))))


def freeze_sources(args,manifest):
    root=args.root.resolve();repo=args.legacy_repo.resolve();out=root/'audit/vector4_v1/source_snapshot'/('attempt_'+args.attempt)
    out.mkdir(parents=True,exist_ok=False)
    paths=[root/'src/run_vector_factorial.py',root/'src/dualmetric_core.py',
           root/'configs/vector4_v1.json',root/'audit/VECTOR_CONTRACT.md']
    paths += [repo/path for path in (
        'research/research_max_struct_r13/src/vector_adapter.py',
        'research/research_max_struct_r13/src/physics_inputs.py',
        'research/research_max_struct_r8/src/vector_pr3d.py',
        'research/research_max_struct_r8/src/vector_cie3d.py',
        'research/research_max_struct_r1/src/vector_backend.py',
        'research/research_max_struct_r1/src/a9_green.py')]
    entries=[]
    for path in paths:
        rel=path.relative_to(repo);dest=out/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
        entries.append(dict(source=str(path),snapshot=str(dest),bytes=dest.stat().st_size,mtime_ns=path.stat().st_mtime_ns))
    dump(out/'PROVENANCE.json',dict(status='FROZEN_BEFORE_CHILD_EXECUTION',files=entries,
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
        source_change_status=subprocess.check_output(['git','status','--short','--',*[str(path.relative_to(repo)) for path in paths]],cwd=repo,text=True),
        truth_access='NONE',provider='unchanged project default; no provider launches'))
    return out


def run_all(args,manifest):
    from concurrent.futures import ThreadPoolExecutor,as_completed
    controller_cpu=time.process_time();controller_wall=time.perf_counter()
    all_children_before=resource.getrusage(resource.RUSAGE_CHILDREN)
    root=args.root.resolve();queue=root/'queues/vector4_v1'/('attempt_'+args.attempt)
    queue.mkdir(parents=True,exist_ok=False);source_snapshot=freeze_sources(args,manifest)
    before=resource.getrusage(resource.RUSAGE_CHILDREN);started=time.perf_counter();jobs=[]
    dump(queue/'QUEUE_STARTED.json',dict(tasks=len(manifest['tasks']),workers=args.workers,
        interpreter=sys.executable,host=sys.platform,source_snapshot=str(source_snapshot),truth_access='NONE'))
    def child(task):
        command=[sys.executable,str(root/'src/run_vector_factorial.py'),'--root',str(root),
            '--legacy-repo',str(args.legacy_repo.resolve()),'--manifest',str(args.manifest.resolve()),
            '--task-id',task['task_id'],'--attempt',args.attempt]
        start=time.perf_counter();log=queue/(task['task_id']+'.log')
        with log.open('w') as stream:
            result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,check=False)
        path=root/'RESULTS'/task['campaign']/task['task_id']/('attempt_'+args.attempt)/'RECEIPT.json'
        receipt=json.loads(path.read_text()) if path.exists() else None
        return dict(task_id=task['task_id'],returncode=result.returncode,wall_seconds=time.perf_counter()-start,
            receipt=str(path),child_status=None if receipt is None else receipt['status'],
            child_reported_CPU_seconds=None if receipt is None else receipt['whole_process_CPU_seconds'],
            log=str(log))
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures={executor.submit(child,task):task for task in manifest['tasks']}
        for future in as_completed(futures):
            task=futures[future]
            try:row=future.result()
            except Exception as exc:row=dict(task_id=task['task_id'],returncode=None,child_status='CONTROLLER_FAILURE',error=repr(exc))
            jobs.append(row);dump(queue/'PROGRESS.json',dict(terminal=len(jobs),expected=len(manifest['tasks']),jobs=jobs))
            print(json.dumps(row),flush=True)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    receipt=dict(status='ALL_EXPECTED_CHILDREN_TERMINAL',expected=len(manifest['tasks']),terminal=len(jobs),jobs=jobs,
        successful=sum(row.get('returncode')==0 and row.get('child_status')=='COMPLETE' for row in jobs),
        failed_or_missing=sum(not(row.get('returncode')==0 and row.get('child_status')=='COMPLETE') for row in jobs),
        workers=args.workers,BLAS_threads=1,host=sys.platform,interpreter=sys.executable,
        whole_children_OS_CPU_seconds=(after.ru_utime+after.ru_stime)-(before.ru_utime+before.ru_stime),
        all_controller_subprocess_OS_CPU_seconds=(after.ru_utime+after.ru_stime)-(all_children_before.ru_utime+all_children_before.ru_stime),
        controller_CPU_seconds=time.process_time()-controller_cpu,
        controller_wall_including_source_freeze_seconds=time.perf_counter()-controller_wall,
        wall_seconds=time.perf_counter()-started,source_snapshot=str(source_snapshot),truth_access='NONE',scoring_performed=False)
    dump(queue/'RECEIPT.json',receipt)
    inventory=[]
    for task in manifest['tasks']:
        folder=root/'RESULTS'/task['campaign']/task['task_id']/('attempt_'+args.attempt)
        files=[dict(path=str(path.relative_to(root)),bytes=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns)
               for path in sorted(folder.rglob('*')) if path.is_file()] if folder.exists() else []
        inventory.append(dict(task_id=task['task_id'],files=files,missing=not bool(files)))
    seal=root/'predictions/vector4_v1'/('attempt_'+args.attempt);seal.mkdir(parents=True,exist_ok=False)
    dump(seal/'PREDICTION_SEAL.json',dict(status='SEALED_TERMINAL_PREDICTIONS_WITH_ALL_FAILURES_RETAINED',
        expected_tasks=len(manifest['tasks']),terminal_tasks=len(jobs),queue_receipt=str(queue/'RECEIPT.json'),
        successful=receipt['successful'],failed_or_missing=receipt['failed_or_missing'],inventory=inventory,
        truth_access_before_seal='NONE',GT_scoring_authorized_by_this_seal=False,scoring_performed=False,
        source_snapshot=str(source_snapshot),material_scope='4 exposed conditions, 3 independent geometries',
        scientific_acceptance='requires parent review; sealing records terminal artifacts only'))
    print(json.dumps(dict(queue_receipt=str(queue/'RECEIPT.json'),prediction_seal=str(seal/'PREDICTION_SEAL.json'),successful=receipt['successful'],failed=receipt['failed_or_missing'])))
    return 0 if receipt['failed_or_missing']==0 else 1


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True);parser.add_argument('--legacy-repo',type=Path,required=True)
    parser.add_argument('--manifest',type=Path,required=True)
    group=parser.add_mutually_exclusive_group(required=True);group.add_argument('--task-id');group.add_argument('--run-all',action='store_true')
    parser.add_argument('--workers',type=int,choices=(1,2),default=2);parser.add_argument('--attempt',default='01')
    args=parser.parse_args();manifest=json.loads(args.manifest.read_text())
    if args.run_all:return run_all(args,manifest)
    task=next(task for task in manifest['tasks'] if task['task_id']==args.task_id)
    return VectorTask(args,manifest,task).execute()


if __name__=='__main__':sys.exit(main())
