from __future__ import annotations
from pathlib import Path
from collections import Counter,OrderedDict
from contextlib import contextmanager
from dataclasses import dataclass
import fcntl,json,os,sys,time,signal
import numpy as np
from scipy.special import ndtr,ndtri
from p1_material.prior import MaterialPrior

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'research/p3l_r4_conditional_transport'
CAPS={'0':0,'I':24000000,'II':6000000,'III':8000000,'IV':20000000}
STAGE_WALL={'0':7200,'I':86400,'II':43200,'III':86400,'IV':115200}
STAGE_GPU={'0':0,'I':172800,'II':57600,'III':28800,'IV':86400}
TOTAL_RHS=58000000
CPU_CAP=94*3600
GPU_CAP=96*3600

def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,Path):return str(x)
    if isinstance(x,complex):return [x.real,x.imag]
    if isinstance(x,dict):return {str(k):plain(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [plain(v) for v in x]
    return x

def write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(plain(obj),ensure_ascii=False,allow_nan=False,indent=2)+'\n');temp.replace(path)

_ARRAY_IO_PROFILE=dict(stored_array_writes=0,stored_array_bytes=0,stored_array_write_seconds=0.)

def array_io_profile():return dict(_ARRAY_IO_PROFILE)

def save(path,**arrays):
    codec=os.environ.get('R4_ARRAY_CODEC','compressed')
    if codec not in ('compressed','stored'):raise ValueError('R4_ARRAY_CODEC must be compressed or stored')
    io_start=time.perf_counter() if codec=='stored' else None
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.tmp')
    with temp.open('wb') as f:
        (np.savez if codec=='stored' else np.savez_compressed)(f,**arrays);f.flush()
    temp.replace(path)
    if io_start is not None:
        _ARRAY_IO_PROFILE['stored_array_writes']+=1
        _ARRAY_IO_PROFILE['stored_array_bytes']+=path.stat().st_size
        _ARRAY_IO_PROFILE['stored_array_write_seconds']+=time.perf_counter()-io_start

def pack(z):return np.concatenate((z.real.T,z.imag.T),axis=1).ravel()

class BudgetStop(RuntimeError):pass

class BoundedTransfers(OrderedDict):
    """Bound native cache lifetime without changing a paid Point's T/J.

    The upstream ordinary dict strongly owns AnchorState and hence its dense
    LU/matrix. Keeping every historical receiver transfer caused campaign
    memory growth. Eviction removes only a cache reference; Point.T/J remain
    valid and any later genuinely missing solve is billed normally.
    """
    def __init__(self,*args,limit=4,**kwargs):
        self.limit=limit;super().__init__(*args,**kwargs)
    def __setitem__(self,key,value):
        super().__setitem__(key,value);self.move_to_end(key)
        while len(self)>self.limit:self.popitem(last=False)

class Book:
    """Atomic campaign reservations. JSONL is authoritative, snapshot recoverable.

    A forward reserves its eight RHS and LU before the unchanged A30 adapter
    starts. The adapter's subsequent full_forward callback consumes that credit,
    so rejected, interrupted and failed attempts cannot disappear or double bill.
    """
    def __init__(self,stage,job,max_wall=7200,*,analysis_only=False,gpu=False,kind=None):
        if stage not in CAPS:raise ValueError(stage)
        if (OUT/'CAMPAIGN_CLOSED.json').exists():raise BudgetStop('R4 closed')
        self.analysis_only=bool(analysis_only);self.gpu=bool(gpu);self.kind=kind or ("analysis" if analysis_only else "physics")
        if stage=="0" and (not analysis_only or gpu):raise BudgetStop("stage0 forbids new physics and GPU")
        # A zero-physics evaluator may read frozen image arrays in its own
        # process. It cannot construct Target or reserve physical actions.
        # Physics/training workers retain the exclusive cache-memory lock.
        self.workerlock=None
        if self.kind in ('physics','training','training_cpu'):
            self.workerlock=(OUT/(self.kind+'.lock')).open('a')
            try:fcntl.flock(self.workerlock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError:
                self.workerlock.close();raise BudgetStop('memory policy: one campaign worker at a time')
        auth=json.loads((OUT/'configs/AUTHORIZATION.json').read_text())
        if auth['campaign_id']!='P3LR4-20261025' or not auth['explicit_user_authorization']:raise BudgetStop('authorization')
        self.stage,self.job=stage,job;self.counts=Counter();self.credit=0;self.context={};self.finished=False
        self.path=OUT/'RESOURCE_LEDGER.jsonl';self.snapshot=OUT/'RESOURCE_LEDGER.json'
        self.jobdir=OUT/stage/'jobs'/job;self.jobdir.mkdir(parents=True,exist_ok=True)
        self.joblock=(self.jobdir/'job.lock').open('a');fcntl.flock(self.joblock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if (self.jobdir/'CLOSED.json').exists() or (self.jobdir/'START.json').exists():raise BudgetStop('existing job: reconcile receipt, never reset or automatically replay')
        self.control=(OUT/'budget.lock').open('a');self.max_wall=float(max_wall);self.start=time.perf_counter()
        self._update('JOB_START',reserve_wall=self.max_wall,analysis_only=self.analysis_only,gpu_job=self.gpu,kind=self.kind)
        write(self.jobdir/'START.json',dict(pid=os.getpid(),job=job,stage=stage,start_unix=time.time(),max_wall=self.max_wall,kind=self.kind,gpu=self.gpu))
        import resource
        self._r0=resource.getrusage(resource.RUSAGE_SELF);self._c0=resource.getrusage(resource.RUSAGE_CHILDREN)
        # Preserve the actual new-stage source even during uncommitted research
        # iteration. Git revision alone cannot identify a working-tree fit.
        import shutil,subprocess
        source_dir=self.jobdir/'source_snapshot';source_dir.mkdir()
        files=[]
        for source in (ROOT/'src/p3r4').glob('*.py'):
            shutil.copyfile(source,source_dir/source.name)
            files.append(dict(path=str(source.relative_to(ROOT)),bytes=source.stat().st_size))
        # Capture entry points as well as the imported physics/neural modules.
        # Older jobs retain their original narrower snapshots honestly.
        tools_dir=source_dir/'tools';tools_dir.mkdir()
        entrypoints=set((ROOT/'tools').glob('*r4*.py'))|{ROOT/'tools/train_operator.py'}
        for source in sorted(entrypoints):
            if source.is_file():
                shutil.copyfile(source,tools_dir/source.name)
                files.append(dict(path=str(source.relative_to(ROOT)),bytes=source.stat().st_size))
        write(self.jobdir/'SOURCE_REVISION.json',dict(
            git_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),array_codec=os.environ.get('R4_ARRAY_CODEC','compressed'),
            argv=sys.argv,source_snapshot=str(source_dir),source_files=files,
            historical_sources='unchanged P1/p3r2 at6eadb469; A30 at717d3d26',checksum_checks=False))
        signal.signal(signal.SIGALRM,lambda a,b:(_ for _ in ()).throw(BudgetStop('job controlled-wall quota')))
        signal.setitimer(signal.ITIMER_REAL,self.max_wall)

    def _latest(self):
        if not self.path.exists():return dict(sequence=0,stages={s:dict(RHS=0,states=0,LU=0,controlled_wall_seconds=0.,physics_wall_seconds=0.,GPU_seconds=0.) for s in CAPS},active_wall={})
        with self.path.open('rb') as f:
            f.seek(0,2);size=f.tell();f.seek(max(0,size-50000));lines=f.read().splitlines()
        return json.loads(lines[-1])['after']

    def _update(self,event,*,rhs=0,states=0,lu=0,controlled=0.,physics=0.,gpu=0.,reserve_wall=0.,release_wall=False,**meta):
        if self.analysis_only and (rhs or states or lu or physics):
            raise BudgetStop('analysis-only receipt forbids physics, LU and RHS')
        fcntl.flock(self.control,fcntl.LOCK_EX)
        try:
            d=self._latest();v=d['stages'][self.stage]
            active_meta=d.setdefault('active_meta',{})
            active_stage=sum(w for k,w in d['active_wall'].items() if active_meta.get(k,{}).get('stage')==self.stage)
            active_gpu=sum(w for k,w in d['active_wall'].items() if active_meta.get(k,{}).get('stage')==self.stage and active_meta.get(k,{}).get('gpu'))
            if reserve_wall and (v['controlled_wall_seconds']+active_stage+reserve_wall>STAGE_WALL[self.stage] or (self.gpu and v['GPU_seconds']+active_gpu+reserve_wall>STAGE_GPU[self.stage])):raise BudgetStop('stage wall/GPU reserve cap')
            if v['controlled_wall_seconds']+controlled>STAGE_WALL[self.stage]+.1 or v['GPU_seconds']+gpu>STAGE_GPU[self.stage]+.1:raise BudgetStop('stage wall/GPU hard cap')
            if v['RHS']+rhs>CAPS[self.stage] or sum(x['RHS'] for x in d['stages'].values())+rhs>TOTAL_RHS:raise BudgetStop('RHS cap')
            used=sum(x['controlled_wall_seconds'] for x in d['stages'].values())
            if reserve_wall and used+sum(d['active_wall'].values())+reserve_wall>CPU_CAP:raise BudgetStop('CPU reservation cap')
            if used+controlled>CPU_CAP or sum(x['GPU_seconds'] for x in d['stages'].values())+gpu>GPU_CAP:raise BudgetStop('CPU/GPU cap')
            if reserve_wall:
                if self.job in d['active_wall']:raise BudgetStop('unreconciled previous job; no automatic replay')
                d['active_wall'][self.job]=reserve_wall
                active_meta[self.job]=dict(stage=self.stage,gpu=self.gpu,pid=os.getpid())
            if release_wall:
                d['active_wall'].pop(self.job,None);active_meta.pop(self.job,None)
            for k,a in [('RHS',rhs),('states',states),('LU',lu),('controlled_wall_seconds',controlled),('physics_wall_seconds',physics),('GPU_seconds',gpu)]:v[k]+=a
            d['sequence']+=1
            row=dict(sequence=d['sequence'],event=event,stage=self.stage,job=self.job,utc_unix=time.time(),RHS=rhs,states=states,LU=lu,controlled_wall_seconds=controlled,physics_wall_seconds=physics,GPU_seconds=gpu,context=self.context,**plain(meta),after=d)
            with self.path.open('a') as f:f.write(json.dumps(row,allow_nan=False,separators=(',',':'))+'\n');f.flush();os.fsync(f.fileno())
            write(self.snapshot,d)
            return d
        finally:fcntl.flock(self.control,fcntl.LOCK_UN)

    def reserve_state(self,material,*,resolution,purpose):
        if resolution!=24:raise ValueError('only original TM24 authorized')
        self._update('FULL_FORWARD_ATTEMPT',states=1,rhs=8,lu=1,purpose=purpose);self.credit+=8

    def reserve_rhs(self,count,role):
        if role=='full_forward' and self.credit>=count:self.credit-=count;return
        self._update('SOLVER_RHS',rhs=int(count),role=role)

    @contextmanager
    def span(self,label,**counts):
        start=time.perf_counter();self.counts.update(counts)
        try:yield
        finally:self.counts[label+'_seconds']+=time.perf_counter()-start

    @contextmanager
    def physical(self,label):
        start=time.perf_counter()
        try:yield
        finally:self._update('PHYSICAL_OPERATION',physics=time.perf_counter()-start,label=label)

    def close(self,status='COMPLETE',error=None):
        if self.finished:return
        signal.setitimer(signal.ITIMER_REAL,0);self.finished=True
        seconds=time.perf_counter()-self.start
        import resource
        r=resource.getrusage(resource.RUSAGE_SELF);c=resource.getrusage(resource.RUSAGE_CHILDREN)
        profile=dict(CPU_seconds=r.ru_utime+r.ru_stime-self._r0.ru_utime-self._r0.ru_stime+c.ru_utime+c.ru_stime-self._c0.ru_utime-self._c0.ru_stime,RSS_peak_bytes=max(r.ru_maxrss,c.ru_maxrss),GPU_allocation_wall_seconds=seconds if self.gpu else 0.)
        d=self._update('JOB_END',controlled=seconds,gpu=seconds if self.gpu else 0.,release_wall=True,profile=profile,status=status,error=error,counts=dict(self.counts),unused_forward_credit=self.credit)
        write(self.jobdir/'CLOSED.json',dict(status=status,error=error,controlled_wall_seconds=seconds,profile=profile,counters=dict(self.counts),stage_cumulative=d['stages'][self.stage]))
        fcntl.flock(self.joblock,fcntl.LOCK_UN);self.joblock.close();self.control.close()
        if self.workerlock is not None:
            fcntl.flock(self.workerlock,fcntl.LOCK_UN);self.workerlock.close()

@dataclass
class Point:
    q:np.ndarray
    state:object
    key:str
    T:np.ndarray|None=None
    J:np.ndarray|None=None
    gradient:np.ndarray|None=None
    gradient_observation:np.ndarray|None=None

class Target:
    def __init__(self,book,*,keep_fields=True):
        if getattr(book,'analysis_only',False):
            raise BudgetStop('analysis-only receipt cannot construct a physics Target')
        self.book=book;self.keep_fields=keep_fields
        self.cfg=json.loads((OUT/'configs/TARGET.json').read_text())
        points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
        self.prior=MaterialPrior(self.cfg,points);self.sigma=self.cfg['noise']['sigma'];self.model=None;self.evaluations=0

    def backend(self):
        if self.model is None:
            p=Path(os.environ.get('B1_P3L_A30_ROOT',str(ROOT.parents[1]/'A30/a30-maxwell-material-manifold')));sys.path.insert(0,str(p/'src'))
            from a30.physics import FrozenMaxwellAdapter
            with self.book.span('paid_background_array_setup'):
                self.model=FrozenMaxwellAdapter.from_background_npz(OUT/'inputs/background/weak.npz',self.book,purpose_prefix='P3LR4_'+self.book.stage)
                self.model._transfers=BoundedTransfers(self.model._transfers,limit=2)
        return self.model

    def evaluate(self,q,purpose='target'):
        chi,_=self.prior.decode(q);self.evaluations+=1;key=f'{self.book.job}_{self.evaluations:07d}'
        folder=self.book.jobdir/'endpoints';self.book.context=dict(self.book.context,endpoint=key,purpose=purpose)
        save(self.book.jobdir/'attempts'/f'{key}.npz',q=q,purpose=np.array(purpose))
        try:
            with self.book.physical('full_forward'):
                state=self.backend().full_state(chi,purpose=purpose)
                if state.residual>1e-10:raise RuntimeError('backward residual contract')
            arrays=dict(q=q,chi=chi,prediction=state.data,backward_residual=state.residual)
            if self.keep_fields:arrays.update(current=state.current,field=state.field)
            save(folder/f'{key}.npz',**arrays)
            return Point(np.asarray(q).copy(),state,key)
        except BaseException as e:
            write(self.book.jobdir/'failures'/f'{key}.json',dict(endpoint=key,error=repr(e),purpose=purpose));raise

    def phi(self,p,y):return float(np.vdot(p.state.data-y,p.state.data-y).real/self.sigma**2)
    def white(self,a):return pack(a)*np.sqrt(2)/self.sigma
    def objective(self,p,y):return self.phi(p,y)+.5*float(p.q[:64]@p.q[:64])

    def jacobian(self,p):
        if p.J is None:
            with self.book.physical('receiver_J_32_RHS'):
                if p.T is None:p.T=self.backend().receiver_transfer(p.state)
                C=np.einsum('qp,ps->qsp',p.T,p.state.field).reshape(256,576)
                a=C@self.prior.derivative(p.q)
                p.J=np.column_stack([self.white(a[:,k].reshape(32,8)) for k in range(66)])
        return p.J

    def phi_gradient(self,p,y):
        y=np.asarray(y)
        if y.shape!=(32,8) or not np.isfinite(y).all():raise ValueError('gradient observation layout')
        if p.gradient is None or p.gradient_observation is None or not np.array_equal(p.gradient_observation,y):
            model=self.backend();res=p.state.data-y
            if p.J is not None:
                p.gradient=p.J.T@self.white(res);self.book.counts['cached_J_gradient_actions']+=1
            else:
                with self.book.physical('source_adjoint_gradient_8_RHS'):
                    a=model.solve_anchor(model.S.conj().T@res,p.state,adjoint=True,role='likelihood_gradient')
                    h=np.sum(np.conj(a)*p.state.field,axis=1)
                    p.gradient=2*np.real(self.prior.derivative(p.q).T@h)/self.sigma**2
            p.gradient_observation=y.copy()
            self.book.counts['gradient_arrays_persisted']+=1
            save(self.book.jobdir/'derivatives'/f'{p.key}_{self.book.counts["gradient_arrays_persisted"]:07d}.npz',q=p.q,gradient_Phi=p.gradient,observation=y)
        return p.gradient

    def persist_cache(self,p,path):
        save(path,q=p.q,chi=p.state.chi,current=p.state.current,field=p.state.field,prediction=p.state.data,LU=p.state.factor[0],pivots=p.state.factor[1],backward_residual=p.state.residual,receiver_T=np.zeros((0,0),complex) if p.T is None else p.T)

    def load_cache(self,path,owner):
        """Paid full-state cache; arithmetic audits add no LU-solver RHS."""
        path=Path(path);model=self.backend()
        with np.load(path,allow_pickle=False) as z:a={k:z[k].copy() for k in z.files}
        from a30.physics import _cpu
        chi=self.prior.decode(a['q'])[0]
        if not np.array_equal(chi,a['chi']):raise ValueError('cache material differs from target decoder')
        L=np.eye(576)-chi[:,None]*model.G;rhs=chi[:,None]*model.incident
        errors=[np.linalg.norm(L@a['current']-rhs)/max(1e-30,np.linalg.norm(rhs)),np.linalg.norm(model.S@a['current']-a['prediction'])/max(1e-30,np.linalg.norm(a['prediction'])),np.linalg.norm(model.incident+model.G@a['current']-a['field'])/max(1e-30,np.linalg.norm(a['field']))]
        if max(errors)>1e-10:raise ValueError('cache Maxwell equations')
        native=_cpu.FullState(chi,a['current'],a['field'],a['prediction'],L,(a['LU'],a['pivots']),float(a['backward_residual']))
        p=Point(a['q'],model._export_state(native,cached=True,path=path),path.stem)
        if a.get('receiver_T',np.zeros((0,0))).size:
            T=a['receiver_T'];error=np.linalg.norm(T@L-model.S)/max(1e-30,np.linalg.norm(model.S))
            if error>1e-10:raise ValueError('cache receiver transfer equation')
            p.T=T;model._transfers[p.state]=T
        self.book.counts['paid_cache_loads']+=1
        self.book._update('CACHE_REUSE',owner=owner,path=str(path),equation_error=max(errors),cold_forward_RHS=8,cold_J_RHS=32 if p.T is not None else 0)
        return p

    def observation(self,parent):
        with np.load(OUT/f'inputs/online/{parent}.npz',allow_pickle=False) as a:
            if set(a.files)!={'observation','target_id','sigma'} or str(a['target_id'])!=self.cfg['target_id'] or float(a['sigma'])!=self.sigma:raise ValueError('observation metadata')
            return a['observation'].copy()

def to_probit(q,*,initialization=False):
    """Exact interior transform; bounded inward move only when explicitly initializing."""
    x=np.array(q,dtype=float,copy=True)
    if x.shape!=(66,) or np.any(x[64:]<0) or np.any(x[64:]>1):raise ValueError('probit input')
    if initialization:x[64:]=np.clip(x[64:],1e-4,1-1e-4)
    elif np.any((x[64:]==0)|(x[64:]==1)):raise ValueError('boundary has infinite probit coordinate')
    x[64:]=ndtri(x[64:]);return x

def from_probit(x):
    q=np.array(x,dtype=float,copy=True);q[64:]=ndtr(q[64:]);return q

def probit_derivative(x):return np.r_[np.ones(64),np.exp(-.5*np.asarray(x)[64:]**2)/np.sqrt(2*np.pi)]
