"""Exact-ROI possibility-cloud engine. No solve without root protocol+enable.
Online functions accept observation-only files; generation and truth scoring are
separate executables. Material support and likelihood are original P1.
"""
from __future__ import annotations
import os
for _n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[_n]='1'
import argparse, copy, fcntl, hashlib, json, resource, sys, time, uuid, warnings
from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy.linalg import lu_factor, lu_solve, null_space, eigh
from scipy.optimize import least_squares, brentq
from scipy.special import ndtr
from scipy.stats import chi2
P=Path(__file__).resolve().parents[1]
B1=P.parents[2]
R4=B1/'b1-p3l-r4-conditional-transport'
OUT=P/'results/possibility_volume'
CONFIG=P/'config/possibility_volume_protocol.json'
ENABLE=P/'config/possibility_volume_enable.json'
SETUP=P/'results/stable_roi/SETUP.npz'
STATS_CONFIG=P/'config/possibility_statistics_frozen.json'
SIGMA=.00119211798514519
NS=[20261030,900901]
ARMS=['FIXED_Q8_RAY','DYNAMIC_Q8_RAY','RANDOM66_FULL']
sys.path.insert(0,str(R4/'src'))
from p1_material.prior import MaterialPrior

class BudgetStop(RuntimeError):pass

def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,Path):return str(x)
    if isinstance(x,dict):return {str(k):plain(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [plain(v) for v in x]
    return x

def write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.'+uuid.uuid4().hex+'.tmp')
    with tmp.open('w') as f:
        json.dump(plain(obj),f,ensure_ascii=False,allow_nan=False,indent=2)
        f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)

def save(path,**arrays):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.'+uuid.uuid4().hex+'.tmp')
    with tmp.open('wb') as f:np.savez_compressed(f,**arrays);f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)

def append(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('a') as f:
        f.write(json.dumps(plain(obj),allow_nan=False,separators=(',',':'))+'\n')
        f.flush();os.fsync(f.fileno())

def load_npz(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def pack(z):return np.concatenate((z.real.T,z.imag.T),axis=1).ravel()
def white(z):return pack(z)*np.sqrt(2)/SIGMA
def rmat(z):return np.concatenate((np.asarray(z).real,np.asarray(z).imag),axis=0)/24
def q_from_x(x):
    q=np.asarray(x,float).copy();q[64:]=ndtr(q[64:]);return q
def decode(prior,x):return prior.decode(q_from_x(x))[0]
def derivative(prior,x):
    return prior.derivative(q_from_x(x))*np.r_[np.ones(64),np.exp(-.5*np.asarray(x)[64:]**2)/np.sqrt(2*np.pi)][None,:]

def file_identity(p):
    p=Path(p);st=p.stat()
    return dict(path=str(p.resolve()),bytes=st.st_size,mtime_ns=st.st_mtime_ns,blake2b=hashlib.blake2b(p.read_bytes()).hexdigest())

def sources():
    paths=[Path(__file__),Path(__file__).with_name('possibility_volume_generate.py'),
           R4/'src/p1_material/prior.py',R4/'research/p3l_r4_conditional_transport/configs/TARGET.json',
           P/'src/stable_roi_operator.py',SETUP]
    if STATS_CONFIG.exists():paths.append(STATS_CONFIG)
    return [file_identity(p) for p in paths]

def draft():
    rows=[]
    for i in range(32):
        rows.append(dict(parent=f'R9V_FRESH_{i:03d}',index=i,
                         prior_seed=NS+[101,i],noise_seed=NS+[102,i],
                         weak_seed_prefix=NS+[103,i],random_seed_prefix=NS+[104,i]))
    c=dict(id='R9_POSSIBILITY_VOLUME_V1',statistics_config_id='R9_POSSIBILITY_VOLUME_STATISTICS_V1',approved_by_root=False,created_utc=time.time(),
       output=str(OUT),private_output=str(P/'private/possibility_volume'),seed_namespace='20261030_R9_POSSIBILITY_VOLUME_V1',
       seed_scope='Unique declared SeedSequence entropy tuples by role and parent; no claim of mathematically disjoint RNG value ranges',
       parents=rows,arms=ARMS,K_total_per_arm=16,radii=[.02,.05,.10,.15],
       pair_count_per_radius=2,antithetic_signs=[1,-1],
       proposal_order='pair0 all radii +/- then pair1 all radii +/-',
       independent_radius_starts=True,weak_rank=8,weak_prior_weight=1.,
       model_prior_weight=1.,anchor_prior_weight=1.,latent_radius=4.,
       anchor=dict(cap_states=64,ftol=1e-9,xtol=1e-9,gtol=1e-7,method='trf',
                   tr_solver='exact',x_scale='jac',start='zeros66',
                   reference='R7 run_e2.direct_anchor; reserve final state slot for best LU restoration'),
       candidate=dict(max_forward_trials=6,max_J_states=4,trust_initial=.5,
                      material_step_fraction=.5,max_decoder_backtracks=9,
                      finite_progress_fraction=.5,radius_relative_tolerance=.02),
       primary=dict(noise_ball=float(chi2.ppf(.95,512)/2),legal=True,
                    displacement_radius_required=True,finite_progress_required=True,
                    selection='minimum true dataPhi among visited legal radius/progress states, no truth'),
       legacy_secondary=dict(DeltaPhi_max=2.,DeltaPrior_max=10.),
       model_caveat='GN model retains MAP prior term1; primary qualification and endpoint selection use absolute true dataPhi',
       dynamic_ray='current weak projection Qwk Qwk.T (x-x0); preserve this ray only within current frame step; no cross-frame transport',
       random_control='own independent paired nativeGaussian66 directions; ordinary full effective66 signed-row corrector is a stronger geometry class',
       always_shared_anchor=True,max_unique_per_primary_cloud=17,
       caps=dict(reduced_RHS256=400000,reduced_LU256=12000,full_RHS576=0,full_LU576=0,
                 controlled_wall_seconds=7200,peak_memory_bytes=8000000000),
       resources=dict(single_worker=True,BLAS_threads=1,live_LU=1,dtype='float64/complex128'),
       shared_setup=dict(path=str(SETUP),historical_study_full_RHS576=2184,
                         already_paid_setup_full_RHS576=264,new_setup_charge=0,
                         identity='exact fixed-host nonuniform ROI all P1 support; cache hit, not a new reduced-model approximation'),
       full576_audits='NOT_RUN; none authorized by this protocol',
       fallback='ABORT_RECORD_FAILURE; no unbudgeted full576 fallback',
       scoring='separate executable after all32 prediction seals; posterior accuracy not a primary requirement',
       weak_line_sign_policy='Unconstrained signed amplitude; historical *_RAY ids mean weak line. Signed material progress remains required.',
       primary_cloud_storage='allselected_with_flags; not qualified-only and not posterior samples',
       statistics_config_identity=file_identity(STATS_CONFIG) if STATS_CONFIG.exists() else dict(path=str(STATS_CONFIG),status='WAITING_ROOT_FROZEN_CONFIG'),
       source_files=sources())
    return c

def estimate(c):
    n=len(c['parents']);k=c['K_total_per_arm'];a=len(c['arms'])
    f=c['candidate']['max_forward_trials'];j=c['candidate']['max_J_states'];ac=c['anchor']['cap_states']
    candidates=n*k*a
    return dict(parents=n,proposals=candidates,per_proposal_reduced_RHS256=8*f+32*j,
                candidate_reduced_RHS256=candidates*(8*f+32*j),
                anchor_reduced_RHS256=n*ac*40,generation_reduced_RHS256=n*8,
                worst_total_reduced_RHS256=candidates*(8*f+32*j)+n*ac*40+n*8,
                candidate_LU256=candidates*f,anchor_LU256=n*ac,generation_LU256=n,
                worst_total_LU256=candidates*f+n*ac+n,
                new_full_RHS576=0,setup_full_RHS576_recharged=0,
                historical_setup_already_charged=264,historical_setup_study_total=2184,
                note='Reservations include rejected/failed work; exact ROI RHS256 are not equated to full576 RHS')

def read_protocol(path=CONFIG,require_approved=True):
    c=json.loads(Path(path).read_text())
    if require_approved and c.get('approved_by_root') is not True:raise RuntimeError('root protocol approval missing')
    if len(c['parents'])!=32 or c['K_total_per_arm']!=16 or c['weak_rank']!=8:raise RuntimeError('protocol dimension mismatch')
    if c['weak_prior_weight']!=1 or c['model_prior_weight']!=1:raise RuntimeError('no unapproved weak/model policy refit')
    if c['candidate']['max_forward_trials']!=6 or c['candidate']['max_J_states']!=4:raise RuntimeError('trial contract')
    if c['caps']['full_RHS576'] or c['caps']['full_LU576']:raise RuntimeError('no full576 work authorized')
    return c

class Book:
    """One global exclusive worker; append-only reservation and actual-attempt accounting."""
    COUNTS=('reduced_RHS256','reduced_LU256','full_RHS576','full_LU576',
            'actual_forward_RHS256','actual_J_RHS256','actual_LU256',
            'completed_forward_RHS256','completed_J_RHS256','failed_actions',
            'states','array_save_attempts','array_save_failures','array_saves','array_bytes','controlled_wall_seconds','CPU_seconds','physics_wall_seconds','save_wall_seconds')
    def __init__(self,c,job,role):
        self.started=time.perf_counter();self.cpu0=time.process_time()
        self.c=c;self.job=job;self.role=role;self.finished=False;self.context={}
        if not ENABLE.exists():raise RuntimeError('physics not enabled')
        en=json.loads(ENABLE.read_text())
        if en.get('enabled') is not True or en.get('protocol_id')!=c['id']:raise RuntimeError('explicit root enable missing/mismatched')
        required={str(SETUP.resolve()),str(STATS_CONFIG.resolve())}
        if not required.issubset({v['path'] for v in c['source_files']}):raise RuntimeError('setup and frozen statistics identities required before physics')
        for src in c['source_files']:
            if hashlib.blake2b(Path(src['path']).read_bytes()).hexdigest()!=src['blake2b']:raise RuntimeError('frozen source changed: '+src['path'])
        self.out=Path(c['output']);self.out.mkdir(parents=True,exist_ok=True)
        self.lock=(P/'shared/physics_worker.lock').open('a');fcntl.flock(self.lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        self.path=self.out/'RESOURCE_LEDGER.jsonl';self.snap=self.out/'RESOURCE_SNAPSHOT.json'
        self.jobdir=self.out/'jobs'/job
        if self.jobdir.exists():self.lock.close();raise RuntimeError('no implicit job replay')
        self.total=self.replay()
        if self.total['active_jobs']:self.lock.close();raise RuntimeError('unreconciled interrupted job; preserve reservations')
        self.counts=Counter()
        self.jobdir.mkdir(parents=True)
        self.event('JOB_BEGIN',pid=os.getpid())
        self.event('EXACT_ROI_CACHE_REUSE',path=c['shared_setup']['path'],
                   already_paid_setup_full_RHS576=264,historical_study_full_RHS576=2184,new_charge=0)
    def replay(self):
        t=dict(sequence=0,active_jobs={},**{k:0 for k in self.COUNTS})
        if not self.path.exists():return t
        raw=self.path.read_bytes()
        if not raw.endswith(b'\n'):raise RuntimeError('incomplete ledger; no silent repair')
        seq=0
        for line in raw.splitlines():
            row=json.loads(line);seq+=1
            if row['after']['sequence']!=seq:raise RuntimeError('ledger sequence')
            t=row['after']
        if self.snap.exists() and json.loads(self.snap.read_text())!=t:raise RuntimeError('ledger/snapshot mismatch')
        return t
    def event(self,event,**kw):
        t=copy.deepcopy(self.total)
        for k in self.COUNTS:
            if k in kw:t[k]+=kw[k]
        for k in ('reduced_RHS256','reduced_LU256','full_RHS576','full_LU576'):
            if t[k]>self.c['caps'][k]:raise BudgetStop(k+' hard cap')
        if event=='JOB_BEGIN':t['active_jobs'][self.job]=dict(pid=kw['pid'],role=self.role)
        elif event=='JOB_END':t['active_jobs'].pop(self.job)
        t['sequence']+=1
        append(self.path,dict(event=event,utc=time.time(),job=self.job,role=self.role,context=self.context,**kw,after=t))
        write(self.snap,t);self.total=t
    def check(self):
        elapsed=time.perf_counter()-self.started
        if self.total['controlled_wall_seconds']+elapsed>self.c['caps']['controlled_wall_seconds']:raise BudgetStop('controlled wall cap')
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
        if rss>self.c['caps']['peak_memory_bytes']:raise BudgetStop('memory cap')
    def reserve_forward(self):
        self.check();self.event('FORWARD_RESERVED',reduced_RHS256=8,reduced_LU256=1,states=1)
    def reserve_J(self):
        self.check();self.event('J_RESERVED',reduced_RHS256=32)
    def save(self,path,**data):
        st=time.perf_counter();self.event('ARRAY_SAVE_ATTEMPT',array_save_attempts=1)
        try:save(path,**data)
        except BaseException as exc:
            self.event('ARRAY_SAVE_FAILURE',array_save_failures=1,save_wall_seconds=time.perf_counter()-st,error=repr(exc));raise
        self.event('ARRAY_SAVE',array_saves=1,array_bytes=Path(path).stat().st_size,save_wall_seconds=time.perf_counter()-st)
    def close(self,error=None):
        if self.finished:return
        self.event('JOB_END',status='FAILED' if error else 'COMPLETE',error=None if error is None else repr(error),
                   controlled_wall_seconds=time.perf_counter()-self.started,CPU_seconds=time.process_time()-self.cpu0,
                   peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
        self.finished=True;self.lock.close()
    def __enter__(self):return self
    def __exit__(self,ty,err,tb):self.close(err)

@dataclass
class Point:
    x:np.ndarray
    q:np.ndarray
    material:np.ndarray
    prediction:np.ndarray
    reduced_current:np.ndarray
    ROI_field:np.ndarray
    key:str
    backward_residual:float
    J:np.ndarray|None=None

class ExactROI:
    """One 256LU owner. Full66 material Jacobian from exactly32 reduced adjoint RHS."""
    def __init__(self,book):
        self.book=book;self.a=load_npz(book.c['shared_setup']['path'])
        cfg=json.loads((R4/'research/p3l_r4_conditional_transport/configs/TARGET.json').read_text())
        with np.load(R4/'research/p3l_r4_conditional_transport/inputs/background/weak.npz') as z:points=z['points'].copy()
        self.prior=MaterialPrior(cfg,points);self.idx=self.a['ROI'].astype(int)
        if not np.array_equal(self.idx,np.where(self.prior.mask)[0]) or len(self.idx)!=256:raise RuntimeError('exact ROI support changed')
        if cfg['noise']['sigma']!=SIGMA or cfg['gain']!=1:raise RuntimeError('frozen target changed')
        for k,shape in [('C',(256,256)),('ER',(256,8)),('O',(32,256)),('Fh',(32,8)),('B',(576,256)),('Eh',(576,8))]:
            if self.a[k].shape!=shape:raise RuntimeError('setup shape '+k)
        self.live=None;self.factor=None;self.serial=0
        from threadpoolctl import threadpool_limits
        self.threads=threadpool_limits(limits=1,user_api='blas')
    def drop_lu(self):self.factor=None;self.live=None
    def evaluate(self,x,purpose):
        x=np.asarray(x,dtype=np.float64)
        if x.shape!=(66,) or not np.isfinite(x).all():raise ValueError('finite x66')
        q=q_from_x(x);material=self.prior.decode(q)[0];self.serial+=1
        key=f'{self.book.job}_{self.serial:07d}';self.book.context.update(endpoint=key,purpose=purpose)
        self.book.save(self.book.jobdir/'attempts'/f'{key}.npz',x=x,q=q,material=material)
        self.drop_lu();self.book.reserve_forward();st=time.perf_counter();physics_ended=False
        try:
            delta=material[self.idx]-self.prior.host
            L=np.eye(256,dtype=np.complex128)-delta[:,None]*self.a['C']
            self.book.event('LU_ATTEMPT',actual_LU256=1)
            with warnings.catch_warnings():
                warnings.simplefilter('error');self.factor=lu_factor(L,check_finite=False)
            rhs=delta[:,None]*self.a['ER']
            self.book.event('FORWARD_SOLVE_ATTEMPT',actual_forward_RHS256=8)
            j=lu_solve(self.factor,rhs,check_finite=False)
            self.book.event('FORWARD_SOLVE_COMPLETED',completed_forward_RHS256=8)
            field=self.a['ER']+self.a['C']@j;prediction=self.a['Fh']+self.a['O']@j
            residual=float(np.linalg.norm(L@j-rhs)/max(np.linalg.norm(rhs),1e-300))
            if not np.isfinite(residual) or residual>1e-10:raise RuntimeError('reduced physical residual')
            p=Point(x.copy(),q.copy(),material,prediction,j,field,key,residual);self.live=p
            self.book.event('FORWARD_COMPLETED',physics_wall_seconds=time.perf_counter()-st,backward_residual=residual);physics_ended=True
            self.persist(p);return p
        except BaseException as exc:
            self.book.event('PHYSICS_FAILURE' if not physics_ended else 'PERSIST_AFTER_PHYSICS_FAILURE',failed_actions=1,physics_wall_seconds=0. if physics_ended else time.perf_counter()-st,error=repr(exc))
            write(self.book.jobdir/'failures'/f'{key}.json',dict(error=repr(exc),endpoint=key,purpose=purpose))
            self.drop_lu();raise
    def jacobian(self,p):
        if p.J is not None:return p.J
        if p is not self.live or self.factor is None:raise RuntimeError('released LU: explicit paid forward restoration required')
        self.book.reserve_J();st=time.perf_counter();physics_ended=False
        try:
            self.book.event('J_SOLVE_ATTEMPT',actual_J_RHS256=32)
            adj=lu_solve(self.factor,self.a['O'].conj().T,trans=2,check_finite=False)
            transfer=adj.conj().T
            self.book.event('J_SOLVE_COMPLETED',completed_J_RHS256=32)
            D=derivative(self.prior,p.x)[self.idx]
            dc=np.einsum('ri,is,ik->rsk',transfer,p.ROI_field,D,optimize=True)
            p.J=np.column_stack([white(dc[:,:,k]) for k in range(66)])
            if p.J.shape!=(512,66) or not np.isfinite(p.J).all():raise RuntimeError('finite full66 reduced J')
            self.book.event('J_COMPLETED',physics_wall_seconds=time.perf_counter()-st);physics_ended=True
            self.persist(p);return p.J
        except BaseException as exc:
            self.book.event('PHYSICS_FAILURE' if not physics_ended else 'PERSIST_AFTER_PHYSICS_FAILURE',failed_actions=1,physics_wall_seconds=0. if physics_ended else time.perf_counter()-st,error=repr(exc));raise
    def persist(self,p,path=None,full_fields=False):
        d=dict(x=p.x,q=p.q,material=p.material,prediction=p.prediction,
               reduced_current=p.reduced_current,ROI_field=p.ROI_field,
               backward_residual=p.backward_residual,exact_ROI=np.array(True))
        if p.J is not None:d['J']=p.J
        if full_fields:
            field=self.a['Eh']+self.a['B']@p.reduced_current
            d.update(field=field,current=p.material[:,None]*field)
        self.book.save(path or self.book.jobdir/'states'/f'{p.key}.npz',**d)
    def close(self):self.drop_lu()
    def phi(self,p,y):return float(np.vdot(p.prediction-y,p.prediction-y).real/SIGMA**2)
    def objective(self,p,y,weight=1.):return self.phi(p,y)+.5*weight*float(p.x@p.x)

def frame(prior,x,J,k=8):
    M=rmat(derivative(prior,x));U,s,Vh=np.linalg.svd(M,full_matrices=False)
    keep=s>max(s[0]*1e-8,1e-10)
    if not keep.any():
        return dict(M=M,T=np.zeros((66,0)),K=np.zeros((512,0)),Qw=np.zeros((66,0)),weaklift=np.zeros((66,0)),eigenvalues=np.array([]),material_singular=s,rank=0)
    T=Vh[keep].T/s[keep][None,:];K=J@T
    ev,W=eigh(K.T@K+T.T@T);weaklift=T@W[:,:min(k,len(ev))];Qw,_=np.linalg.qr(weaklift)
    return dict(M=M,T=T,K=K,Qw=Qw,weaklift=weaklift,eigenvalues=ev,material_singular=s,rank=int(keep.sum()))

def match_radius(prior,x0,v,radius,max_latent=4.):
    m0=decode(prior,x0)
    def f(h):return np.linalg.norm(decode(prior,x0+h*v)-m0)/24-radius
    lo=0.;hi=max(radius/8,1e-8)
    for _ in range(100):
        if np.linalg.norm(hi*v)>max_latent:
            hi=max_latent/max(np.linalg.norm(v),1e-300)
            if f(hi)<0:return None,dict(status='RADIUS_BRACKET_NOT_FOUND_WITHIN_LATENT_CAP',max_RMS=f(hi)+radius)
        if f(hi)>=0:
            h=brentq(f,lo,hi,xtol=1e-14,rtol=1e-12)
            return x0+h*v,dict(status='MATCHED',h=h,latent_step=float(np.linalg.norm(h*v)))
        lo=hi;hi*=1.35
    return None,dict(status='DISTANCE_BRACKET_NOT_FOUND')

def metrics(ph,p,base,y,unit,radius,c):
    dm=rmat(p.material-base.material);rms=float(np.linalg.norm(dm));phi=ph.phi(p,y)
    geometric=abs(rms-radius)<=.02*radius and float(unit@dm)>=.5*radius
    legal=ph.prior.legal(p.q);noise=phi<=c['primary']['noise_ball']
    dp=.5*float(p.x@p.x-base.x@base.x);dphi=phi-ph.phi(base,y)
    return dict(endpoint=p.key,Phi=phi,DeltaPhi=dphi,prior_cost=.5*float(p.x@p.x),
                DeltaPrior=dp,legacy_likelihood_competitive=dphi<=2,legacy_prior_reasonable=dp<=10,
                material_RMS=rms,target_RMS=radius,signed_progress=float(unit@dm),
                radius_progress_qualified=geometric,legal=legal,noise_ball=noise,
                primary_qualified=bool(geometric and legal and noise),
                whitened_pair_L2=float(np.linalg.norm(white(p.prediction-base.prediction))),
                original_prior_support=True,GT_used=False)

def substep(ph,base,p,v0,radius,arm,g0,c,trust):
    begin=time.perf_counter();begin_cpu=time.process_time()
    M=rmat(derivative(ph.prior,p.x));Qw=None;geometry_wall=0.;dynamic_ray_status=None
    if arm=='RANDOM66_FULL':allowed=np.eye(66)
    else:
        if arm=='FIXED_Q8_RAY':Qw=g0['Qw'];ray=v0/np.linalg.norm(v0)
        else:
            st=time.perf_counter();g=frame(ph.prior,p.x,p.J,8);geometry_wall=time.perf_counter()-st
            if g['rank']<8:raise RuntimeError('DYNAMIC_WEAK_RANK_BELOW8')
            Qw=g['Qw']
            weak=Qw@(Qw.T@(p.x-base.x));nn=np.linalg.norm(weak)
            if nn<1e-12:
                weak=Qw@(Qw.T@v0);nn=np.linalg.norm(weak);dynamic_ray_status='ANCHOR_DIRECTION_PROJECTION_FALLBACK'
            if nn<1e-12:raise RuntimeError('DYNAMIC_SELECTED_WEAK_RAY_UNDEFINED')
            ray=weak/nn
        Qs=null_space(Qw.T,rcond=1e-10)
        allowed=np.column_stack((Qs,ray))
        allowed,_=np.linalg.qr(allowed)
    _,sv,Vh=np.linalg.svd(M@allowed,full_matrices=False)
    keep=sv>max(sv[0]*1e-8,1e-10);effective=allowed@Vh[keep].T
    unit=g0['M']@v0;unit/=np.linalg.norm(unit)
    Z=effective@null_space((unit@M@effective)[None,:],rcond=1e-10)
    rr=white(p.prediction-ph.y);A=np.vstack((p.J@Z,Z));res=np.r_[rr,p.x]
    step=Z@(-np.linalg.lstsq(A,res,rcond=1e-10)[0])
    step*=min(1.,trust/max(np.linalg.norm(step),1e-300),.5*radius/max(np.linalg.norm(M@step),1e-300))
    return step,dict(allowed_dimension=allowed.shape[1],GN_dimension=Z.shape[1],
                     dynamic_geometry_wall_seconds=geometry_wall,dynamic_ray_status=dynamic_ray_status,
                     geometry_class='full66 ordinary' if arm=='RANDOM66_FULL' else 'weak-complement plus selected weak ray',
                     dynamic_geometry_new_RHS=0,GN_geometry_model_wall_seconds=time.perf_counter()-begin,
                     GN_geometry_model_CPU_seconds=time.process_time()-begin_cpu),M

def ray_metadata(g0,x,x0,v):
    w0=g0['Qw'].T@v;w=g0['Qw'].T@(x-x0)
    a=float(w@w0/max(w0@w0,1e-300))
    return dict(fixed_anchor_weakline_relative_defect=float(np.linalg.norm(w-a*w0)/max(np.linalg.norm(w),1e-30)),
                anchor_weak_projection_scalar=a,anchor_weak_projection_sign=int(np.sign(a)))

def candidate(ph,base,y,v,radius,arm,g0,dest,c):
    dest.mkdir(parents=True,exist_ok=True);start=time.perf_counter();cpu_start=time.process_time();before=ph.book.total.copy()
    visited=[];history=[];fcount=jcount=0;radius_wall=0.;radius_calls=0
    if v is None:
        row=dict(arm=arm,radius=radius,status='INSUFFICIENT_MATERIAL_OR_WEAK_RANK',primary_qualified=False,
                 forward_trials=0,J_states=0,GT_used=False,cost={k:0 for k in ['reduced_RHS256','reduced_LU256','actual_forward_RHS256','actual_J_RHS256']},wall_seconds=time.perf_counter()-start,selected_path=None)
        write(dest/'RESULT.json',row);return row
    unit=g0['M']@v;unit/=np.linalg.norm(unit)
    radius_start=time.perf_counter();xp,matched=match_radius(ph.prior,base.x,v,radius,c['latent_radius']);radius_wall+=time.perf_counter()-radius_start;radius_calls+=1
    if xp is None:
        row=dict(arm=arm,radius=radius,status=matched['status'],initial_match=matched,
                 primary_qualified=False,forward_trials=0,J_states=0,GT_used=False,
                 cost={k:0 for k in ['reduced_RHS256','reduced_LU256','actual_forward_RHS256','actual_J_RHS256']},wall_seconds=time.perf_counter()-start,selected_path=None)
        write(dest/'RESULT.json',row);return row
    def trial(x,label):
        nonlocal fcount
        if fcount>=6:raise BudgetStop('proposal forward trial cap')
        fcount+=1;p=ph.evaluate(x,label)
        m=metrics(ph,p,base,y,unit,radius,c);m.update(trial=fcount,role='RAW' if fcount==1 else 'CORRECTION_TRIAL',
                                                  accepted=False,has_J=False,x_path=str(ph.book.jobdir/'states'/f'{p.key}.npz'))
        visited.append(m);return p,m
    current,raw=trial(xp,arm+'_raw');ph.jacobian(current);jcount=1;raw['has_J']=True;raw['accepted']=True
    selected=current if raw['radius_progress_qualified'] and raw['legal'] else None
    best=raw if selected is not None else None;trust=.5
    for step_index in range(3):
        if fcount>=6 or jcount>=4:break
        try:step,extra,M=substep(ph,base,current,v,radius,arm,g0,c,trust)
        except RuntimeError as exc:
            history.append(dict(step_index=step_index,stop='CORRECTION_GEOMETRY_FAILURE',error=str(exc)));break
        hrow=dict(step_index=step_index,from_endpoint=current.key,trust=trust,**extra);history.append(hrow)
        accepted=False
        for ls in range(c['candidate']['max_decoder_backtracks']):
            if fcount>=6:break
            proposal=current.x+step*(.5**ls)
            radius_start=time.perf_counter();xx,ma=match_radius(ph.prior,base.x,proposal-base.x,radius,c['latent_radius']);radius_wall+=time.perf_counter()-radius_start;radius_calls+=1
            if xx is None:continue
            e=xx-current.x
            if np.linalg.norm(e)>trust*(1+1e-10) or np.linalg.norm(M@e)>.5*radius*(1+1e-10):continue
            dm=rmat(decode(ph.prior,xx)-base.material)
            if abs(np.linalg.norm(dm)-radius)>.02*radius or unit@dm<.5*radius:continue
            rr=white(current.prediction-y);aug=np.r_[rr,current.x]
            pred=.5*float(aug@aug-np.linalg.norm(aug+np.r_[current.J@e,e])**2)
            if pred<=0 or np.linalg.norm(e)<1e-9:continue
            pp,mm=trial(xx,arm+f'_step{step_index}_line{ls}')
            # Primary selection is true dataPhi, even though the model is MAP.
            if mm['radius_progress_qualified'] and mm['legal'] and (best is None or mm['Phi']<best['Phi']):
                ph.jacobian(pp);jcount+=1;mm['has_J']=True;mm['accepted']=True
                current_metric_Phi=ph.phi(current,y);actual=ph.objective(current,y)-ph.objective(pp,y);ratio=actual/pred
                if ratio<.25:trust*=.5
                elif ratio>.75:trust=min(1.,trust*1.5)
                hrow.update(accepted_endpoint=pp.key,line_search=ls,predicted_MAP_decrease=pred,
                            actual_MAP_decrease=actual,ratio=ratio,dataPhi_decrease=current_metric_Phi-mm['Phi'])
                selected=current=pp;best=mm;accepted=True;break
            mm['rejection']='NO_TRUE_DATAPHI_IMPROVEMENT';trust*=.5
        if not accepted:hrow['stop']='NO_ACCEPTED_DATA_DESCENT_WITHIN_BUDGET';break
    if selected is not None:ph.persist(selected,dest/'SELECTED.npz',full_fields=True)
    rawpoint_path=ph.book.jobdir/'states'/f"{raw['endpoint']}.npz"
    # Raw proposal and every visited true state remain available regardless of final selection.
    row=dict(raw if best is None else best,arm=arm,radius=radius,status='EXECUTED' if selected is not None else 'NO_RADIUS_PROGRESS_ENDPOINT',selected_path=None if selected is None else str(dest/'SELECTED.npz'),forward_trials=fcount,J_states=jcount,
             raw_endpoint=raw['endpoint'],raw_primary_qualified=raw['primary_qualified'],raw_noise_qualified=raw['noise_ball'] and raw['legal'],raw_path=str(rawpoint_path),
             selected_is_raw=selected is not None and selected.key==raw['endpoint'],
             weak_line_metadata=None if selected is None else ray_metadata(g0,selected.x,base.x,v),
             decoder_radius_wall_seconds=radius_wall,decoder_radius_calls=radius_calls,
             GN_geometry_model_wall_seconds=sum(h.get('GN_geometry_model_wall_seconds',0.) for h in history),
             cost={k:ph.book.total[k]-before[k] for k in ['reduced_RHS256','reduced_LU256','actual_forward_RHS256','actual_J_RHS256']},
             wall_seconds=time.perf_counter()-start,visited_count=len(visited))
    write(dest/'VISITED.json',dict(rows=visited,history=history,no_posterior_weights=True))
    ph.book.save(dest/'VISITED_COORDINATES.npz',x=np.array([load_npz(vv['x_path'])['x'] for vv in visited]),
         Phi=np.array([vv['Phi'] for vv in visited]),primary_qualified=np.array([vv['primary_qualified'] for vv in visited]),
         accepted=np.array([vv['accepted'] for vv in visited]),has_J=np.array([vv['has_J'] for vv in visited]))
    row['wall_seconds']=time.perf_counter()-start;row['CPU_seconds']=time.process_time()-cpu_start
    row['cost'].update({k:ph.book.total[k]-before[k] for k in ['physics_wall_seconds','save_wall_seconds','array_save_attempts','array_saves','array_bytes','failed_actions']})
    row['cost'].update(wall_seconds=row['wall_seconds'],CPU_seconds=row['CPU_seconds'],decoder_radius_wall_seconds=radius_wall,
                       GN_geometry_model_wall_seconds=row['GN_geometry_model_wall_seconds'])
    write(dest/'RESULT.json',row);return row

def online_observation(c,parent):
    path=Path(c['output'])/'inputs/online'/f'{parent}.npz'
    a=load_npz(path)
    if set(a)!={'observation','sigma','parent'} or str(a['parent'])!=parent or float(a['sigma'])!=SIGMA:raise RuntimeError('observation-only input contract')
    y=a['observation']
    if y.shape!=(32,8) or not np.isfinite(y).all():raise RuntimeError('finite observation layout')
    return y

def anchor(ph,y,parent,dest,c):
    dest.mkdir(parents=True,exist_ok=True);last=best=None;best_energy=np.inf;fresh=0;cap=64;history=[]
    def paid(x,final=False):
        nonlocal fresh
        if fresh>=cap-(0 if final else 1):raise RuntimeError('ANCHOR_STRICT_STATE_CAP')
        fresh+=1;return ph.evaluate(x,parent+'_direct_probit64')
    def get(x):
        nonlocal last,best,best_energy
        for p in (last,best):
            if p is not None and np.array_equal(p.x,x):return p
        p=paid(x);obj=ph.objective(p,y);last=p
        history.append(dict(attempt=fresh,endpoint=p.key,Phi=ph.phi(p,y),objective=obj,prior=.5*float(x@x)))
        append(dest/'ANCHOR_HISTORY.jsonl',history[-1])
        if obj<best_energy:best=p;best_energy=obj
        return p
    def fun(x):return np.r_[white(get(x).prediction-y),x]
    def jac(x):
        p=get(x)
        if p.J is None and p is not ph.live:p=paid(x)
        return np.vstack((ph.jacobian(p),np.eye(66)))
    status={}
    try:
        r=least_squares(fun,np.zeros(66),jac=jac,method='trf',tr_solver='exact',x_scale='jac',
                        ftol=1e-9,xtol=1e-9,gtol=1e-7,max_nfev=cap)
        status=dict(status=int(r.status),success=bool(r.success),message=str(r.message),nfev=int(r.nfev),njev=int(r.njev),optimality=float(r.optimality))
    except RuntimeError as e:
        if str(e)!='ANCHOR_STRICT_STATE_CAP':raise
        status=dict(status='STRICT_STATE_CAP',success=False)
    if best is None:raise RuntimeError('no finite anchor')
    old=best
    if old is not ph.live:best=paid(old.x,final=True)
    ph.jacobian(best);ph.persist(best,dest/'ANCHOR.npz',full_fields=True)
    gradient=best.J.T@white(best.prediction-y)+best.x
    write(dest/'ANCHOR_STATUS.json',dict(optimizer=status,fresh_state_attempts=fresh,Phi=ph.phi(best,y),
          objective=ph.objective(best,y),gradient_inf=float(max(abs(gradient))),
          noise_ball=ph.phi(best,y)<=c['primary']['noise_ball'],legal=ph.prior.legal(best.q),GT_used=False,
          mode='PROBIT_COORDINATE_POINT_NOT_MATERIAL_MAP',backend='EXACT_NONUNIFORM_ROI256'))
    return best

def load_anchor(path):
    a=load_npz(path)
    return Point(a['x'],a['q'],a['material'],a['prediction'],a['reduced_current'],a['ROI_field'],
                 str(path),float(a['backward_residual']),a['J'])

def parent_run(c,parent):
    root=Path(c['output']);dest=root/'parents'/parent
    if (dest/'COMPLETE.json').exists():raise RuntimeError('completed parent no replay')
    y=online_observation(c,parent);entry=next(r for r in c['parents'] if r['parent']==parent)
    with Book(c,parent,'ONLINE_ANCHOR_AND_CLOUD') as book:
        setup_clock=time.perf_counter();setup_cpu=time.process_time();ph=ExactROI(book);ph.y=y
        runtime_setup=dict(wall_seconds=time.perf_counter()-setup_clock,CPU_seconds=time.process_time()-setup_cpu,setup_cache_hit=True,new_full_RHS576=0)
        try:
            anchor_begin=book.total.copy();anchor_clock=time.perf_counter();anchor_cpu=time.process_time()
            base=anchor(ph,y,parent,dest/'anchor',c)
            anchor_cost={k:book.total[k]-anchor_begin[k] for k in ['reduced_RHS256','reduced_LU256','actual_forward_RHS256','actual_J_RHS256']}
            anchor_wall=time.perf_counter()-anchor_clock;anchor_cost.update(wall_seconds=anchor_wall,CPU_seconds=time.process_time()-anchor_cpu)
            geometry_clock=time.perf_counter();geometry_cpu=time.process_time();g0=frame(ph.prior,base.x,base.J,8);shared_geometry_wall=time.perf_counter()-geometry_clock
            fixed_frame_setup=dict(wall_seconds=shared_geometry_wall,CPU_seconds=time.process_time()-geometry_cpu,common_initial_basis_for_B_and_C=True)
            book.save(dest/'ANCHOR_FRAME.npz',**g0)
            rows=[];proposals=[]
            for pair in range(2):
                for ri,radius in enumerate(c['radii']):
                    weak=np.random.default_rng(np.random.SeedSequence(entry['weak_seed_prefix']+[ri,pair])).normal(size=8)
                    fixed=None
                    if g0['rank']>=8:
                        fixed=g0['weaklift']@weak;fixed/=np.linalg.norm(g0['M']@fixed)
                    rv=np.random.default_rng(np.random.SeedSequence(entry['random_seed_prefix']+[ri,pair])).normal(size=66)
                    rv_norm=np.linalg.norm(g0['M']@rv)
                    rv=rv/rv_norm if rv_norm>1e-12 else None
                    for sign in [1,-1]:
                        pid=f'r{ri}_p{pair}_{"plus" if sign>0 else "minus"}'
                        slot=pair*8+ri*2+(0 if sign>0 else 1)
                        proposals.append(dict(id=pid,slot=slot,radius=radius,pair=pair,sign=sign,weak_seed=entry['weak_seed_prefix']+[ri,pair],
                                              random_seed=entry['random_seed_prefix']+[ri,pair]))
                        order=ARMS[(entry['index']+ri+pair)%3:]+ARMS[:(entry['index']+ri+pair)%3]
                        for arm in order:
                            book.context.update(parent=parent,arm=arm,proposal=pid,radius=radius)
                            direction=rv if arm=='RANDOM66_FULL' else fixed
                            v=None if direction is None else sign*direction
                            proposal_clock=time.perf_counter();proposal_cpu=time.process_time()
                            row=candidate(ph,base,y,v,radius,arm,g0,dest/arm/pid,c)
                            row.update(parent=parent,proposal=pid,slot=slot,pair=pair,sign=sign,arm_order=order,
                                       timing=dict(window_id=c['id']+':'+parent+':'+pid,arm_order_index=order.index(arm),arm_order=order),GT_used=False)
                            write(dest/arm/pid/'RESULT.json',row)
                            row.update(wall_seconds=time.perf_counter()-proposal_clock,CPU_seconds=time.process_time()-proposal_cpu)
                            row['cost'].update(wall_seconds=row['wall_seconds'],CPU_seconds=row['CPU_seconds'])
                            write(dest/arm/pid/'RESULT.json',row);rows.append(row)
                            print('VOLUME_ENDPOINT',parent,arm,pid,row['status'],row.get('Phi'),book.total['reduced_RHS256'],flush=True)
            cloud={}
            for arm in ARMS:
                ar=[r for r in rows if r['arm']==arm]
                selected=[load_npz(dest/arm/r['proposal']/'SELECTED.npz')['x'] for r in ar if r['status']=='EXECUTED']
                if selected:book.save(dest/arm/'PRIMARY_CLOUD.npz',x=np.vstack([base.x,np.array(selected)]),
                                     anchor_included=np.array(True),primary_qualified=np.r_[ph.prior.legal(base.q) and ph.phi(base,y)<=c['primary']['noise_ball'],np.array([r['primary_qualified'] for r in ar if r['status']=='EXECUTED'])],posterior_sample=np.array(False))
                cloud[arm]=dict(proposals=16,available=sum(r['status']=='EXECUTED' for r in ar),
                                primary_qualified=sum(r.get('primary_qualified',False) for r in ar),
                                anchor_noise_qualified=ph.prior.legal(base.q) and ph.phi(base,y)<=c['primary']['noise_ball'],anchor_legal=ph.prior.legal(base.q),
                                max_unique_with_shared_anchor=17,posterior_claim=False)
            parent_wall=time.perf_counter()-book.started;parent_cpu=time.process_time()-book.cpu0
            common_overhead=dict(wall_seconds=max(0.,parent_wall-sum(r['wall_seconds'] for r in rows)-anchor_wall-fixed_frame_setup['wall_seconds']-runtime_setup['wall_seconds']),
                                 CPU_seconds=max(0.,parent_cpu-sum(r['CPU_seconds'] for r in rows)-anchor_cost['CPU_seconds']-fixed_frame_setup['CPU_seconds']-runtime_setup['CPU_seconds']))
            write(dest/'COMPLETE.json',dict(parent=parent,rows=rows,cloud=cloud,proposal_manifest=proposals,
                  parent_end_to_end_wall_seconds=parent_wall,parent_CPU_seconds=parent_cpu,common_unallocated_overhead=common_overhead,
                  primary_cloud_storage='allselected_with_flags',
                  data_only_selection=True,MAP_model_prior_weight=1.,truth_read=False,
                  anchor_path=str(dest/'anchor/ANCHOR.npz'),anchor_cost=anchor_cost,anchor_wall_seconds=anchor_wall,
                  shared_geometry_wall_seconds=shared_geometry_wall,fixed_frame_setup=fixed_frame_setup,runtime_setup=runtime_setup,shared_setup=c['shared_setup'],
                  same_window_order='rotated cyclic arm order per parent/radius/pair; exact order retained on each row',
                  live_cumulative_cost=book.total,exact_ROI_support='all frozen P1 materials'))
        finally:ph.close()

def seal(c):
    root=Path(c['output']);records=[]
    for p in c['parents']:
        dest=root/'parents'/p['parent'];done=dest/'COMPLETE.json'
        if not done.exists():raise RuntimeError('all32 parent predictions required before scoring')
        a=json.loads(done.read_text())
        if len(a['rows'])!=48:raise RuntimeError('all3arms*16 records required including unavailable/failures')
        files=[]
        for path in sorted(dest.rglob('*')):
            if path.is_file():files.append(dict(path=str(path),bytes=path.stat().st_size))
        jobdir=root/'jobs'/p['parent']
        if not jobdir.exists():raise RuntimeError('missing online parent job evidence')
        for path in sorted(jobdir.rglob('*')):
            if path.is_file():files.append(dict(path=str(path),bytes=path.stat().st_size))
        records.append(dict(parent=p['parent'],files=files,cloud=a['cloud']))
    write(root/'PREDICTIONS_SEALED.json',dict(protocol_id=c['id'],utc=time.time(),status='ALL32_ONLINE_PREDICTIONS_SEALED',
          parents=records,truth_read=False,source_files=c['source_files'],statistics_config_id=c['statistics_config_id'],posterior_claim=False,
          raw_and_visited_preserved=True,scoring_may_start_only_after_this_seal=True))

def main():
    p=argparse.ArgumentParser();p.add_argument('command',choices=['draft','plan','parent','all','seal'])
    p.add_argument('--config',default=str(CONFIG));p.add_argument('--parent')
    a=p.parse_args()
    if a.command=='draft':
        c=draft();path=P/'config/possibility_volume_protocol_DRAFT.json'
        if path.exists():raise RuntimeError('draft exists; root must version')
        write(path,c);print(json.dumps(estimate(c),indent=2));return
    if a.command=='plan':
        path=Path(a.config);c=json.loads(path.read_text()) if path.exists() else draft()
        print(json.dumps(estimate(c),indent=2));return
    c=read_protocol(a.config)
    if a.command=='parent':parent_run(c,a.parent)
    elif a.command=='all':
        for r in c['parents']:parent_run(c,r['parent'])
    else:seal(c)
if __name__=='__main__':main()
