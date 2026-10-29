"""R9 E4 actual TM24 fixed-amplitude numerical transfer; historical inputs read-only.

Formula and source prior inherited R6/R8P. All full complex columns charged
before execution; reduced work separately counted. No rigorous rounding or
Schur coefficient enclosure is asserted by this numerical implementation.
"""
from __future__ import annotations
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[k]='1'
import sys,time,json,resource,platform,contextlib
sys.dont_write_bytecode=True
from pathlib import Path
from collections import Counter
import numpy as np
from scipy.linalg import lu_factor,lu_solve,schur,solve_triangular
from scipy.special import ndtr,betaincinv
from scipy.stats import norm
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss
P=Path(__file__).resolve().parents[1]
ROOT=P.parents[4]
R4=ROOT/'Gaussian/B1/b1-p3l-r4-conditional-transport'
R6=ROOT/'Gaussian/B1/b1-p3l-r6-task-optimal/research/p3l_r6_task_optimal_maxwell'
RES=P/'results/probability'; AUD=P/'audit/probability'
SIGMA=.00119211798514519; HOST=.1+.04j
Z0=ndtr(4)-ndtr(-4); A=1/1.3
PH=(ndtr(-A)-ndtr(-4))/Z0
ULO,UHI=ndtr(-4),ndtr(4); RLO,RHI=.8+ULO,.8+UHI; ILO,IHI=.03+.15*ULO,.03+.15*UHI

def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,complex):return [x.real,x.imag]
    if isinstance(x,dict):return {str(k):plain(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [plain(v) for v in x]
    if isinstance(x,Path):return str(x)
    return x

def write(path,data):
    path=Path(path);path.parent.mkdir(exist_ok=True,parents=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(plain(data),indent=2,ensure_ascii=False,allow_nan=False)+'\n');tmp.replace(path)

def pack(x):return np.concatenate((x.real.T,x.imag.T),axis=1).ravel()
def sat(x):
    t=np.clip((np.asarray(x)+1)/2,0,1);return t**3*(10+t*(-15+6*t))
def gl(n,a,b):
    x,w=leggauss(n);return a+(b-a)*(x+1)/2,w*(b-a)/2

def fibre(v,n=96):
    R,I=v.real,v.imag;lo=max(R/RHI,I/IHI);hi=min(R/RLO,I/ILO,1.)
    if not 0<lo<hi:return np.zeros(2)
    zl,zh=(2*betaincinv(3,3,np.array([lo,hi]))-1)/1.3
    z,w=gl(n,zl,zh);q=sat(1.3*z)
    density=float(np.sum(w*norm.pdf(z)/Z0/q**2))/((RHI-RLO)*(IHI-ILO))
    if zh<=0:phase=0.
    elif zl>=0:phase=density
    else:
        z,w=gl(n,0.,zh);q=sat(1.3*z);phase=float(np.sum(w*norm.pdf(z)/Z0/q**2))/((RHI-RLO)*(IHI-ILO))
    return np.array([density,phase])

def prior(v,n=96):
    transition,phase=fibre(v,n)
    pure=PH/((RHI-RLO)*(IHI-ILO)) if RLO<v.real<RHI and ILO<v.imag<IHI else 0.
    return transition+pure,np.array([0.,transition,pure]),phase+pure

def legal(v):
    if v==0:return True
    return 0<max(v.real/RHI,v.imag/IHI)<min(v.real/RLO,v.imag/ILO,1.)

class Book:
    def __init__(self,job,cap=4000):
        self.job=job;self.start=time.perf_counter();self.cpu=time.process_time();self.cap=cap
        self.full_rhs=0;self.previous_CPU=0.;self.counts=Counter();self.timings=Counter()
        self.path=AUD/'ledger.jsonl'
        if self.path.exists():
            for line in self.path.read_text().splitlines():
                row=json.loads(line)
                if row['event']=='RESERVE_FULL':self.full_rhs+=row['full_complex_RHS']
                if row['event']=='JOB_END':self.previous_CPU+=row.get('CPU',0.)
        self.event('JOB_BEGIN',total_previous_full_RHS=self.full_rhs)
    def event(self,event,**kw):
        row=dict(event=event,job=self.job,utc=time.time(),**plain(kw))
        with self.path.open('a') as f:f.write(json.dumps(row,ensure_ascii=False,allow_nan=False)+'\n');f.flush()
    def reserve(self,n,lu,label):
        if self.full_rhs+n>self.cap:raise RuntimeError('BUDGET_EXHAUSTED')
        self.event('RESERVE_FULL',full_complex_RHS=n,full_LU=lu,label=label,cumulative_full_RHS=self.full_rhs+n)
        self.full_rhs+=n;self.counts['full_complex_RHS']+=n;self.counts['full_LU']+=lu
    @contextlib.contextmanager
    def span(self,label,**counts):
        if self.previous_CPU+time.process_time()-self.cpu>3600:raise RuntimeError('REDUCED_CPU_CAP_EXHAUSTED')
        t=time.perf_counter();c=time.process_time();self.counts.update(counts)
        try:yield
        finally:self.timings[label+'_wall']+=time.perf_counter()-t;self.timings[label+'_CPU']+=time.process_time()-c
    def close(self,status,error=None):
        self.event('JOB_END',status=status,error=error,total_full_RHS=self.full_rhs,counts=self.counts,timings=self.timings,wall=time.perf_counter()-self.start,CPU=time.process_time()-self.cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write(RES/(self.job+'_receipt.json'),dict(status=status,error=error,total_full_RHS=self.full_rhs,counts=self.counts,timings=self.timings,wall=time.perf_counter()-self.start,CPU=time.process_time()-self.cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))

def load_inputs():
    cfg=json.loads((R4/'research/p3l_r4_conditional_transport/configs/TARGET.json').read_text())
    geom=json.loads((R4/'research/p1_material_target_20261018/inputs/geometry.json').read_text())
    points=np.asarray(geom['points']);mask=np.linalg.norm(points,axis=1)<=cfg['roi_radius']
    with np.load(R4/'research/p3l_r4_conditional_transport/inputs/background/weak.npz') as z:
        G,S,incident=z['G'],z['S'],z['incident']
    with np.load(R6/'track_a/OBSERVATIONS.npz') as z:
        names=list(z['names']);obs={name:z['observation'][names.index(name)] for name in ('REGULAR','PURE_OVERLAP')};sigmas={name:float(z['sigma'][names.index(name)]) for name in obs}
    assert G.shape==(576,576) and S.shape==(32,576) and incident.shape==(576,8)
    assert int(mask.sum())==256 and np.allclose(cfg['host_chi'],[.1,.04],rtol=0,atol=0)
    assert all(v==SIGMA for v in sigmas.values())
    assert cfg['packing']=='source-major Re32 then Im32'
    return G,S,incident,mask,obs

class Physics:
    def __init__(self,book):
        self.book=book;self.G,self.S,self.incident,self.mask,self.obs=load_inputs();self.ix=np.where(self.mask)[0]
    def setup(self):
        b=self.book;b.reserve(264,1,'nonzero_host_Eh8_and_ROI_basis256')
        with b.span('setup_full_host',host_source_RHS=8,host_basis_RHS=256):
            Ah=np.eye(576,dtype=complex)-self.G*HOST
            factor=lu_factor(Ah,check_finite=False)
            Eh=lu_solve(factor,self.incident,check_finite=False)
            B=lu_solve(factor,self.G[:,self.ix],check_finite=False)
            self.Fh=self.S@(HOST*Eh);C=B[self.ix];ER=Eh[self.ix]
            Op=self.S[:,self.ix]+HOST*(self.S@B)
            host_res=float(np.linalg.norm(Ah@Eh-self.incident)/np.linalg.norm(self.incident))
            basis_res=float(np.linalg.norm(Ah@B-self.G[:,self.ix])/np.linalg.norm(self.G[:,self.ix]))
            del factor,Ah
        with b.span('Schur256_and_transforms',Schur256=1):
            self.T,Q=schur(C,output='complex',check_finite=False)
            self.O=Op@Q;self.ER=Q.conj().T@ER
            schur_res=float(np.linalg.norm(C-Q@self.T@Q.conj().T)/np.linalg.norm(C))
            unitary_res=float(np.linalg.norm(Q.conj().T@Q-np.eye(256),ord=2))
        with b.span('setup_diagnostics'):
            self.Cnorm=float(np.linalg.norm(C,ord=2));self.absT=np.abs(self.T)
            diagnostics=dict(host_residual=host_res,basis_residual=basis_res,Schur_relative_residual=schur_res,Schur_unitary_residual=unitary_res,C_2norm=self.Cnorm,Schur_diagonal_maxabs=float(np.max(np.abs(np.diag(self.T)))),matrix_order=256,qualification='NUMERICAL_PHYSICS; no interval rounding/coefficient enclosure')
        np.savez_compressed(RES/'transfer_setup.npz',T=self.T,O=self.O,ER=self.ER,Fh=self.Fh)
        write(RES/'setup_diagnostics.json',diagnostics)
        return diagnostics
    def load_setup(self):
        with np.load(RES/'transfer_setup.npz') as z:
            self.T,self.O,self.ER,self.Fh=(z[k] for k in ('T','O','ER','Fh'))
        self.absT=np.abs(self.T)
    def reduced(self,v):
        with self.book.span('reduced_triangular256',reduced_nodes=1,reduced_complex_RHS256=8,readout32x256_by256x8=1):
            x=solve_triangular(np.eye(256,dtype=complex)-v*self.T,self.ER,check_finite=False)
            return self.Fh+self.O@(v*x)
    def direct(self,v):
        self.book.reserve(8,1,'direct_TM24')
        with self.book.span('direct_TM24',direct_nodes=1):
            chi=np.full(576,HOST,complex);chi[self.mask]+=v
            # Original A24 current equation; one complete LU alive only here.
            L=np.eye(576,dtype=complex)-chi[:,None]*self.G
            fac=lu_factor(L,check_finite=False);rhs=chi[:,None]*self.incident
            current=lu_solve(fac,rhs,check_finite=False)
            residual=float(np.linalg.norm(L@current-rhs)/np.linalg.norm(rhs))
            out=self.S@current
            del fac,L,current
            return out,residual
    def majorant_derivative(self,a,hr,hi):
        """R8P triangular bound in ordinary float; diagnostic, not certificate."""
        radius=float(np.hypot(hr,hi));vmax=abs(a)+radius
        diagonal=np.diag(self.T)
        d=np.maximum(0,np.abs(1-a*diagonal)-radius*np.abs(diagonal))
        if np.any(d<=0):return None
        x=np.abs(self.ER).copy()
        for _ in range(2):
            y=np.zeros_like(x)
            for i in range(255,-1,-1):y[i]=(x[i]+vmax*(self.absT[i,i+1:]@y[i+1:]))/d[i]
            x=y
        output=np.abs(self.O)@x
        return dict(radius=radius,white_output_radius=float(np.sqrt(2)/SIGMA*radius*np.linalg.norm(output)),minimum_diagonal_margin=float(d.min()),max_majorant=float(np.max(x)))

def zero_audit():
    G,S,inc,mask,obs=load_inputs()
    sys.path.insert(0,str(R4/'src'))
    from p1_material.prior import MaterialPrior
    cfg=json.loads((R4/'research/p3l_r4_conditional_transport/configs/TARGET.json').read_text());points=json.loads((R4/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'];p=MaterialPrior(cfg,points)
    rows=[]
    for z0,u,v in [(-2,.2,.7),(0,.4,.2),(.5,.8,.6),(2,.1,.9)]:
        q=np.zeros(66);q[0]=z0;q[64:]=[u,v];chi,a=p.decode(q);V=sat(1.3*z0)*complex(.8+u,.03+.15*v)
        rows.append(dict(q0=z0,material_formula_max_error=float(np.max(abs(chi-(HOST+mask*V)))),passive=bool(np.min(chi.imag)>=0),ROI=int(mask.sum())))
    source_paths=[R4/'src/p1_material/prior.py',R4/'research/p3l_r4_conditional_transport/configs/TARGET.json',R4/'research/p3l_r4_conditional_transport/inputs/background/weak.npz',R6/'track_a/OBSERVATIONS.npz',R6/'track_a/LOCAL_REFINEMENT_RESULTS.json',R6/'track_a/stratified.py',R6/'track_a/local_material_reference.py']
    report=dict(status='PASS',new_full_complex_RHS=0,formula_replays=rows,input_files=[dict(path=str(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in source_paths],normalizer_Z0=Z0,host_atom_prior=PH,pure_layer_prior=PH,transition_prior=1-2*PH,noise=dict(sigma_complex=SIGMA,real_and_imag_variance=SIGMA**2/2),historical_reference_qualification='LOCAL_GH15_21_ONLY; global_tail_certificate=false; not promoted',note='No file hashes added; input paths, size and array contract checks recorded.')
    assert all(r['material_formula_max_error']<2e-15 and r['passive'] for r in rows)
    write(AUD/'zero_RHS_contract_audit.json',report)


def probe():
    if (RES/'probe_receipt.json').exists():raise RuntimeError('probe already executed; append-only no replay')
    b=Book('probe');status='FAILED';err=None
    try:
        ph=Physics(b);diag=ph.setup();cfg=json.loads((P/'config/probability_frozen.json').read_text());rows=[]
        for vr,vi in cfg['validation_amplitudes']:
            v=complex(vr,vi);red=ph.reduced(v);direct,res=ph.direct(v)
            white=float(np.linalg.norm(pack(red-direct))*np.sqrt(2)/SIGMA)
            rows.append(dict(V=[vr,vi],white_difference=white,backward_residual=res,legal_in_R6_support=legal(v)))
        diag.update(max_white_difference=max(x['white_difference'] for x in rows),max_direct_residual=max(x['backward_residual'] for x in rows))
        assert diag['max_white_difference']<=1e-8 and diag['max_direct_residual']<=1e-10
        # Compare old independently paid stored GH node phi, without charging old physics as new.
        hist=[]
        for name in ph.obs:
            with np.load(R6/f'track_a/raw/{name}_material_GH21.npz') as z:
                for j in [0,len(z['V'])//2,len(z['V'])-1]:
                    v=complex(*z['V'][j]);f=ph.reduced(v);phi=float(np.sum(abs(f-ph.obs[name])**2)/SIGMA**2)
                    hist.append(dict(case=name,node=j,V=[v.real,v.imag],old_phi=float(z['phi'][j]),new_phi=phi,phi_absolute_difference=abs(phi-float(z['phi'][j]))))
        write(RES/'physical_validation.json',dict(status='NUMERICAL_PHYSICS',diagnostics=diag,validation=rows,historical_paid_node_replay=hist,selected_cell_majorants=[dict(case=name,cell=ph.majorant_derivative(complex(*np.load(R6/f'track_a/raw/{name}_material_GH21.npz')['mode']),.001,.001)) for name in ph.obs],strict_physics_enclosure=None))
        status='PASS_NUMERICAL_PHYSICS';print(json.dumps(plain(diag)),flush=True)
    except BaseException as e:err=repr(e);write(AUD/'probe_failure.json',dict(error=err));raise
    finally:b.close(status,err)

if __name__=='__main__':
    if sys.argv[1]=='audit':zero_audit()
    elif sys.argv[1]=='probe':probe()
