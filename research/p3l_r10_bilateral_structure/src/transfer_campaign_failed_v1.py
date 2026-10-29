"""R10 fixed-shape exact transfer and conditional property integrals.
Original P1 inputs read only. Every attempted complex solve is reserved first.
"""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[key]='1'
import sys,time,json,csv,resource,fcntl,argparse
from pathlib import Path
import numpy as np
from scipy.linalg import schur,lu_factor,lu_solve,solve_triangular,get_lapack_funcs
from scipy.optimize import least_squares
from scipy.special import ndtr,logsumexp
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss
P=Path(__file__).absolute().parents[1]; B1=P.parents[2]
R9=B1/'b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace'
R4=B1/'b1-p3l-r4-conditional-transport'; V=R9/'results/possibility_volume'
OUT=P/'results/transfer'; SIGMA=.00119211798514519
LO=np.array([.8,.03]);HI=np.array([1.8,.18]);HOST=.1+.04j
sys.path.insert(0,str(R4/'src'))
from p1_material.prior import MaterialPrior

def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,complex):return [x.real,x.imag]
    if isinstance(x,Path):return str(x)
    if isinstance(x,dict):return {k:plain(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [plain(v) for v in x]
    return x

def write(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(plain(x),indent=2,allow_nan=False)+'\n')
def pack(z):return np.concatenate((z.real.T,z.imag.T),axis=1).ravel()
def white(z):return pack(z)*np.sqrt(2)/SIGMA

def freeze():
    dest=P/'config/transfer_frozen.json'
    if dest.exists():raise RuntimeError('freeze exists')
    cfg=dict(id='R10_SHAPE_AMPLITUDE_CONDITIONAL_V1',frozen_utc=time.time(),parents=[f'R9V_FRESH_{i:03d}' for i in range(8)],
        conditional_parents=[f'R9V_FRESH_{i:03d}' for i in range(3)],amplitude_rectangle=[LO,HI],uniform_density=1/.15,
        timing_repeats=30,amplitude_schedule='u_r=(r+.5)/30, u_i=frac((r+.5)*golden_ratio); same nodes both methods, rotate execution order',
        M=[1,4,16,64,256,1024],RHS256_cap=50000,LU256_cap=2000,RHS576_cap=0,wall_seconds_cap=900,RSS_bytes_cap=6000000000,
        BLAS_threads=1,host_setup_cache=str(R9/'results/stable_roi/SETUP.npz'),historical_setup='264 RHS576 and 1 LU576 in R9; no new full solve',
        shape='fixed original Gaussian64 of exposed anchors; nonnegative quintic occupancy; P=sqrt(a)',
        algorithm='complex Schur of P C P, exact triangular solves with complex amplitude; no Hermitian eigensolver',
        optimizer='bounded least_squares original amplitudes with exact analytic holomorphic derivative, max100 state evaluations/parent',
        integration='bounded original-coordinate composite GL; core box mode +/-8 marginal Hessian SD clipped to full prior; eight tail rectangles all included; refine core8/12/20 tails2/3/5; independently rotated Hessian GH8/12, GH18 if need and within budget',
        event_thresholds=[1.3,.105],event='original material ROI mean above host + threshold*ROI fraction mean, equivalent delta coordinate exceeds fixed threshold',
        numeric_tolerances=dict(relative_evidence=2e-4,amplitude_mean=2e-5,event_probability=5e-4),
        qualification='numerical finite fixed-shape conditional only; no rigorous tails/rounding enclosure, no Full66 posterior; axis full rectangle coverage alone not certified',
        no_GT=True,max_conditional_nodes_per_parent=1500)
    write(dest,cfg);print('FROZEN',dest)

class Book:
    def __init__(self,c):
        self.c=c;self.start=time.perf_counter();self.cpu=time.process_time();self.rhs=self.lu=0;self.seq=0
        self.path=OUT/'COST_LEDGER.jsonl'
        if self.path.exists():raise RuntimeError('no automatic replay')
        self.event('BEGIN',config=c['id'],pid=os.getpid())
    def event(self,event,**kw):
        self.seq+=1
        row=dict(sequence=self.seq,event=event,utc=time.time(),RHS256=self.rhs,LU256=self.lu,RHS576=0,**kw)
        with self.path.open('a') as f:f.write(json.dumps(plain(row),allow_nan=False)+'\n');f.flush()
    def reserve(self,rhs,lu,label):
        if self.rhs+rhs>self.c['RHS256_cap'] or self.lu+lu>self.c['LU256_cap']:raise RuntimeError('SOLVE_BUDGET_STOP')
        if time.perf_counter()-self.start>self.c['wall_seconds_cap']:raise RuntimeError('WALL_BUDGET_STOP')
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>self.c['RSS_bytes_cap']:raise RuntimeError('RSS_BUDGET_STOP')
        self.rhs+=rhs;self.lu+=lu;self.event('RESERVE_ATTEMPT',rhs=rhs,lu=lu,label=label)
    def close(self,status,error=None):
        row=dict(status=status,error=error,RHS256=self.rhs,LU256=self.lu,RHS576=0,wall_seconds=time.perf_counter()-self.start,CPU_seconds=time.process_time()-self.cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        self.event('END',**{k:v for k,v in row.items() if k not in ('RHS256','LU256','RHS576')});write(OUT/'RECEIPT.json',row)

class Shape:
    def __init__(self,parent,setup,prior,book):
        self.parent=parent;self.book=book;self.Fh=setup['Fh'];self.C=setup['C'];self.ER=setup['ER'];self.O=setup['O'];self.roi=setup['ROI']
        with np.load(V/'parents'/parent/'anchor/ANCHOR.npz') as z:self.x=z['x'].copy();self.old_prediction=z['prediction'].copy()
        q=self.x.copy();q[64:]=ndtr(q[64:]);self.anchor_delta=complex(.8+q[64],.03+.15*q[65]);self.material,self.meta=prior.decode(q);self.a=self.meta['fraction'][self.roi];self.amean=float(self.a.mean())
        with np.load(V/'inputs/online'/f'{parent}.npz') as z:
            assert set(z.files)=={'observation','sigma','parent'} and float(z['sigma'])==SIGMA
            self.y=z['observation'].copy()
        t=time.perf_counter();p=np.sqrt(self.a);self.B=p[:,None]*self.C*p[None,:];self.L=self.O*p[None,:];self.R=p[:,None]*self.ER
        self.formation_wall=time.perf_counter()-t
        t=time.perf_counter();self.T,self.Q=schur(self.B,output='complex',check_finite=False);self.LQ=self.L@self.Q;self.QR=self.Q.conj().T@self.R
        self.schur_wall=time.perf_counter()-t;self.I=np.eye(256,dtype=complex)
        self.setup=dict(parent=parent,formation_wall=self.formation_wall,Schur_transforms_wall=self.schur_wall,order=256,
            Schur_relative_residual=float(np.linalg.norm(self.B-self.Q@self.T@self.Q.conj().T)/max(np.linalg.norm(self.B),1e-300)),
            nonHermitian_relative=float(np.linalg.norm(self.B-self.B.conj().T)/max(np.linalg.norm(self.B),1e-300)),
            min_fraction=float(self.a.min()),max_fraction=float(self.a.max()),fraction_ROI_mean=self.amean)
        np.savez_compressed(OUT/f'{parent}_setup.npz',B=self.B,L=self.L,R=self.R,T=self.T,Q=self.Q,LQ=self.LQ,QR=self.QR,x=self.x,fraction=self.meta['fraction'],observation=self.y,Fh=self.Fh)
        book.event('SHAPE_SETUP',**self.setup)
    def forward(self,d,method='schur',derivative=False,diagnostic=False):
        self.book.reserve(16 if derivative else 8,int(method=='direct'),f'{self.parent}:{method}:derivative={derivative}')
        t=time.perf_counter()
        if method=='direct':
            A=self.I-d*self.B;fac=lu_factor(A,check_finite=False);X=lu_solve(fac,self.R,check_finite=False);F=self.Fh+d*(self.L@X)
            if derivative:Y=lu_solve(fac,X,check_finite=False);DF=self.L@Y
        else:
            A=self.I-d*self.T;X=solve_triangular(A,self.QR,check_finite=False);F=self.Fh+d*(self.LQ@X)
            if derivative:Y=solve_triangular(A,X,check_finite=False);DF=self.LQ@Y
        elapsed=time.perf_counter()-t
        diag={}
        if diagnostic:
            rhs=self.R if method=='direct' else self.QR
            diag['relative_residual']=float(np.linalg.norm(A@X-rhs)/max(np.linalg.norm(rhs),1e-300))
            if method=='direct':
                gecon=get_lapack_funcs('gecon',(A,));rcond,info=gecon(fac[0],np.linalg.norm(A,1));diag['condition1_estimate']=float(1/rcond);diag['gecon_info']=int(info)
        return F,(DF if derivative else None),elapsed,diag

def timing(shape,c):
    rows=[]
    for r in range(c['timing_repeats']):
        u=np.array([(r+.5)/30,((r+.5)*((np.sqrt(5)-1)/2))%1]);d=complex(*(LO+(HI-LO)*u));out={}
        for method in (['direct','schur'] if r%2==0 else ['schur','direct']):out[method]=shape.forward(d,method,diagnostic=True)
        ef=float(np.linalg.norm(white(out['direct'][0]-out['schur'][0])))
        row=dict(parent=shape.parent,repeat=r,delta_real=d.real,delta_imag=d.imag,LU_kernel_wall=out['direct'][2],Schur_node_wall=out['schur'][2],white_port_error=ef,**out['direct'][3],Schur_relative_residual=out['schur'][3]['relative_residual']);rows.append(row)
    np.savez_compressed(OUT/f'{shape.parent}_timing_predictions.npz',delta=np.array([complex(r['delta_real'],r['delta_imag']) for r in rows]),last_direct=out['direct'][0],last_Schur=out['schur'][0])
    return rows

def integrate(shape,c):
    saved=[];cache={};nopt=0
    def evalopt(v):
        nonlocal nopt
        key=tuple(v)
        if key not in cache:
            nopt+=1
            if nopt>100:raise RuntimeError('OPTIMIZER_STATE_CAP')
            f,df,_,_=shape.forward(complex(*v),derivative=True);rr=white(f-shape.y);jj=np.column_stack((white(df),white(1j*df)));cache[key]=(rr,jj)
        return cache[key]
    x0=np.array([shape.anchor_delta.real,shape.anchor_delta.imag]);opt=least_squares(lambda v:evalopt(v)[0],x0,jac=lambda v:evalopt(v)[1],bounds=(LO,HI),x_scale='jac',ftol=1e-11,xtol=1e-11,gtol=1e-9,max_nfev=80)
    mode=opt.x;rr,jj=evalopt(mode);shift=float(.5*rr@rr);H=jj.T@jj
    # Full objective Hessian by central differences of analytic gradient, without GT.
    steps=np.maximum(1e-7,1e-3/np.sqrt(np.diag(H)))
    HH=np.zeros((2,2))
    for j in range(2):
        plus=mode.copy();minus=mode.copy();plus[j]=min(HI[j]-1e-12,mode[j]+steps[j]);minus[j]=max(LO[j]+1e-12,mode[j]-steps[j]);rp,jp=evalopt(plus);rm,jm=evalopt(minus);HH[:,j]=(jp.T@rp-jm.T@rm)/(plus[j]-minus[j])
    HH=(HH+HH.T)/2;ee,QQ=np.linalg.eigh(HH);hessian_used='full objective' if np.min(ee)>0 else 'Gauss Newton fallback';H=HH if np.min(ee)>0 else H;cov=np.linalg.inv(H);sd=np.sqrt(np.diag(cov));cl=np.maximum(LO,mode-8*sd);ch=np.minimum(HI,mode+8*sd)
    count=0
    def node(v,w,method,level,region):
        nonlocal count
        count+=1
        if count>c['max_conditional_nodes_per_parent']:raise RuntimeError('CONDITIONAL_NODE_CAP')
        f,_,_,_=shape.forward(complex(*v));phi=float(.5*np.linalg.norm(white(f-shape.y))**2)
        saved.append((v[0],v[1],w,phi,method,level,region));return phi
    def moments(nodes,weights,phis,region):
        v=np.asarray(nodes);w=np.asarray(weights);phi=np.asarray(phis);lw=np.log(w)+shift-phi-np.log(.15);z=logsumexp(lw);pn=np.exp(lw-z)
        mean=pn@v;var=pn@((v-mean)**2);events=pn@np.column_stack((v[:,0]>1.3,v[:,1]>.105))
        return dict(log_evidence=float(z-shift),scaled_evidence=float(np.exp(z)),amplitude_mean=mean,amplitude_variance=var,event_probability=events,material_ROI_mean=[HOST.real+shape.amean*mean[0],HOST.imag+shape.amean*mean[1]],tail_probability=float(np.sum(pn[np.asarray(region)!='core'])))
    axis=[]
    for lev,(nc,nt) in enumerate([(8,2),(12,3),(20,5)]):
        xs=[LO[0],cl[0],ch[0],HI[0]];ys=[LO[1],cl[1],ch[1],HI[1]];nodes=[];weights=[];phis=[];regions=[]
        for ix in range(3):
            for iy in range(3):
                if xs[ix+1]<=xs[ix] or ys[iy+1]<=ys[iy]:continue
                reg='core' if ix==1 and iy==1 else 'tail';n=nc if reg=='core' else nt;x,w=leggauss(n);xv=xs[ix]+(xs[ix+1]-xs[ix])*(x+1)/2;yv=ys[iy]+(ys[iy+1]-ys[iy])*(x+1)/2
                for i in range(n):
                    for j in range(n):
                        v=np.array([xv[i],yv[j]]);ww=w[i]*w[j]*(xs[ix+1]-xs[ix])*(ys[iy+1]-ys[iy])/4;nodes.append(v);weights.append(ww);regions.append(reg);phis.append(node(v,ww,'original_GL',lev,reg))
        row=moments(nodes,weights,phis,regions);row.update(method='original_GL',core_order=nc,tail_order=nt,node_count=len(nodes));axis.append(row)
    # Independent rotated Hessian coordinates, Gaussian-weighted rule; illegal nodes zero.
    eig,Q=np.linalg.eigh(H);transform=Q@np.diag(1/np.sqrt(eig));gh=[]
    for n in [8,12]:
        zz,ww=hermgauss(n);nodes=[];weights=[];phis=[];regions=[]
        for i in range(n):
            for j in range(n):
                z=np.sqrt(2)*np.array([zz[i],zz[j]]);v=mode+transform@z
                if np.any(v<LO) or np.any(v>HI):continue
                w=ww[i]*ww[j]*2*np.linalg.det(transform)*np.exp(float(.5*z@z));nodes.append(v);weights.append(w);regions.append('core');phis.append(node(v,w,'rotated_GH',n,'local_transformed'))
        row=moments(nodes,weights,phis,regions);row.update(method='rotated_GH',order=n,node_count=len(nodes));gh.append(row)
    fine=axis[-1];coarse=axis[-2];independent=gh[-1]
    diff=dict(relative_evidence_axis=float(abs(np.expm1(coarse['log_evidence']-fine['log_evidence']))),relative_evidence_cross=float(abs(np.expm1(independent['log_evidence']-fine['log_evidence']))),amplitude_mean_axis=float(np.max(abs(coarse['amplitude_mean']-fine['amplitude_mean']))),amplitude_mean_cross=float(np.max(abs(independent['amplitude_mean']-fine['amplitude_mean']))),event_axis=float(np.max(abs(coarse['event_probability']-fine['event_probability']))),event_cross=float(np.max(abs(independent['event_probability']-fine['event_probability']))))
    tol=c['numeric_tolerances'];passed=max(diff['relative_evidence_axis'],diff['relative_evidence_cross'])<tol['relative_evidence'] and max(diff['amplitude_mean_axis'],diff['amplitude_mean_cross'])<tol['amplitude_mean'] and max(diff['event_axis'],diff['event_cross'])<tol['event_probability']
    result=dict(parent=shape.parent,status='NUMERICAL_CONDITIONAL_REFERENCE_PASS' if passed else 'REFERENCE_UNQUALIFIED',scope='original Uniform2 conditional on exposed fixedGaussian64 shape, not Full66 marginal',mode=mode,profile_phi=shift,profile_material_ROI_mean=[HOST.real+shape.amean*mode[0],HOST.imag+shape.amean*mode[1]],optimizer_success=bool(opt.success),optimizer_states=nopt,Hessian=H,Hessian_type=hessian_used,marginal_SD=sd,core_box=[cl,ch],axis_references=axis,independent_references=gh,differences=diff,conditional_node_count=count,strict_certificate=False,full_support_numeric='original-coordinate rules cover complete rectangle; tail convergence numerical only',tail_certificate=False,prior_density=1/.15,likelihood='exp(-||F-y||^2/sigma^2)',event_thresholds=[1.3,.105],GT_access=False)
    write(OUT/f'{shape.parent}_conditional.json',result)
    with (OUT/f'{shape.parent}_integration_nodes.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['delta_real','delta_imag','integration_weight','phi','method','level','region']);writer.writerows(saved)
    np.savez_compressed(OUT/f'{shape.parent}_conditional_arrays.npz',mode=mode,Hessian=H,amplitude_mean=fine['amplitude_mean'],amplitude_variance=fine['amplitude_variance'],material_mean=HOST+shape.meta['fraction']*(fine['amplitude_mean'][0]+1j*fine['amplitude_mean'][1]),event_probability=fine['event_probability'])
    return result

def run():
    c=json.loads((P/'config/transfer_frozen.json').read_text());OUT.mkdir(parents=True,exist_ok=True)
    lock=(OUT/'physics_worker.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);book=Book(c);status='FAILED';error=None
    try:
        t=time.perf_counter()
        with np.load(c['host_setup_cache']) as z:setup={k:z[k].copy() for k in z.files}
        cfg=json.loads((R4/'research/p3l_r4_conditional_transport/configs/TARGET.json').read_text())
        with np.load(R4/'research/p3l_r4_conditional_transport/inputs/background/weak.npz') as z:points=z['points'].copy()
        prior=MaterialPrior(cfg,points);load_wall=time.perf_counter()-t
        book.event('HISTORICAL_CACHE_REUSE',path=c['host_setup_cache'],load_wall=load_wall,historical_RHS576=264,historical_LU576=1,new_RHS576=0)
        timing_rows=[];conditional=[];setups=[]
        for parent in c['parents']:
            shape=Shape(parent,setup,prior,book);setups.append(shape.setup);timing_rows+=timing(shape,c)
            if parent in c['conditional_parents']:conditional.append(integrate(shape,c))
            write(OUT/'PARTIAL.json',dict(completed_parents=[s['parent'] for s in setups],conditionals=conditional,RHS256=book.rhs))
        with (OUT/'TIMINGS.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(timing_rows[0]));w.writeheader();w.writerows(timing_rows)
        groups=[]
        for s in setups:
            rows=[r for r in timing_rows if r['parent']==s['parent']];d=np.array([r['LU_kernel_wall'] for r in rows]);q=np.array([r['Schur_node_wall'] for r in rows]);fixed=s['Schur_transforms_wall'];delta=float(np.median(d)-np.median(q));be=int(np.floor(fixed/delta)+1) if delta>0 else None
            groups.append(dict(parent=s['parent'],LU_p50=float(np.median(d)),LU_range=[float(d.min()),float(d.max())],Schur_node_p50=float(np.median(q)),Schur_node_range=[float(q.min()),float(q.max())],Schur_setup_wall=fixed,break_even_M_p50=be,M_comparison=[dict(M=m,direct_p50_total=s['formation_wall']+m*float(np.median(d)),Schur_p50_total=s['formation_wall']+fixed+m*float(np.median(q)),Schur_cold_amortized_pernode=(s['formation_wall']+fixed)/m+float(np.median(q))) for m in c['M']]))
        result=dict(status='COMPLETE',timing_groups=groups,setups=setups,conditional_results=conditional,cache_load_wall=load_wall,max_white_port_error=max(r['white_port_error'] for r in timing_rows),max_condition1_estimate=max(r['condition1_estimate'] for r in timing_rows),empirical_reduced_order='NOT_RUN; exact256 transfer already changes amortized query cost; no uncharged lower-order scan',certified_reduced_order=None,RHS256=book.rhs,LU256=book.lu,RHS576=0,scope='8 exposed DEV shapes, 30 legitimate amplitudes per shape; 3 fixed-shape originalproperty conditional integrations; no population/general posterior claims')
        write(OUT/'RESULTS.json',result);status='COMPLETE';print(json.dumps(plain(dict(status=status,RHS256=book.rhs,LU256=book.lu,max_white_port_error=result['max_white_port_error'],conditional_statuses=[r['status'] for r in conditional]))),flush=True)
    except BaseException as e:
        error=repr(e);write(OUT/'FAILURE.json',dict(error=error));raise
    finally:book.close(status,error);lock.close()

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('command',choices=['freeze','run']);ns=a.parse_args();freeze() if ns.command=='freeze' else run()
