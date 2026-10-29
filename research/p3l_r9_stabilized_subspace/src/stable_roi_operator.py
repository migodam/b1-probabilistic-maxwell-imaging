"""Exact fixed-host ROI block elimination; ordinary strong solver control.

All P1 nonuniform ROI materials remain supported. This is standard Schur /
Woodbury algebra, not a novel physics identity or a posterior claim.
Only run after freeze, and while holding the shared physical worker lock.
"""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[key]='1'
import argparse, fcntl, hashlib, json, resource, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.special import ndtr

P=Path(__file__).resolve().parents[1]
B1=P.parents[2]
R4=B1/'b1-p3l-r4-conditional-transport'
H7=B1/'b1-p3l-r7-material-fibers/research/p3l_r7_material_solution_fibers'
OUT=P/'results/stable_roi'
SIGMA=.00119211798514519
sys.path.insert(0,str(R4/'src'))
from p1_material.prior import MaterialPrior

def write(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def append(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('a') as f:
        f.write(json.dumps(obj,ensure_ascii=False,allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
def pack(a):return np.concatenate((a.real.T,a.imag.T),axis=1).ravel()
def white(a):return pack(a)*np.sqrt(2)/SIGMA
def inputs():
    with np.load(R4/'research/p3l_r4_conditional_transport/inputs/background/weak.npz') as z:
        bg={k:z[k].copy() for k in z.files}
    cfg=json.loads((R4/'research/p3l_r4_conditional_transport/configs/TARGET.json').read_text())
    return bg,cfg
def context_paths():
    paths=[]
    for i in range(8):
        parent=f'R7_LOCKED_{i:03d}';hp=H7/'E2/parents'/parent
        paths.append((parent,'anchor',hp/'ANCHOR.npz'))
        paths.append((parent,'FullJ_endpoint',hp/'B_FULLJ/RMS_0.0200/ENDPOINT.npz'))
    return paths

def freeze():
    OUT.mkdir(parents=True,exist_ok=True)
    # Use only paths resolved from actual R7 records; no GT files are inspected.
    resolved=[]
    for parent,role,path in context_paths():
        if not path.exists():
            candidates=sorted((H7/'E2/parents'/parent).glob('**/ENDPOINT.npz'))
            candidates=[q for q in candidates if 'B_FULLJ' in str(q) and ('0.020' in str(q) or '0.02' in str(q))]
            if len(candidates)!=1:raise RuntimeError((parent,role,[str(x) for x in candidates]))
            path=candidates[0]
        resolved.append(dict(parent=parent,role=role,path=str(path),blake2b=hashlib.blake2b(path.read_bytes()).hexdigest()))
    c=dict(name='FIXED_HOST_EXACT_ROI_STRONG_SOLVER_CONTROL',created=time.time(),
           purpose='after development identify cost location; all-arms-available exact solver control, not an advantage assigned to frozen Q',
           history='16 already exposed R7 material states, 8 independent parent identities but DEV only',
           contexts=resolved,repeats=3,full_RHS_cap=2184,full_setup_RHS=264,
           direct_per_evaluation_RHS=40,reduced_per_evaluation_RHS256=40,
           expected_new_full_RHS=264+16*3*40,
           original_dtype='float64/complex128',max_wall_seconds=300,BLAS_threads=1,
           quality_white_prediction_max=1e-8,quality_relative_J_max=1e-9,
           cost_window='kernel setup+forward+J; equal IO recorded outside kernels; no historical timing ratios',
           operators='general nonuniform diagonal ROI increment; full nonlinear exact block elimination',
           physics_identity='P1 nonzero host outside and inside ROI; source-major observations',
           success='>=20% cold setup-inclusive same-quality wall saving; DEV finite benchmark only')
    path=P/'config/stable_roi_freeze.json'
    if path.exists():raise RuntimeError('freeze exists')
    write(path,c);print('FROZEN',len(resolved),c['expected_new_full_RHS'])

class Ledger:
    def __init__(self):
        self.n=0;self.full=0;self.lu_full=0;self.small=0;self.lu_small=0
        self.path=OUT/'COST_LEDGER.jsonl'
        if self.path.exists():raise RuntimeError('no automatic replay')
    def row(self,event,full_RHS=0,full_LU=0,reduced_RHS256=0,reduced_LU256=0,**kw):
        if self.full+full_RHS>2184:raise RuntimeError('physical cap')
        self.n+=1;self.full+=full_RHS;self.lu_full+=full_LU;self.small+=reduced_RHS256;self.lu_small+=reduced_LU256
        append(self.path,dict(sequence=self.n,utc=time.time(),event=event,full_RHS=full_RHS,full_LU=full_LU,
             reduced_RHS256=reduced_RHS256,reduced_LU256=reduced_LU256,**kw,
             totals=dict(full_RHS=self.full,full_LU=self.lu_full,reduced_RHS256=self.small,reduced_LU256=self.lu_small)))

def run():
    c=json.loads((P/'config/stable_roi_freeze.json').read_text())
    if not (P/'config/stable_roi_physics_enable.json').exists():raise RuntimeError('not enabled')
    lock=(P/'shared/physics_worker.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    book=Ledger();started=time.perf_counter();cpu=time.process_time();rows=[]
    book.row('JOB_BEGIN',pid=os.getpid(),contract=c['name'])
    try:
        bg,cfg=inputs();G,S,Einc=bg['G'],bg['S'],bg['incident']
        prior=MaterialPrior(cfg['prior'] if 'prior' in cfg else cfg,bg['points'])
        idx=np.where(prior.mask)[0];host=prior.host
        book.row('SETUP_RESERVE',full_RHS=264,full_LU=1)
        t=time.perf_counter();Ah=np.eye(576)-host*G;fac=lu_factor(Ah,check_finite=False)
        Eh=lu_solve(fac,Einc,check_finite=False);B=lu_solve(fac,G[:,idx],check_finite=False)
        C=B[idx];ER=Eh[idx];O=S[:,idx]+host*S@B;Fh=S@(host*Eh)
        setup=time.perf_counter()-t
        residuals=dict(host=float(np.linalg.norm(Ah@Eh-Einc)/np.linalg.norm(Einc)),basis=float(np.linalg.norm(Ah@B-G[:,idx])/np.linalg.norm(G[:,idx])))
        del fac,Ah
        np.savez_compressed(OUT/'SETUP.npz',G=G,S=S,Einc=Einc,Eh=Eh,B=B,C=C,ER=ER,O=O,Fh=Fh,ROI=idx)
        book.row('SETUP_END',wall_seconds=setup,residuals=residuals)
        for repeat in range(c['repeats']):
            for ci,item in enumerate(c['contexts']):
                with np.load(item['path']) as z:a={k:z[k].copy() for k in z.files}
                x=a['x'];q=x.copy();q[64:]=ndtr(x[64:]);m=prior.decode(q)[0]
                dm=prior.derivative(q)*np.r_[np.ones(64),np.exp(-x[64:]**2/2)/np.sqrt(2*np.pi)][None,:]
                dr=dm[idx];delta=m[idx]-host
                out={};times={};resid={}
                order=['direct','reduced'] if (repeat+ci)%2==0 else ['reduced','direct']
                for method in order:
                    book.row('EVALUATION_RESERVE',full_RHS=40 if method=='direct' else 0,full_LU=int(method=='direct'),
                             reduced_RHS256=40 if method=='reduced' else 0,reduced_LU256=int(method=='reduced'),
                             parent=item['parent'],role=item['role'],repeat=repeat,method=method)
                    t=time.perf_counter()
                    if method=='direct':
                        L=np.eye(576)-m[:,None]*G;fac=lu_factor(L,check_finite=False)
                        j=lu_solve(fac,m[:,None]*Einc,check_finite=False)
                        field=Einc+G@j;f=S@j
                        transfer=lu_solve(fac,S.T,trans=1,check_finite=False).T
                        deriv=np.einsum('ri,is,ik->rsk',transfer,field,dm,optimize=True)
                        rr=np.linalg.norm(L@j-m[:,None]*Einc)/np.linalg.norm(m[:,None]*Einc)
                    else:
                        L=np.eye(256)-delta[:,None]*C;fac=lu_factor(L,check_finite=False)
                        j=lu_solve(fac,delta[:,None]*ER,check_finite=False)
                        field=ER+C@j;f=Fh+O@j
                        transfer=lu_solve(fac,O.T,trans=1,check_finite=False).T
                        deriv=np.einsum('ri,is,ik->rsk',transfer,field,dr,optimize=True)
                        rr=np.linalg.norm(L@j-delta[:,None]*ER)/max(np.linalg.norm(delta[:,None]*ER),1e-300)
                    J=np.column_stack([white(deriv[:,:,k]) for k in range(66)])
                    elapsed=time.perf_counter()-t;del fac,L
                    out[method]=(f,J);times[method]=elapsed;resid[method]=float(rr)
                    book.row('EVALUATION_END',method=method,parent=item['parent'],role=item['role'],repeat=repeat,wall_seconds=elapsed,backward_residual=float(rr))
                ef=float(np.linalg.norm(white(out['direct'][0]-out['reduced'][0])))
                ej=float(np.linalg.norm(out['direct'][1]-out['reduced'][1])/np.linalg.norm(out['direct'][1]))
                reference=float(np.linalg.norm(white(out['direct'][0]-a['prediction'])))
                row=dict(parent=item['parent'],role=item['role'],repeat=repeat,wall=times,backward=resid,white_prediction_difference=ef,relative_J_difference=ej,
                         historical_white_prediction_difference=reference,pass_quality=ef<=c['quality_white_prediction_max'] and ej<=c['quality_relative_J_max'])
                rows.append(row)
                t=time.perf_counter();np.savez_compressed(OUT/f"{item['parent']}_{item['role']}_r{repeat}.npz",x=x,material=m,prediction_direct=out['direct'][0],prediction_reduced=out['reduced'][0],J_direct=out['direct'][1],J_reduced=out['reduced'][1]);savewall=time.perf_counter()-t
                book.row('SAVE',wall_seconds=savewall);write(OUT/'PARTIAL.json',dict(rows=rows))
        td=sum(r['wall']['direct'] for r in rows);tr=sum(r['wall']['reduced'] for r in rows)
        summary=dict(status='COMPLETE',quality_pass=all(r['pass_quality'] for r in rows),rows=rows,setup_wall=setup,
                     pooled_direct_wall=td,pooled_reduced_wall=tr,setup_inclusive_saving=1-(setup+tr)/td,
                     marginal_saving=1-tr/td,full_RHS=book.full,reduced_RHS256=book.small,full_LU=book.lu_full,reduced_LU256=book.lu_small,
                     peak_RSS=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                     scope='8 exposed DEV parents, 16 fixed material states, 3 rotated timing repeats; ordinary exact block elimination, no population speed claim',
                     setup_inclusive_count_basis='one acquisition/background setup across all candidates; also report per observation amortization')
        write(OUT/'RESULTS.json',summary);print(json.dumps({k:v for k,v in summary.items() if k!='rows'}))
        book.row('JOB_END',status='COMPLETE',wall_seconds=time.perf_counter()-started,CPU_seconds=time.process_time()-cpu)
    except BaseException as exc:
        book.row('JOB_END',status='FAILED',error=repr(exc),wall_seconds=time.perf_counter()-started,CPU_seconds=time.process_time()-cpu)
        write(OUT/'FAILURE.json',dict(error=repr(exc)));raise
    finally:lock.close()

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('command',choices=['freeze','run']);ns=a.parse_args()
    freeze() if ns.command=='freeze' else run()
