"""Predeclared block-Krylov transfer comparison, no new materials or GT reads.
Numerical coercivity+full residual gates; not interval-certified physics.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[k]='1'
from pathlib import Path
import json,time,resource,csv,argparse
import numpy as np
from scipy.linalg import schur,solve_triangular,svd,eigh
P=Path(__file__).resolve().parents[1]; OUT=P/'results/reduced_response'; CFG=P/'config/reduced_response_frozen.json'
SIGMA=.00119211798514519

def write(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def freeze():
    if CFG.exists():raise RuntimeError('already frozen')
    write(CFG,dict(id='R10_BILATERAL_COERCIVE_BLOCK_KRYLOV_V1',utc=time.time(),parents=[f'R9V_FRESH_{i:03d}' for i in range(8)],orders=[32,64,96,128],amplitudes='same30 predeclared original-prior amplitude nodes as transfer benchmark',RHS_all_orders_conservative_cap=10000,LU256_cap=0,wall_cap_seconds=180,memory_cap_bytes=2000000000,BLAS_threads=1,white_port_bound_gate=.01,likelihood_potential_bound_gate=.1,rank_relative_cutoff=1e-12,scope='8 exposed DEV anchor shapes, no new parent or GT use; exact fullSchur reference at all30nodes; full residual and numerical coercivity before accept; fallback costs charged via shared reference and separately allocated per-method cost',cost='lower-order RHS counted separately by dimension and conservatively against10000reserve; setups+residual norm+errorbound+fallback included; no continuum/roundingcertificate',selection='first among32,64,96,128 passing both bounds; compare fixedr and progressive adaptive rank, no empirical-error gate'))

class Book:
    def __init__(self,c):
        OUT.mkdir(parents=True,exist_ok=True);self.c=c;self.f=(OUT/'COST_LEDGER.jsonl').open('x');self.start=time.perf_counter();self.cpu=time.process_time();self.n=0;self.rhs={};self.schur=0;self.matvec=0
    def event(self,event,**kw):
        self.n+=1;self.f.write(json.dumps(dict(sequence=self.n,utc=time.time(),event=event,**kw))+'\n');self.f.flush()
    def reserve(self,order,label):
        if sum(self.rhs.values())+8>self.c['RHS_all_orders_conservative_cap']:raise RuntimeError('RHS cap')
        if time.perf_counter()-self.start>self.c['wall_cap_seconds']:raise RuntimeError('wallcap')
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>self.c['memory_cap_bytes']:raise RuntimeError('RSS cap')
        self.rhs[order]=self.rhs.get(order,0)+8;self.event('TRIANGULAR_RHS_ATTEMPT',order=order,RHS=8,label=label)
    def close(self,status,error=None):
        d=dict(status=status,error=error,RHS_by_order=self.rhs,RHS_conservative_total=sum(self.rhs.values()),RHS256=self.rhs.get(256,0),LU256=0,Schur_decompositions=self.schur,Green_matrix_columns_applied=self.matvec,wall_seconds=time.perf_counter()-self.start,CPU_seconds=time.process_time()-self.cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        self.event('JOB_END',**d);write(OUT/'RECEIPT.json',d);self.f.close()

def basis(B,R,rmax,tol,book):
    U,s,_=svd(R,full_matrices=False,check_finite=False);keep=s>tol*s[0];V=U[:,keep];block=V.copy(); widths=[V.shape[1]]
    while V.shape[1]<rmax and block.shape[1]:
        book.matvec+=block.shape[1];W=B@block
        for _ in range(2):W-=V@(V.conj().T@W)
        U,s,_=svd(W,full_matrices=False,check_finite=False)
        if len(s)==0 or s[0]<=tol*np.linalg.norm(B):break
        keep=s>tol*s[0];block=U[:,keep][:,:rmax-V.shape[1]]
        V=np.column_stack((V,block));widths.append(block.shape[1])
    return V,widths

def run():
    c=json.loads(CFG.read_text());b=Book(c);rows=[];setups=[];references=[]
    try:
        for parent in c['parents']:
            a=dict(np.load(P/'results/transfer'/f'{parent}_setup.npz')); ds=np.load(P/'results/transfer'/f'{parent}_timing_predictions.npz')['delta'];B,L,R=a['B'],a['L'],a['R'];Fh=a['Fh'];y=a['observation']
            t=time.perf_counter();T,Q=schur(B,output='complex',check_finite=False);LQ=L@Q;QR=Q.conj().T@R;full_setup=time.perf_counter()-t;b.schur+=1
            t=time.perf_counter();beta=float(eigh((B-B.conj().T)/(2j),eigvals_only=True,subset_by_index=(0,0),check_finite=False)[0]);Lnorm=float(np.linalg.norm(L,'fro'));energy_setup=time.perf_counter()-t
            t=time.perf_counter();V,widths=basis(B,R,max(c['orders']),c['rank_relative_cutoff'],b);basis_setup=time.perf_counter()-t
            models={};cumulative=energy_setup+basis_setup
            if V.shape[1]<min(c['orders']):raise RuntimeError('KRYLOV_RANK_BELOW32: explicit failed probe; do not claim reduced realization')
            for order in c['orders']:
                if V.shape[1]<order:continue
                t=time.perf_counter();Vr=V[:,:order];b.matvec+=order;BV=B@Vr;Hr=Vr.conj().T@BV;Tr,Qr=schur(Hr,output='complex',check_finite=False);b.schur+=1;VrQ=Vr@Qr;BVQ=BV@Qr;LVQ=L@VrQ;RQ=Qr.conj().T@(Vr.conj().T@R);ss=time.perf_counter()-t
                models[order]=(Tr,VrQ,BVQ,LVQ,RQ);cumulative+=ss
                setups.append(dict(parent=parent,order=order,energy_setup=energy_setup,basis_to128_setup=basis_setup,reduced_model_setup=ss,full_Schur_setup=full_setup,beta=beta,bmin=.03/(1.8**2+.03**2)+beta,L_frobenius=Lnorm,basis_effective_rank=int(V.shape[1]),cumulative_model_setup=cumulative))
            Fs=[]; approximations={r:[] for r in models}
            for i,delta in enumerate(ds):
                b.reserve(256,f'{parent}:{i}:full');t=time.perf_counter();uf=solve_triangular(np.eye(256)-delta*T,QR,check_finite=False);f=Fh+delta*(LQ@uf);tf=time.perf_counter()-t;Fs.append(f);references.append(dict(parent=parent,node=i,wall_seconds=tf))
                for order in (list(models) if i%2==0 else list(models)[::-1]):
                    Tr,VrQ,BVQ,LVQ,RQ=models[order];b.reserve(order,f'{parent}:{i}:rank{order}');t=time.perf_counter();u=solve_triangular(np.eye(order)/delta-Tr,RQ,check_finite=False);fa=Fh+LVQ@u;solve_time=time.perf_counter()-t
                    t=time.perf_counter();res=R-(VrQ@u)/delta+BVQ@u;bn=float(delta.imag/abs(delta)**2+beta);bound=Lnorm*np.linalg.norm(res)/bn if bn>0 else 1e100;ew=float(np.sqrt(2)/SIGMA*bound);rw=float(np.sqrt(2)/SIGMA*np.linalg.norm(fa-y));potential=rw*ew+.5*ew*ew;gate_time=time.perf_counter()-t
                    ae=float(np.sqrt(2)/SIGMA*np.linalg.norm(f-fa));qualified=bool(bn>0 and ew<=c['white_port_bound_gate'] and potential<=c['likelihood_potential_bound_gate'])
                    rows.append(dict(parent=parent,node=i,order=order,delta_real=float(delta.real),delta_imag=float(delta.imag),numerical_coercivity=bn,residual_Fnorm=float(np.linalg.norm(res)),white_port_bound=ew,actual_white_port_error=ae,bound_dominates_sampled_error=bool(ae<=ew*(1+1e-8)+1e-10),potential_bound=potential,qualified=qualified,solve_seconds=solve_time,gate_seconds=gate_time,online_seconds=solve_time+gate_time,full_Schur_seconds=tf,with_fallback_seconds=solve_time+gate_time+(0 if qualified else tf)));approximations[order].append(fa)
            np.savez_compressed(OUT/f'{parent}.npz',delta=ds,full_response=np.stack(Fs),**{f'approx_r{r}':np.stack(approximations[r]) for r in models},V128=V)
        for name,data in [('NODE_RESULTS',rows),('SETUP_COSTS',setups),('REFERENCE_COSTS',references)]:
            with (OUT/f'{name}.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
        summary={}
        for r in c['orders']:
            rr=[x for x in rows if x['order']==r];ss=[x for x in setups if x['order']==r]
            if not rr:
                summary[str(r)]=dict(status='NOT_AVAILABLE_KRYLOV_BREAKDOWN',nodes=0);continue
            fallback_parents={x['parent'] for x in rr if not x['qualified']}
            fallback_setup=sum(x['full_Schur_setup'] for x in ss if x['parent'] in fallback_parents)
            summary[str(r)]=dict(fallback_parent_count=len(fallback_parents),fallback_Schur_setup_seconds=fallback_setup,nodes=len(rr),qualified=sum(x['qualified'] for x in rr),max_actual_white_error=max(x['actual_white_port_error'] for x in rr),max_white_bound=max(x['white_port_bound'] for x in rr),sampled_bounds_all_hold=all(x['bound_dominates_sampled_error'] for x in rr),online_seconds=sum(x['online_seconds'] for x in rr),online_with_fallback_seconds=sum(x['with_fallback_seconds'] for x in rr),full_reference_seconds=sum(x['full_Schur_seconds'] for x in rr),setup_conservative_seconds=sum(x['energy_setup']+x['basis_to128_setup']+x['reduced_model_setup'] for x in ss),full_setup_seconds=sum(x['full_Schur_setup'] for x in ss))
        adaptive=[]
        for parent in c['parents']:
            for i in range(30):
                rr=sorted([x for x in rows if x['parent']==parent and x['node']==i],key=lambda x:x['order']);cost=0.;chosen=256
                for x in rr:
                    cost+=x['online_seconds']
                    if x['qualified']:chosen=x['order'];break
                if chosen==256:cost+=rr[0]['full_Schur_seconds']
                adaptive.append(dict(parent=parent,node=i,chosen_order=chosen,online_seconds=cost,exact_reference_seconds=rr[0]['full_Schur_seconds']))
        cold=[]
        for r in c['orders']:
            q=summary[str(r)]
            if q['nodes']==0:continue
            nparent=len(c['parents']);nn=q['nodes'];setup=q['setup_conservative_seconds']+q['fallback_Schur_setup_seconds']
            for M in (1,4,16,30,64,256,1024):
                reduced=setup/nparent+M*q['online_with_fallback_seconds']/nn;full=q['full_setup_seconds']/nparent+M*q['full_reference_seconds']/nn
                cold.append(dict(order=r,M=M,reduced_plus_gate_plus_fallback_seconds=reduced,full_Schur_seconds=full,saving=1-reduced/full,scope='measured mean node cost amortization,not actual new M queries'))
        with (OUT/'COLD_AMORTIZATION.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(cold[0]));w.writeheader();w.writerows(cold)
        write(OUT/'RESULTS.json',dict(cold_amortization=cold,protocol_id=c['id'],fixed_order=summary,adaptive=adaptive,qualification='NUMERICAL_DISCRETE_RESIDUAL_BOUND; not outward-rounded certificate, no unseen-shape or all-amplitude proof from sampled outcomes',source='all8 exposed fixed shapes, same30 priorrectangle nodes; noGT',model='full256 reference retains fullnonlinearVIE',compression='Krylov/Galerkin and scalaramplitude reuse known methods, no new minimality theorem'))
        b.close('COMPLETE')
    except BaseException as e:b.close('FAILED',repr(e));raise

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('action',choices=['freeze','run']);n=a.parse_args();globals()[n.action]()
