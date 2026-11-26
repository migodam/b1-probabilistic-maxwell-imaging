"""Data-only, strictly checked inner ROI peak intervals in the fixed port model."""
import argparse, json, os, sys, time
from pathlib import Path
for k in ('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.linalg import lu_factor,lu_solve
from scipy.sparse.linalg import LinearOperator,gmres
from flint import arb,acb,ctx
from threadpoolctl import threadpool_limits
from verified_wave import P, VerifiedWave, up, fro, opnorm, rd, arb_upper
ctx.prec=128;U=4.5
def cball(z,e=0):
    return acb(arb(float(z.real),e),arb(float(z.imag),e))
def evaluate_curve(f,y,r,eF,e0,eE,col,eC,h,eh,base,high):
    t=arb(float(high))-arb(float(base)); z=1-t*cball(h,eh)
    if not abs(z)>0:return None
    q0=arb(0);q1=arb(0);cross=arb(0)
    for j in range(f.shape[1]):
        ej=cball(e0[j],eE)
        for i in range(f.shape[0]):
            d=cball(f[i,j],eF)-cball(y[i,j]); new=z*d+t*cball(col[i],eC)*ej
            q0+=abs(d)**2;q1+=abs(new)**2;cross+=(d.conjugate()*new).real
    r2=arb(float(r))**2
    result=dict(zero_or_base_ratio_upper=arb_upper(q0.nonnegative_part().sqrt()/r),end_ratio_upper=arb_upper(q1.nonnegative_part().sqrt()/(abs(z)*r)))
    # A real quadratic on [0,1] in the Bernstein basis. A finite subdivision
    # proves the full path including near-interior extrema and denominators.
    A=q0-r2;B=cross-r2*z.real;C=q1-r2*abs(z)**2
    parts=[(A,B,C)];worst=None
    for depth in range(5):
        worst=max(max(a.upper(),b.upper(),c.upper()) for a,b,c in parts)
        if worst<0:break
        nxt=[]
        for a,b,c in parts:
            d=(a+b)/2; e=(b+c)/2; mid=(d+e)/2;nxt.extend([(a,d,mid),(mid,e,c)])
        parts=nxt
    result.update(full_segment_strict=bool(worst<0),quadratic_ratio_upper=arb_upper(worst/r2),subdivisions=2**depth)
    return result
def run(case,roi,root='profiles',replay_accepted=False):
    src=P/'results'/root/case/f'roi_{roi:02}.npz'
    if not src.is_file():return
    out=src.with_name(src.stem+'_inner.json')
    previous=json.loads(out.read_text()) if out.is_file() else None
    if previous and not replay_accepted:return
    if replay_accepted and (not previous or not previous.get('selected')):return
    selected=previous['selected'] if replay_accepted else None
    start=time.perf_counter();cpu=time.process_time()
    with np.load(src) as a:x=a['material'].ravel();mask=a['roi'].ravel();j=a['current']
    with np.load(P/'ports/public/cases'/f'{case}.npz') as public:y=public['Y'];r=float(public['radius']);opid=public['operator_id'].item()
    with np.load(P/'ports/public/operators'/f'{opid}.npz') as a:
        g=VerifiedWave(a);dense=a['G'] if 'G' in a else None
    u=x.max(); active=np.flatnonzero(mask); cmax=float(x[mask].max())
    inv,info=g.inverse(x)
    b=x[:,None]*g.B;ej=g.solve_error(x,j,b,up(rd(1)*u*fro(g.B)),inv)
    f=g.S@j; ef=up(g.snorm*ej+rd(g.n)*g.snorm*fro(j))
    E=g.B+g.apply(j,True)
    eE=up(g.gnorm*ej+g.full_green_err*fro(j)+rd(1)*(fro(g.B)+g.gnorm*fro(j)))
    ratio=up((fro(f-y)+ef+rd(2)*(fro(f)+fro(y)))/r)
    if ratio>1:
        print(case,roi,'base replay failed',ratio,flush=True)
    if dense is not None:
        A=np.eye(g.n)-x[:,None]*dense;lu=lu_factor(A,check_finite=False)
        # An adjoint identifies weakly read cells; select before reading labels.
        rx=lu_solve(lu,g.S.conj().T,trans=2,check_finite=False).conj().T
    else:
        mv=LinearOperator((g.n,g.n),matvec=lambda v:v-x*g.apply(v),dtype=complex)
        # Strict correctness uses the residual, independently of GMRES success.
        # The finite native prescreen uses field magnitude and uniform cell coverage.
        rx=None
    if rx is not None: score=np.linalg.norm(rx[:,active],axis=0)*np.linalg.norm(E[active],axis=1)
    else:score=np.linalg.norm(E[active],axis=1)
    ntry=min(20,len(active))
    ids=np.unique(np.r_[active[np.argsort(score,kind='stable')[:ntry]],active[np.linspace(0,len(active)-1,ntry,dtype=int)]])
    delta=np.eye(g.n,dtype=float)[:,ids] if dense is not None else np.zeros((g.n,len(ids)))
    if selected is not None:
        ids=np.array([selected['pixel']],dtype=int)
        delta=np.zeros((g.n,1));delta[ids[0],0]=1
    if dense is not None:kern=lu_solve(lu,delta,check_finite=False)
    else:
        delta[ids,np.arange(len(ids))]=1
        kern=np.column_stack([gmres(mv,col,rtol=3e-13,restart=80,maxiter=50)[0] for col in delta.T])
    gc=g.apply(kern,True); fc=g.S@kern
    best=None; rows=[]; seed=ratio<=1
    # q=4.5 is tested first; if it fails, attempt a literal feasible improvement.
    limits=np.unique(np.r_[U, np.linspace(cmax,U,7)])[::-1]
    if selected is not None:limits=[selected['cell_to']]
    for high in limits:
        for k,i in enumerate(ids):
            if high<x[i]:continue
            z=kern[:,k:k+1]; defect=z[:,0]-x*gc[:,k]-delta[:,k]
            ez=up(inv*(fro(defect)+u*g.full_green_err*fro(z)+rd(4)*(fro(z)+u*fro(gc[:,k])+1)))
            ec=up(g.snorm*ez+rd(g.n)*g.snorm*fro(z))
            eh=up(g.gnorm*ez+g.full_green_err*fro(z))
            rec=evaluate_curve(f,y,r,ef,E[i],eE,fc[:,k],ec,gc[i,k],eh,x[i],high)
            if rec is None:continue
            rec.update(pixel=int(i),peak_from=cmax,peak_to=float(max(cmax,high)),cell_from=float(x[i]),cell_to=float(high))
            rows.append(rec)
            if rec['end_ratio_upper']<1 and rec['full_segment_strict'] and seed:
                best=rec;break
        if best is not None:break
    meta=dict(case=case,roi=roi,mask_cells=len(active),strict_feasible_base=seed,base_residual_ratio_upper=ratio,
              attainable_peak_from=cmax if seed else None,attainable_peak_to=best['peak_to'] if best else (cmax if seed else None),
              exact_full_prior_peak_range=bool(best and cmax==0 and best['peak_to']==U),selected=best,
              proof=info,precision=128,represented_operator='stored binary64 spectrum + exact low-rank loading + stored B,S,Y',
              selected_cell_count=len(ids),tested=len(rows),y_error=ef,field_error=eE,
              seconds_wall=time.perf_counter()-start,seconds_cpu=time.process_time()-cpu)
    if selected is not None:
        if best is None or best['peak_from']!=selected['peak_from'] or best['peak_to']!=selected['peak_to']:
            raise ArithmeticError('Independent replay could not recover the original feasible segment')
        keys=('selected_cell_count','tested','seconds_wall','seconds_cpu')
        meta['accepted_point_replay']={k:meta[k] for k in keys}
        for k in keys:meta[k]=previous[k]
    out.write_text(json.dumps(meta,indent=2,allow_nan=False))
    print(case,roi,ratio,meta['attainable_peak_from'],meta['attainable_peak_to'],len(rows),round(meta['seconds_wall'],2),flush=True)
if __name__=='__main__':
    pa=argparse.ArgumentParser();pa.add_argument('cases',nargs='+');pa.add_argument('--kind',default='profiles');pa.add_argument('--roi',type=int,nargs='*');pa.add_argument('--replay-accepted',action='store_true');args=pa.parse_args()
    with threadpool_limits(1):
        for case in args.cases:
            for pth in sorted((P/'results'/args.kind/case).glob('roi_??.npz')):
                ri=int(pth.stem[-2:])
                if args.roi is None or ri in args.roi:run(case,ri,args.kind,args.replay_accepted)
