"""Exploratory R12 Maxwell normal-memory mechanism probe, not blind confirmation.

Uses the exact published projectors and already saved physical starts.
Predicts the first constitutive projection switch along an AP normal ray,
then tests real DR continuation and the limitations of that prediction.
No GT enters initialization, projection, thresholds or candidate selection.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[key]='1'
import sys,json,time,argparse,importlib.util
from pathlib import Path
import numpy as np
from scipy.stats import chi2
from scipy.linalg import lu_factor,lu_solve

ROOT=Path('/workspace/scratch/759c14241356')
REPO=ROOT/'research_repo'
BASE=REPO/'research/research_max_struct_r12/supplement_paired'
path=BASE/'runtime/field_current.py'
spec=importlib.util.spec_from_file_location('published_field_current',path)
fc=importlib.util.module_from_spec(spec);sys.modules[spec.name]=fc;spec.loader.exec_module(fc)

def pnorm(X,Q): return float(np.hypot(np.linalg.norm(X),np.linalg.norm(Q)))

def thresholds(cx,cq,dx,dq,m,scale=6.,lo=.02,hi=20.):
    theta=np.arctan(m/scale);s=np.sin(theta)[:,None];c=np.cos(theta)[:,None]
    v=s*cx+c*cq; w=c*dx-s*dq
    a=np.sum(abs(v)**2,axis=1); b=np.sum(abs(w)**2,axis=1)
    cross=np.real(np.sum(v.conj()*w,axis=1))
    tang=s*dx+c*dq
    lower=m<=lo+1e-9; upper=m>=hi-1e-9; interior=~(lower|upper)
    crit=np.full(len(m),np.inf)
    ok=b>1e-30
    crit[ok&interior]=np.sqrt(a[ok&interior]/b[ok&interior])
    for active,other in ((lower,np.arctan(hi/scale)),(upper,np.arctan(lo/scale))):
        use=active&ok; cc=cross[use]/np.tan(other-theta[use]);rad=np.sqrt(cc*cc+a[use]*b[use])
        crit[use]=(rad-cc)/b[use]
    return dict(tcrit=crit,a=a,b=b,cross=cross,interior=interior,lower=lower,upper=upper,
                tangent_normal_residual=np.linalg.norm(tang,axis=1))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--case',required=True);ap.add_argument('--steps',type=int,default=3000)
    ap.add_argument('--follow',type=int,default=500);args=ap.parse_args()
    out=ROOT/'work/normal_memory'/args.case;out.mkdir(parents=True,exist_ok=True)
    cpu=time.process_time();wall=time.perf_counter()
    pub=dict(np.load(BASE/'inputs'/(args.case+'.npz')))
    initial=dict(np.load(BASE/'RESULTS'/(args.case+'__ap')/'initial.npz'))
    G,E,S,Y=pub['G'],pub['E0'],pub['S'],pub['Y']
    delta=float(pub['sigma'])*np.sqrt(chi2.ppf(.99,2*Y.size)/2)
    P=fc.DenseFieldCurrentProjectors(G,E,S,Y,delta,scale=6.,lower=.02,upper=20.,grouped_ports=True)
    X=initial['j']/6.;Q=initial['Q'].copy();history=[]
    for k in range(args.steps):
        cx,cq,m=P.material(X,Q);ax,aq=P.affine_ball(cx,cq)
        if (k+1)%100==0:
            history.append(dict(iteration=k+1,AP_step=pnorm(ax-X,aq-Q),gap=pnorm(ax-cx,aq-cq)))
        X,Q=ax,aq
    # a=(X,Q) belongs to A; c=PC(a) exactly, d=a-c exactly C-normal.
    cx,cq,m=P.material(X,Q);dx,dq=X-cx,Q-cq
    ax2,aq2=P.affine_ball(cx,cq)
    gap=pnorm(dx,dq); pair_error=pnorm(ax2-X,aq2-Q)
    th=thresholds(cx,cq,dx,dq,m)
    crit=th['tcrit']; icrit=int(np.argmin(crit)); tcrit=float(crit[icrit])
    # Below and above the predicted event; this tests observed nonzero material
    # switches and finite DR timing, not an algebra error-floor identity.
    ray=[]
    factors=[0.,.5,1.,max(1.,.95*tcrit),max(1.,.999*tcrit),1.001*tcrit,1.05*tcrit]
    for t in sorted(set(factors)):
        rx,rq,rm=P.material(cx+t*dx,cq+t*dq)
        ray.append(dict(t=t,material_RMS=float(np.sqrt(np.mean((rm-m)**2))),
                        max_material_change=float(np.max(abs(rm-m))),
                        first_cell_material=float(rm[icrit]),shadow_motion=pnorm(rx-cx,rq-cq)))
    # Actual ordinary DR continuation from the AP state versus one proposed
    # event step. Both use the same projectors and full physical audit cadence.
    follow={}; final_arrays={};physics={'LU':0,'RHS':0}
    def score(mm):
        L=np.eye(len(G))-mm[:,None]*G;lu=lu_factor(L,check_finite=False)
        jj=lu_solve(lu,mm[:,None]*E,check_finite=False)
        physics['LU']+=1;physics['RHS']+=E.shape[1]
        rr=S@jj-Y
        return dict(residual_over_delta=float(np.linalg.norm(rr)/delta),material=mm.copy())
    # No inherited prior is used for a new solver comparison; all audits below
    # are mechanism trajectories and no best endpoint is chosen using GT.
    init_score=score(m)
    for arm,t0 in [('ordinary',1.),('event',max(1.,1.001*tcrit))]:
        zx,zq=cx+t0*dx,cq+t0*dq;trace=[];first_change=None
        for k in range(args.follow):
            px,pq,pm=P.material(zx,zq)
            motion=float(np.max(abs(pm-m)))
            if first_change is None and motion>1e-3: first_change=k
            bx,bq=P.affine_ball(2*px-zx,2*pq-zq)
            if k<10 or k%25==0 or k==args.follow-1:
                row={'iteration':k,'material_maxchange':motion,'shadow_gap':pnorm(bx-px,bq-pq),
                     'normal_ray_error':pnorm(zx-(cx+(t0+k)*dx),zq-(cq+(t0+k)*dq))}
                if k%50==0 or k==args.follow-1:
                    sc=score(pm);row.update(residual_over_delta=sc['residual_over_delta'])
                trace.append(row)
            zx,zq=zx+bx-px,zq+bq-pq
        follow[arm]={'t0':t0,'first_material_change_step':first_change,'trace':trace}
        final_arrays[arm+'_m']=pm;final_arrays[arm+'_X']=zx;final_arrays[arm+'_Q']=zq
    np.savez_compressed(out/'states.npz',cX=cx,cQ=cq,dX=dx,dQ=dq,m=m,**th,**final_arrays)
    res=dict(case=args.case,AP_iterations=args.steps,follow_iterations=args.follow,scale=6.,delta=delta,
             gap=gap,AP_pair_error=pair_error,relative_pair_error=pair_error/max(gap,1e-300),
             first_cell=icrit,first_cell_initial_m=float(m[icrit]),first_tcrit=tcrit,
             threshold_quantiles=np.quantile(crit[np.isfinite(crit)],[0,.1,.5,.9,1]).tolist(),
             ray=ray,AP_history=history,follow=follow,
             initial_AP_shadow_residual_over_delta=init_score['residual_over_delta'],
             physical_cost=physics,projection_cost=P.costs,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
             scope='exposed-state mechanism study; not a matched-cost algorithm superiority test; no GT used')
    (out/'result.json').write_text(json.dumps(res,indent=2)+'\n')
    print(json.dumps({k:res[k]for k in ['case','gap','relative_pair_error','first_tcrit','first_cell','initial_AP_shadow_residual_over_delta','cpu_seconds','wall_seconds']},indent=2),flush=True)

if __name__=='__main__': main()
