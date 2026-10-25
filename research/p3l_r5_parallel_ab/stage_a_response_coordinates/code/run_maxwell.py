"""A_v1: counted original-P1 two-dimensional cubature and response integration.
No full66 posterior claims. Full likelihood is evaluated at every paid endpoint.
"""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,time,gc,traceback,argparse
from collections import OrderedDict
import numpy as np
from scipy.special import ndtr,ndtri,logsumexp
from scipy.optimize import least_squares
from scipy.stats import multivariate_normal
from numpy.polynomial.hermite import hermgauss
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT))
from research.p3l_r5_parallel_ab.shared.physics_adapter import Target,save,write,SIGMA,BudgetStop
OUT=Path(__file__).resolve().parents[1]
PAIRS=[(0,1),(0,64),(12,27)]
LOGPRIOR_NORM=np.log(2*np.pi)+2*np.log(ndtr(4)-ndtr(-4))
def plain(v):
    if isinstance(v,np.ndarray):return plain(v.tolist())
    if isinstance(v,np.generic):return plain(v.item())
    if isinstance(v,complex):return [v.real,v.imag]
    if isinstance(v,dict):return {str(k):plain(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [plain(x) for x in v]
    return v

def report(path,d):write(path,plain(d))
def quadnodes(n):
    t,w=hermgauss(n);i,j=np.meshgrid(np.arange(n),np.arange(n),indexing='ij')
    return np.sqrt(2)*np.c_[t[i.ravel()],t[j.ravel()]],(w[i]*w[j]/np.pi).ravel()
def gaussian_log(x,m,C):
    x=np.asarray(x);v=x-m;sg,ld=np.linalg.slogdet(C)
    if sg<=0:raise ValueError('SPD covariance required')
    return -.5*np.einsum('...i,ij,...j->...',v,np.linalg.inv(C),v)-np.log(2*np.pi)-.5*ld
class Model:
    def __init__(self,target,pair):
        self.t=target;self.pair=pair;self.V=np.eye(66)[:,pair];self.cache={};self.last=None;self.calls=0;self.jcalls=0
        # Explicit restart reuse of only own paid compact endpoints. No LU assumed.
        for path in sorted((target.jobdir/'endpoints').glob('*.npz')):
            with np.load(path,allow_pickle=False) as a:
                x=a['x'];z=x[list(pair)];inactive=x.copy();inactive[list(pair)]=0
                if np.any(inactive!=0):continue
                key=tuple(np.round(z,14));self.cache[key]=dict(z=z.copy(),f=target.white(a['prediction'])*SIGMA,chi=a['chi'].copy(),key=path.stem)
        for path in sorted((target.jobdir/'derivatives').glob('*.npz')):
            with np.load(path,allow_pickle=False) as a:
                if 'J' not in a or a['J'].shape!=(512,2) or not np.array_equal(a['V'],self.V):continue
                z=a['x'][list(pair)];key=tuple(np.round(z,14))
                if key in self.cache:self.cache[key]['J']=a['J']*SIGMA
    def full(self,z,purpose):
        z=np.asarray(z,float);key=tuple(np.round(z,14))
        if key in self.cache:return self.cache[key]
        p=self.t.evaluate(self.V@z,purpose=purpose);self.last=(key,p);self.calls+=1
        r=dict(z=z.copy(),f=self.t.white(p.prediction)*SIGMA,chi=p.chi.copy(),key=p.key)
        self.cache[key]=r
        return r
    def jac(self,z,purpose,V=None):
        z=np.asarray(z,float);key=tuple(np.round(z,14));r=self.full(z,purpose)
        if V is None and 'J' in r:return r['J']
        if self.last is None or self.last[0]!=key:
            # Re-evaluation counted honestly; compact response cache is not a free LU.
            p=self.t.evaluate(self.V@z,purpose=purpose+'_state_reacquire');self.last=(key,p);self.calls+=1
        else:p=self.last[1]
        J=self.t.jacobian(p,self.V if V is None else self.V@V,purpose=purpose)*SIGMA;self.jcalls+=1
        if V is None:r['J']=J
        return J
    def close(self):self.last=None;gc.collect()

def p_log(z,f,d,sig):return -.5*float(z@z)-LOGPRIOR_NORM-.5*float((f-d)@(f-d))/sig**2

def discover(model,grid,d,sig):
    costs=np.array([p_log(r['z'],r['f'],d,sig) for r in grid]);order=np.argsort(costs)[::-1]
    seeds=[]
    for k in order:
        z=grid[k]['z']
        if all(np.linalg.norm(z-q)>.8 for q in seeds):seeds.append(z)
        if len(seeds)==3:break
    modes=[];attempts=[]
    for ii,z0 in enumerate(seeds):
        memo={}
        def fun(z):
            r=model.full(z,'mode_residual');return np.r_[z,(r['f']-d)/sig]
        def jac(z):return np.r_[np.eye(2),model.jac(z,'mode_direct_two_JVP')/sig]
        sol=least_squares(fun,z0,jac=jac,bounds=(-3.999999,3.999999),ftol=1e-9,xtol=1e-10,gtol=1e-7,max_nfev=24)
        z=sol.x;r=model.full(z,'mode_final');J=model.jac(z,'mode_final_J')/sig;H=np.eye(2)+J.T@J;C=np.linalg.inv(H)
        value=p_log(z,r['f'],d,sig)
        attempts.append(dict(seed=z0,mode=z,log_target=value,status=int(sol.status),optimality=sol.optimality,nfev=sol.nfev))
        if not any(np.linalg.norm(z-m['z'])<1e-3 for m in modes):modes.append(dict(z=z.copy(),C=C,J=J,r=r,logp=value,optimality=sol.optimality))
    best=max(m['logp'] for m in modes);modes=[m for m in modes if m['logp']>=best-20]
    modes=sorted(modes,key=lambda m:-m['logp'])[:2]
    return modes,attempts

def original(model,modes,d,sig,n):
    comps=[dict(z=np.zeros(2),C=np.eye(2))]+modes;zq,wq=quadnodes(n);nodes=[];lweights=[];lpvals=[];components=[]
    for ci,c in enumerate(comps):
        L=np.linalg.cholesky(c['C'])
        for zz,ww in zip(zq,wq):
            x=c['z']+L@zz
            if np.any(abs(x)>4):continue
            r=model.full(x,'original_cubature_n'+str(n));lp=p_log(x,r['f'],d,sig)
            lg=logsumexp([gaussian_log(x,k['z'],k['C']) for k in comps])-np.log(len(comps))
            nodes.append(r);lweights.append(np.log(ww)-np.log(len(comps))+lp-lg);lpvals.append(lp);components.append(ci)
    lw=np.array(lweights);logZ=logsumexp(lw);w=np.exp(lw-logZ)
    return summarize(nodes,w,d,sig,modes,logZ),dict(nodes=nodes,w=w,lp=np.array(lpvals),logZ=logZ,lw=lw,components=np.array(components))

def summarize(nodes,w,d,sig,modes,logZ):
    X=np.array([r['z'] for r in nodes]);M=np.array([r['chi'] for r in nodes]);F=np.array([r['f'] for r in nodes]);mean=w@X;mm=w@M;fm=w@F
    cov=(X-mean).T@(w[:,None]*(X-mean));mv=w@(abs(M-mm)**2)
    dv=F-d;ph=.5*np.einsum('ij,ij->i',dv,dv)/sig**2
    # Predeclared branch diagnostic is sign split; modes are numerical, not proven branches.
    branch=np.array([w[X[:,0]<0].sum(),w[X[:,0]>=0].sum()])
    return dict(logZ=float(logZ),mean=mean,cov=cov,material_mean=mm,material_sd=np.sqrt(mv),material_total_SD=float(np.sqrt(mv.mean())),response_mean=fm,response_variance_trace=float(w@np.sum((F-fm)**2,axis=1)/sig**2),mean_phi=float(w@ph),sign_branch_mass=branch,node_count=len(nodes),quadrature_relative_ESS=float(1/np.sum(w*w)/len(w)),max_normalized_node_weight=float(w.max()))

def difference(a,b):
    sd=max(b['material_total_SD'],1e-8)
    return dict(relative_Z_change=float(abs(np.expm1(a['logZ']-b['logZ']))),latent_mean_difference=float(np.linalg.norm(a['mean']-b['mean'])),material_mean_RMS=float(np.sqrt(np.mean(abs(a['material_mean']-b['material_mean'])**2))),material_mean_error_over_SD=float(np.sqrt(np.mean(abs(a['material_mean']-b['material_mean'])**2))/sd),sign_branch_mass_max_difference=float(np.max(abs(a['sign_branch_mass']-b['sign_branch_mass']))))

def response(model,modes,d,sig,n):
    """One weak latent + one actual selected response per chart.
    Inverses stay in a fixed radius8 GN ellipsoid and have original-box support.
    Partition is Gaussian weights against identity component, which is integrated
    independently below. This is numerical atlas evidence, no injectivity proof.
    """
    zq,wq=quadnodes(n);nodes=[];lw=[];records=[];ablw={'no_determinant':[],'no_normal_residual':[],'standard_response_gaussian':[]};fail=0;maxerr=0.
    comps=[dict(z=np.zeros(2),C=np.eye(2))]+modes
    # Known POU uses strictly positive Gaussian components; fixed chart support
    # is ball q<64 with compact bump. Identity handles uncovered remainder.
    def bumps(x):
        out=[]
        for m in modes:
            v=x-m['z'];r2=float(v@np.linalg.solve(m['C'],v))/64
            out.append(np.exp(-1/(1-r2)+1) if r2<1 else 0.)
        return np.asarray(out)
    def rho(x,c):
        b=bumps(x);return b[c]/(1+b.sum())
    for ci,m in enumerate(modes):
        # W is strongest GN direction, V its fixed orthogonal complement.
        _,_,vh=np.linalg.svd(m['J'],full_matrices=False);W=vh[0];V=np.array([-W[1],W[0]])
        U=m['J']@W;U/=np.linalg.norm(U)
        J0=m['J']*sig;B0=float(U@(J0@W));a_sd=float(np.sqrt(V@m['C']@V));s0=float(U@(m['r']['f']-d)/sig)
        for t,ww in zip(zq,wq):
            a=a_sd*t[0];s=t[1] # selected residual is centered at zero, prior center preserved
            x=m['z']+V*a+W*(sig*(s-s0)/B0)
            ok=False;its=0
            for its in range(12):
                if np.any(abs(x)>4.05):break
                v=x-m['z'];r2=float(v@np.linalg.solve(m['C'],v))
                if r2>=70:break
                r=model.full(x,'response_inverse_n'+str(n));err=float(U@(r['f']-d)/sig-s)
                if abs(err)<2e-8:ok=True;break
                step=sig*err/B0
                if abs(step)>2:break
                x=x-W*step
            if not ok or np.any(abs(x)>4):fail+=1;continue
            rr=rho(x,ci)
            if rr<=0:continue
            jb=model.jac(x,'response_determinant_one_JVP',V=W[:,None])[:,0];B=float(U@jb)
            if not np.isfinite(B) or abs(B)<1e-12:fail+=1;continue
            maxerr=max(maxerr,abs(err));pprior=-.5*float(x@x)-LOGPRIOR_NORM
            resid=(r['f']-d)/sig;normal2=float(resid@resid-s*s)
            # da proposal N(0,a_sd²), ds proposal N(0,1); Lebesgue factors retained.
            logq_a=-.5*(a/a_sd)**2-.5*np.log(2*np.pi)-np.log(a_sd)
            logq_s=-.5*s*s-.5*np.log(2*np.pi)
            common=np.log(ww)+np.log(rr)+pprior-.5*s*s-logq_a-logq_s
            determinant=np.log(sig)-np.log(abs(B));normal=-.5*max(0,normal2)
            nodes.append(r);lw.append(common+determinant+normal)
            ablw['no_determinant'].append(common+normal)
            ablw['no_normal_residual'].append(common+determinant)
            ablw['standard_response_gaussian'].append(np.log(ww)+np.log(rr))
            records.append(dict(chart=ci,a=a,s=s,B=B,normal_residual2=normal2,log_remaining_weight=pprior+normal-np.log(abs(B))+np.log(rr),inverse_error=err,iterations=its+1,log_unnormalized_quad_mass=common+determinant+normal))
    # Identity POU component on actual original-coordinate quadrature. Its
    # residual and prior are exact at nodes; omitted inverse mass is not hidden.
    ident_summary,ident=original(model,modes,d,sig,n if n>=9 else 9)
    for r,l in zip(ident['nodes'],ident['lw']):
        r0=1/(1+bumps(r['z']).sum());nodes.append(r);lw.append(l+np.log(r0))
        for name in ablw:ablw[name].append(l+np.log(r0))
    if not lw:raise RuntimeError('no response nodes')
    lw=np.array(lw);logZ=logsumexp(lw);w=np.exp(lw-logZ)
    out=summarize(nodes,w,d,sig,modes,logZ)
    ab={k:summarize(nodes,np.exp(np.array(v)-logsumexp(v)),d,sig,modes,float(logsumexp(v))) for k,v in ablw.items()}
    diag=dict(inverse_failures=fail,inverse_attempts=len(modes)*n*n,max_whitened_inverse_error=maxerr,records=records,identity_mass_fraction=float(np.exp(logsumexp(lw[-len(ident['nodes']):])-logZ)),chart_count=len(modes),reliability='NUMERICAL_NOT_CERTIFIED',global_curvature_bound='NOT_AVAILABLE',global_Delta='NOT_AVAILABLE')
    return out,ab,diag,dict(nodes=nodes,w=w,lp=np.array([p_log(r['z'],r['f'],d,sig) for r in nodes]),logZ=logZ)

def baseline_metrics(model,m,d,sig,ref,refraw):
    Cgn=m['C'];h=np.maximum(.02*np.sqrt(np.diag(Cgn)),1e-6);x=m['z'];f=lambda z:-p_log(z,model.full(z,'observed_Hessian_Laplace')['f'],d,sig)
    f0=f(x);H=np.empty((2,2));
    for k in range(2):
        v=np.eye(2)[k]*h[k];H[k,k]=(f(x+v)-2*f0+f(x-v))/h[k]**2
    u=np.array([h[0],0]);v=np.array([0,h[1]])
    H[0,1]=H[1,0]=(f(x+u+v)-f(x+u-v)-f(x-u+v)+f(x-u-v))/(4*h[0]*h[1])
    covs={'GN_LIS_fullrank':Cgn};status={'RTO':'NOT_RUN: nonlinear branch-normalized RTO density and success normalization not established','ordinary_Spline':'NOT_RUN: no trained normalized spline baseline in this stage','manifold_lifting':'NOT_RUN: full constrained auxiliary-state implementation absent'}
    if np.linalg.eigvalsh(H).min()>0:covs['Laplace_observed_Hessian']=np.linalg.inv(H)
    else:status['Laplace_observed_Hessian']='NOT_RUN: observed Hessian not SPD'
    X=np.array([r['z'] for r in refraw['nodes']]);w=refraw['w'];lp=refraw['lp']-refraw['logZ'];ans={}
    for name,C in covs.items():
        norm=multivariate_normal.cdf([4,4],mean=x,cov=C,lower_limit=[-4,-4],rng=np.random.default_rng(820))
        lq=gaussian_log(X,x,C)-np.log(norm);rat=np.exp(np.clip(lq-lp,-745,690))
        # Reverse integral is evaluated on same independent original cubature.
        tv=.5*float(w@abs(1-rat));kp=float(w@(lp-lq));kq=float(w@(rat*(lq-lp)));qmass=float(w@rat)
        mean=np.einsum('n,nj->j',w*rat,np.array([r['chi'] for r in refraw['nodes']]))/qmass
        ans[name]=dict(TV=tv,KL_p_q=kp,KL_q_p=kq,q_integral_check=qmass,material_mean_error_over_SD=float(np.sqrt(np.mean(abs(mean-ref['material_mean'])**2))/max(ref['material_total_SD'],1e-8)),covariance=C,normalization=norm,reliability='NUMERICAL_NOT_CERTIFIED')
    return ans,status,H

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--space',type=int,default=-1);args=parser.parse_args()
    t=Target(stage='A',job='A1_response_v1',output_dir=OUT/'physics/A1_response_v1')
    rng=np.random.default_rng(20261026051);noise_rng=np.random.default_rng(20261026052)
    parents=ndtri(ndtr(-4)+(ndtr(4)-ndtr(-4))*rng.random((3,2,2)))
    innovations=noise_rng.normal(size=(3,2,512))
    save(OUT/'datasets/PARENTS_AND_NOISE.npz',parents=parents,innovations=innovations,coordinate_pairs=np.array(PAIRS),noise_multipliers=np.array([1,.5,.25]))
    allrows=[];start=time.time()
    for si,pair in enumerate(PAIRS):
        if args.space>=0 and si!=args.space:continue
        model=Model(t,pair);grid=[]
        for a in np.linspace(-3.999999,3.999999,13):
            for b in np.linspace(-3.999999,3.999999,13):grid.append(model.full(np.array([a,b]),'global_branch_search_space'+str(si)))
        report(OUT/f'physical/SPACE_{si}_GRID.json',dict(space=si,pair=pair,nodes=len(grid),cost=t.book.totals()))
        for pi in range(2):
            # Truth used exclusively to generate each declared observation.
            truth=model.full(parents[si,pi],'prior_predictive_parent');fstar=truth['f'].copy()
            for ni,mult in enumerate([1.,.5,.25]):
                name=f'S{si}_P{pi}_N{ni}';sig=SIGMA*mult;d=fstar+sig*innovations[si,pi]
                save(OUT/f'datasets/{name}.npz',d=d,sigma=sig,parent=np.array(f'S{si}_P{pi}'),coordinate_pair=np.array(pair))
                before=t.book.totals();now=time.time();modes,attempts=discover(model,grid,d,sig)
                report(OUT/f'physical/{name}_MODES.json',dict(attempts=attempts,modes=[{k:v for k,v in m.items() if k not in ('J','r')} for m in modes],cost=t.book.totals()))
                originals={};last=None;raw=None
                for n in [9,15,21]:
                    ss,rr=original(model,modes,d,sig,n);originals[str(n)]=ss
                    report(OUT/f'physical/{name}_ORIGINAL_{n}.json',ss)
                    if last is not None:
                        dd=difference(ss,last)
                        if n>=15 and dd['material_mean_error_over_SD']<.002 and dd['relative_Z_change']<.002:raw=rr;last=ss;break
                    last=ss;raw=rr
                ref=last;responses={};ablate={};diags={};previous=None
                for n in [5,9]:
                    ss,ab,dg,rraw=response(model,modes,d,sig,n);responses[str(n)]=ss;ablate[str(n)]=ab;diags[str(n)]=dg
                    report(OUT/f'physical/{name}_RESPONSE_{n}.json',dict(summary=ss,ablations=ab,diagnostics=dg))
                    previous=ss
                bmetrics,bstatus,H=baseline_metrics(model,modes[0],d,sig,ref,raw)
                conv_orig=difference(list(originals.values())[-1],list(originals.values())[-2]);conv_resp=difference(responses['9'],responses['5']);cross=difference(responses['9'],ref)
                row=dict(observation=name,space=si,parent=pi,noise_multiplier=mult,pair=pair,original_convergence=conv_orig,response_convergence=conv_resp,cross_coordinate_agreement=cross,baselines=bmetrics,baseline_status=bstatus,mode_count=len(modes),original_orders=list(originals),response_orders=[5,9],inverse_failures=diags['9']['inverse_failures'],identity_mass_fraction=diags['9']['identity_mass_fraction'],cost_before=before,cost_after=t.book.totals(),wall_seconds=time.time()-now,status='COMPLETED_NUMERICAL_NOT_CERTIFIED',full66='NOT_RUN')
                row['gates']=dict(reference_convergence=bool(conv_orig['material_mean_error_over_SD']<=.01 and conv_orig['relative_Z_change']<=.01),response_mean=bool(cross['material_mean_error_over_SD']<=.01),response_normalization=bool(cross['relative_Z_change']<=.01),coverage='NOT_CERTIFIED',TV='NOT_CERTIFIED',cost_saving='NOT_ESTABLISHED',scale_up='HOLD')
                report(OUT/f'physical/{name}_RESULT.json',row);allrows.append(row)
                report(OUT/'RUN_STATUS.json',dict(status='RUNNING',completed_observations=len(allrows),results=allrows,cost=t.book.totals(),elapsed_seconds=time.time()-start))
                print(json.dumps(plain(dict(observation=name,completed=len(allrows),cross=cross,cost=t.book.totals()))),flush=True)
                model.close()
        del model;gc.collect()
    report(OUT/'RUN_STATUS.json',dict(status='COMPLETE_NUMERICAL_NOT_CERTIFIED',completed_observations=len(allrows),results=allrows,cost=t.book.totals(),elapsed_seconds=time.time()-start));print('COMPLETE',flush=True)
if __name__=='__main__':
    try:main()
    except BaseException as e:
        report(OUT/'RUN_FAILURE.json',dict(error=repr(e),traceback=traceback.format_exc(),status='FAILED_OR_BUDGET_STOP'))
        raise
