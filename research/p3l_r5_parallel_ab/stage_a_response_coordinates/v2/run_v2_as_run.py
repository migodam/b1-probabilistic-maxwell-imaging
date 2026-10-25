"""Version2: bounded conditional integration, disjoint response branches, local QR-RTO."""
from run_maxwell import *
from scipy.optimize import brentq
from scipy.special import log_ndtr
from numpy.polynomial.legendre import leggauss
V2=OUT/'v2';V2.mkdir(exist_ok=True)

def boxline(m,V,W,a):
    z=m+V*a;lo=-np.inf;hi=np.inf
    for j in range(2):
        if abs(W[j])<1e-14:
            if abs(z[j])>4:return None
        else:
            q=sorted([(-4-z[j])/W[j],(4-z[j])/W[j]]);lo=max(lo,q[0]);hi=min(hi,q[1])
    return (lo,hi) if lo<hi else None

def geometry(m):
    _,_,vh=np.linalg.svd(m['J'],full_matrices=False);W=vh[0];V=np.array([-W[1],W[0]]);U=m['J']@W;U/=np.linalg.norm(U)
    return V,W,U
class Conditional:
    def __init__(self,model,m,d,sig):self.model=model;self.m=m;self.d=d;self.sig=sig;self.V,self.W,self.U=geometry(m);self.lines={};self.root_calls=0;self.root_fail=0
    def line(self,a):
        key=round(float(a),14)
        if key in self.lines:return self.lines[key]
        limits=boxline(self.m['z'],self.V,self.W,a)
        if limits is None:return None
        lo,hi=limits;bs=np.linspace(lo,hi,9);rs=[self.model.full(self.m['z']+self.V*a+self.W*b,'v2_branch_scan') for b in bs]
        gs=np.array([self.U@(r['f']-self.d)/self.sig for r in rs]);diff=np.diff(gs)
        # Whole interval admitted only after monotone finite probe; otherwise
        # each disjoint subinterval is a separate numerical branch candidate.
        mono=bool(np.all(diff>0) or np.all(diff<0));indices=[(0,8)] if mono else [(i,i+1) for i in range(8)]
        line=dict(a=a,bounds=(lo,hi),bs=bs,rs=rs,gs=gs,branches=[(bs[i],bs[j],gs[i],gs[j]) for i,j in indices],probe_monotone=mono)
        self.lines[key]=line;return line
    def inverse(self,line,s,branch):
        blo,bhi,glo,ghi=branch;self.root_calls+=1
        if not min(glo,ghi)-1e-9<=s<=max(glo,ghi)+1e-9:return None
        def fun(b):
            r=self.model.full(self.m['z']+self.V*line['a']+self.W*b,'v2_bounded_inverse');return float(self.U@(r['f']-self.d)/self.sig-s)
        try:
            b=brentq(fun,blo,bhi,xtol=2e-12,rtol=4*np.finfo(float).eps,maxiter=35)
            r=self.model.full(self.m['z']+self.V*line['a']+self.W*b,'v2_inverse_endpoint');err=fun(b)
            if abs(err)>2e-7:self.root_fail+=1;return None
            return b,r,err
        except (ValueError,RuntimeError):self.root_fail+=1;return None

def outer_nodes(c,n):
    t,w=leggauss(n);v=c.V;dot=v@c.m['z'];radius=4*np.abs(v).sum();lo=-radius-dot;hi=radius-dot
    return (lo+hi)/2+(hi-lo)/2*t,w*(hi-lo)/2

def nested_original(c,na,nb):
    aa,wa=outer_nodes(c,na);t,wb=hermgauss(nb);t=t*np.sqrt(2);wb=wb/np.sqrt(np.pi)
    nodes=[];lw=[];line_diagnostics=[]
    for a,wwa in zip(aa,wa):
        line=c.line(a)
        if line is None:continue
        centers=[]
        for branch in line['branches']:
            if min(branch[2],branch[3])<=0<=max(branch[2],branch[3]):
                ans=c.inverse(line,0.,branch)
                if ans is not None:centers.append(ans[0])
        if not centers:
            # No selected-response zero: retain original-coordinate mass rather
            # than deleting this entire outer node.
            centers=[float(line['bs'][np.argmin([np.linalg.norm(r['f']-c.d) for r in line['rs']])])]
        components=[]
        for b in centers:
            h=1e-3*max(np.sqrt(c.W@c.m['C']@c.W),1e-4);blo,bhi=line['bounds'];h=min(h,(bhi-blo)/100)
            bl=max(blo,b-h);bh=min(bhi,b+h)
            rlo=c.model.full(c.m['z']+c.V*a+c.W*bl,'v2_original_width_FD');rhi=c.model.full(c.m['z']+c.V*a+c.W*bh,'v2_original_width_FD')
            J=(rhi['f']-rlo['f'])/(bh-bl)/c.sig;sd=1/np.sqrt(1+J@J);components.append((b,float(sd)))
        for b,sd in components:
            for tt,wwb in zip(t,wb):
                bb=b+sd*tt
                if not line['bounds'][0]<=bb<=line['bounds'][1]:continue
                x=c.m['z']+c.V*a+c.W*bb;r=c.model.full(x,'v2_original_conditional_GH')
                lq=logsumexp([-.5*((bb-bc)/sdc)**2-np.log(sdc)-.5*np.log(2*np.pi) for bc,sdc in components])-np.log(len(components))
                nodes.append(r);lw.append(np.log(wwa)+np.log(wwb)-np.log(len(components))+p_log(x,r['f'],c.d,c.sig)-lq)
        line_diagnostics.append(dict(a=a,branches=len(line['branches']),centers=centers,probe_monotone=line['probe_monotone']))
    lw=np.array(lw);lz=logsumexp(lw);w=np.exp(lw-lz);s=summarize(nodes,w,c.d,c.sig,[c.m],lz)
    return s,dict(nodes=nodes,w=w,lp=np.array([p_log(r['z'],r['f'],c.d,c.sig) for r in nodes]),lw=lw,logZ=lz),line_diagnostics

def trunc_normal_nodes(lo,hi,n):
    # Stable Gaussian measure on a finite selected-response interval.
    lo,hi=sorted([lo,hi]);t,w=leggauss(n);u=(t+1)/2;w=w/2
    if hi<=0:
        lhi=log_ndtr(hi);llo=log_ndtr(lo);lm=lhi+np.log1p(-np.exp(min(-1e-16,llo-lhi)))
        mass=np.exp(lm)
        if mass==0:return [],[],lm
        p=ndtr(lo)+u*mass;s=ndtri(np.clip(p,np.finfo(float).tiny,1-np.finfo(float).eps))
    elif lo>=0:
        lhi=log_ndtr(-lo);llo=log_ndtr(-hi);lm=lhi+np.log1p(-np.exp(min(-1e-16,llo-lhi)))
        mass=np.exp(lm)
        if mass==0:return [],[],lm
        p=ndtr(-hi)+u*mass;s=-ndtri(np.clip(p,np.finfo(float).tiny,1-np.finfo(float).eps))
    else:
        mass=ndtr(hi)-ndtr(lo);lm=np.log(mass);s=ndtri(ndtr(lo)+u*mass)
    return s,w,lm

def response_bounded(c,na,ns):
    aa,wa=outer_nodes(c,na);nodes=[];lw=[];no_det=[];no_normal=[];rec=[];fail=0;underflow=0
    for a,wwa in zip(aa,wa):
        line=c.line(a)
        if line is None:continue
        for bi,branch in enumerate(line['branches']):
            ss,ws,lmass=trunc_normal_nodes(branch[2],branch[3],ns)
            if len(ss)==0:underflow+=1;continue
            for s,wws in zip(ss,ws):
                ans=c.inverse(line,float(s),branch)
                if ans is None:fail+=1;continue
                b,r,err=ans;Jb=c.model.jac(r['z'],'v2_response_actual_determinant_JVP',V=c.W[:,None])[:,0];B=float(c.U@Jb)
                if abs(B)<1e-12:fail+=1;continue
                res=(r['f']-c.d)/c.sig;normal2=max(0,float(res@res-s*s));prior=-.5*float(r['z']@r['z'])-LOGPRIOR_NORM
                base=np.log(wwa)+np.log(wws)+lmass+.5*np.log(2*np.pi)+prior
                det=np.log(c.sig)-np.log(abs(B));norm=-.5*normal2
                nodes.append(r);lw.append(base+det+norm);no_det.append(base+norm);no_normal.append(base+det)
                rec.append(dict(a=a,s=s,b=b,branch=bi,B=B,normal_residual2=normal2,inverse_error=err,log_remaining_weight=prior+norm-np.log(abs(B)),log_quadrature_mass=base+det+norm))
    if not nodes:raise RuntimeError('bounded response integration has no nodes')
    lw=np.array(lw);lz=logsumexp(lw);w=np.exp(lw-lz);summary=summarize(nodes,w,c.d,c.sig,[c.m],lz)
    absumm={}
    for name,ls in [('no_determinant',no_det),('no_normal_residual',no_normal)]:
        zz=logsumexp(ls);absumm[name]=summarize(nodes,np.exp(np.array(ls)-zz),c.d,c.sig,[c.m],zz)
    weak=[]
    for a in aa:
        ix=np.array([i for i,r in enumerate(rec) if r['a']==a]);
        if not len(ix):continue
        wi=w[ix];mass=wi.sum();strong=np.array([rec[i]['s'] for i in ix]);sm=float(wi@strong/mass)
        weak.append(dict(a=a,mass=float(mass),strong_mean=sm,strong_variance=float(wi@((strong-sm)**2)/mass)))
    return summary,absumm,dict(records=rec,weak_marginal=weak,inverse_failures=fail,underflow_response_intervals=underflow,reliability='NUMERICAL_NOT_CERTIFIED',coverage='branch intervals partition declared box; monotonicity only finite-probed, unresolved internal folds remain possible',global_Delta='NOT_AVAILABLE',global_curvature='NOT_AVAILABLE')

def rto(c,n=8):
    J=c.m['J'];Q,R=np.linalg.qr(np.r_[np.eye(2),J],mode='reduced');Qx=Q[:2];Qf=Q[2:];x0=c.m['z'];r0=c.m['r'];D0=Qx.T+Qf.T@J;g0=Qx.T@x0+Qf.T@((r0['f']-c.d)/c.sig)
    rng=np.random.default_rng(20261026057);eta=rng.normal(size=(n,2));rows=[]
    for e in eta:
        x=x0+np.linalg.solve(D0,e-g0);good=False
        for it in range(15):
            if np.any(abs(x)>4):break
            rr=c.model.full(x,'QR_RTO_full_forward_inverse');g=Qx.T@x+Qf.T@((rr['f']-c.d)/c.sig);er=g-e
            if np.linalg.norm(er)<2e-7:good=True;break
            x=x-np.linalg.solve(D0,er)
        if not good:rows.append(dict(success=False,iterations=it+1));continue
        jj=c.model.jac(x,'QR_RTO_density_two_JVP')/c.sig;D=Qx.T+Qf.T@jj;sg,ld=np.linalg.slogdet(D)
        lq=-.5*float(e@e)-np.log(2*np.pi)+ld;lp=p_log(x,rr['f'],c.d,c.sig)
        rows.append(dict(success=True,x=x,eta=e,iterations=it+1,inverse_error=float(np.linalg.norm(er)),determinant_sign=float(sg),ideal_log_q=float(lq),log_IS_weight=float(lp-lq),chi=rr['chi']))
    suc=[r for r in rows if r['success']];ans=dict(attempted=n,succeeded=len(suc),rows=[{k:v for k,v in r.items() if k!='chi'} for r in rows],status='LOCAL_NUMERICAL_RTO_ATTEMPT',density_qualification='NOT_GLOBAL: fixed QR projected augmented map genuinely implemented; finite inverse success set/branch multiplicity normalization unproved; ideal map logq not a qualified unconditional law')
    if suc:
        l=np.array([r['log_IS_weight'] for r in suc]);w=np.exp(l-logsumexp(l));ans['local_self_normalized_IS_ESS']=float(1/(w@w));ans['local_material_mean']=sum(ww*r['chi'] for ww,r in zip(w,suc))
    return ans

def main():
    t=Target('A','A1_response_v2',output_dir=OUT/'physics/A1_response_v2');rng=np.random.default_rng(20261026051);noise_rng=np.random.default_rng(20261026052)
    parents=ndtri(ndtr(-4)+(ndtr(4)-ndtr(-4))*rng.random((3,2,2)));innovations=noise_rng.normal(size=(3,2,512))
    models=[];obs=[];results=[];start=time.time()
    for si,pair in enumerate(PAIRS):
        model=Model(t,pair)
        # Explicit compact endpoint reuse from paid v1, no uncounted LU reuse.
        old=t.jobdir;t.jobdir=OUT/'physics/A1_response_v1';legacy=Model(t,pair);t.jobdir=old;model.cache.update(legacy.cache);del legacy;gc.collect();models.append(model)
        for pi in range(2):
            rr=model.full(parents[si,pi],'v2_prior_predictive_parent')
            for ni,mult in enumerate([1.,.5,.25]):
                name=f'S{si}_P{pi}_N{ni}';sig=SIGMA*mult;d=rr['f']+sig*innovations[si,pi];save(V2/f'datasets/{name}.npz',d=d,sigma=sig,coordinate_pair=np.array(pair),parent=np.array(f'S{si}_P{pi}'))
                obs.append(dict(name=name,si=si,pi=pi,ni=ni,mult=mult,sig=sig,d=d))
    report(V2/'OBSERVATION_STATUS.json',dict(observations_generated=18,parents=6,paired_noise=True,cost=t.book.totals()))
    # Coupled weak-volume progression first; all18 remain predeclared.
    obs.sort(key=lambda o:({1:0,0:1,2:2}[o['si']],o['pi'],o['ni']))
    for o in obs:
        model=models[o['si']];name=o['name'];d=o['d'];sig=o['sig'];before=t.book.totals();nc=model.calls;now=time.time()
        try:
            old=OUT/f'v1_failed_preserved/physical/{name}_MODES.json'
            if old.exists():
                dd=json.loads(old.read_text());m=dd['modes'][0];m['z']=np.array(m['z']);m['C']=np.array(m['C']);m['r']=model.full(m['z'],'v2_paid_mode_response');m['J']=model.jac(m['z'],'v2_mode_J')/sig;attempts=dd['attempts']
            else:
                grid=[]
                for a in np.linspace(-3.999999,3.999999,9):
                    for b in np.linspace(-3.999999,3.999999,9):grid.append(model.full(np.array([a,b]),'v2_global_seed_grid'))
                modes,attempts=discover(model,grid,d,sig);m=modes[0]
            c=Conditional(model,m,d,sig);orig={};raw=None
            for na,nb in [(12,5),(24,9)]:
                ss,raw,dg=nested_original(c,na,nb);orig[f'{na}x{nb}']=ss;report(V2/f'{name}_ORIGINAL_{na}x{nb}.json',dict(summary=ss,line_diagnostics=dg,cost=t.book.totals()))
            ref=orig['24x9'];oconv=difference(ref,orig['12x5']);responses={};ab={};diags={}
            for na,ns in [(12,5),(12,9)]:
                ss,aa,dg=response_bounded(c,na,ns);responses[str(ns)]=ss;ab[str(ns)]=aa;diags[str(ns)]=dg;report(V2/f'{name}_RESPONSE_{na}x{ns}.json',dict(summary=ss,ablations=aa,diagnostics=dg,cost=t.book.totals()))
            rconv=difference(responses['9'],responses['5']);cross=difference(responses['9'],ref)
            bmet,bst,H=baseline_metrics(model,m,d,sig,ref,raw)
            rt=rto(c,8);report(V2/f'{name}_QR_RTO.json',rt)
            row=dict(observation=name,space=o['si'],parent=o['pi'],noise_multiplier=o['mult'],status='COMPLETED_NUMERICAL_NOT_CERTIFIED',original_convergence=oconv,response_convergence=rconv,cross_coordinate_agreement=cross,baseline_metrics=bmet,baseline_status=bst,RTO_status=rt['status'],RTO_successes=rt['succeeded'],inverse_failures=diags['9']['inverse_failures'],cost_before=before,cost_after=t.book.totals(),new_forward_calls=model.calls-nc,wall_seconds=time.time()-now,reference_qualification='UNQUALIFIED unless convergence passed; no rigorous enclosure',coverage='UNQUALIFIED finite monotonicity probes do not certify global injectivity',scale_up='HOLD')
            if rt.get('local_material_mean') is not None:row['RTO_local_mean_error_over_SD']=float(np.sqrt(np.mean(abs(rt['local_material_mean']-ref['material_mean'])**2))/max(ref['material_total_SD'],1e-8))
            row['numerical_reference_convergence_pass']=bool(oconv['material_mean_error_over_SD']<=.01 and oconv['relative_Z_change']<=.01)
            report(V2/f'{name}_RESULT.json',row);results.append(row);report(V2/'RUN_STATUS.json',dict(status='RUNNING',completed=len(results),results=results,cost=t.book.totals()));print(json.dumps(plain(row)),flush=True);model.close()
        except BudgetStop as e:
            results.append(dict(observation=name,status='BUDGET_STOP',error=str(e)));break
    completed={r['observation'] for r in results}
    for o in obs:
        if o['name'] not in completed:results.append(dict(observation=o['name'],status='NOT_COMPLETED_WITHIN_ALLOCATED_BUDGET'))
    report(V2/'RUN_STATUS.json',dict(status='CLOSED_NUMERICAL_NOT_CERTIFIED',results=results,cost=t.book.totals(),wall_seconds=time.time()-start));print('CLOSED',flush=True)
if __name__=='__main__':
    try:main()
    except BaseException as e:report(V2/'RUN_FAILURE.json',dict(error=repr(e),traceback=traceback.format_exc()));raise
