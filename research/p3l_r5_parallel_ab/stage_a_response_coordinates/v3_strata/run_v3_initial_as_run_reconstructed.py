"""Mixed-rank coupled (0,64) submodel: exact rank0 host, rank1 upper, rank2 middle."""
from run_v2 import *
from scipy.special import logsumexp
S=OUT/'v3_strata';S.mkdir(exist_ok=True)
CUT=10/13;Z1=ndtr(4)-ndtr(-4);PMASS=(ndtr(-CUT)-ndtr(-4))/Z1
PHI=lambda z:np.exp(-.5*z*z)/np.sqrt(2*np.pi)
HOSTMEAN=(PHI(-4)-PHI(-CUT))/(ndtr(-CUT)-ndtr(-4))
HOSTVAR=1+((-4)*PHI(-4)-(-CUT)*PHI(-CUT))/(ndtr(-CUT)-ndtr(-4))-HOSTMEAN**2
PVAR=1-8*PHI(4)/Z1
class MatModel(Model):
    def __init__(self,t,pair):
        super().__init__(t,pair);self.materials={r['chi'].tobytes():r['f'] for r in self.cache.values()}
    def full(self,z,purpose):
        z=np.asarray(z,float);key=tuple(np.round(z,14))
        if key in self.cache:return self.cache[key]
        chi,_=self.t.decoder(self.V@z);mk=chi.tobytes()
        if mk in self.materials:
            r=dict(z=z.copy(),f=self.materials[mk],chi=chi,key='exact_decoder_material_reuse');self.cache[key]=r;return r
        r=super().full(z,purpose);self.materials[mk]=r['f'];return r

def loadmodel(t,pair):
    m=MatModel(t,pair);job=t.jobdir
    for name in ['A1_response_v1','A1_response_v2']:
        t.jobdir=OUT/'physics'/name;old=Model(t,pair);m.cache.update(old.cache);del old
    t.jobdir=job;m.materials={r['chi'].tobytes():r['f'] for r in m.cache.values()};return m

def grid(model,n):
    x,w=leggauss(n);nodes=[];weights=[]
    for a,wa in zip(CUT*x,CUT*w):
        for b,wb in zip(4*x,4*w):nodes.append(model.full(np.array([a,b]),'v3_transition_original_GL'));weights.append(wa*wb)
    return nodes,np.array(weights)

def pieces(model,m,d,sig,n,gnodes,gweights):
    # Lower host strip has exact positive prior mass and constant material/F.
    host=model.full(np.array([HOSTMEAN,0.]),'v3_exact_host_stratum');lh=np.log(PMASS)-.5*float((host['f']-d)@(host['f']-d))/sig**2
    nodes=[host];lw=[lh];labels=['host']
    # Upper strip: field is gauge; only property remains. Exact field tail mass.
    tt,ww=hermgauss(n);tt=np.sqrt(2)*tt;ww=ww/np.sqrt(np.pi)
    bmean=float(m['z'][1]);bsd=max(float(np.sqrt(m['C'][1,1])),1e-5)
    for mu,sd in [(0.,1.),(bmean,bsd)]:
        for z,w in zip(tt,ww):
            b=mu+sd*z
            if abs(b)>4:continue
            r=model.full(np.array([-HOSTMEAN,b]),'v3_pure_inclusion_property_GH')
            lqb=logsumexp([-.5*b*b-.5*np.log(2*np.pi),-.5*((b-bmean)/bsd)**2-np.log(bsd)-.5*np.log(2*np.pi)])-np.log(2)
            lp=-.5*b*b-.5*np.log(2*np.pi)-np.log(Z1)-.5*float((r['f']-d)@(r['f']-d))/sig**2
            nodes.append(r);lw.append(np.log(PMASS)+np.log(w)-np.log(2)+lp-lqb);labels.append('pure')
    # Transition partition: original Lebesgue GL + adaptive original latent GH.
    C=m['C'];mu=m['z'];L=np.linalg.cholesky(C)
    for r,w in zip(gnodes,gweights):
        x=r['z'];lqm=gaussian_log(x,mu,C);lqp=-.5*x@x-np.log(2*np.pi);lrho=lqp-logsumexp([lqp,lqm])
        nodes.append(r);lw.append(np.log(w)+p_log(x,r['f'],d,sig)+lrho);labels.append('transition')
    zq,wq=quadnodes(n)
    for z,w in zip(zq,wq):
        x=mu+L@z
        if not (-CUT<x[0]<CUT and abs(x[1])<4):continue
        r=model.full(x,'v3_transition_mode_adapted_original_GH');lqm=gaussian_log(x,mu,C);lqp=-.5*x@x-np.log(2*np.pi);lrho=lqm-logsumexp([lqp,lqm])
        nodes.append(r);lw.append(np.log(w)+p_log(x,r['f'],d,sig)+lrho-lqm);labels.append('transition')
    return nodes,np.array(lw),labels

def summary_parts(model,m,d,sig,n,g):
    nodes,lw,labels=pieces(model,m,d,sig,n,*g);lz=logsumexp(lw);w=np.exp(lw-lz);s=summarize(nodes,w,d,sig,[m],lz);mass={k:float(sum(ww for ww,ll in zip(w,labels) if ll==k)) for k in ['host','transition','pure']}
    s['stratum_mass']=mass;s['cov']+=mass['host']*np.diag([HOSTVAR,PVAR])+mass['pure']*np.diag([HOSTVAR,0.])
    return s,dict(nodes=nodes,w=w,lw=lw,lp=np.array([p_log(r['z'],r['f'],d,sig) for r in nodes]),logZ=lz,labels=labels)

def response_transition(model,m,d,sig,hostobs,na,ns):
    # Physical axes are fixed: a=x0, b=x64. Host/pure strips handled separately.
    W=np.array([0.,1.]);probe=model.full(np.array([0.,0.]),'v3_public_transition_chart_anchor');J=model.jac(np.zeros(2),'v3_property_chart_basis');U=J[:,1];U/=np.linalg.norm(U)
    if hostobs:
        a,w=leggauss(na);aa=CUT*a;wa=CUT*w
    else:
        a,w=hermgauss(na);sd=np.sqrt(m['C'][0,0]);aa=m['z'][0]+np.sqrt(2)*sd*a;wa=w/np.sqrt(np.pi);wa=wa*np.exp(.5*((aa-m['z'][0])/sd)**2)*np.sqrt(2*np.pi)*sd
    nodes=[];lw=[];rec=[];fail=0;skip=0
    for x0,wwa in zip(aa,wa):
        if not -CUT<x0<CUT:continue
        blo,bhi=-4.,4.;rl=model.full(np.array([x0,blo]),'v3_response_branch_left');rh=model.full(np.array([x0,bhi]),'v3_response_branch_right');gl=float(U@(rl['f']-d)/sig);gh=float(U@(rh['f']-d)/sig)
        ss,ws,lmass=trunc_normal_nodes(gl,gh,ns)
        if not len(ss):skip+=1;continue
        for s,ww in zip(ss,ws):
            def fun(b):return float(U@(model.full(np.array([x0,b]),'v3_property_Brent_inverse')['f']-d)/sig-s)
            try:b=brentq(fun,blo,bhi,xtol=2e-12,maxiter=40)
            except BudgetStop:raise
            except (ValueError,RuntimeError):fail+=1;continue
            r=model.full(np.array([x0,b]),'v3_response_endpoint');jb=model.jac(r['z'],'v3_response_det_JVP',V=W[:,None])[:,0];B=float(U@jb)
            if abs(B)<1e-12:fail+=1;continue
            res=(r['f']-d)/sig;normal2=max(0,float(res@res-s*s));prior=-.5*r['z']@r['z']-LOGPRIOR_NORM
            ll=np.log(wwa)+np.log(ww)+lmass+.5*np.log(2*np.pi)+prior-.5*normal2+np.log(sig)-np.log(abs(B))
            nodes.append(r);lw.append(ll);rec.append(dict(a=x0,b=b,s=s,B=B,normal_residual2=normal2,log_mass=ll))
    return nodes,np.array(lw),dict(records=rec,inverse_failures=fail,underflow_intervals=skip,qualification='NUMERICAL_MONOTONE_PROPERTY_BRANCH_NOT_GLOBAL_CERTIFICATE')

def main():
    t=Target('A','A1_strata_v3',output_dir=OUT/'physics/A1_strata_v3');model=loadmodel(t,(0,64));grids={24:grid(model,24),48:grid(model,48)};rows=[]
    for pi in range(2):
        for ni,mult in enumerate([1,.5,.25]):
            name=f'S1_P{pi}_N{ni}'
            with np.load(OUT/f'v2/datasets/{name}.npz') as z:d=z['d'];sig=float(z['sigma'])
            host=model.full(np.array([HOSTMEAN,0.]),'v3_host_observation_diagnostic');hostphi=.5*np.sum((host['f']-d)**2)/sig**2;hostobs=hostphi<512
            coarse=grids[24][0];modes,attempts=discover(model,coarse,d,sig);m=modes[0]
            report(S/f'{name}_MODES.json',dict(attempts=attempts,z=m['z'],C=m['C'],host_phi=hostphi))
            orig={};rr={}
            for ng,n in [(24,15),(48,21)]:
                ss,raw=summary_parts(model,m,d,sig,n,grids[ng]);orig[str(ng)]=ss;rr[str(ng)]=raw;report(S/f'{name}_ORIGINAL_{ng}.json',ss)
            ref=orig['48'];conv=difference(ref,orig['24']);resp={}
            for na,ns in [(12,5),(24 if hostobs else 12,9)]:
                nn,ll,dg=response_transition(model,m,d,sig,hostobs,na,ns)
                raw=rr['48'];baseidx=[i for i,k in enumerate(raw['labels']) if k!='transition'];bn=[raw['nodes'][i] for i in baseidx];bl=raw['lw'][baseidx];alln=bn+nn;alll=np.r_[bl,ll];lz=logsumexp(alll);w=np.exp(alll-lz);ss=summarize(alln,w,d,sig,[m],lz);ss['host_mass']=float(np.exp(bl[0]-lz));resp[str(ns)]=ss;report(S/f'{name}_RESPONSE_{ns}.json',dict(summary=ss,diagnostics=dg))
            rt=rto(Conditional(model,m,d,sig),8);report(S/f'{name}_QR_RTO.json',rt)
            row=dict(observation=name,host_phi=hostphi,host_prior_mass=PMASS,original_convergence=conv,response_convergence=difference(resp['9'],resp['5']),cross_coordinate_agreement=difference(resp['9'],ref),stratum_masses=ref['stratum_mass'],RTO_successes=rt['succeeded'],RTO_status=rt['status'],cost=t.book.totals(),qualification='NUMERICAL_NOT_CERTIFIED',scale_up='HOLD')
            rows.append(row);report(S/f'{name}_RESULT.json',row);report(S/'RUN_STATUS.json',dict(status='RUNNING',results=rows,cost=t.book.totals()));print(json.dumps(plain(row)),flush=True);model.close()
    report(S/'RUN_STATUS.json',dict(status='CLOSED_NUMERICAL_NOT_CERTIFIED',results=rows,cost=t.book.totals()))
if __name__=='__main__':
    try:main()
    except BaseException as e:report(S/'RUN_FAILURE.json',dict(error=repr(e),traceback=traceback.format_exc()));raise
