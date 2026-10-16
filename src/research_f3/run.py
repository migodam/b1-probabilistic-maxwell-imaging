"""One bounded F3 experiment. GT is excluded from every online code path."""
import json,sys,time,traceback,os
import numpy as np
from scipy.optimize import lsq_linear, minimize, brentq
from scipy.linalg import subspace_angles
from scipy.ndimage import gaussian_filter
from .common import ROOT,OUT,OLD,Book,BudgetStop,whiten,dump,mass,rms,TAU
from .prior import QuotientPrior,SCALE,LIMIT,reflected,finite_number
from prior_constrained import prior
sys.path.insert(0,str(prior.A30/'src'))
from a30.physics import FrozenMaxwellAdapter,_cpu

def response_jac(T,state):return np.einsum('qp,ps->qsp',T,state.field).reshape(256,576)

def geometry(P,q,C,sigma):
    D=P.derivative(q);DM=np.r_[D.real,D.imag]/24;Q,s,V=np.linalg.svd(DM,full_matrices=False);rank=int(sum(s>s[0]*1e-7));E=V[:rank].T/s[:rank]
    Z=C@D;A=np.column_stack([whiten(Z[:,i].reshape(32,8),sigma) for i in range(7)])
    U,ss,VT=np.linalg.svd(A@E,full_matrices=False);v=E@VT[-1];w=D@v
    if (q[6]<1e-7 and v[6]<0) or (q[6]>1-1e-7 and v[6]>0):v=-v;w=-w
    elif 1e-7<=q[6]<=1-1e-7 and w[np.argmax(abs(w))].real<0:v=-v;w=-w
    return dict(D=D,DM=DM,Q=Q[:,:rank],s=s,rank=rank,E=E,A=A,response_s=ss,v=v,w=w)

def kkt(q,A,D,r):
    g=A.T@r;pg=g.copy();exact=g.copy()
    for tol,vec in [(1e-8,pg),(0.,exact)]:vec[((q<=tol)&(g>=0))|((q>=1-tol)&(g<=0))]=0
    fit=lsq_linear(A,-r,bounds=(-q,1-q),tol=1e-13,lsq_solver='exact',max_iter=1000)
    delta=fit.x;pred=max(0,float((r@r-np.linalg.norm(r+A@delta)**2)/2));step=rms(D@delta)
    return dict(gradient=g,projected_gradient_inf=float(np.max(abs(pg))),exact_projected_gradient_inf=float(np.max(abs(exact))),
        scaled_gradient_inf=float(np.max(abs(g*np.where(g<0,1-q,q)))),material_GN_step=step,predicted_Phi_decrease=pred,
        qualified=bool(np.max(abs(pg))<=.001 and step<=SCALE/1000 and pred<=1e-6)),delta

class Run:
    def __init__(self,recover=False):
        self.recover=recover;self.cfg=json.loads((OUT/'configs/PROTOCOL.json').read_text());self.book=Book('F3_RECOVERY_1' if recover else 'F3_PILOT');self.records=json.loads((OUT/'results/RECORDS.json').read_text()) if recover else [];self.models={};self.pool={};self.jacs={};self.anchors=[];self.paths={}
    def model(self,bg):
        if bg not in self.models:
            path=OUT/f'inputs/background/{bg}.npz'
            if not path.exists():path=prior.INPUTS/'background'/f'{bg}.npz'
            self.models[bg]=FrozenMaxwellAdapter.from_background_npz(path,self.book,purpose_prefix='RESEARCH_F3')
        return self.models[bg]
    def load(self,row):
        case=row['case'];meta=json.loads((OLD/f'inputs/online/{case}.json').read_text());P=QuotientPrior(meta)
        with np.load(OLD/f'inputs/online/{case}.npz') as z:
            if set(z.files)!={'observation','sigma'}:raise RuntimeError('online input contract')
            y=z['observation'];sigma=float(z['sigma'])
        q=P.from_eight(np.array(row['final_p']));model=self.model(meta['background_id']);path=OLD/f'results/F1/{case}/state_11.npz'
        with np.load(path) as z:chi=z['chi'];current=z['current'];field=z['field'];pred=z['prediction'];factor=(z['LU'],z['pivots']);res=float(z['backward_residual'])
        if np.max(abs(P.decode(q)[0]-chi))>1e-12:raise RuntimeError('quotient material image mismatch')
        L=np.eye(576)-chi[:,None]*model.G
        if np.linalg.norm(L@current-chi[:,None]*model.incident)/np.linalg.norm(chi[:,None]*model.incident)>1e-10:raise RuntimeError('legacy state equation')
        ss=model._export_state(_cpu.FullState(chi,current,field,pred,L,factor,res),cached=True,path=path)
        jacpath=path.with_name(path.stem+'_JAC.npz')
        with np.load(jacpath) as z:T=z['transfer'];C=z['C']
        self.pool[(case,q.tobytes())]=ss;self.jacs[(case,q.tobytes())]=(T,C)
        self.paths[(case,q.tobytes())]=str(path)
        self.book.append(dict(event='LEGACY_CACHE_LOAD',case=case,source=str(path),legacy_state_cost=1,legacy_RHS_cost=40,new_states=0,new_RHS=0,material_roundoff_max=float(np.max(abs(P.decode(q)[0]-chi)))))
        return dict(case=case,P=P,q=q,state=ss,y=y,sigma=sigma,model=model,reference_Phi=float(np.linalg.norm(whiten(pred-y,sigma))**2/2))
    def full(self,ctx,q,key,phase,method):
        cache=(ctx['case'],q.tobytes())
        if cache in self.pool:
            self.book.append(dict(event='SAME_Q_CACHE_HIT',phase=phase,method=method,key=key,new_states=0,new_RHS=0));return self.pool[cache]
        x,row=ctx['P'].decode(q);cert=prior.verify_exact_endpoint(x,row,background_id=ctx['P'].bg)
        if cert['prior_status']!='P0_EXACT':raise RuntimeError('outside original P0')
        self.book.context=dict(case=ctx['case'],phase=phase,method=method,endpoint=key)
        ss=ctx['model'].full_state(x,purpose=key)
        path=OUT/f'results/endpoints/{key}.npz';path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('wb') as f:
            np.savez_compressed(f,q=q,chi=ss.chi,current=ss.current,field=ss.field,prediction=ss.data,LU=ss.factor[0],pivots=ss.factor[1],backward_residual=ss.residual)
            f.flush();os.fsync(f.fileno())
        self.pool[cache]=ss;self.paths[cache]=str(path.relative_to(ROOT));dump(path.with_suffix('.json'),dict(case=ctx['case'],phase=phase,method=method,GT_used=False,prior_status='P0_EXACT',backward_residual=ss.residual));return ss
    def jac(self,ctx,q,ss,phase,method):
        k=(ctx['case'],q.tobytes())
        if k not in self.jacs:
            self.book.context=dict(case=ctx['case'],phase=phase,method=method,operation='receiver_Jacobian')
            T=ctx['model'].receiver_transfer(ss);C=response_jac(T,ss);self.jacs[k]=(T,C)
            p=OUT/f'results/jacobians/{ctx["case"]}_{len(self.jacs):02d}.npz';p.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(p,q=q,T=T,C=C,D=ctx['P'].derivative(q))
        return self.jacs[k]
    def phi(self,ctx,ss):return float(np.linalg.norm(whiten(ss.data-ctx['y'],ctx['sigma']))**2/2)
    def record(self,ctx,q,ss,method,phase,baseq,base,**extra):
        for existing in self.records:
            if existing['case']==ctx['case'] and existing['method']==method and existing['phase']==phase and np.array_equal(existing['q'],q):return existing
        delta=ss.chi-base.chi;lp=ctx['P'].logprior(q);phi=self.phi(ctx,ss);x0=ctx['P'].decode(baseq)[0]
        D0=ctx['P'].derivative(baseq);D=ctx['P'].derivative(q);Q0=np.linalg.svd(np.r_[D0.real,D0.imag]/24,full_matrices=False)[0][:,:7];Q=np.linalg.svd(np.r_[D.real,D.imag]/24,full_matrices=False)[0][:,:7]
        dr=rms(delta);smooth=rms(gaussian_filter(delta.real.reshape(24,24),1)+1j*gaussian_filter(delta.imag.reshape(24,24),1));cells=int(np.sum(abs(delta)>=.012))
        row=dict(case=ctx['case'],method=method,phase=phase,q=q,array=self.paths[(ctx['case'],q.tobytes())],material_RMS=dr,material_MAE=float(np.mean(abs(delta))),
            Phi=phi,delta_Phi=phi-ctx['reference_Phi'],reference='original Research F final material, frozen before F3',
            data_feasible=bool(2*phi<=TAU**2),likelihood_competitive=bool(phi-ctx['reference_Phi']<=2),
            log_prior=finite_number(lp),log_posterior_unnormalized=finite_number(lp-phi),prior_density_zero=not np.isfinite(lp),
            prior_status='P0_EXACT',posterior_supported_positive_density=bool(np.isfinite(lp)),
            task_displacement=bool(dr>=SCALE),structure_check=bool(cells>=4 and smooth>=.002),changed_cells=cells,smoothed_RMS=smooth,
            data_drift=float(np.linalg.norm(whiten(ss.data-base.data,ctx['sigma']))),
            tangent_largest_angle_degrees=float(np.degrees(np.max(subspace_angles(Q0,Q)))),
            rho=(2*q[6]-1)*LIMIT,distance_to_quotient_boundary=float(min(q[6],1-q[6])),backward_residual=ss.residual,
            posterior_samples=False,posterior_chain_proposal=method=='MH',probability_mass=None,**extra)
        self.records.append(row);dump(OUT/'results/RECORDS.json',self.records);return row
    def optimize(self,ctx):
        q=ctx['q'].copy();ss=ctx['state'];hist=[];used=0;stop='ALLOCATION'
        for iteration in range(9):
            T,C=self.jac(ctx,q,ss,'F1.5','quotient_projected_GN');geo=geometry(ctx['P'],q,C,ctx['sigma']);r=whiten(ss.data-ctx['y'],ctx['sigma']);info,delta=kkt(q,geo['A'],geo['D'],r)
            hist.append(dict(iteration=iteration,q=q,Phi=self.phi(ctx,ss),residual=float(np.linalg.norm(r)),geometry_rank=geo['rank'],response_s=geo['response_s'],**info))
            dump(OUT/f'results/F15/{ctx["case"]}_HISTORY.json',hist)
            if info['qualified']:stop='NEW_F3_STATIONARITY_MET';break
            if used>=8:break
            delta*=min(1,.1/max(np.linalg.norm(delta,np.inf),1e-300));accepted=False
            for alpha in [1.,.5,.25]:
                if used>=8:break
                proposal=np.clip(q+alpha*delta,0,1);proposal[proposal<1e-8]=0;proposal[proposal>1-1e-8]=1
                ep=self.full(ctx,proposal,f'{ctx["case"]}_F15_{used}','F1.5','quotient_projected_GN');used+=1
                val=self.phi(ctx,ep);hist[-1].setdefault('trials',[]).append(dict(alpha=alpha,Phi=val,actual_decrease=self.phi(ctx,ss)-val))
                if val<=self.phi(ctx,ss)+1e-12:q=proposal;ss=ep;accepted=True;break
            if not accepted:stop='NO_FEASIBLE_DECREASE_WITH_FROZEN_BACKTRACKING';break
        # Final diagnostics use the existing accepted endpoint; no new state.
        T,C=self.jac(ctx,q,ss,'F1.5','final_stationarity');geo=geometry(ctx['P'],q,C,ctx['sigma']);info,_=kkt(q,geo['A'],geo['D'],whiten(ss.data-ctx['y'],ctx['sigma']))
        info['stationarity_only_pass']=info['qualified']
        info['qualified']=bool(info['qualified'] and self.phi(ctx,ss)<=TAU**2/2 and geo['rank']==7 and
            np.sum(geo['response_s']>geo['response_s'][0]*1e-7)==7 and ss.residual<1e-10)
        result=dict(case=ctx['case'],q=q,states=used,stop=stop,Phi=self.phi(ctx,ss),residual=float(np.sqrt(2*self.phi(ctx,ss))),old_gate2_remains='FAIL',response_s=geo['response_s'],material_s=geo['s'],**info)
        dump(OUT/f'results/F15/{ctx["case"]}_FINAL.json',result);ctx.update(explore_q=q,explore_state=ss,stationary=info['qualified']);self.anchors.append(result)
        print(json.dumps(dict(phase='F1.5',case=ctx['case'],qualified=info['qualified'],gradient=info['projected_gradient_inf'],states=used,cost=self.book.totals())),flush=True)
    def ray(self,P,origin,direction,base_chi,target):
        norm=max(np.linalg.norm(direction),1e-300);direction=direction/norm
        limits=[(1-origin[i])/direction[i] if direction[i]>0 else -origin[i]/direction[i] for i in range(7) if abs(direction[i])>1e-16]
        cap=min(limits);cap=max(0,cap)*(1-1e-10)
        if cap<=0:return origin.copy(),False
        def dist(a):return rms(P.decode(np.clip(origin+a*direction,0,1))[0]-base_chi)-target
        prev=0
        for end in np.linspace(0,cap,65)[1:]:
            if dist(end)>=0:
                root=brentq(dist,prev,end,xtol=1e-14);return np.clip(origin+root*direction,0,1),True
            prev=end
        return np.clip(origin+cap*direction,0,1),False
    def retract(self,P,qseed,base_chi,w,target,goal=None):
        goal=base_chi+target*w if goal is None else goal
        def fun(q):d=P.decode(q)[0]-goal;return float(rms(d)**2)
        def jac(q):d=P.decode(q)[0]-goal;return 2*np.real(P.derivative(q).conj().T@d)/576
        def sphere(q):return rms(P.decode(q)[0]-base_chi)**2-target**2
        def spherejac(q):return 2*np.real(P.derivative(q).conj().T@(P.decode(q)[0]-base_chi))/576
        def progress(q):return float(np.real(np.vdot(w,P.decode(q)[0]-base_chi))/576-.5*target)
        def progressjac(q):return np.real(P.derivative(q).conj().T@w)/576
        fit=minimize(fun,qseed,jac=jac,bounds=[(0,1)]*7,method='SLSQP',constraints=[{'type':'eq','fun':sphere,'jac':spherejac},{'type':'ineq','fun':progress,'jac':progressjac}],options=dict(ftol=1e-14,maxiter=100))
        return np.clip(fit.x,0,1),dict(decoder_only_status=int(fit.status),decoder_only_message=fit.message,progress=progress(fit.x)+.5*target,sphere_error=sphere(fit.x))
    def candidates(self,ctx):
        P=ctx['P'];q0=ctx['explore_q'];base=ctx['explore_state'];T,C=self.jac(ctx,q0,base,'F2','shared_geometry');g0=geometry(P,q0,C,ctx['sigma']);w=g0['w'];v=g0['v'];native={};targets=[SCALE*(1+1e-6),2*SCALE*(1+1e-6)]
        np.savez_compressed(OUT/f'results/{ctx["case"]}_EXPLORATION_GEOMETRY.npz',q0=q0,chi0=base.chi,prediction0=base.data,**g0)
        # Four stationary arms. They have equal endpoint budgets, explicit RHS
        # costs; ordinary native prefixes below provide the matched-RHS control.
        if ctx['stationary']:
            for method in ['A_frozen','B_updated','C_native','D_intrinsic_feedback']:
                q=q0.copy();ss=base;geo=g0;CM=None;II=None;previous_error=None
                if method=='D_intrinsic_feedback':
                    eps=.0005/max(np.linalg.norm(v),1.);central=True
                    if np.any(q0-eps*v<0) or np.any(q0-eps*v>1):central=False
                    for reduction in range(20):
                        tests=[q0+eps*v,q0-eps*v] if central else [q0+eps*v,q0+2*eps*v]
                        if all(np.all(t>=0) and np.all(t<=1) for t in tests):break
                        eps*=.5
                    else:
                        self.book.append(dict(event='CURVATURE_BRANCH_NOT_RUN',case=ctx['case'],reason='no legal stencil'));continue
                    center=P.decode(q0)[0]
                    def d2_at(h):
                        if central:return (P.decode(q0+h*v)[0]-2*center+P.decode(q0-h*v)[0])/h**2
                        return (P.decode(q0+2*h*v)[0]-2*P.decode(q0+h*v)[0]+center)/h**2
                    coarse=d2_at(eps);fine=d2_at(eps/2)
                    d2=(4*fine-coarse)/3 if central else 2*fine-coarse
                    self.book.context=dict(case=ctx['case'],phase='F2_C1',method=method,operation='paid_Maxwell_directional_Hessian')
                    current=ctx['model'].solve_anchor(w[:,None]*base.field,base,role='F3_Hessian_current')
                    cm_complex=T@(2*w[:,None]*ctx['model']._cpu.green(current));CM=whiten(cm_complex,ctx['sigma'])
                    dm=mass(d2);normal=dm-g0['Q']@(g0['Q'].T@dm);II=24*(normal[:576]+1j*normal[576:])
                    cp=whiten((C@II).reshape(32,8),ctx['sigma'])
                    np.savez_compressed(OUT/f'results/{ctx["case"]}_INTRINSIC.npz',w=w,v=v,CM=CM,II=II,Cprior=cp,Cintrinsic=CM+cp,stencil_step=eps,
                        central_stencil=central,D2_coarse=coarse,D2_fine=fine,D2_extrapolated=d2)
                for j,target in enumerate(targets):
                    completed=next((r for r in self.records if r['case']==ctx['case'] and r['method']==method and r['phase']=='F2_C1' and abs(r.get('target_RMS',-1)-target)<1e-14),None)
                    if self.recover and completed is not None:
                        q=np.array(completed['q']);ss=self.pool[(ctx['case'],q.tobytes())]
                        if method=='C_native':native[j]=(q,ss,completed)
                        continue
                    if method=='B_updated' and j:
                        _,Ci=self.jac(ctx,q,ss,'F2_C1',method);geo=geometry(P,q,Ci,ctx['sigma'])
                        if np.real(np.vdot(geo['w'],w))<0:geo['v']=-geo['v'];geo['w']=-geo['w']
                        proposal,matched=self.ray(P,q,geo['v'],base.chi,target)
                    else:proposal,matched=self.ray(P,q0,v,base.chi,target)
                    extra=dict(target_RMS=target,ray_target_reached=matched)
                    if method in ['C_native','D_intrinsic_feedback']:
                        goal=base.chi+target*w
                        if method=='D_intrinsic_feedback':
                            goal=goal+.5*target*target*II
                            A=g0['A'];D=g0['D'];H=g0['DM'].T@g0['DM'];constraint=(H@v).reshape(1,7)
                            _,_,Vh=np.linalg.svd(constraint,full_matrices=True);N=Vh[1:].T
                            predicted=target*(A@v)+.5*target*target*(CM+cp)
                            if j:
                                _,Ci=self.jac(ctx,q,ss,'F2_C1',method);current_geo=geometry(P,q,Ci,ctx['sigma']);A=current_geo['A'];D=current_geo['D']
                                # Measured finite response at previous endpoint
                                # is fed forward without re-solving a resolvent.
                                predicted=predicted+previous_error
                            delta=-N@np.linalg.pinv(A@N,rcond=1e-8)@predicted
                            correction=D@delta;correction*=min(1,target/max(rms(correction),1e-300));goal=goal+correction
                            extra['normal_correction_RMS']=rms(correction)
                        proposal,retract_info=self.retract(P,proposal,base.chi,w,target,goal);extra.update(retract_info)
                    ss=self.full(ctx,proposal,f'{ctx["case"]}_{method}_{j}','F2_C1',method)
                    rec=self.record(ctx,proposal,ss,method,'F2_C1',q0,base,**extra)
                    linear=whiten((C@(ss.chi-base.chi)).reshape(32,8),ctx['sigma']);actual=whiten(ss.data-base.data,ctx['sigma'])
                    rec['actual_material_J_prediction_error']=float(np.linalg.norm(actual-linear));rec['target_match']=bool(abs(rec['material_RMS']-target)<=.02*target)
                    if method=='C_native':native[j]=(proposal,ss,rec)
                    if method=='D_intrinsic_feedback':
                        prediction=target*(g0['A']@v)+.5*target*target*(CM+cp)
                        previous_error=actual-prediction;rec['intrinsic_prediction_error']=float(np.linalg.norm(previous_error));rec['finite_measured_feedback_used']=bool(j)
                        np.savez_compressed(OUT/f'results/{ctx["case"]}_TRANSPORT_PAIR_{j}.npz',base_chi=base.chi,candidate_chi=ss.chi,
                            q0=q0,q1=proposal,direction=w,target_RMS=target,normal_correction=correction,II=II,CM=CM,
                            actual_response=actual,prediction=prediction,finite_prediction_error=previous_error)
                    q=proposal
        # Same-anchor ordinary material-native proposals: eight fixed directions
        # and target radii, no GT and no physics projection masquerading as MH.
        rng=np.random.default_rng(26101762);basis=g0['E'];Aq=g0['A']@basis
        _,s,Vt=np.linalg.svd(Aq,full_matrices=False)
        for i in range(8):
            level=0 if i<4 else 1;target=targets[level]
            if i in [0,4] and level in native:
                proposal,ss,_=native[level];rec=self.record(ctx,proposal,ss,'native_random','F2_C2',q0,base,target_RMS=target,reused_C1=True,cold_RHS=8)
                rec['target_match']=bool(abs(rec['material_RMS']-target)<=.02*target)
                continue
            coefficients=rng.normal(size=7)/np.maximum(s,1e-12);direction=basis@(Vt.T@coefficients)
            ww=P.derivative(q0)@direction;direction/=max(rms(ww),1e-300);ww/=max(rms(ww),1e-300)
            if (q0[6]<=1e-7 and direction[6]<0) or (q0[6]>=1-1e-7 and direction[6]>0):direction=-direction;ww=-ww
            qp,matched=self.ray(P,q0,direction,base.chi,target);proposal,extra=self.retract(P,qp,base.chi,ww,target)
            ss=self.full(ctx,proposal,f'{ctx["case"]}_native_random_{i}','F2_C2','native_random')
            rec=self.record(ctx,proposal,ss,'native_random','F2_C2',q0,base,target_RMS=target,ray_target_reached=matched,cold_RHS=8,**extra)
            rec['target_match']=bool(abs(rec['material_RMS']-target)<=.02*target)
        # Fixed diagonal reflected Gaussian: symmetric in the seven q coords.
        # No state-dependent covariance, selection, truncation or clipping law.
        raw=1/np.maximum(np.linalg.norm(g0['A'],axis=0),1.);norm=np.sqrt(np.sum(np.linalg.norm(g0['DM'],axis=0)**2*raw**2));sd=raw*SCALE/max(norm,1e-300);sd*=min(1,.15/max(sd))
        q=q0.copy();ss=base;lp=P.logprior(q)-self.phi(ctx,ss);chain=[];rng=np.random.default_rng(26101771)
        for i in range(8):
            qp=reflected(q+sd*rng.normal(size=7));ep=self.full(ctx,qp,f'{ctx["case"]}_MH_{i}','F2_C2','MH');epost=P.logprior(qp)-self.phi(ctx,ep)
            if not np.isfinite(epost):loga=-np.inf
            elif not np.isfinite(lp):loga=0.
            else:loga=min(0,epost-lp)
            accept=bool(np.log(rng.random())<loga)
            row=self.record(ctx,qp,ep,'MH','F2_C2',q0,base,MH_accepted=accept,log_acceptance=finite_number(loga),proposal_sd=sd,
                boundary_initialization_extension=bool(not np.isfinite(lp)),proposal_ratio_log=0.)
            if accept:q=qp;ss=ep;lp=epost
            chain.append(dict(step=i,accepted=accept,q=q.copy(),Phi=self.phi(ctx,ss),log_target=finite_number(lp),initial_boundary_zero_density=not np.isfinite(lp),material_RMS_from_initial=rms(ss.chi-base.chi)))
            dump(OUT/f'results/{ctx["case"]}_MH_CHAIN.json',dict(steps=chain,posterior_target='exp(-Phi)*p_rho*2L wrt dq',finite_chain_accuracy='UNASSESSED_8_STEPS',burnin_not_removed=True,proposal='fixed diagonal reflected Gaussian',sd=sd))
        dump(OUT/'results/RECORDS.json',self.records)
    def run(self):
        error=None;contexts=[]
        try:
            oldrows=json.loads((OLD/'results/F1/SUMMARY.json').read_text())
            for row in oldrows:
                ctx=self.load(row)
                if self.recover:
                    # Reload exact paid states and paid receiver transfers.
                    # No LU application, no repeated full state, no new budget.
                    case=ctx['case'];model=ctx['model']
                    for path in sorted((OUT/'results/endpoints').glob(case+'*.npz')):
                        with np.load(path) as z:
                            qq=z['q'];chi=z['chi'];native=_cpu.FullState(chi,z['current'],z['field'],z['prediction'],np.eye(576)-chi[:,None]*model.G,(z['LU'],z['pivots']),float(z['backward_residual']))
                        self.pool[(case,qq.tobytes())]=model._export_state(native,cached=True,path=path);self.paths[(case,qq.tobytes())]=str(path.relative_to(ROOT))
                        self.book.append(dict(event='PAID_STATE_RECOVERY',source=str(path.relative_to(ROOT)),new_states=0,new_RHS=0))
                    for path in sorted((OUT/'results/jacobians').glob(case+'*.npz')):
                        with np.load(path) as z:self.jacs[(case,z['q'].tobytes())]=(z['T'],z['C'])
                    result=json.loads((OUT/f'results/F15/{case}_FINAL.json').read_text());qq=np.array(result['q']);ctx.update(explore_q=qq,explore_state=self.pool[(case,qq.tobytes())],stationary=result['qualified']);self.anchors.append(result)
                else:self.optimize(ctx)
                contexts.append(ctx)
            dump(OUT/'results/F15/SUMMARY.json',self.anchors)
            for ctx in contexts:self.candidates(ctx);print(json.dumps(dict(case=ctx['case'],phase='F2_complete',cost=self.book.totals())),flush=True)
        except BaseException as e:
            error=repr(e);dump(OUT/('FAILURE_RECOVERY.json' if self.recover else 'FAILURE.json'),dict(error=error,traceback=traceback.format_exc(),cost=self.book.totals()));raise
        finally:
            self.book.end(error);dump(OUT/('RUN_FINAL_CLOSED.json' if self.recover else 'RUN_CLOSED.json'),dict(status='CLOSED',error=error,cost=self.book.totals(),no_more_physics=True));dump(OUT/'results/RECORDS.json',self.records)

if __name__=='__main__':Run(recover='--recover' in sys.argv).run()
