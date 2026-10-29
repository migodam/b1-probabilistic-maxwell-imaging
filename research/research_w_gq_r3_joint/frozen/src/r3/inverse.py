"""Measurement-only, scaled bounded trust-region Gauss-Newton baseline.

The proper prior is independent truncated Normal(3,1.5) wall eps and
Normal(.06,.15) wall loss, and uniform over all frozen target bounds.
Constants independent of z are dropped. No templates, finite target labels,
truth values, scoring data, shape choice by truth, or output image clipping.
"""
import json,time
from pathlib import Path
import numpy as np
from scipy.optimize import lsq_linear
from .io import Trace,save_endpoint
from .accounting import plain,BudgetExceeded

def whiten(v,sigma):
    a=np.asarray(v);return np.r_[a.real.ravel(),a.imag.ravel()]/sigma

def prior(model,z):
    n=(len(z)-7)//2;r=np.r_[(z[:n]-3.)/1.5,(z[n:2*n]-.06)/.15];J=np.zeros((2*n,len(z)));J[np.arange(2*n),np.arange(2*n)]=np.r_[np.full(n,1/1.5),np.full(n,1/.15)];return r,J

def residual(model,y,z,pred,with_jac=False):
    r=whiten(pred.data-y,model.config['sigma']);rp,jp=prior(model,z);rr=np.r_[r,rp]
    if not with_jac:return rr
    jc=pred.jacobian;J=np.vstack((jc.real,jc.imag))/model.config['sigma'];return rr,np.vstack((J,jp))

def objective(model,y,z,pred):
    rr=residual(model,y,z,pred);return .5*float(rr@rr)

def optimize(model,y,z0,trace,label,*,max_evaluations=100,max_accepted=75,active=None,gtol=2e-4,ftol=1e-8,jacobian_mode=True):
    spec=model.spec;scale=np.asarray(spec['scale']);lo=np.asarray(spec['lower']);hi=np.asarray(spec['upper']);z=np.asarray(z0,float).copy();active=np.arange(len(z)) if active is None else np.asarray(active,int);n=len(active)
    # Stage-specific fixed values may include contrast=0 for wall-only initialization.
    p=model.forward(z,jacobian=jacobian_mode,label=label+' base');r,J=residual(model,y,z,p,True);f=.5*float(r@r);cid=trace.record(z,p,f,stage=label,kind='base');best=(f,z.copy(),p,cid);accepted=0;evals=1;delta=1.;damping=1e-3;small=0;stop='MAX_EVALUATIONS';history=[]
    while evals<max_evaluations and accepted<max_accepted:
        model.book.check();Js=J[:,active]*scale[active];g=Js.T@r;H=Js.T@Js;D=np.sqrt(np.maximum(np.diag(H),1e-6));qlo=np.maximum((lo[active]-z[active])/scale[active],-delta);qhi=np.minimum((hi[active]-z[active])/scale[active],delta)
        projected=g.copy();projected[(qlo>=-1e-10)&(g>0)]=0.;projected[(qhi<=1e-10)&(g<0)]=0.;pg=float(np.max(np.abs(projected)/np.maximum(D,1.)))
        if pg<gtol:stop='PROJECTED_GRADIENT';break
        try:
            sub=lsq_linear(np.vstack((Js,np.sqrt(damping)*np.diag(D))),np.r_[-r,np.zeros(n)],bounds=(qlo,qhi),method='trf',tol=1e-10,max_iter=100)
            step=sub.x
            if jacobian_mode=='ADAPTIVE_CURRENT_CHART':
                direction=np.zeros(len(z));direction[active]=scale[active]*step;audits,fallback=model.audit_chart_step(z,direction,p,tolerance=.1);trace.event(stage=label,event='CHART_STEP_AUDIT',audits=audits,fallback=fallback is not None)
                if fallback is not None:
                    p=fallback;r,J=residual(model,y,z,p,True);Js=J[:,active]*scale[active];g=Js.T@r;H=Js.T@Js;D=np.sqrt(np.maximum(np.diag(H),1e-6));sub=lsq_linear(np.vstack((Js,np.sqrt(damping)*np.diag(D))),np.r_[-r,np.zeros(n)],bounds=(qlo,qhi),method='trf',tol=1e-10,max_iter=100);step=sub.x
            zp=z.copy();zp[active]+=scale[active]*step
            # Numerical bound roundoff is handled only within 1e-12; this is
            # optimizer feasibility, never post-hoc reconstruction clipping.
            if np.any(zp[active]<lo[active]-1e-12) or np.any(zp[active]>hi[active]+1e-12):raise RuntimeError('infeasible TR subproblem')
            zp[active]=np.minimum(np.maximum(zp[active],lo[active]),hi[active])
            pp=model.forward(zp,label=label+f' trial{evals}');rr=residual(model,y,zp,pp);ff=.5*float(rr@rr);evals+=1;tid=trace.record(zp,pp,ff,stage=label,kind='trial',parent_candidate=cid)
        except BudgetExceeded:raise
        except Exception as e:
            model.book.failure('R3',label,e);trace.event(stage=label,event='FAILED_TRIAL',error=str(e));evals+=1;damping=min(damping*10,1e12);delta*=.5
            if delta<1e-6:stop='FAILED_TRUST_REGION';break
            continue
        predicted=-float(g@step+.5*step@H@step);actual=f-ff;rho=actual/predicted if predicted>0 else -1.;ok=actual>0 and rho>.03
        trace.event(stage=label,event='TR_DECISION',candidate=tid,accepted=ok,objective_before=f,objective_trial=ff,actual_decrease=actual,predicted_decrease=predicted,rho=rho,scaled_step_norm=float(np.linalg.norm(step)),delta=delta,damping=damping,projected_gradient=pg)
        if ok:
            fprev=f;z,p,f,cid=zp,pp,ff,tid;accepted+=1;history.append(dict(candidate=cid,objective=f,pg=pg));print(label,'accepted',accepted,'eval',evals,'obj',round(f,6),'pg',round(pg,5),flush=True)
            if f<best[0]:best=(f,z.copy(),p,cid)
            small=small+1 if actual<=ftol*(1+fprev) else 0
            if small>=3:stop='RELATIVE_OBJECTIVE';break
            if rho>.75:damping=max(damping/3,1e-9);delta=min(delta*1.6,3.)
            elif rho<.25:damping=min(damping*3,1e8);delta=max(delta*.65,1e-5)
            # Evaluating same state reuses its LU/primal and pays only shared adjoints.
            p=model.forward(z,jacobian=jacobian_mode,label=label+' accepted Jacobian');r,J=residual(model,y,z,p,True)
            np.savez_compressed(trace.out/f'jacobian_{cid:05d}.npz',z=z,data=p.data,jacobian=p.jacobian,jacobian_mode=str(jacobian_mode))
        else:
            damping=min(damping*6,1e12);delta*=.5
            if delta<1e-6:stop='TRUST_REGION_STALLED';break
    if accepted>=max_accepted:stop='MAX_ACCEPTED'
    # Final exact full-J diagnostic is shared by every method and never selects
    # branches or changes optimizer steps. Approximate projected gradients can
    # otherwise make a fixed-chart stop look like true stationarity.
    final_approx=model.last_chart_prediction.jacobian.copy() if model.last_chart_prediction is not None and np.array_equal(model.last_chart_z,best[1]) and jacobian_mode is not True else None
    finalp=model.forward(best[1],jacobian=True,label=label+' final exact stationarity audit');fr,fJ=residual(model,y,best[1],finalp,True);fJs=fJ[:,active]*scale[active];fg=fJs.T@fr;fD=np.sqrt(np.maximum(np.sum(fJs*fJs,axis=0),1e-6));atlo=best[1][active]<=lo[active]+1e-9;athi=best[1][active]>=hi[active]-1e-9;fg[atlo&(fg>0)]=0.;fg[athi&(fg<0)]=0.;finalpg=float(np.max(np.abs(fg)/np.maximum(fD,1.)))
    final_arrays=dict(z=best[1],data=finalp.data,exact_jacobian=finalp.jacobian,residual=finalp.data-y,parameter_scales=scale,prior_residual=prior(model,best[1])[0],prior_jacobian=prior(model,best[1])[1])
    if final_approx is not None:final_arrays['approximate_jacobian']=final_approx
    np.savez_compressed(trace.out/(label+'_FINAL_JACOBIAN.npz'),**final_arrays);best=(best[0],best[1],finalp,best[3])
    summary=dict(stage=label,objective=best[0],selected_candidate=best[3],accepted_steps=accepted,evaluations=evals,stop=stop,final_exact_projected_gradient=finalpg,final_approximate_jacobian_saved=final_approx is not None,history=history)
    trace.event(event='START_FINISHED',**summary);return best[1],best[2],summary

def invert(model,y,out,method='FULL_GN',options=None):
    options={} if options is None else options;out=Path(out);trace=Trace(out);y=np.asarray(y,complex)
    if y.shape!=(len(model.config['frequencies']),4,16):raise ValueError('measurement array must be (3,4,16)')
    if method not in ['FULL_GN','FIXED_CURRENT_CHART','ADAPTIVE_CURRENT_CHART']:raise ValueError('Unknown method')
    model.observation=y.copy();model.audit_dir=out/'jacobian_audits'
    if method!='FULL_GN':model.initialize_charts(rank=16)
    (out/'ONLINE_INPUT_CONTRACT.json').write_text(json.dumps(dict(config=model.config,spec=model.spec,prior='Independent truncated Normal wall eps mean3 sd1.5, wall loss mean.06 sd.15; independent uniform target bounds',method=method,options=options,input='measurement y only'),default=plain,indent=2))
    z=np.asarray(model.spec['initial'],float);nwall=len(z)-7;zw=z.copy();zw[-3: -1]=0.
    wallz,_,wallrec=optimize(model,y,zw,trace,'wall_only_initialization',max_evaluations=options.get('wall_evaluations',35),max_accepted=options.get('wall_accepted',25),active=np.arange(nwall),gtol=1e-3)
    candidates=[]
    # 18 public continuous seeds. No old finite12 labels or stored PDE templates.
    for x in [-.6,0.,.6]:
        for yy in [1.8,2.1,2.4]:
            for rx,ry,ang in [(.12,.12,0.),(.15,.09,.35)]:
                zi=z.copy();zi[:nwall]=wallz[:nwall];zi[-7:]=[x,yy,rx,ry,3.,.06,ang];pi=model.forward(zi,label='measurement-only coarse seed');fi=objective(model,y,zi,pi);ci=trace.record(zi,pi,fi,stage='coarse_initialization',kind='coarse');candidates.append((fi,zi,pi,ci))
    candidates.sort(key=lambda r:r[0]);seeds=[]
    for row in candidates:
        if not seeds or all(np.linalg.norm((row[1][-7:-5]-s[1][-7:-5])/np.array([.5,.35]))>.35 for s in seeds):seeds.append(row)
        if len(seeds)>=options.get('starts',4):break
    if len(seeds)<options.get('starts',4):
        for row in candidates:
            if all(row[3]!=s[3] for s in seeds):seeds.append(row)
            if len(seeds)>=options.get('starts',4):break
    np.savez_compressed(out/'INITIALIZATION.npz',wall_z=wallz,seeds=np.asarray([s[1] for s in seeds]),coarse_objectives=np.asarray([s[0] for s in seeds]))
    endpoints=[]
    for k,row in enumerate(seeds):
        try:
            zz,pp,rec=optimize(model,y,row[1],trace,f'joint_start_{k}',max_evaluations=options.get('max_evaluations',120),max_accepted=options.get('max_accepted',100),jacobian_mode=True if method=='FULL_GN' else method);endpoints.append((rec['objective'],zz,pp,rec));np.savez_compressed(out/f'START_{k}_ENDPOINT.npz',z=zz,data=pp.data);np.savez_compressed(out/f'START_{k}_EXACT_JACOBIAN.npz',z=zz,data=pp.data,exact_jacobian=pp.jacobian,residual=pp.data-y,parameter_scales=np.asarray(model.spec['scale']),prior_residual=prior(model,zz)[0],prior_jacobian=prior(model,zz)[1])
        except BudgetExceeded as e:
            trace.event(event='DEADLINE_STOP',error=str(e));break
    if not endpoints:raise RuntimeError('No joint endpoint')
    endpoints.sort(key=lambda r:r[0]);f,zz,pp,rec=endpoints[0];summary=dict(method=method,objective=f,wall_initialization=wallrec,selected_start=rec['stage'],starts=[e[3] for e in endpoints],data_chi2=float(np.sum(np.abs(pp.data-y)**2)/model.config['sigma']**2),parameters=dict(zip(model.spec['names'],zz.tolist())),selection='minimum full measurement objective including frozen prior',converged=rec['stop'] in ['PROJECTED_GRADIENT','RELATIVE_OBJECTIVE'] and rec['final_exact_projected_gradient']<.005)
    with np.load(out/(rec['stage']+'_FINAL_JACOBIAN.npz')) as final_arrays:
        if 'approximate_jacobian' in final_arrays:np.savez_compressed(out/'endpoint_approximate_jacobian.npz',**{k:final_arrays[k] for k in final_arrays.files})
    save_endpoint(model,zz,pp,out,summary);return summary
