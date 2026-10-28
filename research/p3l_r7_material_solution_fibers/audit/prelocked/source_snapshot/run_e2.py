"""The sole R7 physical runner. Generation, online work and truth scoring separate.

No online function opens offline_truth. No object is dropped for bad convergence.
All trial states are durable before any objective-based selection.
"""
from __future__ import annotations
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[k]='1'
import argparse, hashlib, json, sys, time, traceback
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2
from physics import OUT,REPO,Book,Target,write,save,totals,SIGMA,probit_derivative
from material_jets import (decode,material_derivative,real_material,decoder_coefficients,
    geometry,jets,match_distance,polynomial,tasks,feasible_progress,gn_step)


def stamp_files(paths):
    return [dict(path=str(p),bytes=p.stat().st_size,blake2b=hashlib.blake2b(p.read_bytes()).hexdigest()) for p in paths]


def counters(book):
    return dict(RHS=int(book.counts['actual_forward_RHS']+book.counts['actual_derivative_RHS']),
                forward_RHS=int(book.counts['actual_forward_RHS']),
                derivative_RHS=int(book.counts['actual_derivative_RHS']))


def diff_counts(after,before):return {k:after[k]-before[k] for k in after}


def append(path,row):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('a') as f:f.write(json.dumps(row,allow_nan=False,separators=(',',':'))+'\n');f.flush()


def saved_geometry(target,p,path,native_J=False):
    with np.load(path,allow_pickle=False) as f:
        q=f['q'] if 'q' in f else f['q0']
        if not np.array_equal(q,p.q):raise ValueError('cached geometry coordinate mismatch')
        T=f['receiver_T'].copy()
        source_J=f['J'].copy()
        source_prediction=f['prediction'].copy()
    model=target.backend()
    L=np.eye(576)-p.state.chi[:,None]*model.G
    err=float(np.linalg.norm(T@L-model.S)/np.linalg.norm(model.S))
    prediction_error=float(np.linalg.norm(source_prediction-p.state.data))
    if err>1e-10 or prediction_error>1e-10:raise ValueError('cached transfer identity')
    p.T=T
    J=target.jacobian(p) # receiver solve is bypassed only for validated SAME material T.
    old=source_J*probit_derivative(p.x)[None,:] if native_J else source_J
    je=float(np.linalg.norm(J-old)/max(np.linalg.norm(J),1e-300))
    if je>1e-10:raise ValueError('cached chain-rule identity')
    target.book._event('HISTORICAL_GEOMETRY_REUSE',source=str(path),transfer_error=err,
                       prediction_error=prediction_error,J_relative_error=je,cold_derivative_RHS=32)


def generate():
    manifest=json.loads((OUT/'configs/E2_PARENT_MANIFEST.json').read_text())
    with Book('E2','GENERATE8',role='DATA_GENERATION',max_wall=300) as book:
        target=Target(book)
        try:
            for row in manifest['locked']:
                parent=row['parent']
                destination=OUT/'raw/online'/f'{parent}.npz'
                if destination.exists():raise FileExistsError('no implicit data replay')
                rng=np.random.default_rng(np.random.SeedSequence(row['prior_seed_sequence']))
                x=rng.normal(size=66)
                p=target.evaluate(x,parent+'_truth_generation')
                nr=np.random.default_rng(np.random.SeedSequence(row['noise_seed_sequence']))
                noise=SIGMA/np.sqrt(2)*(nr.normal(size=(32,8))+1j*nr.normal(size=(32,8)))
                save(OUT/'raw/offline_truth'/f'{parent}.npz',x=x,q=p.q,material=p.state.chi,
                     prediction=p.state.data,noise=noise,observation=p.state.data+noise,
                     sigma=np.array(SIGMA),parent=np.array(parent))
                save(destination,observation=p.state.data+noise,sigma=np.array(SIGMA),parent=np.array(parent))
                print('GENERATED',parent,counters(book),flush=True)
            write(OUT/'raw/GENERATION_SEALED.json',dict(parents=[r['parent'] for r in manifest['locked']],
                  online_fields=['observation','sigma','parent'],cost=counters(book),new_truth_not_read_by_online=True))
        finally:target.close()


def direct_anchor(target,y,parent,folder,cap=64):
    folder.mkdir(parents=True,exist_ok=True)
    last=best=None
    best_energy=np.inf
    fresh=0
    def get(x):
        nonlocal last,best,best_energy,fresh
        for p in (last,best):
            if p is not None and np.array_equal(p.x,x):return p
        if fresh>=cap:raise RuntimeError('ANCHOR_STRICT_STATE_CAP')
        fresh+=1
        p=target.evaluate(x,parent+'_direct_probit64')
        energy=target.objective(p,y)
        append(folder/'ANCHOR_HISTORY.jsonl',dict(attempt=fresh,endpoint=p.key,
          phi=target.phi(p,y),objective=energy,prior=.5*float(x@x),x_norm=float(np.linalg.norm(x))))
        last=p
        if energy<best_energy:best=p;best_energy=energy
        return p
    def fun(x):return np.r_[target.white(get(x).state.data-y),x]
    def jac(x):
        p=get(x)
        if p.J is None and p is not target._live_point:
            # Rare rejected best-without-J: paid restoration, never an unowned LU.
            p=target.evaluate(x,parent+'_restore_anchor_J')
        return np.vstack((target.jacobian(p),np.eye(66)))
    status={}
    try:
        r=least_squares(fun,np.zeros(66),jac=jac,method='trf',tr_solver='exact',x_scale='jac',
                        ftol=1e-9,xtol=1e-9,gtol=1e-7,max_nfev=cap)
        status=dict(status=int(r.status),success=bool(r.success),message=str(r.message),
                    nfev=int(r.nfev),njev=int(r.njev),optimality=float(r.optimality))
    except RuntimeError as e:
        if str(e)!='ANCHOR_STRICT_STATE_CAP':raise
        status=dict(status='STRICT_STATE_CAP',success=False)
    if best is None:raise RuntimeError('no anchor state')
    old=best
    if old is not target._live_point:
        best=target.evaluate(old.x,parent+'_restore_best_live_LU')
        if old.T is not None:best.T=old.T.copy()
    target.jacobian(best)
    target.persist_cache(best,folder/'ANCHOR.npz')
    residual=target.white(best.state.data-y)
    gradient=best.J.T@residual+best.x
    write(folder/'ANCHOR_STATUS.json',dict(optimizer=status,fresh_state_attempts=fresh,
         phi=target.phi(best,y),objective=target.objective(best,y),gradient_inf=float(max(abs(gradient))),
         noise_ball=bool(target.phi(best,y)<=chi2.ppf(.95,512)/2),
         GT_used=False,mode='PROBIT_COORDINATE_POINT_NOT_MATERIAL_MAP',cost=counters(target.book)))
    return best


def load_parent(target,parent,folder):
    if parent.startswith('R6_CDEV_'):
        registry=json.loads((OUT/'audit/history/DEV4_RAW_ANCHOR_REGISTRY.json').read_text())['rows']
        row=next(r for r in registry if r['parent']==parent)
        with np.load(row['observation'],allow_pickle=False) as f:y=f['observation'].copy()
        with np.load(row['context'],allow_pickle=False) as f:x=f['x'].copy()
        p=target.evaluate(x,parent+'_historical_center_live_LU')
        saved_geometry(target,p,Path(row['physical_cache']),native_J=True)
        target.persist_cache(p,folder/'ANCHOR.npz')
        write(folder/'ANCHOR_STATUS.json',dict(status='PAID_HISTORICAL_CENTER_REUSED',
             role='EXPOSED_DEVELOPMENT_NOT_NEW_INFERENCE',phi=target.phi(p,y),new_cost=counters(target.book)))
        return p,y
    # This function deliberately cannot derive or open the offline-truth path.
    with np.load(OUT/'raw/online'/f'{parent}.npz',allow_pickle=False) as f:
        if set(f.files)!={'observation','sigma','parent'}:raise ValueError('online input leakage')
        if float(f['sigma'])!=SIGMA:raise ValueError('sigma changed')
        y=f['observation'].copy()
    return direct_anchor(target,y,parent,folder),y


def metrics(target,p,base,y,unit,distance):
    dm=real_material(p.state.chi-base.state.chi)
    pair=target.white(p.state.data-base.state.data)
    phi=target.phi(p,y);phi0=target.phi(base,y)
    info=decode(target.prior,p.x)[1]
    base_info=decode(target.prior,base.x)[1]
    def stratum(a):
        s=a['field']/target.prior.temperature
        return np.where(s<=-1,0,np.where(s>=1,2,1))
    seam_crossings=int(np.sum((stratum(info)!=stratum(base_info))&target.prior.mask))
    rms=float(np.linalg.norm(dm))
    progress=float(unit@dm)
    prior=.5*float(p.x@p.x);prior0=.5*float(base.x@base.x)
    return dict(endpoint=p.key,material_RMS=rms,target_RMS=float(distance),
       distance_matched=bool(abs(rms-distance)<=.02*distance),signed_progress=progress,
       progress_pass=bool(progress>=.5*distance),Phi=phi,DeltaPhi=phi-phi0,
       likelihood_competitive=bool(phi-phi0<=2),noise_ball=bool(phi<=chi2.ppf(.95,512)/2),
       whitened_pair_L2=float(np.linalg.norm(pair)),group_A_pair_L2=float(np.linalg.norm(pair[:256])),
       group_B_pair_L2=float(np.linalg.norm(pair[256:])),prior_cost=prior,DeltaPrior=prior-prior0,
       reasonable_prior_heuristic=bool(prior-prior0<=10),legal=bool(target.prior.legal(p.q)),
       source_order='source-major Re32,Im32;A first256 real components',
       transition_cells=int(np.sum((info['fraction']>0)&(info['fraction']<1)&target.prior.mask)),
       host_cells=int(np.sum((info['fraction']==0)&target.prior.mask)),
       pure_cells=int(np.sum((info['fraction']==1)&target.prior.mask)),
       seam_crossings=seam_crossings,
       tasks=tasks(target.prior,p.state.chi).tolist(),GT_used=False)


def correct(target,base,y,xstart,v,distance,folder,*,frozen=False,max_states=3,rowsel=slice(None)):
    unit=real_material(material_derivative(target.prior,base.x)@v)
    unit/=np.linalg.norm(unit)
    x=xstart.copy();best=None;best_obj=np.inf;trust=.5
    previous=None;history=[]
    for i in range(max_states):
        p=target.evaluate(x,folder.name+f'_trial{i}')
        rr=target.white(p.state.data-y)[rowsel]
        obj=.5*float(rr@rr+p.x@p.x)
        feasible=feasible_progress(target.prior,base.x,x,unit,distance)
        rho=None
        if previous is not None:
            prev_obj,pred=previous
            rho=(prev_obj-obj)/max(pred,1e-300) if pred>0 else -1.
            if rho<.25:trust*=.5
            elif rho>.75:trust=min(1.,trust*1.5)
        row=dict(iteration=i,endpoint=p.key,objective=obj,Phi=target.phi(p,y),
                 feasible=bool(feasible),trust_radius=trust,actual_predicted_ratio=rho)
        history.append(row)
        if feasible and obj<best_obj:best=p;best_obj=obj
        if i==max_states-1:break
        current=best if best is not None else p
        if frozen:J=base.J
        else:
            if current.J is None and current is not target._live_point:
                # Rejected trial has no need of a fresh derivative; its predecessor
                # already owned a J when it proposed this trial.
                raise RuntimeError('missing paid best Jacobian')
            J=target.jacobian(current)
        r=target.white(current.state.data-y)[rowsel]
        J=J[rowsel]
        step,_=gn_step(target.prior,base.x,current.x,v,J,r,trust=trust,
                       material_limit=.5*distance)
        proposal=current.x+step
        xn,match=match_distance(target.prior,base.x,
                lambda h:base.x+h*(proposal-base.x),distance)
        row['projection']=match
        if xn is None or not feasible_progress(target.prior,base.x,xn,unit,distance):
            row['stop']='NO_FEASIBLE_RETRACTED_STEP';break
        actual_step=xn-current.x
        augmented=np.r_[r,current.x]
        model_change=np.r_[J@actual_step,actual_step]
        pred=.5*float(augmented@augmented-np.linalg.norm(augmented+model_change)**2)
        row.update(predicted_reduction=pred,step_norm=float(np.linalg.norm(actual_step)),
                   projection_induced_step_change=float(np.linalg.norm(actual_step-step)))
        if pred<=0 or np.linalg.norm(actual_step)<1e-9:
            row['stop']='NONPOSITIVE_PREDICTED_REDUCTION';break
        previous=(.5*float(r@r+current.x@current.x),pred)
        x=xn
    if best is None:
        write(folder/'CORRECTION.json',dict(status='NO_FEASIBLE_ENDPOINT',history=history));return None
    target.persist_cache(best,folder/'ENDPOINT.npz')
    write(folder/'CORRECTION.json',dict(status='BUDGETED_CORRECTOR_NOT_CONVERGED_PROFILE',
          history=history,selected_endpoint=best.key,selected_objective=best_obj))
    return best


def physical_probe():
    folder=OUT/'E2/probe';folder.mkdir(parents=True,exist_ok=True)
    with Book('E2','FIRST_PHYSICAL_PROBE',parent='R6_CDEV_0000',max_wall=600) as book:
        target=Target(book)
        try:
            base,y=load_parent(target,'R6_CDEV_0000',folder)
            g=geometry(target.prior,base.x,base.J,2026102810)
            j=jets(target,base,g['v'],g['strong_U'])
            save(folder/'COMPENSATED_JETS.npz',**{k:v for k,v in j.items() if isinstance(v,np.ndarray)})
            write(folder/'COMMON_COMPENSATION.json',j['groups'])
            x=base.x.copy();prediction=base.state.data.copy();J=base.J.copy();D=g['D']
            # Independent source-adjoint inner-product check, paid even though J exists.
            rng=np.random.default_rng(2026102812)
            raw=rng.normal(size=(32,8))+1j*rng.normal(size=(32,8))
            adj=target.solve_anchor(target.backend().S.conj().T@raw,base,adjoint=True,role='QA_source_adjoint')
            h=np.sum(np.conj(adj)*base.state.field,axis=1)
            grad=2*np.real(D.T@h)/SIGMA**2
            adjoint_rel=float(np.linalg.norm(grad-J.T@target.white(raw))/np.linalg.norm(grad))
            FD=[]
            for name,v in [('weak',g['v']),('random',g['random_v']),('property',np.eye(66)[64])]:
                v=v/max(np.linalg.norm(v),1e-300)
                prediction_J=J@v
                for hh in (2e-3,1e-3,5e-4):
                    pp=target.evaluate(x+hh*v,'QA_'+name+'_plus')
                    yp=pp.state.data.copy()
                    pm=target.evaluate(x-hh*v,'QA_'+name+'_minus')
                    fd=target.white((yp-pm.state.data)/(2*hh))
                    FD.append(dict(direction=name,h=hh,relative_error=float(np.linalg.norm(fd-prediction_J)/max(np.linalg.norm(prediction_J),1e-300))))
            # Use finite positive amplitudes before declaring the expensive branch ready.
            trials=[]
            for radius in (.0001,.0002):
                xp,m=match_distance(target.prior,x,lambda hh:polynomial(x,j['v'],hh,j['x2'],j['x3']),radius)
                if xp is None:trials.append(dict(radius=radius,status='INFEASIBLE'));continue
                pp=target.evaluate(xp,'QA_jet_endpoint')
                actual=target.white(pp.state.data-prediction)
                hh=m['h'];pred=hh*j['c1']+hh**2*j['c2']+hh**3*j['c3']
                trials.append(dict(radius=radius,h=hh,status='EXECUTED',relative_error=float(np.linalg.norm(actual-pred)/max(np.linalg.norm(actual),1e-300))))
            passed=(max(j['checks'].values())<2e-4 and adjoint_rel<2e-4 and
              all(min(r['relative_error'] for r in FD if r['direction']==d)<2e-4 for d in ('weak','random','property')))
            receipt=dict(status='PASS' if passed else 'FAIL_PAUSE',jet_algebra=j['checks'],
               adjoint_relative_error=adjoint_rel,finite_differences=FD,finite_jet_trials=trials,
               actual_cost=counters(book),background_audit=dict(target.backend().background_audit),
               precondition='No locked physics until this real derivative gate passes')
            write(folder/'PROBE_RESULT.json',receipt)
            print(json.dumps(receipt,indent=2),flush=True)
        finally:target.close()


def process_parent(parent):
    if json.loads((OUT/'E2/probe/PROBE_RESULT.json').read_text())['status']!='PASS':
        raise RuntimeError('physical derivative gate closed')
    folder=OUT/'E2/parents'/parent
    if (folder/'COMPLETE.json').exists():raise FileExistsError('completed parent: no implicit replay')
    folder.mkdir(parents=True,exist_ok=True)
    with Book('E2',parent+'_REV1',parent=parent,max_wall=1200) as book:
        target=Target(book)
        try:
            base,y=load_parent(target,parent,folder)
            save(folder/'OBSERVATION.npz',observation=y,sigma=SIGMA)
            start=time.perf_counter();anchor_cost=counters(book)
            idx=int(parent.rsplit('_',1)[1])
            g=geometry(target.prior,base.x,target.jacobian(base),np.random.SeedSequence([2026102810,idx,9 if parent.startswith('R7') else 8]))
            save(folder/'BASE_GEOMETRY.npz',**g)
            j=jets(target,base,g['v'],g['strong_U'])
            save(folder/'COMPENSATED_JETS.npz',**{k:v for k,v in j.items() if isinstance(v,np.ndarray)})
            write(folder/'JET_QA.json',j['checks']);write(folder/'COMMON_COMPENSATION.json',j['groups'])
            jet_cost=diff_counts(counters(book),anchor_cost)
            if max(j['checks'].values())>2e-4:raise RuntimeError('jet algebra gate')
            ja=jets(target,base,g['v_group_A'],g['strong_U'])
            save(folder/'GROUP_A_COMPENSATED_JETS.npz',**{k:v for k,v in ja.items() if isinstance(v,np.ndarray)})
            write(folder/'GROUP_A_COMMON_COMPENSATION.json',ja['groups'])
            rows=[];endpoints={};rawjet_h=None
            for distance in (.0025,.005,.01,.02):
                # Cyclic order uses parent ID only; no method receives selected warm starts.
                order=['A_FROZEN','B_FULLJ','C_JETS','D_RANDOM']
                order=order[idx%4:]+order[:idx%4]
                for arm in order:
                    dest=folder/arm/f'RMS_{distance:.4f}';dest.mkdir(parents=True,exist_ok=True)
                    begin=counters(book);clock=time.perf_counter()
                    v=g['random_v'] if arm=='D_RANDOM' else g['v']
                    curve=(lambda hh:polynomial(base.x,v,hh,j['x2'],j['x3'])) if arm=='C_JETS' else (lambda hh:base.x+hh*v)
                    xp,match=match_distance(target.prior,base.x,curve,distance)
                    if xp is None:
                        row=dict(parent=parent,arm=arm,radius=distance,status=match['status'],cost=diff_counts(counters(book),begin),wall_seconds=time.perf_counter()-clock)
                        rows.append(row);write(dest/'RESULT.json',row);continue
                    if arm=='C_JETS' and distance==.02:rawjet_h=match['h']
                    p=correct(target,base,y,xp,v,distance,dest,frozen=(arm=='A_FROZEN'))
                    if p is None:
                        row=dict(parent=parent,arm=arm,radius=distance,status='CORRECTION_INFEASIBLE')
                    else:
                        row=metrics(target,p,base,y,real_material(material_derivative(target.prior,base.x)@v),distance)
                        row.update(parent=parent,arm=arm,radius=distance,status='EXECUTED',initial_distance_match=match)
                        if distance==.02:endpoints[arm]=p.x.copy()
                    row.update(cost=diff_counts(counters(book),begin),wall_seconds=time.perf_counter()-clock)
                    if arm=='C_JETS':row['shared_jet_setup_cost']=jet_cost
                    rows.append(row);write(dest/'RESULT.json',row)
                    print('ENDPOINT',parent,arm,distance,row.get('Phi'),row.get('material_RMS'),counters(book)['RHS'],flush=True)
            # Predetermined paths. A chord failure does not prove disconnectedness.
            for name in ('B_FULLJ_CHORD','C_RAW_CUBIC'):
                if name=='B_FULLJ_CHORD' and 'B_FULLJ' not in endpoints:continue
                if name=='C_RAW_CUBIC' and rawjet_h is None:continue
                pathdir=folder/'paths'/name;pathdir.mkdir(parents=True,exist_ok=True)
                pr=[]
                for t in (0,.25,.5,.75,1):
                    xx=(base.x+t*(endpoints['B_FULLJ']-base.x)) if name=='B_FULLJ_CHORD' else polynomial(base.x,g['v'],t*rawjet_h,j['x2'],j['x3'])
                    p=target.evaluate(xx,name+str(t));target.jacobian(p)
                    gg=geometry(target.prior,p.x,p.J,2026102813)
                    save(pathdir/f'NODE_{t:.2f}.npz',x=p.x,material=p.state.chi,prediction=p.state.data,
                         current=p.state.current,field=p.state.field,J=p.J,residual_white=target.white(p.state.data-y),
                         material_singular=gg['singular_material'],response_singular=gg['generalized_joint_singular'],
                         material_frame=gg['Q'],weak_material=gg['weak_material'])
                    row=metrics(target,p,base,y,g['weak_material'],float(np.linalg.norm(real_material(p.state.chi-base.state.chi))))
                    row.update(t=t,path=name,material_effective_rank=gg['rank'])
                    if name=='C_RAW_CUBIC':
                        hh=t*rawjet_h;pred=hh*j['c1']+hh*hh*j['c2']+hh**3*j['c3']
                        actual=target.white(p.state.data-base.state.data)
                        row.update(jet_prediction_error=float(np.linalg.norm(actual-pred)),
                          jet_relative_error=float(np.linalg.norm(actual-pred)/max(np.linalg.norm(actual),1e-20)),
                          material_speed_at_origin=float(np.linalg.norm(g['weak_material'])),h=hh)
                    pr.append(row)
                write(pathdir/'PATH.json',dict(nodes=pr,sampled_path_qualified=all(a['noise_ball'] and a['likelihood_competitive'] for a in pr),
                          exact_fiber_claim=False,qualification='five nodes only; no continuous-between-node certificate'))
            # Source-conditional compensation endpoints using common lower-order path.
            group_rows=[]
            for group in ('A','B','JOINT'):
                z2=np.asarray(j['groups'][group]['z2']);cx2=g['strong_U']@z2
                xp,ma=match_distance(target.prior,base.x,lambda hh:polynomial(base.x,g['v'],hh,cx2),.02)
                if xp is None:group_rows.append(dict(group=group,status='INFEASIBLE'));continue
                p=target.evaluate(xp,'COMMON_GROUP_'+group)
                row=metrics(target,p,base,y,g['weak_material'],.02)
                hh=ma['h'];actual=target.white(p.state.data-base.state.data)
                pp=hh*j['c1']+hh*hh*(j['b2']+j['A']@z2)
                row.update(group=group,status='EXECUTED',h=hh,predicted_group_A=float(np.linalg.norm(pp[:256])),
                           predicted_group_B=float(np.linalg.norm(pp[256:])),prediction_error=float(np.linalg.norm(actual-pp)))
                group_rows.append(row)
            write(folder/'GROUP_FINITE_ENDPOINTS.json',group_rows)
            # A deliberately single-group weak material direction, same U and
            # complete eight-source output. Extra B channels are NOT a new layout.
            group_corrected=[]
            for group,sel in [('A',slice(0,256)),('B',slice(256,512)),('JOINT',slice(None))]:
                dest=folder/'group_A_weak'/group;dest.mkdir(parents=True,exist_ok=True)
                xp,ma=match_distance(target.prior,base.x,lambda hh:base.x+hh*g['v_group_A'],.02)
                if xp is None:group_corrected.append(dict(group=group,status='INFEASIBLE'));continue
                begin=counters(book)
                pp=correct(target,base,y,xp,g['v_group_A'],.02,dest,rowsel=sel)
                if pp is None:group_corrected.append(dict(group=group,status='CORRECTION_INFEASIBLE'));continue
                row=metrics(target,pp,base,y,real_material(g['D']@g['v_group_A']),.02)
                row.update(group=group,status='EXECUTED',cost=diff_counts(counters(book),begin))
                group_corrected.append(row)
            write(folder/'GROUP_A_WEAK_CORRECTED.json',group_corrected)
            # Decoder-only curvature control: same legal decoder and material RMS.
            xp,ma=match_distance(target.prior,base.x,
                 lambda hh:polynomial(base.x,g['v'],hh,j['decoder_only_x2'],j['decoder_only_x3']),.02)
            if xp is not None:
                pp=target.evaluate(xp,'DECODER_ONLY_CURVATURE_CONTROL')
                target.persist_cache(pp,folder/'DECODER_ONLY_CURVATURE_CONTROL.npz')
                write(folder/'DECODER_ONLY_CURVATURE_CONTROL.json',metrics(target,pp,base,y,g['weak_material'],.02))
            write(folder/'COMPLETE.json',dict(parent=parent,anchor_cost=anchor_cost,jet_cost=jet_cost,
               cost=counters(book),controlled_main_seconds=time.perf_counter()-start,endpoint_results=rows,
               online_GT_used=False,new_layout_tested=False,posterior_claim=False))
            print('PARENT_COMPLETE',parent,counters(book),flush=True)
        finally:target.close()


def seal():
    parents=json.loads((OUT/'configs/E2_PARENT_MANIFEST.json').read_text())['locked']
    paths=[]
    for r in parents:
        p=OUT/'E2/parents'/r['parent']/'COMPLETE.json'
        if not p.exists():raise RuntimeError('all eight locked parents required before seal')
        paths.extend(sorted(p.parent.rglob('*.json')))
        paths.extend(sorted(p.parent.rglob('*.npz')))
    write(OUT/'E2/PREDICTIONS_SEALED.json',dict(utc=time.time(),parents=[r['parent'] for r in parents],
       files=stamp_files(paths),no_GT_reader_in_online=True,ledger=totals(),status='SEALED_BEFORE_OFFLINE_SCORING'))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['probe','generate','parent','seal']);parser.add_argument('--parent')
    args=parser.parse_args()
    if args.command=='probe':physical_probe()
    elif args.command=='generate':generate()
    elif args.command=='parent':process_parent(args.parent)
    else:seal()


if __name__=='__main__':main()
