"""Real Maxwell static information/finite-prior/finite-task stress panel.

Known-support fixtures test mechanisms; they are not unknown-wall recovery.
All posterior reference atoms have an explicit finite proper counting prior.
"""
from pathlib import Path
from dataclasses import replace
import gc,json,time
import numpy as np
from bgga.accounting import BudgetBook,_default
from bgga.fixed_physics import FixedEvaluator
from bgga.stress_metrics import stress_metrics,finite_pair_diagnostics
from wslam.forward import Evaluator
from wslam.world import World
from wslam.physics import pack,packj

ROOT=Path(__file__).resolve().parents[1]
SCALES=np.array([1.,.1,.3,.2,1.,.1,.04,.005,.005])

def normweights(log):
    x=log-np.max(log);p=np.exp(x);return p/p.sum()

def fixtures():
    answer=[]
    for j,e in enumerate([2.5,3.75,5.]):
        answer.append(dict(id=f'C1_{j}',family='C1',physical_parent=j,
            support=np.ones((2,2),bool),epsinf=e,ell=.08,minimal=False))
    for family in ['C2','C3']:
        for j,col in enumerate([8,11,14]):
            s=np.zeros((6,24),bool);s[2,col]=True
            answer.append(dict(id=f'{family}_{j}',family=family,physical_parent=3+j,
                support=s,epsinf=4.,ell=.08,minimal=family=='C3'))
    for j,s in enumerate([[[1,1],[1,1]],[[1,0],[1,1]],[[1,1],[0,0]]]):
        answer.append(dict(id=f'C4_{j}',family='C4',physical_parent=6+j,
            support=np.array(s,bool),epsinf=3.75,ell=.125,minimal=True))
    return answer

def finite_prior(fixed,w,case,cache):
    key=case['physical_parent']
    if key not in cache:
        atoms=[];means=[];prior=[]
        walls=[(0.,0.,1/3)]+[(de,dl,1/6) for de in [-.1,.1] for dl in [-.01,.01]]
        poses=[(0.,0.,.5),(-.04,-.02,.25),(.04,.02,.25)]
        for de,dl,pw in walls:
            for dx,dy,pp in poses:
                ww=replace(w,epsinf=float(w.epsinf)+de,ell=float(w.ell)+dl,pose=w.pose+np.array([dx,dy]))
                data=fixed.evaluate(ww,label=f'{case["id"]} bounded finite15 atom').data.copy()
                atoms.append([ww.epsinf,ww.ell,*ww.pose]);means.append(data);prior.append(pw*pp)
        cache[key]=dict(atoms=np.array(atoms),means=np.array(means),prior=np.array(prior))
    return cache[key]

def posterior_summary(atoms,prior,means,y,select):
    predictions=np.array([pack(x.ravel()[select]) for x in means])
    phi=.5*np.sum((predictions-y)**2,axis=1)
    weights=normweights(np.log(prior)-phi)
    mean=weights@atoms;D=atoms-mean;cov=(D.T*weights)@D
    pm=prior@atoms;Pd=atoms-pm;pc=(Pd.T*prior)@Pd
    std=np.array([.1,.01,.04,.02])
    ratio=np.trace(cov[:2,:2]/std[:2,None]/std[None,:2])/np.trace(pc[:2,:2]/std[:2,None]/std[None,:2])
    return dict(weights=weights,mean=mean,covariance=cov,prior_covariance=pc,
        wall_variance_fraction_of_prior=float(ratio),prior_id='BGGA_G0_LOCAL_FINITE15_V1',
        labels='exact defined finite prior; known support/intrinsic/instrument, static target; not arbitrary-wall posterior')

def mean_pair(fixed,w,base_data,B,sel,error):
    wall=replace(w,epsinf=float(w.epsinf)+.1)
    y=pack(base_data.ravel()[sel]);bp=packj(B[sel,:2])
    # Only current data choose the best endpoint; no material truth score.
    u,s,vt=np.linalg.svd(bp,full_matrices=False)
    threshold=float(error);keep=s>threshold
    pose=w.pose.copy();best=None;trials=[];fixed_probe=None
    for k in range(3):
        state=fixed.evaluate(replace(wall,pose=pose),label='declared finite compensation endpoint')
        data=state.data.copy();r=pack(data.ravel()[sel])-y
        if k==0:fixed_probe=data.copy()
        value=float(np.linalg.norm(r));trials.append(dict(trial=k,pose=pose.tolist(),data_distance=value))
        if best is None or value<best[0]:best=(value,replace(wall,pose=pose.copy()),data)
        if not np.any(keep):break
        step=-vt[keep].T@((u[:,keep].T@r)/s[keep])
        pose=np.clip(pose+step,np.array([-1.,1.6]),np.array([1.,2.6]))
    return best[1],best[2],fixed_probe,dict(trials=trials,best_data_distance=best[0],
        frozen_pose_J_rank=int(keep.sum()),threshold=threshold,
        optimality='feasible finite pair only; not global or converged profile optimum')

def main():
    book=BudgetBook();book.authorization('G0')
    out=ROOT/'results/g0_stress';out.mkdir(parents=True,exist_ok=True)
    plan=json.loads((ROOT/'G0_STRESS_PLAN.json').read_text())
    (out/'PLAN.json').write_text(json.dumps(plan,indent=2,ensure_ascii=False)+'\n')
    fixed=FixedEvaluator(book);ref=Evaluator(book,'G0',h=.17,n=384,grid_refinement=3)
    common={};prior_cache={}
    rows=json.loads((out/'RESULTS.json').read_text()) if (out/'RESULTS.json').exists() else []
    for c in fixtures():
        case=c['id'];key=c['physical_parent'];select=np.array([8]) if c['minimal'] else np.arange(32)
        if c['family']=='C4' and not book.declared_completion:
            book.event('G0','predeclared_C4_required_verification_closeout',
                fixture_ids=['C4_0','C4_1','C4_2'],scope='already planned finite pairs, finite priors and all25 future endpoint checks; no new campaign')
            book.declared_completion=True
        w=World(c['support'],c['epsinf'],c['ell'],(.2,2.))
        existing=out/(case+'.npz');saved=None
        if existing.exists():
            with np.load(existing,allow_pickle=False) as z:saved={k:z[k].copy() for k in z.files}
            js=np.column_stack((saved['A'],saved['B'][:,:5]))
            common[key]=dict(data=saved['data'],A=saved['A'],B=saved['B'],Jref=saved['J_reference'],
                error=float(np.linalg.norm(packj(js-saved['J_reference'])*SCALES[:7],2)),
                reference_data=saved['reference_data'],white_mesh=float(np.linalg.norm(pack(saved['data']-saved['reference_data']))),
                duality=dict(status='passed code check; original scalar not preserved in pre-serialization failure'))
            prior_cache[key]=dict(atoms=saved['prior_atoms'],prior=saved['prior_weights'],means=saved['prior_means'])
            oldrow=next((r for r in rows if r['case']==case),None)
            if oldrow is not None and c['family']!='C4':continue
            if oldrow is not None and oldrow['future25_status']=='EXECUTED_ALL25_AND_INDEPENDENT_ARGMAX':continue
        if key not in common:
            st=fixed.evaluate(w,label=case);A,B=fixed.derivatives(st)
            data=st.data.copy();J=np.column_stack((A,B[:,:5]))
            rr=ref.evaluate(w,jacobian=True);Jr=rr.jacobian.copy();rd=rr.data.copy();ref.clear()
            e=float(np.linalg.norm(packj(J-Jr)*SCALES[:7],2))
            common[key]=dict(data=data,A=A,B=B,Jref=Jr,error=e,
                reference_data=rd,white_mesh=float(np.linalg.norm(pack(data-rd))),
                duality=st.derivative_audit)
        raw=common[key];A,B=raw['A'],raw['B'];data=raw['data'];error=raw['error']
        As=packj(A[select])*SCALES[:2];Bs=packj(B[select])*SCALES[2:]
        # K is the full-response current declared target j0. It is NOT the
        # derivative of the whole25-probe max; that family is checked finitely.
        KA=packj(A if c['family']!='C3' else A[select])*SCALES[:2]
        K=np.column_stack((KA,np.zeros((len(KA),7))))
        diagnostics={}
        for name,n in [('pose_only',2),('pose_intrinsic',5),('pose_intrinsic_instrument',7)]:
            kk=np.column_stack((KA,np.zeros((len(KA),n))))
            diagnostics[name]=stress_metrics(As,Bs[:,:n],kk,prior_precision=np.eye(2+n),operator_error=error)
        rng=np.random.default_rng(27102700+key)
        noise=rng.normal(size=data.shape)+1j*rng.normal(size=data.shape)
        y=pack((data+.01*noise).ravel()[select])
        grid=finite_prior(fixed,w,c,prior_cache)
        post=posterior_summary(grid['atoms'],grid['prior'],grid['means'],y,select)
        if saved is None:
            endpoint,endpoint_data,fixed_probe,correction=mean_pair(fixed,w,data,B,select,error)
            # Each finite endpoint receives its own independent fitted-grid check.
            rr=ref.evaluate(endpoint);endpoint_ref=rr.data.copy();ref.clear()
            rr=ref.evaluate(replace(endpoint,pose=w.pose));fixed_probe_ref=rr.data.copy();ref.clear()
        else:
            endpoint_data=saved['endpoint_data'];endpoint_ref=saved['endpoint_reference']
            fixed_probe=saved['same_fixed_probe'];fixed_probe_ref=saved['same_fixed_probe_reference']
            if 'endpoint_pose' in saved:pose=saved['endpoint_pose']
            else:
                # Failed C1_0 had durable raw arrays and both fitted endpoints;
                # recover only its own endpoint pose from the paid event record.
                paid=[r for r in book.records() if r.get('label')=='evaluation_complete'
                    and r.get('world',{}).get('epsinf')==float(w.epsinf)+.1
                    and r.get('world',{}).get('support')==w.support.astype(int).tolist()]
                if len(paid)<2:raise ValueError('Cannot recover missing endpoint metadata')
                pose=np.array(paid[-2]['world']['target_pose'])
            endpoint=replace(w,epsinf=float(w.epsinf)+.1,pose=pose)
            correction=dict(status='raw saved endpoint reused, no new PDE',
                best_data_distance=float(np.linalg.norm(pack((endpoint_data-data).ravel()[select]))),
                original_trial_history='NOT_PRESERVED_BEFORE_SERIALIZATION_FAILURE' if 'endpoint_pose' not in saved else 'previous metadata',
                optimality='feasible finite pair only, not profile optimum')
            if oldrow is not None:
                correction=oldrow['corrector'];raw['duality']=oldrow['duality']
        pair_num=pack((endpoint_data-data).ravel()[select])
        pair_ref=pack((endpoint_ref-raw['reference_data']).ravel()[select])
        task_num=pack((fixed_probe-data).ravel()[select] if c['family']=='C3' else (fixed_probe-data))
        task_ref=pack((fixed_probe_ref-raw['reference_data']).ravel()[select] if c['family']=='C3' else (fixed_probe_ref-raw['reference_data']))
        pair=finite_pair_diagnostics(pair_num,task_num,
            current_error=float(np.linalg.norm(pair_num-pair_ref)),
            task_error=float(np.linalg.norm(task_num-task_ref)),
            task_norm='euclidean',task_contract_id='same fixed declared target j0 in both walls; full25 family audited separately')
        path=out/(case+'.npz')
        np.savez_compressed(path,A=A,B=B,J_reference=raw['Jref'],A_scaled=As,B_scaled=Bs,K_j0=K,
            data=data,reference_data=raw['reference_data'],endpoint_data=endpoint_data,endpoint_reference=endpoint_ref,
            same_fixed_probe=fixed_probe,same_fixed_probe_reference=fixed_probe_ref,
            endpoint_pose=endpoint.pose,
            prior_atoms=grid['atoms'],prior_weights=grid['prior'],prior_means=grid['means'],posterior_weights=post['weights'],
            posterior_covariance=post['covariance'],observation=y,selected_complex_indices=select)
        row=dict(case=case,family=c['family'],physical_parent=key,
            known_support_fixture=True,independent_of_class2_world=c['family']!='C3',
            measurement='TX0 RX8' if c['minimal'] else '2TX all16RX',
            empirical_operator_error=error,error_label='measured scaled J difference, not a uniform bound',
            white_mesh_difference=raw['white_mesh'],duality=raw['duality'],diagnostics=diagnostics,
            finite_posterior={k:v for k,v in post.items() if k!='weights'},
            finite_pair=pair,corrector=correction,endpoint=endpoint.metadata(),
            future25_status='PREDECLARED_CLOSEOUT_PENDING' if c['family']=='C4' else 'not this fixture task',
            category='intended only; inspect measured absolute/spectral/prior/task quantities')
        rows=[r for r in rows if r['case']!=case];rows.append(row)
        (out/'RESULTS.json').write_text(json.dumps(rows,default=_default,indent=2,allow_nan=False)+'\n')
        print(json.dumps(dict(case=case,white_mesh=raw['white_mesh'],operator_error=error,
            wall_variance_fraction=post['wall_variance_fraction_of_prior'],joint_pair_distance=float(np.linalg.norm(pair_num)))),flush=True)
        if c['family']=='C4':
            # This verification was in PLAN before any stress output. At90%
            # only this already-declared endpoint audit may finish, no new task.
            fixed.declared_completion=True
            nodes=np.sort(np.cos(np.pi*np.arange(5)/4));future=[];rawpairs=[]
            checkpoint=out/(case+'_future_partial.npz')
            if checkpoint.exists():
                with np.load(checkpoint,allow_pickle=False) as z:
                    rawpairs=[(tuple(p),a.copy(),b.copy()) for p,a,b in zip(z['poses'],z['first'],z['second'])]
                future=[float(np.linalg.norm(pack(b-a))) for _,a,b in rawpairs]
            probe_index=0
            for yy in 2.1+.5*nodes:
                for xx in nodes:
                    if probe_index<len(rawpairs):probe_index+=1;continue
                    pos=(float(xx),float(yy))
                    f1=fixed.evaluate(replace(w,pose=pos),label=case+' declared future probe first').data.copy()
                    f2=fixed.evaluate(replace(endpoint,pose=pos),label=case+' declared future probe second').data.copy()
                    future.append(float(np.linalg.norm(pack(f2-f1))));rawpairs.append((pos,f1,f2))
                    np.savez_compressed(checkpoint,poses=np.array([p[0] for p in rawpairs]),
                        first=np.array([p[1] for p in rawpairs]),second=np.array([p[2] for p in rawpairs]))
                    probe_index+=1
            maximum=int(np.argmax(future));pos,f1,f2=rawpairs[maximum]
            # Matched independent argmax validation, not borrowed from j0.
            q1=ref.evaluate(replace(w,pose=pos)).data.copy();ref.clear()
            q2=ref.evaluate(replace(endpoint,pose=pos)).data.copy();ref.clear()
            err1=float(np.linalg.norm(pack(q1-f1)));err2=float(np.linalg.norm(pack(q2-f2)))
            row['future25']=dict(distances=future,max=float(max(future)),argmax=maximum,position=pos,
                independent_argmax_distance=float(np.linalg.norm(pack(q2-q1))),
                endpoint_numerical_differences=[err1,err2],
                empirical_margin_lower=float(max(future)-err1-err2),
                norm='max25 whitened total-response norm, not concatenation',
                continuum_certificate=False)
            row['future25_status']='EXECUTED_ALL25_AND_INDEPENDENT_ARGMAX'
            np.savez_compressed(out/(case+'_future25.npz'),poses=np.array([p[0] for p in rawpairs]),
                first=np.array([p[1] for p in rawpairs]),second=np.array([p[2] for p in rawpairs]),
                independent_first=q1,independent_second=q2)
            (out/'RESULTS.json').write_text(json.dumps(rows,default=_default,indent=2)+'\n')
            fixed.declared_completion=False
        gc.collect()
    (out/'COMPLETE.json').write_text(json.dumps(dict(fixtures=12,unique_physical_parents=9,
        status='COMPLETE_AS_DECLARED_NOT_ALL_INTENDED_CATEGORIES_ASSERTED',new_scientific_thresholds=False),indent=2)+'\n')

if __name__=='__main__':main()
