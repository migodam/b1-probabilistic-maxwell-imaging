"""Online proposal/correction: observation, support and discovered history only.

Proxy optimizers never accept truth, teacher or reference. Formal archive
entries require independently charged original Mie endpoints. No probabilities.
"""
import time
import numpy as np
from scipy.optimize import least_squares
from b1_joint.material import prior_sample,geometry
from .target import BOUNDS,TAU,to_p,profile
from .material_metrics import distance,structural_difference,nested_select

EDIT_BOXES=[(0,0,0,0),(.05,.13,0,0),(.14,.21,0,0),(.22,.2925518253046813,0,0),
            (0,0,.025,.06),(0,0,.065,.11702073012187252),(.14,.2925518253046813,.025,.06),(.05,.21,.065,.11702073012187252)]

def correct(target,p,y,sigma,box=None,max_nfev=25):
    lo=BOUNDS[:,0].copy();hi=BOUNDS[:,1].copy()
    if box is not None:lo[:2]=[box[0],box[2]];hi[:2]=[box[1],box[3]]
    free=np.flatnonzero(hi>lo);fixed=np.clip(np.array(p,float),lo,hi)
    start=fixed[free];width=(hi-lo)[free]
    def residual(q):
        z=fixed.copy();z[free]=lo[free]+q*width
        return profile(target.predict_proxy(z)[0],y,sigma)[1]
    t=time.perf_counter();before=target.proxy.evaluations
    fit=least_squares(residual,(start-lo[free])/width,bounds=(np.zeros(len(free)),np.ones(len(free))),max_nfev=max_nfev,ftol=1e-7,xtol=1e-7,gtol=1e-7)
    result=fixed.copy();result[free]=lo[free]+fit.x*width
    return result,dict(proxy_residual=float(np.linalg.norm(fit.fun)),nfev=fit.nfev,proxy_evaluations=target.proxy.evaluations-before,wall_seconds=time.perf_counter()-t,converged=bool(fit.success),protected_geometry=box)

def basic_reconstruction(target,y,sigma,rng):
    p=np.array([to_p(x) for x in prior_sample(rng,96)])
    r=np.linalg.norm(profile(target.predict_proxy(p),y,sigma)[1],axis=1)
    trials=[correct(target,z,y,sigma) for z in p[np.argsort(r)[:4]]]
    return min(trials,key=lambda a:a[1]['proxy_residual'])

def endpoint(target,p,y,sigma,kind,history):
    f,_,_,info=target.forward(p,kind=kind)
    g,r=profile(f,y,sigma);norm=float(np.linalg.norm(r));valid=norm<=TAU and info.get('Qabs',-1)>=-1e-10 and info.get('boundary_backward_residual',1)<1e-10
    return dict(p=np.asarray(p).tolist(),gain=float(g),prediction=f.tolist(),fitted_prediction=(float(g)*f).tolist(),residual=norm,valid=bool(valid),boundary_residual=info.get('boundary_backward_residual'),Qabs=info.get('Qabs'),history=history,posterior_weight=None,physics='original full vector-Mie n14',cost=dict(new_material_frequency_attempts=info.get('new_physics_attempts',1),boundary_RHS=info.get('dense_boundary_solves',0)))

def dedup_archive(archive,c):
    if not c['valid']:return 'DATA_INCOMPATIBLE'
    if any(distance(c['p'],a['p'])<.03 for a in archive):return 'DUPLICATE_MATERIAL'
    archive.append(c);return 'NEW_MATERIAL'

def search(target,y,sigma,base,arm,n,rng,*,revision=1,mdn_draws=None):
    archive=[base.copy()] if base['valid'] else [];attempts=[];history=[]
    successes=np.zeros(len(EDIT_BOXES));uses=np.zeros(len(EDIT_BOXES));near=np.zeros(len(EDIT_BOXES));before=target.proxy.evaluations;t0=time.perf_counter()
    population=np.array([to_p(x) for x in prior_sample(rng,max(256,n))])
    pop_r=np.linalg.norm(profile(target.predict_proxy(population),y,sigma)[1],axis=1)
    for t in range(n):
        target.budget.check();edit_id=None;box=None;parent=archive[t%len(archive)] if archive else base
        if arm in ('finite_static','adaptive'):
            if t%5==4: # fixed 20% global, non-greedy prior exploration
                p=to_p(prior_sample(rng,1)[0]);source='global-prior'
            else:
                if arm=='finite_static' or revision==0:edit_id=t%len(EDIT_BOXES)
                elif t<(2 if revision>=3 else len(EDIT_BOXES)):edit_id=t
                else:
                    # Reward only charged, valid, genuinely new material. Failed
                    # categories and repeated families reduce the next score.
                    ucb=(successes+.5)/(uses+1)+.25*np.sqrt(np.log(t+1)/(uses+1))
                    if revision>=2:
                        # Two cheap probes per region prevent wasting full solves
                        # on already impossible annular edits. No hard rejection.
                        forecasts=[]
                        for b in EDIT_BOXES:
                            zz=np.array(parent['p']);zz[:2]=[(b[0]+b[1])/2,(b[2]+b[3])/2]
                            zz,_=correct(target,zz,y,sigma,b,max_nfev=8)
                            rr=float(np.linalg.norm(profile(target.predict_proxy(zz)[0],y,sigma)[1]))
                            novel=min([distance(zz,a['p']) for a in archive] or [1.])
                            forecasts.append(novel/(1+max(rr-TAU,0)**2))
                        ucb+=np.asarray(forecasts)
                    edit_id=int(np.argmax(ucb))
                box=EDIT_BOXES[edit_id];p=np.array(parent['p']);p[:2]=[rng.uniform(box[0],box[1]),rng.uniform(box[2],box[3])];source='protected-finite-edit'
        elif arm=='multistart':p=to_p(prior_sample(rng,1)[0]);source='independent-prior-restart'
        elif arm=='population':
            # Residual-ranked tempered population with material repulsion;
            # not Bayesian SMC and not a posterior sample.
            if archive:
                novelty=np.array([min(distance(z,a['p']) for a in archive) for z in population])
                score=pop_r-.35*np.log(novelty+1e-3)
            else:score=pop_r
            j=int(np.argmin(score));p=population[j].copy();population[j]=to_p(prior_sample(rng,1)[0]);pop_r[j]=np.linalg.norm(profile(target.predict_proxy(population[j])[0],y,sigma)[1]);source='tempered-repulsive-population'
        elif arm=='mdn':p=to_p(mdn_draws[t]);source='frozen-C1-full-distill-proposal'
        else:raise ValueError(arm)
        raw=target.predict_proxy(p)[0];rawres=float(np.linalg.norm(profile(raw,y,sigma)[1]))
        corrected,receipt=correct(target,p,y,sigma,box)
        target.budget.context=dict(arm=arm,step=t)
        c=endpoint(target,corrected,y,sigma,'SEARCH_ENDPOINT',dict(parent_p=parent['p'],proposal_p=p.tolist(),source=source,edit_id=edit_id,corrector=receipt,raw_proxy_residual=rawres,step=t,revision=revision,prior_support_unchanged=True))
        status=dedup_archive(archive,c)
        if edit_id is not None:
            uses[edit_id]+=1;successes[edit_id]+=status=='NEW_MATERIAL'
        c['discovery_status']=status;attempts.append(c)
        history.append(dict(step=t+1,distinct_worlds=len(archive),structural_worlds=structure_count(archive),status=status,cumulative_endpoint_calls=t+1,wall_seconds=time.perf_counter()-t0))
    selected=nested_select(archive)
    return dict(arm=arm,archive=archive,attempts=attempts,selected_indices=selected,selected_online_index=int(np.argmin([c['residual'] for c in archive])) if archive else None,history=history,endpoint_calls=n,wall_seconds=time.perf_counter()-t0,proxy_evaluations=target.proxy.evaluations-before,edit_successes=successes.tolist(),edit_uses=uses.tolist(),weight_reliability='candidate_only',weights=None)

def structure_count(archive):
    reps=[]
    for c in archive:
        if all(structural_difference(c['p'],a['p']) for a in reps):reps.append(c)
    return len(reps)
