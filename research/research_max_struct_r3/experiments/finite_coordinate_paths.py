"""Two formula-first finite material path controls; one disclosed R2 case.

Single-cell data stationary roots + boundaries + zero, and nearest-neighbor
two-cell swaps. Exact prior is evaluated on every candidate; this is not an
exact coordinate optimizer of the prior-augmented objective. No GT or grid.
"""
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
import json,sys,time,traceback
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];R2=REPO/'research/research_max_struct_r2'
sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_fbg_r3_experiments.src.full_maxwell_corrector import corrector
from research.research_max_struct_r2.experiments.run_blind import dump,compact,metric

CASE='obs_b20d5a55677d4eca';LOW=.02;HIGH=20.

def real_roots(coeff):
    c=np.trim_zeros(np.asarray(coeff,float),'f')
    if len(c)<2:return []
    roots=np.roots(c/max(abs(c)))
    return [float(x.real)for x in roots if abs(x.imag)<=128*np.finfo(float).eps*max(1.,abs(x.real))]

def resolvent(model,m,cost):
    state=model.state(m)
    with model._physics.operation('FINITE_FULL_RESOLVENT',rhs=model.n,resolvent_RHS=model.n):
        R=la.lu_solve(state['lu'],np.eye(model.n),check_finite=False)
    C=model.G@R;O=model.S@R
    cost['full_inverse_anchors']+=1;cost['G_resolvent_columns']+=model.n;cost['S_resolvent_columns']+=model.n
    return state,C,O

def onecell_scan(model,m,Y,state,C,O,cost):
    r=state['Y']-Y;Q=state['Q'];scale=float(np.vdot(Y,Y).real)
    prior=model.reg.evaluate(m)[0];base=.5*float(np.vdot(r,r).real)/scale+prior
    A=np.einsum('rs,ri,is->i',r.conj(),O,Q,optimize=True)
    B=np.sum(abs(O)**2,axis=0)*np.sum(abs(Q)**2,axis=1)
    c=np.diag(C);a=c.real;b=abs(c)**2;u=A.real;v=B-2*np.real(c.conj()*A)
    cost['all_cell_response_contractions']+=1
    best=None;failures=[];candidates=[]
    for i in range(model.n):
        lo,hi=LOW-m[i],HIGH-m[i]
        roots=real_roots([-u[i]*b[i]-v[i]*a[i],v[i],u[i]])
        cost['real_quadratic_root_problems']+=1
        ds=[(lo,'lower_bound'),(hi,'upper_bound'),(0.,'zero')]+[(d,'data_stationary_root')for d in roots if lo<=d<=hi]
        seen=set()
        for d,kind in ds:
            trial=m.copy();trial[i]=np.clip(m[i]+d,LOW,HIGH);d=float(trial[i]-m[i])
            if d in seen:continue
            seen.add(d);cost['one_cell_candidates']+=1
            den=1-c[i]*d;D=float(abs(den)**2)
            if D==0 or not np.isfinite(D):
                failures.append(dict(cell=i,d=d,error='singular_or_nonfinite_onecell_denominator'));continue
            t=d/den;delta=np.outer(O[:,i],Q[i])*t
            val=.5*float(np.vdot(r+delta,r+delta).real)/scale+float(model.reg.evaluate(trial)[0])
            cost['full_prior_candidate_evaluations']+=1
            candidates.append((i,d,val))
            if not np.isfinite(val):failures.append(dict(cell=i,d=d,error='nonfinite_finite_response_objective'));continue
            if d!=0 and val<base and (best is None or val<best['predicted_Phi']):
                minpoint=np.clip(a[i]/b[i],lo,hi)if b[i]>0 else lo
                minD=float(abs(1-c[i]*minpoint)**2)
                best=dict(coordinate=int(i),d=d,candidate_kind=kind,real_data_stationary_roots=roots,
                    legal_delta_interval=[float(lo),float(hi)],candidate_denominator_squared=D,
                    legal_interval_min_denominator_squared=minD,min_denominator_location=float(minpoint),
                    data_rational_coefficients=dict(u=float(u[i]),v=float(v[i]),a=float(a[i]),b=float(b[i])),
                    predicted_Phi=float(val),material_distance_l2=float(abs(d)),trial=trial)
    return best,base,failures,np.asarray(candidates,float)

def fixed_edges(shape):
    ny,nx=map(int,shape);edges=[]
    for y in range(ny):
        for x in range(nx):
            i=y*nx+x
            if x+1<nx:edges.append((i,i+1))
            if y+1<ny:edges.append((i,i+nx))
    return edges

def determinant_path_min(Csmall,d):
    M=d[:,None]*Csmall
    poly=np.array([1.,-np.trace(M),np.linalg.det(M)],complex)
    mod=np.polynomial.polynomial.polymul(poly,poly.conj()).real
    derivative=np.arange(1,len(mod))*mod[1:]
    roots=real_roots(derivative[::-1]);ts=[0.,1.]+[t for t in roots if 0<=t<=1]
    vals=[float(abs(np.polynomial.polynomial.polyval(t,poly))**2)for t in ts]
    k=int(np.argmin(vals))
    return dict(path_min_determinant_squared=vals[k],path_min_location=float(ts[k]),
                path_stationary_roots=roots,qualification='floating polynomial minimum, not interval enclosure')

def swap_scan(model,m,Y,state,C,O,edges,cost):
    r=state['Y']-Y;Q=state['Q'];scale=float(np.vdot(Y,Y).real)
    base=.5*float(np.vdot(r,r).real)/scale+float(model.reg.evaluate(m)[0])
    best=None;failures=[];candidates=[]
    for index,(i,j)in enumerate(edges):
        cost['swap_candidates']+=1
        d=np.array([m[j]-m[i],m[i]-m[j]])
        if np.all(d==0):cost['identity_swaps']+=1;candidates.append((i,j,base));continue
        cells=np.array([i,j]);small=C[np.ix_(cells,cells)];K=np.eye(2)-d[:,None]*small
        cost['two_cell_small_solves']+=1;cost['two_cell_small_RHS']+=model.ns
        try:z=la.solve(K,d[:,None]*Q[cells],assume_a='gen',check_finite=True)
        except Exception as exc:failures.append(dict(edge=index,cells=[i,j],error=repr(exc)));continue
        delta=O[:,cells]@z;trial=m.copy();trial[i],trial[j]=m[j],m[i]
        val=.5*float(np.vdot(r+delta,r+delta).real)/scale+float(model.reg.evaluate(trial)[0])
        cost['full_prior_candidate_evaluations']+=1;candidates.append((i,j,val))
        if not np.isfinite(val):failures.append(dict(edge=index,cells=[i,j],error='nonfinite_objective'));continue
        if val<base and(best is None or val<best['predicted_Phi']):
            best=dict(edge_index=index,coordinates=[i,j],d=d.tolist(),material_distance_l2=float(la.norm(d)),
                predicted_Phi=float(val),candidate_determinant_squared=float(abs(la.det(K))**2),
                prior_preserved='same functional, re-evaluated after swap; material histogram exactly preserved by proposal',trial=trial)
    if best is not None:
        cells=np.array(best['coordinates']);best.update(determinant_path_min(C[np.ix_(cells,cells)],np.array(best['d'])))
        cost['selected_swap_path_polynomial_problems']+=1
    return best,base,failures,np.asarray(candidates,float)

def run(arm):
    directory=ROOT/'results/development'/arm;directory.mkdir(parents=True,exist_ok=True)
    path=directory/'RECEIPT.json'
    if path.exists():raise FileExistsError('Preserve bounded control disposition')
    start=time.perf_counter();cpu=time.process_time();model=None;cost=Counter();history=[]
    rec=dict(status='RUNNING',run_id=arm,case_id=CASE,classification='FORMULA_DRIVEN_POST_HOC_DEVELOPMENT',
        started_utc=datetime.now(timezone.utc).isoformat(),truth_access='NONE',code='experiments/finite_coordinate_paths.py',
        seed='same scaled_original_restart stationaryendpoint',max_accepted_moves=10,prior='unchanged historical fullprior',
        material_box=[LOW,HIGH],corrector=dict(method='originalphysical L-BFGS-B',maxiter=1000,maxfun=4000),
        selection='strictly lower actual originalPhi; initialincumbent retained; deterministiccell/edgetiebreak',
        threads=1,device='CPU',hash_checks=False)
    dump(path,rec)
    try:
        with np.load(R2/'data_public'/(CASE+'.npz'),allow_pickle=False)as f:public={k:f[k]for k in f.files}
        with np.load(ROOT/'results/development/scaled_original_restart/prediction.npz',allow_pickle=False)as f:m=f['m'].copy()
        initial_m=m.copy();model=MaxwellModel(public);Y=public['Y'];initial=metric(model,Y,m);incumbents=[initial]
        edges=fixed_edges(model.shape)
        rec['candidate_dictionary']='576physicalcells roots/bounds/0'if arm=='onecell_roots'else '1104fixednearestneighbor undirected4connectededges'
        rec['dictionary_size']=model.n if arm=='onecell_roots'else len(edges)
        for iteration in range(10):
            tick=time.perf_counter();ct=time.process_time();state,C,O=resolvent(model,m,cost)
            if arm=='onecell_roots':best,base,failures,table=onecell_scan(model,m,Y,state,C,O,cost)
            elif arm=='neighbor_swaps':best,base,failures,table=swap_scan(model,m,Y,state,C,O,edges,cost)
            else:raise ValueError(arm)
            np.savez_compressed(directory/f'scan_{iteration:02d}.npz',candidate_table=table,anchor_m=m)
            item=dict(iteration=iteration,anchor_original_Phi=base,failed_candidates=failures,
                      cumulative_auxiliary_cost=dict(cost),accepted=False)
            if best is None:item['stop']='NO_ENUMERATED_IMPROVEMENT'
            else:
                trial=best.pop('trial');actual=metric(model,Y,trial);item.update(best,actual_Phi=actual['objective'],
                    prediction_fullforward_difference=actual['objective']-best['predicted_Phi'])
                cost['best_candidate_full_validations']+=1
                if actual['objective']<base:
                    item['accepted']=True;m=trial;incumbents.append(actual)
                    cost['accepted_moves']+=1
                    np.savez_compressed(directory/f'accepted_{iteration:02d}.npz',m=m)
                else:item['stop']='BEST_FORMULA_CANDIDATE_DID_NOT_STRICTLY_DECREASE_ACTUAL_PHI'
            item.update(wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-ct)
            history.append(item);dump(directory/'GREEDY_HISTORY.json',history)
            print(json.dumps(dict(event='FINITE_GREEDY',arm=arm,iteration=iteration,accepted=item['accepted'],Phi=item.get('actual_Phi',base),stop=item.get('stop'))),flush=True)
            if not item['accepted']:break
        seed=min(incumbents,key=lambda c:c['objective'])
        corr=corrector(model,Y,seed['m'],dict(method='lbfgs',max_nfev=1000,box=(LOW,HIGH)))
        dump(directory/'CORRECTOR.json',compact(corr));incumbents.append(metric(model,Y,corr['m']))
        best=min(incumbents,key=lambda c:c['objective']);state=model.state(best['m']);value,g=model.data_value_gradient(best['m'],Y)
        pg=best['m']-np.clip(best['m']-g,LOW,HIGH)
        np.savez_compressed(directory/'prediction.npz',m=best['m'],Y_prediction=state['Y'],j=state['j'],Q=state['Q'],gradient=g,projected_gradient=pg)
        rec.update(status='SEALED_DEVELOPMENT_PREDICTION',finished_utc=datetime.now(timezone.utc).isoformat(),
            prediction=str((directory/'prediction.npz').relative_to(ROOT)),physical_objective=float(value),
            physical_residual=float(la.norm(state['Y']-Y)/la.norm(Y)),physical_whitened_residual=float(la.norm(state['Y']-Y)/float(public['sigma'])),
            projected_gradient_inf=float(la.norm(pg,np.inf)),projected_gradient_l2=float(la.norm(pg)),
            material_distance_from_seed_l2=float(la.norm(best['m']-initial_m)),corrector_stop=corr['stop'],
            greedy_accepted=int(cost['accepted_moves']),scans=len(history),auxiliary_cost=dict(cost),corrector_cost=corr['cost'],
            final_greedy_Phi=seed['objective'],initial_Phi=initial['objective'],
            scope='finite enumerated material proposal control; knownWoodbury/Schur; no globalcoordinateoptwithprior or generalbasinclaim')
    except Exception as exc:rec.update(status='FAILED_DEVELOPMENT_DISPOSITION',error=repr(exc),traceback=traceback.format_exc(),history=history)
    finally:
        rec.update(physical_model_cost=dict(model.cost)if model else {},physical_operation_wall=dict(model.timings)if model else {},
            physical_operation_cpu=dict(model._physics.cpu_timings)if model else {},wall_seconds=time.perf_counter()-start,
            cpu_seconds=time.process_time()-cpu,auxiliary_cost=dict(cost),FFT=0,GPU=0,
            timer_scope='load,allanchors/full576RHSresolvents,allfinitecandidates/fullprior,fullvalidations,corrector,finalgradient,serialization; receiptwriteexcluded')
        dump(path,rec)
    print(json.dumps({k:rec.get(k)for k in ['run_id','status','greedy_accepted','physical_residual','projected_gradient_inf','wall_seconds','physical_model_cost','error']}),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=1):
        for arm in ('onecell_roots','neighbor_swaps'):run(arm)
