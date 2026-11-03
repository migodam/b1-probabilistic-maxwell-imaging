"""All unordered exact two-cell exchanges: one finite material neighborhood.

No distance/gradient/GT screening, parameter grid or changed prior. Batching
limits temporary memory only. Two-by-two response plus rank-two Gram score
avoids constructing a receiver-by-source response for every pair.
"""
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
import json,sys,time,traceback,resource
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];R2=REPO/'research/research_max_struct_r2'
sys.path.insert(0,str(REPO));sys.path.insert(0,str(ROOT/'experiments'))
from finite_coordinate_paths import resolvent,determinant_path_min
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_fbg_r3_experiments.src.full_maxwell_corrector import corrector
from research.research_max_struct_r2.experiments.run_blind import metric,dump,compact

CASE='obs_b20d5a55677d4eca';BATCH=8192

def prior_cache(m,reg):
    D=reg.D.tocsc();degree=max(np.diff(D.indptr));ids=np.full((len(m),degree),-1,int)
    coeff=np.zeros((len(m),degree))
    for i in range(len(m)):
        a,b=D.indptr[i:i+2];ids[i,:b-a]=D.indices[a:b];coeff[i,:b-a]=D.data[a:b]
    diff=reg.D@m
    return dict(ids=ids,coeff=coeff,diff=diff,smooth=reg.smooth,tv=reg.tv,epsilon=reg.epsilon,edges=reg.edges)

def prior_increment(i,j,delta,cache):
    """Exact full historical prior difference; only incident edges change."""
    ids=cache['ids'];coeff=cache['coeff'];diff=cache['diff']
    ei,ej=ids[i],ids[j];mi,mj=ei>=0,ej>=0
    zi=diff[np.maximum(ei,0)];zj=diff[np.maximum(ej,0)]
    ci=delta[:,None]*coeff[i];cj=-delta[:,None]*coeff[j]
    def energy(x):
        return cache['smooth']*x*x+cache['tv']*(np.sqrt(x*x+cache['epsilon']**2)-cache['epsilon'])
    change=np.sum((energy(zi+ci)-energy(zi))*mi,axis=1)+np.sum((energy(zj+cj)-energy(zj))*mj,axis=1)
    # An edge shared by the two cells must receive BOTH endpoint increments
    # once, rather than two separate one-endpoint energy changes.
    same=(ei[:,:,None]==ej[:,None,:])&mi[:,:,None]&mj[:,None,:]
    z=zi[:,:,None];di=ci[:,:,None];dj=cj[:,None,:]
    correction=energy(z+di+dj)-energy(z+di)-energy(z+dj)+energy(z)
    change+=np.sum(correction*same,axis=(1,2))
    return change/cache['edges']

def subblocks(M,i,j):
    return np.stack([M[i,i],M[i,j],M[j,i],M[j,j]],axis=1).reshape(-1,2,2)

def scan(model,m,Y,state,C,O,pairs,cost):
    tick=time.perf_counter();ct=time.process_time();r=state['Y']-Y;Q=state['Q']
    scale=float(np.vdot(Y,Y).real);prior=float(model.reg.evaluate(m)[0]);base=.5*float(np.vdot(r,r).real)/scale+prior
    A=(O.conj().T@r)@Q.conj().T;BO=O.conj().T@O;BQ=Q@Q.conj().T
    cost['rank2_cross_Gram_builds']+=1;cost['rank2_O_Gram_builds']+=1;cost['rank2_Q_Gram_builds']+=1
    cache=prior_cache(m,model.reg);values=np.full(len(pairs[0]),base);best=None;failures=[];batch_peak=0
    for begin in range(0,len(pairs[0]),BATCH):
        end=min(begin+BATCH,len(pairs[0]));ii=pairs[0][begin:end];jj=pairs[1][begin:end]
        delta=m[jj]-m[ii];active=delta!=0
        cost['all_pair_candidates']+=len(ii);cost['identity_swaps']+=int(np.count_nonzero(~active))
        ids=np.arange(begin,end)[active];i,j,d=ii[active],jj[active],delta[active]
        if len(i)==0:continue
        ca=subblocks(C,i,j);da=np.column_stack([d,-d]);mat=np.eye(2)[None]-da[:,:,None]*ca
        rhs=np.zeros_like(mat);rhs[:,0,0]=d;rhs[:,1,1]=-d
        cost['two_cell_small_solves']+=len(i);cost['two_cell_small_RHS']+=2*len(i);cost['small_solve_batch_calls']+=1
        try:K=np.linalg.solve(mat,rhs)
        except Exception as exc:
            failures.append(dict(batch=[begin,end],error=repr(exc),action='perpairfallback'))
            K=np.full_like(mat,np.nan)
            for z in range(len(i)):
                cost['fallback_small_solves']+=1;cost['fallback_small_RHS']+=2
                try:K[z]=np.linalg.solve(mat[z],rhs[z])
                except Exception as one:failures.append(dict(cells=[int(i[z]),int(j[z])],error=repr(one)))
        ab=subblocks(A,i,j);ob=subblocks(BO,i,j);qb=subblocks(BQ,i,j)
        cross=2*np.real(np.einsum('bij,bij->b',ab.conj(),K))
        norm=np.real(np.trace(K.conj().transpose(0,2,1)@ob@K@qb,axis1=1,axis2=2))
        pd=prior_increment(i,j,d,cache);val=base+.5*(cross+norm)/scale+pd
        cost['rank2_Gram_scores']+=len(i);cost['exact_incident_edge_prior_scores']+=len(i)
        values[ids]=val;valid=np.isfinite(val)
        if np.any(valid):
            candidates=np.where(valid,val,np.inf);z=int(np.argmin(candidates))
            if val[z]<base and(best is None or val[z]<best['predicted_Phi']):
                a,b=int(i[z]),int(j[z]);trial=m.copy();trial[a],trial[b]=m[b],m[a]
                best=dict(pair_index=int(ids[z]),coordinates=[a,b],d=[float(d[z]),float(-d[z])],trial=trial,
                    predicted_Phi=float(val[z]),predicted_prior_change=float(pd[z]),predicted_norm_data_change=float(cross[z]+norm[z]),
                    candidate_determinant_squared=float(abs(np.linalg.det(mat[z]))**2),
                    material_distance_l2=float(np.sqrt(2)*abs(d[z])),
                    grid_manhattan_distance=int(abs(a//24-b//24)+abs(a%24-b%24)))
        if np.any(~valid):failures.append(dict(batch=[begin,end],nonfinite_scores=int(np.count_nonzero(~valid))))
        batch_peak=max(batch_peak,sum(x.nbytes for x in (ca,da,mat,rhs,K,ab,ob,qb,cross,norm,pd,val)))
    if best is not None:
        cells=np.array(best['coordinates']);best.update(determinant_path_min(C[np.ix_(cells,cells)],np.array(best['d'])))
        cost['selected_path_polynomial_problems']+=1
    return best,base,values,dict(wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-ct,
        measured_listed_batch_array_peak_bytes=batch_peak,shared_rank2_Gram_bytes=A.nbytes+BO.nbytes+BQ.nbytes,
        failures=failures,batch_size=BATCH,scope='memoryfigures count listedarrays, not allocatorcompletepeak')

def main():
    directory=ROOT/'results/development/all_pair_swaps';directory.mkdir(parents=True,exist_ok=True)
    path=directory/'RECEIPT.json'
    if path.exists():raise FileExistsError('Preserve one authorized allpair disposition')
    start=time.perf_counter();cpu=time.process_time();model=None;cost=Counter();history=[]
    rec=dict(status='RUNNING',run_id='all_pair_swaps',case_id=CASE,classification='FORMULA_DRIVEN_POST_HOC_DEVELOPMENT',
        code='experiments/all_pair_swaps.py',started_utc=datetime.now(timezone.utc).isoformat(),truth_access='NONE',
        seed='same scaled_original_restart endpoint',dictionary='all165600unordered distinctcellpairs, no distanceorgradientfilter',
        batch_size=BATCH,batch_semantics='memorypartition only; allpairs scored',max_accepted_moves=10,
        prior='exact unchanged historical fullprior, incidentedgeincrement with sharededgecorrection',
        material_box=[.02,20.],corrector=dict(method='originalphysical L-BFGS-B',maxiter=1000,maxfun=4000),
        selection='lowest enumerated originalPhi, physicalfullforwardverificationbeforestrictaccept; retaininitial',
        threads=1,device='CPU',hash_checks=False)
    dump(path,rec)
    try:
        with np.load(R2/'data_public'/(CASE+'.npz'),allow_pickle=False)as f:public={k:f[k]for k in f.files}
        with np.load(ROOT/'results/development/scaled_original_restart/prediction.npz',allow_pickle=False)as f:m=f['m'].copy()
        model=MaxwellModel(public);Y=public['Y'];initial_m=m.copy();initial=metric(model,Y,m);incumbents=[initial]
        pairs=np.triu_indices(model.n,k=1)
        for it in range(10):
            tick=time.perf_counter();ct=time.process_time();state,C,O=resolvent(model,m,cost)
            best,base,values,detail=scan(model,m,Y,state,C,O,pairs,cost)
            np.savez_compressed(directory/f'scan_{it:02d}.npz',anchor_m=m,pair_i=pairs[0],pair_j=pairs[1],candidate_Phi=values)
            row=dict(iteration=it,anchor_Phi=base,candidates=len(values),scan_cost=detail,accepted=False)
            if best is None:row['stop']='NO_ENUMERATED_SWAP_IMPROVEMENT'
            else:
                trial=best.pop('trial');actual=metric(model,Y,trial);cost['best_candidate_full_validations']+=1
                pd_actual=float(model.reg.evaluate(trial)[0]-model.reg.evaluate(m)[0])
                row.update(best,actual_Phi=actual['objective'],prediction_fullforward_difference=actual['objective']-best['predicted_Phi'],
                           local_prior_full_prior_difference=best['predicted_prior_change']-pd_actual)
                if actual['objective']<base:
                    m=trial;row['accepted']=True;cost['accepted_moves']+=1;incumbents.append(actual)
                    np.savez_compressed(directory/f'accepted_{it:02d}.npz',m=m)
                else:row['stop']='BEST_FINITE_SCORE_DID_NOT_DECREASE_FULL_PHI'
            row.update(wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-ct,cumulative_auxiliary_cost=dict(cost))
            history.append(row);dump(directory/'GREEDY_HISTORY.json',history)
            print(json.dumps(dict(event='ALL_PAIR_SCAN',iteration=it,accepted=row['accepted'],Phi=row.get('actual_Phi',base),coordinates=row.get('coordinates'),stop=row.get('stop'))),flush=True)
            if not row['accepted']:break
        seed=min(incumbents,key=lambda c:c['objective'])
        corr=corrector(model,Y,seed['m'],dict(method='lbfgs',max_nfev=1000,box=(.02,20.)));dump(directory/'CORRECTOR.json',compact(corr))
        incumbents.append(metric(model,Y,corr['m']));best=min(incumbents,key=lambda c:c['objective'])
        state=model.state(best['m']);value,g=model.data_value_gradient(best['m'],Y);pg=best['m']-np.clip(best['m']-g,.02,20.)
        np.savez_compressed(directory/'prediction.npz',m=best['m'],Y_prediction=state['Y'],j=state['j'],Q=state['Q'],gradient=g,projected_gradient=pg)
        rec.update(status='SEALED_DEVELOPMENT_PREDICTION',finished_utc=datetime.now(timezone.utc).isoformat(),
            prediction=str((directory/'prediction.npz').relative_to(ROOT)),physical_objective=float(value),
            physical_residual=float(np.linalg.norm(state['Y']-Y)/np.linalg.norm(Y)),
            physical_whitened_residual=float(np.linalg.norm(state['Y']-Y)/float(public['sigma'])),
            projected_gradient_inf=float(np.linalg.norm(pg,np.inf)),projected_gradient_l2=float(np.linalg.norm(pg)),
            accepted_moves=int(cost['accepted_moves']),scans=len(history),initial_Phi=initial['objective'],final_greedy_Phi=seed['objective'],
            material_distance_from_seed_l2=float(np.linalg.norm(best['m']-initial_m)),corrector_stop=corr['stop'],corrector_cost=corr['cost'],
            qualification='best finite enumerated swaps only; swap preserveshistogram beforecorrector, not a globalmaterial or multisteppath guarantee')
    except Exception as exc:rec.update(status='FAILED_DEVELOPMENT_DISPOSITION',error=repr(exc),traceback=traceback.format_exc(),history=history)
    finally:
        peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        rec.update(wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,
            physical_model_cost=dict(model.cost)if model else {},auxiliary_cost=dict(cost),FFT=0,GPU=0,
            process_peak_memory_bytes=int(peak if sys.platform=='darwin' else 1024*peak),
            memory_scope='processlifetimepeak includesimports; listedbatcharrays separatelyrecorded',
            timer_scope='load,fullresolvents,allpairs/priors,fullverification,corrector,finalgradient,serialization; ownreceiptwriteexcluded')
        dump(path,rec)
    print(json.dumps({k:rec.get(k)for k in ['run_id','status','accepted_moves','scans','physical_residual','projected_gradient_inf','wall_seconds','physical_model_cost','auxiliary_cost','error']}),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=1):main()
