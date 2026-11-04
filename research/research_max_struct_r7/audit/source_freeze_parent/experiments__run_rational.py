"""Disclosed failed-anchor finite-action study; GT is absent from this runner."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.linalg import norm
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r7.experiments.run_case import load_public,dump
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r7.src.rational_search import RationalAnchor
from research.research_max_struct_r7.src.benchmark_solver import PaidStateObserver,refine_observed

def run(case):
    start,cpu=time.perf_counter(),time.process_time();out=ROOT/'RESULTS/rational'/case;out.mkdir(parents=True,exist_ok=False)
    tasks=json.loads((REPO/'research/research_max_struct_r6/configs/primary_final_manifest.json').read_text())['tasks']
    task=next(t for t in tasks if t['case_id']==case and t['arm']=='cie_dual')
    public=load_public(task['public_files']);model=MaxwellModel(public);Y=public['Y']
    with np.load(REPO/task['output_dir']/'prediction.npz') as z:m=z['m'].copy()
    mimage=m.reshape(24,24);energy=np.sum(abs(model.E)**2,axis=1)
    dy,dx=np.gradient(mimage);edges=np.hypot(dx,dy).ravel()
    # Choose blocks with visible illumination/material variation + geometric cover.
    blocks=np.arange(576).reshape(12,2,12,2).transpose(0,2,1,3).reshape(-1,4)
    score=(energy*(1+edges))[blocks].sum(axis=1)
    ids=np.unique(np.r_[np.argsort(score)[-24:],np.linspace(0,143,12,dtype=int)])
    pool=np.unique(blocks[ids].ravel());anchor=RationalAnchor(model,Y,m,pool)
    E=anchor.E[pool];O=anchor.O[:,pool];C=anchor.C[pool]
    Z=O.conj().T@anchor.error@E.conj().T/anchor.scale
    gn=(O.conj().T@O)*(E.conj()@E.T)/anchor.scale
    H=gn.real+2*(C*np.conj(Z)).real # symmetrize residual second derivative
    H=(H+H.T)/2
    positions=public.get('coords',public.get('centers'))[pool]
    dist=np.linalg.norm(positions[:,None]-positions[None,:],axis=2)
    normalized=abs(H-gn.real)/np.sqrt(np.maximum(gn.diagonal().real[:,None]*gn.diagonal().real[None,:],1e-30))
    pairs=[];used=set()
    def add(i,j,kind):
        i,j=sorted((int(pool[i]),int(pool[j])))
        if i!=j and (i,j) not in used:pairs.append(([i],[j],kind));used.add((i,j))
    for mask,number,kind in [(dist>=.16,24,'distant_residual_curvature'),((dist>.001)&(dist<=.06),12,'adjacent_control')]:
        candidates=np.argwhere(np.triu(mask,1));order=np.argsort(normalized[mask & np.triu(np.ones_like(mask,dtype=bool),1)])[::-1]
        for row in candidates[order][:number]:add(*row,kind)
    rng=np.random.default_rng(2026110409)
    while sum(p[2]=='random_pair_control' for p in pairs)<12:
        i,j=rng.choice(len(pool),2,replace=False);add(i,j,'random_pair_control')
    # Two 2x2 patches are rank8, not rank2 despite two scalar amplitudes.
    patch_ids=ids[np.argsort(score[ids])[-12:]]
    patch_pairs=[]
    for a in patch_ids:
        far=[b for b in patch_ids if b>a and np.linalg.norm(positions[[np.flatnonzero(pool==blocks[a,0])[0]]][0]-positions[[np.flatnonzero(pool==blocks[b,0])[0]]][0])>=.16]
        for b in far:patch_pairs.append((blocks[a].tolist(),blocks[b].tolist(),'two_nonlocal_rank8_patches'))
    patch_pairs=patch_pairs[:24]
    records=[];maps={};objects=[]
    for index,(a,b,kind) in enumerate(pairs+patch_pairs):
        action=anchor.action([a,b]);result=action.search(points=9)
        grid=result.pop('grid_points');values=result.pop('grid_Phi')
        maps[f'x_{index}']=grid;maps[f'Phi_{index}']=values
        # Restricted line barrier is descriptive only, not a full-space basin theorem.
        path=np.linspace(0,1,101)[:,None]*result['best_x'][None,:]
        pathvalues=np.array([action.evaluate(x) for x in path])
        result.update(action_id=index,kind=kind,queries=action.queries,path_max_above_origin=float(pathvalues.max()-result['origin_Phi']),material_step_norm=float(norm(action.material(result['best_x'])-m)))
        records.append(result);objects.append(action)
    np.savez_compressed(out/'OBJECTIVE_MAPS.npz',**maps)
    dump(out/'SEARCH_RECORDS.json',records)
    # Always verify the top two distinct families plus continued anchor control.
    selected=[min([r for r in records if r['state_rank']==rank],key=lambda r:r['best_Phi']) for rank in (2,8)]
    endpoints=[('anchor_continuation',m,None)]
    for r in selected:endpoints.append(('rank'+str(r['state_rank'])+'_finite',objects[r['action_id']].material(r['best_x']),r))
    # Fixed same-family local exact search comparator, counted full correction.
    bestlocal=min(records,key=lambda r:r['local_origin_Phi'])
    endpoints.append(('local_exact_origin',objects[bestlocal['action_id']].material(bestlocal['local_origin_x']),bestlocal))
    validations=[]
    for label,candidate,source in endpoints:
        here=out/label;here.mkdir()
        def snapshot(name,row):
            np.savez_compressed(here/(name+'.npz'),**{k:row[k] for k in ('m','j','Q','Y_prediction','gradient') if k in row})
        observed=PaidStateObserver(model,Y,snapshot,[100,1000,3000]);observed.budget_cpu_offset=time.process_time()-cpu
        before=dict(model.cost);st=model.state(candidate);actual=model.physical_score(candidate,Y)
        pred=actual if source is None else source['local_origin_Phi'] if label=='local_exact_origin' else source['best_Phi']
        gap=float(abs(actual-pred))
        # New wiring qualification at a finite endpoint, not frozen algebra repetition.
        if gap>1e-9*(1+abs(actual)):raise RuntimeError('Finite endpoint does not match original Maxwell objective')
        refined,corr=refine_observed(observed,Y,candidate,maximum=3000,accepted_checkpoints=(100,1000,3000))
        observed.save_budget_states(here);snapshot('prediction',observed.best)
        dump(here/'CORRECTOR.json',corr)
        row=dict(label=label,source_action=source,finite_endpoint_original_Phi=actual,prediction_error=gap,corrected_original_Phi=observed.best['objective'],physical_relative_residual=observed.best['physical_relative_residual'],cost={k:model.cost.get(k,0)-before.get(k,0) for k in model.cost},stop=corr['stop'])
        dump(here/'RECEIPT.json',row);validations.append(row)
    # Full-forward comparison on six prespecified already-ranked finite proposals;
    # no new endpoint is selected by this engineering audit or by GT.
    t,c=time.perf_counter(),time.process_time();check=[]
    for r in sorted(records,key=lambda z:z['best_Phi'])[:6]:
        cand=objects[r['action_id']].material(r['best_x']);f=model.physical_score(cand,Y)
        check.append(dict(action_id=r['action_id'],actual_Phi=f,predicted_Phi=r['best_Phi']))
    rec=dict(case_id=case,truth_access='NONE',anchor_original_Phi=anchor.base_data+anchor.base_prior,
        candidate_pool_cells=pool.tolist(),actions=len(records),rank2_actions=len(pairs),rank8_actions=len(patch_pairs),small_objective_queries=sum(r['queries'] for r in records),setup_wall_seconds=anchor.setup_wall,setup_cpu_seconds=anchor.setup_cpu,
        validations=validations,full_forward_six_check=check,full_forward_six_cpu_seconds=time.process_time()-c,full_forward_six_wall_seconds=time.perf_counter()-t,physical_counts=dict(model.cost),physical_dimension=576,
        wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,selection_rule='Original normalized Maxwell Phi plus unchanged smooth/TV prior only',status='SEALED_NO_GT')
    dump(out/'RECEIPT.json',rec);print(json.dumps(dict(case_id=case,anchor=rec['anchor_original_Phi'],endpoints=[(r['label'],r['corrected_original_Phi']) for r in validations],cpu=rec['cpu_seconds'])),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=1):
        for case in ['obs_40ac5b79e8a74ac3','obs_82ac794e228243bb','obs_42cb233e4f224bb3']:run(case)
