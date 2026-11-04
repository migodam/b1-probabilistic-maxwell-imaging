"""One structural refinement: rank32 phase/curvature-coherent finite actions.

The rank2/rank8 failure is retained. All proposals use the disclosed anchor,
source/receiver residual, legal box, and full physical Hessian; no GT.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r7.experiments.run_case import load_public,dump
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r7.src.rational_search import RationalAnchor
from research.research_max_struct_r7.src.benchmark_solver import PaidStateObserver,refine_observed

def run(case):
    start,cpu=time.perf_counter(),time.process_time();out=ROOT/'RESULTS/rational_coherent'/case;out.mkdir(parents=True,exist_ok=False)
    tasks=json.loads((REPO/'research/research_max_struct_r6/configs/primary_final_manifest.json').read_text())['tasks']
    task=next(t for t in tasks if t['case_id']==case and t['arm']=='cie_dual')
    public=load_public(task['public_files']);model=MaxwellModel(public);Y=public['Y']
    with np.load(REPO/task['output_dir']/'prediction.npz') as z:m=z['m'].copy()
    anchor=RationalAnchor(model,Y,m,np.arange(model.n));E=anchor.E;O=anchor.O;C=anchor.C
    Z=O.conj().T@anchor.error@E.conj().T/anchor.scale
    GN=((O.conj().T@O)*(E.conj()@E.T)).real/anchor.scale
    T=C*np.conj(Z);H=GN+(T+T.T).real
    reg=model.reg;d=reg.D@m
    w=(2*reg.smooth+reg.tv*reg.epsilon**2/(d*d+reg.epsilon**2)**1.5)/reg.edges
    H+=(reg.D.T@reg.D.multiply(w[:,None])).toarray();H=(H+H.T)/2
    free=(m>.02+1e-5)&(m<20-1e-5);ix=np.flatnonzero(free)
    vals,vecs=eigh(H[np.ix_(ix,ix)],subset_by_index=(0,min(3,len(ix)-1)))
    vectors=[]
    for k in range(vecs.shape[1]):
        h=np.zeros(model.n);h[ix]=vecs[:,k];vectors.append(('weak_physical_Hessian_'+str(k),h))
    g=np.real(np.sum(O.conj()*(anchor.error@E.conj().T),axis=0))/anchor.scale+reg.evaluate(m)[1]
    vectors.append(('negative_physical_gradient',-g))
    # Material interface displacement, inferred only from current reconstruction.
    gy,gx=np.gradient(m.reshape(24,24))
    vectors.extend([('current_material_x_translation',gx.ravel()),('current_material_y_translation',gy.ravel())])
    rng=np.random.default_rng(2026110411)
    vectors.extend([('same_rank_random_'+str(i),rng.normal(size=model.n)) for i in range(3)])
    rows=[];actions=[]
    for name,h in vectors:
        for sign in (1,-1):
            direction=sign*h
            pos=np.flatnonzero((direction>0)&(m<19.9));neg=np.flatnonzero((direction<0)&(m>.12))
            if len(pos)<16 or len(neg)<16:continue
            a=pos[np.argsort(abs(direction[pos]))[-16:]];b=neg[np.argsort(abs(direction[neg]))[-16:]]
            wa=(20.-m[a]);wb=-(m[b]-.02)
            action=anchor.action([a,b],[wa,wb]);action.bounds=[(0.,1.),(0.,1.)]
            result=action.search(points=7);grid=result.pop('grid_points');values=result.pop('grid_Phi')
            np.savez_compressed(out/f'map_{len(rows):02d}.npz',x=grid,Phi=values)
            result.update(name=name,sign=sign,weights=[wa.tolist(),wb.tolist()],action_id=len(rows))
            path=np.linspace(0,1,101)[:,None]*result['best_x'];pathvalues=np.array([action.evaluate(x) for x in path])
            result.update(path_max_above_origin=float(pathvalues.max()-result['origin_Phi']),queries=action.queries)
            rows.append(result);actions.append(action)
    dump(out/'SEARCH_RECORDS.json',rows)
    chosen=min(rows,key=lambda r:r['best_Phi']);candidate=actions[chosen['action_id']].material(chosen['best_x'])
    def snapshot(name,row):np.savez_compressed(out/(name+'.npz'),**{k:row[k] for k in ('m','j','Q','Y_prediction','gradient') if k in row})
    observed=PaidStateObserver(model,Y,snapshot,[100,1000,3000]);observed.budget_cpu_offset=time.process_time()-cpu
    state=model.state(candidate);actual=model.physical_score(candidate,Y);gap=abs(actual-chosen['best_Phi'])
    if gap>1e-9:raise RuntimeError('rank32 reduced endpoint disagrees with original full physics')
    refined,corr=refine_observed(observed,Y,candidate,maximum=3000,accepted_checkpoints=(100,1000,3000))
    observed.save_budget_states(out);snapshot('prediction',observed.best);dump(out/'CORRECTOR.json',corr)
    result=dict(case_id=case,truth_access='NONE',status='SEALED_NO_GT',free_face_dimension=len(ix),lowest_free_physical_Hessian_eigenvalues=vals,original_Phi=anchor.base_data+anchor.base_prior,best_finite_Phi=actual,finite_prediction_error=gap,corrected_Phi=observed.best['objective'],physical_relative_residual=observed.best['physical_relative_residual'],selected=chosen,actions=len(rows),small_objective_queries=sum(r['queries'] for r in rows),physical_counts=dict(model.cost),wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,physical_dimension=576,full_rank_setup='576 selected resolvent columns, explicitly paid; rank32 candidate solves')
    dump(out/'RECEIPT.json',result);print(json.dumps({k:result[k] for k in ('case_id','original_Phi','best_finite_Phi','corrected_Phi','cpu_seconds')}),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=1):
        for case in ['obs_40ac5b79e8a74ac3','obs_82ac794e228243bb','obs_42cb233e4f224bb3']:run(case)
