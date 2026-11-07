"""Fixed public-endpoint segments: finite physical barriers, not a new optimizer."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import sys,json,time,csv
import numpy as np
ROOT=Path(__file__).absolute().parents[1]
sys.path.insert(0,str(ROOT.parents[1]))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel

def load(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def main():
    started=time.perf_counter();out=ROOT/'RESULTS/exposed_endpoint_paths';out.mkdir(parents=True,exist_ok=False)
    inventory=json.loads((ROOT/'delegated/diagnostics/STATE_INVENTORY.json').read_text())
    old=ROOT.parents[2]/'r10_worktree/research/research_max_struct_r10'
    cfg=dict(status='FROZEN_BEFORE_PATH_EVALUATIONS',cases=['r10_twod_0001','r10_twod_0004','r10_twod_0022'],
        methods=['cie_single','matrix_green19'],t=[0,.001,.01,.05,.1,.25,.5,.75,1],
        direction='Larger originalPhi endpoint toward smaller originalPhi endpoint; exact ties keep listed order',
        scope='POSTHOC_EXPOSED_DIAGNOSTIC: exposed objects selected before this diagnostic; no GT input or algorithm selection',
        objective='Inherited full scalar TM Maxwell originalPhi including the same R(m)',
        question='Do a legal straight segment and its initial derivative show a physical barrier between already found endpoints? No assertion about all paths or a global minimum.')
    (out/'FROZEN_CONFIG.json').write_text(json.dumps(cfg,indent=2)+'\n')
    rows=[];summary=[];cost=[]
    for cid in cfg['cases']:
        obj=next(o for o in inventory['objects'] if o['case_id']==cid);public={}
        for p in obj['portable_public_files']:public.update(load(p))
        model=MaxwellModel(public);Y=public['Y'];ends=[]
        for method in cfg['methods']:
            d=old/'RESULTS/twod'/f'{cid}__{method}__noise_ball_99'
            m=load(d/'prediction.npz')['m'];meta=json.loads((d/'prediction.json').read_text())
            ends.append(dict(method=method,m=m,saved_Phi=float(meta['objective']),path=str(d/'prediction.npz')))
        ends=sorted(ends,key=lambda v:-v['saved_Phi']);a,b=ends;direction=b['m']-a['m'];pathrows=[]
        mats=[];predictions=[];gradients=[]
        for t in cfg['t']:
            m=a['m']+t*direction;phi,g=model.data_value_gradient(m,Y)
            state=model.state(m);R=float(model.reg.evaluate(m)[0]);projected=m-np.clip(m-g,.02,20)
            row=dict(case_id=cid,t=t,original_Phi=phi,data_term=phi-R,prior=R,
                directional_derivative=float(g@direction),projected_gradient_inf=float(np.max(abs(projected))),
                lower_active_cells=int(np.sum(m<=.02+1e-8)),upper_active_cells=int(np.sum(m>=20-1e-8)),
                material_min=float(m.min()),material_max=float(m.max()),
                relative_data_residual=float(np.linalg.norm(state['Y']-Y)/np.linalg.norm(Y)))
            pathrows.append(row);rows.append(row);mats.append(m);predictions.append(state['Y']);gradients.append(g)
        np.savez_compressed(out/(cid+'.npz'),t=cfg['t'],m=mats,Y_prediction=predictions,physical_gradient=gradients)
        p0,p1=pathrows[0]['original_Phi'],pathrows[-1]['original_Phi']
        summary.append(dict(case_id=cid,from_method=a['method'],to_method=b['method'],from_path=a['path'],to_path=b['path'],
            initial_Phi=p0,final_Phi=p1,initial_segment_derivative=pathrows[0]['directional_derivative'],
            initial_projected_gradient_inf=pathrows[0]['projected_gradient_inf'],
            initial_lower_active_cells=pathrows[0]['lower_active_cells'],initial_upper_active_cells=pathrows[0]['upper_active_cells'],
            largest_sampled_Phi=max(x['original_Phi'] for x in pathrows),
            sampled_barrier_above_start=max(x['original_Phi'] for x in pathrows)-p0,
            interior_samples_above_start=sum(x['original_Phi']>p0 for x in pathrows[1:-1]),
            scope='A finite sampled segment only; no GT, no search, no certificate of local minimality or all-path barrier'))
        cost.append(dict(case_id=cid,counts=dict(model.cost)))
    for name,data in [('ROWS',rows),('SUMMARY',summary),('COSTS',cost)]:
        (out/(name+'.json')).write_text(json.dumps(data,indent=2)+'\n')
    with (out/'ROWS.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    receipt=dict(status='COMPLETE',truth_access='NONE',new_inverse_runs=0,physical_states=27,
        CPU_seconds_near_exit=time.process_time(),wall_seconds=time.perf_counter()-started,
        cost_scope='Whole child CPU includes imports; external process receipt is authoritative when available')
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
