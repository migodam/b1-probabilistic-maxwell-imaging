"""Meaningful local geometry ablation on already paid R9 material endpoints.

Reads only x,Jx,receiver_T,field from eight exposed DEV endpoints. No truth,
observations, nonlinear Maxwell state solver, LU, or linear system solve.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[k]='1'
import sys,time,json,csv
from pathlib import Path
import numpy as np
from scipy.linalg import subspace_angles
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.p3l_r12_bilateral_fixedpoint.shared.physics import Problem
from research.p3l_r12_bilateral_fixedpoint.agent_a.modes import proposal_basis,material_canonical_basis


def run():
    st=time.perf_counter();cpu=time.process_time()
    problem=Problem();idx=problem.idx;setup=problem.setup
    base=REPO.parent/'b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/results/exploration_rev1/parents'
    rows=[];details=[]
    for i in range(8):
        parent=f'R7_LOCKED_{i:03d}'
        path=base/parent/'B_RAY/RMS_0.0200/ENDPOINT.npz'
        if not path.exists():
            details.append(dict(parent=parent,status='MISSING_CACHED_ENDPOINT',path=str(path)));continue
        # Explicit keys. np.load is lazy: no other array is accessed.
        with np.load(path,allow_pickle=False) as z:
            x=z['x'].copy();J=z['Jx'].copy();L=z['receiver_T'][:,idx].copy();R=z['field'][idx,:].copy()
        _,_,_,D=problem.decode(x,derivative=True)
        Wref,info=proposal_basis(J,D,'full',16,return_info=True)
        B=np.r_[D.real,D.imag]/24.
        Eref=B@Wref
        opt=float(np.linalg.norm(J@Wref,'fro')**2/Wref.shape[1])
        for kind in ('full','fullGram','receiver','source','pair','host'):
            t=time.perf_counter()
            ll,rr=(setup['O'],setup.get('ER',setup.get('E'))) if kind=='host' else (L,R)
            W,inf=proposal_basis(J,D,kind,16,ll,rr,return_info=True)
            sec=time.perf_counter()-t
            angles=subspace_angles(Eref,B@W)
            true=float(np.linalg.norm(J@W,'fro')**2/W.shape[1])
            latent=np.linalg.norm(W,axis=0)
            rows.append(dict(parent=parent,kind=kind,material_rank=inf['material_rank'],weak_rank=W.shape[1],
                             true_current_white_cost_per_unit_materialRMS=true,
                             ratio_to_fullj_weak_cost=true/opt,
                             max_material_principal_angle_deg=float(np.degrees(max(angles))),
                             mean_material_principal_angle_deg=float(np.degrees(np.mean(angles))),
                             largest_unitRMS_latent_norm=float(max(latent)),
                             decomposition_seconds=sec))
            details.append(dict(inf,parent=parent,requested_kind=kind,source_path=str(path),read_keys=['x','Jx','receiver_T','field']))
    out=ROOT/'agent_a'
    with (out/'CACHED_MODE_COMPARISON.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary=[]
    for kind in ('full','fullGram','receiver','source','pair','host'):
        rr=[r for r in rows if r['kind']==kind]
        summary.append(dict(kind=kind,contexts=len(rr),
                            true_cost_ratio_range=[min(r['ratio_to_fullj_weak_cost'] for r in rr),max(r['ratio_to_fullj_weak_cost'] for r in rr)],
                            max_angle_range_deg=[min(r['max_material_principal_angle_deg'] for r in rr),max(r['max_material_principal_angle_deg'] for r in rr)],
                            median_decomposition_seconds=float(np.median([r['decomposition_seconds'] for r in rr]))))
    receipt=dict(status='COMPLETE_ZERO_NEW_PHYSICS',contexts=len(rows)//6,comparisons=len(rows),
                 historical_status='previously exposed DEV endpoints, not fresh confirmation',
                 metrics='local current-J response cost per unit material RMS; material weak-subspace angles; decomposition cost',
                 physical_LU=0,physical_RHS=0,GT_arrays_read=False,observations_read=False,
                 wall_seconds=time.perf_counter()-st,CPU_seconds=time.process_time()-cpu,
                 summary=summary,details=details)
    (out/'CACHED_MODE_COMPARISON.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='details'},indent=2))

if __name__=='__main__':run()
