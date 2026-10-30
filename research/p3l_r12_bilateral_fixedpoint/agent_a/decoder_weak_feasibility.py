"""Original decoder-only finite trust-radius audit; no physics/GT/observations."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import sys,time,json,csv
from pathlib import Path
import numpy as np
from scipy.linalg import svd
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.p3l_r12_bilateral_fixedpoint.shared.physics import Problem
from research.p3l_r12_bilateral_fixedpoint.agent_a.modes import proposal_basis


def run():
    t=time.perf_counter();pr=Problem();idx=pr.idx;rows=[];details=[]
    base=REPO.parent/'b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/results/exploration_rev1/parents'
    for i in range(8):
        parent=f'R7_LOCKED_{i:03d}';p=base/parent/'B_RAY/RMS_0.0200/ENDPOINT.npz'
        with np.load(p,allow_pickle=False) as z:x=z['x'].copy();J=z['Jx'].copy();L=z['receiver_T'][:,idx].copy();R=z['field'][idx,:].copy()
        _,_,m0,D=pr.decode(x,True)
        for kind in ('full','full_prior','host_prior','pair_prior','receiver_prior','source_prior'):
            ll,rr=(pr.setup['O'],pr.setup.get('ER',pr.setup.get('E'))) if kind.startswith('host') else (L,R)
            W,inf=proposal_basis(J,D,kind,16,ll,rr,return_info=True)
            _,s,vh=svd(W,full_matrices=False)
            # All declared basis columns plus most latent-efficient unit-material combination.
            vectors=[W[:,j] for j in range(W.shape[1])]+[W@vh[-1]]
            actual=[];predicted=[];penalties=[]
            for v in vectors:
                v=v/np.linalg.norm(v)
                for sign in (-1.,1.):
                    dx=sign*4.*v
                    _,_,m1,_=pr.decode(x+dx,False)
                    actual.append(float(np.linalg.norm(m1-m0)/24.))
                    predicted.append(float(np.linalg.norm(D@dx)/24.))
                    penalties.append(float(.5*(np.dot(x+dx,x+dx)-np.dot(x,x))))
            rows.append(dict(parent=parent,kind=kind,weak_rank=W.shape[1],
                             least_latent_norm_per_unitRMS_weakspan=float(min(s)),
                             largest_latent_norm_per_unitRMS_weakspan=float(max(s)),
                             local_max_RMS_at_latent_radius4=float(4./min(s)),
                             tested_signed_rays=len(actual),
                             tested_actual_RMS_min=float(min(actual)),tested_actual_RMS_max=float(max(actual)),
                             rays_reaching_RMS_0p01=sum(a>=.01 for a in actual),
                             rays_reaching_RMS_0p02=sum(a>=.02 for a in actual),
                             tested_prior_penalty_delta_min=float(min(penalties)),tested_prior_penalty_delta_max=float(max(penalties))))
    with (ROOT/'agent_a/DECODER_WEAK_FEASIBILITY.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary=[]
    for kind in ('full','full_prior','host_prior','pair_prior','receiver_prior','source_prior'):
        rr=[r for r in rows if r['kind']==kind]
        summary.append(dict(kind=kind,contexts=len(rr),
                            min_latentnorm_per_unitRMS_range=[min(r['least_latent_norm_per_unitRMS_weakspan'] for r in rr),max(r['least_latent_norm_per_unitRMS_weakspan'] for r in rr)],
                            local_max_RMS_at_radius4_range=[min(r['local_max_RMS_at_latent_radius4'] for r in rr),max(r['local_max_RMS_at_latent_radius4'] for r in rr)],
                            tested_finite_max_RMS_range=[min(r['tested_actual_RMS_max'] for r in rr),max(r['tested_actual_RMS_max'] for r in rr)],
                            tested_rays_reaching_0p01=sum(r['rays_reaching_RMS_0p01'] for r in rr),
                            tested_rays_reaching_0p02=sum(r['rays_reaching_RMS_0p02'] for r in rr)))
    result=dict(status='DECODER_ONLY_COMPLETE',physical_RHS=0,physical_LU=0,GT_arrays_read=False,observation_arrays_read=False,
                historical_context='8 exposed R9 endpoints; finite fixed-ray decoder check, not optimization or proof of whole trust-ball infeasibility',
                trust_radius=4.,directions='16 basis columns and most latent-efficient weakspan direction, both signs',
                wall_seconds=time.perf_counter()-t,summary=summary,details=rows)
    (ROOT/'agent_a/DECODER_WEAK_FEASIBILITY.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='details'},indent=2))
if __name__=='__main__':run()
