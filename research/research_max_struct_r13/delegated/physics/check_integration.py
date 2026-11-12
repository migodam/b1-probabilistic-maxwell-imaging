"""Exercise public-only model constructor and mature vector auxiliary seam."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,numpy as np
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3];sys.path.insert(0,str(REPO))
from research.research_max_struct_r13.src.vector_adapter import VectorPublicModel,run_vector_auxiliary
from research.research_max_struct_r13.delegated.physics.generate_inputs import dump

def main():
    path=HERE/'new3d/runner_public/slotted_high_frame.npz'
    with np.load(path,allow_pickle=False) as z:public={k:z[k] for k in z.files}
    m=VectorPublicModel(public);Y=public['Y'];z=public['m_init'];audits=[]
    r=run_vector_auxiliary(m,Y,z,'dual_cie',lambda:None,
        lambda x,s,v,meta:audits.append(dict(**meta,Phi=v)),audit_every=1,outer=2,inner=1)
    D,R=m.data_terms(z,Y);v,g=m.data_value_gradient(z,Y)
    dump(HERE/'INTEGRATION_CHECK.json',dict(status='COMPLETE',public_only=True,
        read_path=str(path),auxiliary=r,paid_audits=audits,
        objective_decomposition_absolute_difference=abs(v-D-R),counts=dict(m.cost),
        scope='two-cycle interface exercise only, no algorithm comparison or scientific acceptance'))
    print('Vector public constructor and two-cycle mature Dual seam complete',flush=True)
if __name__=='__main__':main()
