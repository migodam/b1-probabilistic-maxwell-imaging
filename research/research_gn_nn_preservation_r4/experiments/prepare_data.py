import os,sys,json
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
STUDY=Path(__file__).resolve().parents[1];sys.path.insert(0,str(STUDY.parents[1]))
from research.research_gn_nn_preservation_r4.src.data import prepare_dataset
if __name__=='__main__':
    r=prepare_dataset(STUDY/'data',STUDY.parent/'research_z_alg_r4_diversity/public_data/operators.npz')
    print(json.dumps({k:v for k,v in r.items() if k!='forward_receipts'}))
