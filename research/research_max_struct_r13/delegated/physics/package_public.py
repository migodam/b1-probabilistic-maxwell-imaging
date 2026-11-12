"""Package operator-only metadata with already public observations; no truth."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,numpy as np
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3];sys.path.insert(0,str(REPO))
from research.research_max_struct_r13.src.vector_adapter import DenseVectorModel
from research.research_max_struct_r13.delegated.physics.generate_inputs import dump

def main():
    out=HERE/'new3d/runner_public';out.mkdir(exist_ok=False)
    model=DenseVectorModel();public=model.public_operators();rows=[]
    for path in sorted((HERE/'new3d/data_public').glob('*.npz')):
        with np.load(path,allow_pickle=False) as z:obs={k:z[k] for k in z.files}
        np.savez_compressed(out/path.name,**dict(public,**obs))
        rows.append(dict(case=path.stem,path=str(out/path.name),contains_truth=False,
            G_shape=list(model.G.shape),E_shape=list(model.E.shape),S_shape=list(model.S.shape)))
    dump(out/'PACKAGING_RECEIPT.json',dict(status='COMPLETE',rows=rows,truth_access='NONE'))
if __name__=='__main__':main()
