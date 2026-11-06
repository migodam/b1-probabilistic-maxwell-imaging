"""One public saved-field coupling audit; no truth, target design or solves."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
import time,json,sys,argparse
from pathlib import Path
import numpy as np
from scipy import linalg as la
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.src.vector_grid import VectorGrid

def main():
    p=argparse.ArgumentParser();p.add_argument('--basis',required=True);p.add_argument('--output',required=True)
    p.add_argument('--kind',choices=['spectral','source'],default='spectral');args=p.parse_args()
    wall,cpu=time.perf_counter(),time.process_time();out=ROOT/args.output;out.mkdir(parents=True,exist_ok=False)
    old=REPO.parent/'r10_worktree/research/research_max_struct_r10'
    cfg=json.loads((REPO/'research/research_max_struct_r10/configs/vector_actual_contract.json').read_text())
    model=VectorGrid(12,.6,cfg['acquisition'])
    with np.load(args.basis,allow_pickle=False) as z:U=z['right_vectors'];GU=z['images'];s=z['singular_values']
    B=U@(U.conj().T@GU);rows=[]
    for case in ('v04','v09'):
        for label in ('initial_auxiliary_joint','stage3_joint','prediction'):
            path=old/'RESULTS/vector'/(case+'_GlobalDual')/(label+'.npz')
            with np.load(path,allow_pickle=False) as z:m=z['material'];J=z['J'];Q=z['Q']
            fields=dict(J=J,Q=Q,defect=J-np.repeat(m,3)[:,None]*Q)
            ws=1/np.sum(abs(model.inc)**2,axis=0)
            row=dict(case=case,label=label,source_file=str(path),fields={})
            for key,X in fields.items():
                coordinate=U.conj().T@X;projected=U@coordinate
                green=model.green(X);right=GU@coordinate;gal=B@coordinate
                norm=lambda v:float(np.sum(abs(v)**2*ws))
                den=norm(X);gden=norm(green)
                row['fields'][key]=dict(weighted_current_energy_capture=norm(projected)/den if den else None,
                    right_G_action_relative_error=np.sqrt(norm(green-right)/gden) if gden else None,
                    galerkin_G_action_relative_error=np.sqrt(norm(green-gal)/gden) if gden else None)
            rows.append(row)
    receipt=dict(status='PUBLIC_SAVED_FIELDS_ONLY',basis=str(Path(args.basis).resolve()),rank=U.shape[1],
        M=model.M,N=model.N,basis_arrays_bytes=U.nbytes+GU.nbytes,device_U_B_bytes=U.nbytes+B.nbytes,
        singular_values=s.tolist(),last_to_first=float(s[-1]/s[0]),allstrong_certified=False,
        basis_kind=args.kind,singular_values_are_G_spectrum=args.kind=='spectral',
        interpretation=('Finite returned spectral capture; remainder is uncaptured, not weak or null.' if args.kind=='spectral' else
            'First two public background source fields; listed singular values are seed-matrix spectrum, not Green spectrum. Unknown D_m multiplication not covered.'),
        source_metric='1/sum(abs(Einc)**2,axis=0), all60 sources',rows=rows,counts=dict(model.counts),
        physical_primal_RHS=0,physical_adjoint_RHS=0,GT_access=False,
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:receipt[k] for k in ('rank','last_to_first','CPU_seconds','wall_seconds')}))
    for row in rows:print(json.dumps(dict(case=row['case'],label=row['label'],fields=row['fields'])))

if __name__=='__main__':main()
