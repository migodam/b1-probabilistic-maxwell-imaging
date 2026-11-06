"""Public-operator-only numerical rank of actual projected Fourier spaces.

No material target/GT/forward solve. For orthonormal U_F,V,
rank((I-VV*)U_F)=K-nullity((I-U_F U_F*)V). The latter has only r columns.
This is routine finite subspace bookkeeping, not a new spectral theorem.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import time
T=time.perf_counter();C=time.process_time()
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
import numpy as np
from scipy import fft,linalg as la
from scipy.stats import chi2
from research.research_max_struct_r10.src.vector_grid import VectorGrid


def projected_rank(V,n,mf):
    M,r=V.shape;freq=np.rint(np.fft.fftfreq(n)*n).astype(int)
    mask=np.isin(freq,np.r_[np.arange(mf),np.arange(-mf,0)])
    spatial=mask[:,None,None]&mask[None,:,None]&mask[None,None,:]
    K=3*int(spatial.sum())
    if r:
        hat=fft.fftn(V.reshape((n,)*3+(3,r)),axes=(0,1,2),norm='ortho')
        hat[~spatial]=0
        projected=fft.ifftn(hat,axes=(0,1,2),norm='ortho').reshape(M,r)
        residual=V-projected;values=la.svdvals(residual,check_finite=False)
        cutoff=np.finfo(float).eps*max(M,K,r) # Cancellation is relative to norm(V)=1.
        intersection=r-int(np.sum(values>cutoff))
    else:values=np.empty(0);cutoff=np.finfo(float).eps*max(M,K);intersection=0
    rank=K-intersection
    return dict(mf=mf,nominal_Fourier_columns_per_illumination=K,
        numerical_projected_basis_rank_per_illumination=rank,
        complex_coefficients_stored_all60=K*60,effective_complex_current_DOF_all60=rank*60,
        effective_real_current_DOF_all60=2*rank*60,
        intersection_dimension=intersection,projection_residual_singular_values=values.tolist(),cutoff=cutoff,
        minimum_principal_residual=float(values.min()) if len(values) else None,
        scalar_material_unknowns_shared_all60=1728,
        method='SVD of explicit (I-UF UF*)V via actual native3D FFT mask, then finite rank identity; no G action')


def main():
    output=ROOT/'RESULTS/vector/representation_rank';output.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((ROOT/'configs/vector_actual_contract.json').read_text());model=VectorGrid(12,.6,cfg['acquisition'])
    U,s,Vh=la.svd(model.S,full_matrices=False,check_finite=False)
    supported=int(np.sum(s>np.finfo(float).eps*max(model.S.shape)*s[0]))
    backgrounds={str(mf):projected_rank(Vh[:supported].conj().T,12,mf) for mf in (2,3,4)}
    cases=[]
    for public_id in ('v02_matched12','v09_matched12','v02_cross18to12','v09_cross18to12','v04_matched12','v08_matched12'):
        folder=ROOT/'data_public/vector'/public_id
        if not (folder/'observations.npz').exists():continue
        Y=np.load(folder/'observations.npz',allow_pickle=False)['Y'];meta=json.loads((folder/'PUBLIC.json').read_text())
        coord=U.conj().T@Y;outside=la.norm(Y-U@coord)**2
        tail=outside+np.r_[np.cumsum(np.sum(abs(coord)**2,axis=1)[::-1])[::-1],0.]
        radius=meta['sigma']**2*chi2.ppf(.99,2*Y.size)/2
        eligible=np.flatnonzero(tail[:supported+1]<=radius)
        if not len(eligible):cases.append(dict(public_id=public_id,status='NO_SUPPORTED_PREFIX',forward_RHS=0));continue
        rank=int(eligible[0]);row=dict(public_id=public_id,receiver_retained_rank=rank,
            vector_components_per_illumination=5184,full_complex_current_array_entries_all60=311040,
            material_shared_real_DOF=1728,stages=[projected_rank(Vh[:rank].conj().T,12,mf) for mf in (2,3,4)],
            scope='Public S/Y reconstruction of same native projected space; no material GT. Source-host SVD phase does not define a scalar-material projection.')
        (output/(public_id+'.json')).write_text(json.dumps(row,indent=2)+'\n');cases.append(row)
    receipt=dict(status='PUBLIC_OPERATOR_ONLY_NUMERICAL_RANK_BOOKKEEPING',CPU_seconds=time.process_time()-C,
        wall_seconds=time.perf_counter()-T,platform='Mac',cases=cases,full_supported_receiver_space_rank=supported,
        full_receiver_space_background_stages=backgrounds,
        no_material_GT=True,no_physical_forward_RHS=True,
        full_current_storage_is_not_independent_optimization_DOF=True,
        numerical_rank_not_a_continuum_or_identifiability_certificate=True)
    (output/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(CPU_seconds=receipt['CPU_seconds'],cases=[dict(public_id=x['public_id'],receiver_rank=x.get('receiver_retained_rank'),projected_ranks=[y['numerical_projected_basis_rank_per_illumination'] for y in x.get('stages',[])]) for x in cases])))


if __name__=='__main__':main()
