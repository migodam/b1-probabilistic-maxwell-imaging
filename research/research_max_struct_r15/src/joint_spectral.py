"""S-conditioned G-energy decomposition, fixed tau=.1; public operators only.

Bits in U11/U10/U01/U00 are S-sector then conditional G-energy. This is
not a simultaneous SVD or an exact intersection of global strong spaces.
No forward solve, material reconstruction, truth access or online optimizer.
"""
from __future__ import annotations
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_k]='1'
import argparse,json,time
from pathlib import Path
_MODULE_START_WALL=time.perf_counter()
import numpy as np
from scipy.linalg import svd,eigh

TAU=.1
NUMERICAL_RTOL=1e-12


def conditional_four_spaces(G,S,*,real_material=False,diagnostics=True):
    """Return bases/frames/spectra/metadata for the parent's fixed construction.

    G and S share their column space R^N or C^N. For real material, stack real
    and imaginary rows of each original operator before all SVD/Gram steps.
    Full right frames (N columns) retain exact numerical null-space columns.
    Tiny negative eigenvalues are clipped only in a diagnostic copy; raw
    eigenvalues determine the threshold split and remain in returned arrays.
    """
    cpu=time.process_time();wall=time.perf_counter();cost={}
    G=np.asarray(G);S=np.asarray(S)
    if G.ndim!=2 or S.ndim!=2 or G.shape[1]!=S.shape[1] or G.shape[1]==0:raise ValueError('Finite nonempty G/S matrices must share N columns')
    if not np.isfinite(G).all() or not np.isfinite(S).all():raise ValueError('Operators must be finite')
    original_shapes=dict(G=list(G.shape),S=list(S.shape));n=G.shape[1]
    if real_material:G=np.concatenate((G.real,G.imag),axis=0);S=np.concatenate((S.real,S.imag),axis=0)
    else:G=np.asarray(G,complex);S=np.asarray(S,complex)
    _,sigmaS,vhS=svd(S,full_matrices=True,check_finite=False);cost['full_S_SVD']=1
    _,sigmaG,vhG=svd(G,full_matrices=True,check_finite=False);cost['full_G_SVD']=1
    frameS=vhS.conj().T;frameG=vhG.conj().T
    smax=float(sigmaS[0]) if len(sigmaS) else 0.;gmax=float(sigmaG[0]) if len(sigmaG) else 0.
    # Explicit zero singular values in a full right frame satisfy the declared
    # >=0 strong threshold if the entire operator is zero; no hidden rank pad.
    full_sigmaS=np.zeros(n);full_sigmaS[:len(sigmaS)]=sigmaS
    full_sigmaG=np.zeros(n);full_sigmaG[:len(sigmaG)]=sigmaG
    Sstrong=full_sigmaS>=TAU*smax;Gstrong=full_sigmaG>=TAU*gmax
    Vs=frameS[:,Sstrong];Vw=frameS[:,~Sstrong]
    H=G.conj().T@G;H=(H+H.conj().T)/2;cost['G_Gram_products']=1
    threshold=TAU**2*gmax**2;spaces={};spectra={};clips={};sector_frames={};eighs=0;negative=[]
    for label,V,strong_name,weak_name in [('Sstrong',Vs,'U11','U10'),('Sweak',Vw,'U01','U00')]:
        if V.shape[1]:
            compressed=V.conj().T@H@V;compressed=(compressed+compressed.conj().T)/2
            values,Z=eigh(compressed,check_finite=False);eighs+=1
        else:values=np.zeros(0);Z=np.empty((0,0),dtype=G.dtype)
        # Ascending Hermitian eigenvalues are not singular values. Preserve raw.
        roundoff=128*np.finfo(float).eps*max(1,n)*max(gmax*gmax,np.finfo(float).tiny)
        tiny=(values<0)&(values>=-roundoff);large=values < -roundoff
        diag=values.copy();diag[tiny]=0
        selected=values>=threshold
        sector_frames[label]=V@Z;spaces[strong_name]=V@Z[:,selected];spaces[weak_name]=V@Z[:,~selected]
        spectra[label]=values;clips[label]=diag
        if np.any(large):negative.append(dict(sector=label,count=int(large.sum()),min_eigenvalue=float(values.min()),roundoff_tolerance=roundoff))
    cost['conditional_Hermitian_eigh']=eighs;cost['conditional_Gram_compressions']=eighs
    ortho=complete=None
    if diagnostics:
        allU=np.column_stack([spaces[k] for k in ('U11','U10','U01','U00')]);eye=np.eye(n)
        cost['validation_full_Gram_products']=1;cost['validation_completeness_products']=1
        ortho=float(np.linalg.norm(allU.conj().T@allU-eye,'fro'));complete=float(np.linalg.norm(allU@allU.conj().T-eye,'fro'))
    cross=Vs.conj().T@H@Vw;cost['cross_S_sector_G_Gram_product']=1
    # Small overlap Gram evaluates the same commutator norm; no repeated NxN
    # projector identities are needed inside each online inversion.
    overlap=Vs.conj().T@frameG[:,Gstrong]
    overlap_gram=overlap.conj().T@overlap
    comm2=2*float(np.trace(overlap_gram).real-np.vdot(overlap_gram,overlap_gram).real)
    cost['global_projector_small_overlap_Gram']=1
    nrS=int(np.count_nonzero(sigmaS>NUMERICAL_RTOL*smax)) if smax else 0
    nrG=int(np.count_nonzero(sigmaG>NUMERICAL_RTOL*gmax)) if gmax else 0
    metadata=dict(schema='r15_S_conditioned_G_energy_v1',scalar_domain='real_material' if real_material else 'current_complex',
        input_shapes=original_shapes,effective_shapes=dict(G=list(G.shape),S=list(S.shape)),N=n,tau=TAU,numerical_rank_rtol=NUMERICAL_RTOL,
        bits='S sector first, conditional G energy second: U11=Vs strong, U10=Vs weak, U01=Vw strong, U00=Vw weak',
        definition='Full S right frame then eig(V_sector* G*G V_sector); not simultaneous SVD or exact intersection of global S/G strong spaces',
        real_material_rule='stack(Re operator, Im operator) before all steps; original frames mean pre-conditional full right frames in the selected scalar domain',
        eigen_order='ascending; raw eigenvalues used for >= threshold selection',S_sigma_max=smax,G_sigma_max=gmax,
        S_strong_threshold=TAU*smax,conditional_G_energy_threshold=threshold,
        zero_operator_threshold_policy='Literal declared >= threshold including zero; numerical null ranks separately reported; no synthetic padding',
        ranks={**{k:int(v.shape[1]) for k,v in spaces.items()},'S_strong':int(Sstrong.sum()),'S_weak':int((~Sstrong).sum()),'global_G_strong':int(Gstrong.sum()),'S_numerical_rank':nrS,'S_numerical_null':n-nrS,'G_numerical_rank':nrG,'G_numerical_null':n-nrG},
        rank_sum=sum(v.shape[1] for v in spaces.values()),orthonormality_Frobenius=ortho,completeness_Frobenius=complete,
        cross_S_sector_G_Gram_Frobenius=float(np.linalg.norm(cross,'fro')),
        global_projector_commutator_Frobenius=float(np.sqrt(max(0.,comm2))),
        repeated_basis_identity_validation=bool(diagnostics),
        identity_validation_null_reason=None if diagnostics else 'Previously characterized public geometry; no repeated identity checks in online solve',
        global_projector_noncommutator_norm_definition='Frobenius norm of Ps Pg - Pg Ps for independent global tau-strong projectors',
        conditional_eigen_diagnostics={k:dict(min_eigenvalue=float(v.min()) if len(v) else None,tiny_negative_clipped_for_diagnostics=int(np.count_nonzero(v!=clips[k]))) for k,v in spectra.items()},
        significant_negative_eigenvalue_issues=negative,operation_counts=cost,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        cost_scope='Public operator SVD/Gram/eigh/diagnostics only, no physical LU/RHS; validation products charged',truth_access='NONE')
    return dict(**spaces,Vs=Vs,Vw=Vw,G_right_frame=frameG,S_right_frame=frameS,G_singular_values=sigmaG,S_singular_values=sigmaS,
                Sstrong_conditional_right_frame=sector_frames['Sstrong'],Sweak_conditional_right_frame=sector_frames['Sweak'],
                Sstrong_conditional_eigenvalues=spectra['Sstrong'],Sweak_conditional_eigenvalues=spectra['Sweak'],
                Sstrong_conditional_eigenvalues_diagnostic=clips['Sstrong'],Sweak_conditional_eigenvalues_diagnostic=clips['Sweak'],metadata=metadata)


def joint_spectral(G,S,*,real_material=False):
    """Named entry point; identical construction to conditional_four_spaces."""
    return conditional_four_spaces(G,S,real_material=real_material)


def write_json(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def characterize_public(paths,output,*,real_material=False):
    """Public G/S-only batch; exact array_equal cache, one receipt per input."""
    output=Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True);start=time.perf_counter();cpu=time.process_time();cache=[];executions=[]
    try:
        for index,path in enumerate(paths):
            c=time.process_time();w=time.perf_counter();path=Path(path).resolve()
            row=dict(index=index,public_path=str(path),public_bytes=path.stat().st_size if path.is_file() else None,real_material=real_material,truth_access='NONE',status='RUNNING')
            try:
                if any(x in ('sealed_truth','truth','data_private') for x in path.parts):raise ValueError('Private paths are forbidden')
                with np.load(path,allow_pickle=False) as z:G=z['G'].copy();S=z['S'].copy()
                row.update(keys_read=['G','S'],input_G_shape=list(G.shape),input_S_shape=list(S.shape))
                existing=next((item for item in cache if np.array_equal(item['G'],G) and np.array_equal(item['S'],S)),None)
                if existing is None:
                    result=conditional_four_spaces(G,S,real_material=real_material);group=f'operator_group_{len(cache):02d}';folder=output/group;folder.mkdir()
                    arrays={k:v for k,v in result.items() if isinstance(v,np.ndarray)}
                    np.savez_compressed(folder/'SPACES.npz',**arrays);write_json(folder/'METADATA.json',result['metadata'])
                    item=dict(G=G,S=S,group=group,result=result);cache.append(item);reused=False
                else:item=existing;reused=True
                row.update(status='COMPLETE',operator_group=item['group'],reused_exact_equal_G_S=reused,
                    reused_compute_not_charged_again=reused,characterization=item['result']['metadata'],new_characterization_CPU_seconds=0. if reused else item['result']['metadata']['CPU_seconds'])
            except BaseException as exc:row.update(status='FAILED',error=repr(exc));raise
            finally:
                row.update(CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w);executions.append(row);write_json(output/f'EXECUTION_{index:02d}.json',row)
        return executions
    finally:
        write_json(output/'RECEIPT.json',dict(status='COMPLETE' if len(executions)==len(paths) and all(x['status']=='COMPLETE' for x in executions) else 'FAILED',
            requested_public_paths=[str(Path(x).resolve()) for x in paths],completed=len(executions),operator_groups=len(cache),real_material=real_material,
            executions=[f'EXECUTION_{i:02d}.json' for i in range(len(executions))],CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,
            whole_process_CPU_seconds_at_receipt=time.process_time(),wall_since_before_numpy_scipy_imports=time.perf_counter()-_MODULE_START_WALL,
            truth_access='NONE',physics_calls=0,physical_LU=0,physical_complex_RHS=0,cache_policy='array_equal for complete G and S, no hashes',
            BLAS_threads={k:os.environ.get(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')}))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--public',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--real-material',action='store_true');a=p.parse_args()
    rows=characterize_public(a.public,a.output,real_material=a.real_material)
    print(json.dumps([dict(public=x['public_path'],operator_group=x['operator_group'],reused=x['reused_exact_equal_G_S'],ranks=x['characterization']['ranks'],CPU_seconds=x['CPU_seconds']) for x in rows]))

if __name__=='__main__':main()
