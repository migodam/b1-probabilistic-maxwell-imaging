"""No-GT shifted auxiliary-denominator diagnostics on24 saved joint states."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,csv,time,traceback
from collections import Counter
import numpy as np
from scipy.linalg import lu_factor,lu_solve,svdvals
from threadpoolctl import threadpool_limits
O=Path(__file__).resolve().parent;STUDY=O.parents[1];REPO=STUDY.parents[1];R9=REPO.parent/'r9_worktree'
sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.src.lowrank_denominator import LowRankDenominator
from research.research_max_struct_r4.src.cie_onehop import OneHopAbsorption
CASES=('obs_b20d5a55677d4eca','obs_40ac5b79e8a74ac3')
def js(p):return json.loads(p.read_text())
def load(p):
    with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def stamp(p):
    s=p.stat();return dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
def dump(n,v):(O/n).write_text(json.dumps(v,indent=2,default=lambda v:v.item() if isinstance(v,np.generic) else str(v))+'\n')
def write(n,rows):
    if not rows:return
    with (O/n).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
def norm(x):return float(np.linalg.norm(x))
def energy(x):return float(np.vdot(x,x).real)
def ratio(x,y):return float(x/y) if y else None
def main():
    start=time.perf_counter();cpu=time.process_time();rows=[];missing=[];failures=[];sources=[];cost=Counter();localcost=Counter();smallcost=Counter();arrays={}
    for case in CASES:
        pp=R9/'research/research_max_struct_r9/data_public/development'/f'{case}.npz'
        op=R9/'research/research_max_struct_r9/data_public/development/operators.npz';pub=load(op);pub.update(load(pp));G=pub['G'];n=len(G)
        sources.extend([stamp(pp),stamp(op)])
        green_dir=STUDY/'RESULTS/denominator_dev'/f'{case}__green_strong__galerkin'
        bp=green_dir/'DENOMINATOR_BASIS.npz';basis=load(bp);U,B=basis['U'],basis['B'];sources.append(stamp(bp))
        lowrank=LowRankDenominator(G,U,'galerkin',smallcost)
        # The actual sealed B is retained, after recording its discrepancy from reconstructed public PGP.
        B_source_difference=norm(lowrank.B-B);lowrank.B=B
        local=OneHopAbsorption(G,(24,24),localcost)
        dirs={'R9_global_dual':R9/'research/research_max_struct_r9/RESULTS/development'/f'{case}__global_dual__fixed_001',
              'R9_a_plus':R9/'research/research_max_struct_r9/RESULTS/development'/f'{case}__a_plus__fixed_001',
              'R10_green_strong':green_dir}
        for method,d in dirs.items():
            for st in range(4):
                p=d/f'stage_{st:02d}_joint.npz';mp=p.with_suffix('.json')
                if not p.exists() or not mp.exists():missing.append(dict(case_id=case,method=method,stage=st,path=str(p),status='MISSING_SAVED_STATE'));continue
                state_start=time.perf_counter();state_cpu=time.process_time()
                try:
                    z=load(p);meta=js(mp);sources.extend([stamp(p),stamp(mp)])
                    m,J,Q=z['m'],z['J'],z['Q'];beta=float(meta['beta_material']);e=J-m[:,None]*Q;ns=J.shape[1]
                    a0=1+beta*m;z0=m/a0;T0=U.conj().T@(z0[:,None]*B)
                    ev=np.linalg.eigvals(T0);rho0=float(np.max(abs(ev)));smin=float(svdvals(np.eye(len(T0))-T0)[-1])
                    za=m/(1+(beta-np.diag(G))*m)
                    local_matrices=za[local.indices][:,:,None]*local.goff
                    local_ev=np.linalg.eigvals(local_matrices);rhos=np.max(abs(local_ev),axis=1);maxlocal=float(np.max(rhos))
                    A=np.eye(n,dtype=complex)+beta*np.diag(m)-m[:,None]*G
                    t=time.perf_counter();f=lu_factor(A,check_finite=False);Rfull=beta*lu_solve(f,e,check_finite=False)
                    shifted_wall=time.perf_counter()-t;cost['shifted_auxiliary_576_LU']+=1;cost['shifted_auxiliary_complex_RHS']+=ns
                    RA=local.apply(m,beta,e);RK=beta*lowrank.solve(lowrank.factor(m,beta),e)
                    nf,na,nk=map(norm,(Rfull,RA,RK));ge=G@e;ke=B@(U.conj().T@e)
                    row=dict(case_id=case,method=method,stage=st,mf=meta['basis_info']['mf'],beta_material=beta,beta_current=meta['beta_current'],rank=U.shape[1],
                        truth_access='NONE',status='COMPLETED',projected_global_T0_spectral_radius=rho0,global_projected_Neumann_rho_gt1=bool(rho0>1),
                        min_singular_I_minus_T0=smin,local4x4_max_spectral_radius=maxlocal,local4x4_median_spectral_radius=float(np.median(rhos)),
                        local4x4_blocks_rho_gt1=int(np.sum(rhos>1)),local4x4_argmax_block=int(np.argmax(rhos)),
                        e_norm=norm(e),R_full_norm=nf,R_onehop_norm=na,R_lowrank19_norm=nk,
                        R_full_squared_energy=energy(Rfull),R_onehop_squared_energy=energy(RA),R_lowrank19_squared_energy=energy(RK),
                        onehop_relative_error_to_full=ratio(norm(RA-Rfull),nf),lowrank19_relative_error_to_full=ratio(norm(RK-Rfull),nf),
                        onehop_norm_ratio_to_full=ratio(na,nf),lowrank19_norm_ratio_to_full=ratio(nk,nf),
                        K19_error_on_e_norm=norm(ge-ke),K19_error_on_e_relative_to_Ge=ratio(norm(ge-ke),norm(ge)),
                        shifted_solve_relative_residual=ratio(norm(A@Rfull-beta*e),norm(beta*e)),
                        shifted_auxiliary_576_LU=1,shifted_auxiliary_complex_RHS=ns,physical_LU=0,physical_complex_RHS=0,
                        shifted_factor_solve_wall_seconds=shifted_wall,state_wall_seconds=time.perf_counter()-state_start,state_CPU_seconds=time.process_time()-state_cpu,
                        public_reconstructed_B_minus_saved_B_norm=B_source_difference,source=str(p))
                    rows.append(row);key=f'{case}__{method}__stage{st}'
                    arrays[key+'__global_T0_eigenvalues']=ev;arrays[key+'__local_block_spectral_radii']=rhos
                    print(json.dumps({k:row[k] for k in ('case_id','method','stage','projected_global_T0_spectral_radius','local4x4_max_spectral_radius','onehop_relative_error_to_full','lowrank19_relative_error_to_full')}),flush=True)
                except Exception as exc:failures.append(dict(case_id=case,method=method,stage=st,error=repr(exc),traceback=traceback.format_exc()))
    write('SPECTRUM_AND_RESIDUALS.csv',rows);np.savez_compressed(O/'SPECTRAL_ARRAYS.npz',**arrays)
    dump('MISSING_STATES.json',missing);dump('FAILED_ATTEMPTS.json',failures);dump('SOURCE_PROVENANCE.json',list({r['path']:r for r in sources}.values()))
    dump('EXACT_DEFINITIONS.json',dict(beta='Saved beta_material, not beta_current',global_operator='K=PGP19 from fixed public G right strong basis at tau.1; T0=U*D_(m/(1+beta*m))B; nonzero eigenvalues of D_z0 K equal those of T0',global_limit='Reported rho is projected-global K19, not full dense rho(D_z0 G). min singular(I-T0) is a small-system diagnostic, not min singular(A_G).',local='For each saved2x2spatialblock L_b=D_(m/(1+(beta-Gii)m)) Goff,b; report max of144block spectral radii',full='R_full=beta (I+beta D_m-D_m G)^-1 e, e=J-D_m Q with saved auxiliary J/Q',onehop='R_A=(I+D_z Goff,block) beta D_(1+(beta-Gii)m)^-1 e; exact inherited OneHopAbsorption.apply',lowrank='R_K=beta (I+beta D_m-D_m K)^-1 e; sealed G19 Galerkin basis',error='relative_error_to_full=||R_method-R_full||F/||R_full||F; norm is not squared energy; energy fields are squared',scope='Auxiliary metric comparison on saved states. No GT, no new optimization, no physical forward or recovery certificate.',interpretation='rho(T0)>1 obstructs operator Neumann convergence for projected-global feedback; it does not imply singular(I-T0), failure of shifted solve, or explain every residual difference. Local rho<1 does not certify negligible omitted global coupling.'))
    dump('RECEIPT.json',dict(status='COMPLETE' if not failures and not missing else 'COMPLETE_WITH_GAPS',expected_states=24,completed_states=len(rows),missing_states=len(missing),failed_states=len(failures),truth_access='NONE',new_inverse_runs=0,physical_LU=0,physical_complex_RHS=0,shifted_auxiliary_cost=dict(cost),lowrank_cost=dict(smallcost),onehop_cost=dict(localcost),wall_seconds=time.perf_counter()-start,CPU_seconds=time.process_time()-cpu,GPU=False,BLAS_threads=1,scope='Shifted denominator solves are576-dimensional auxiliary LU, not original Maxwell physical LU. No candidate changes or historical figure modifications.'))
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
