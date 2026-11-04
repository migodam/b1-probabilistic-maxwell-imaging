"""Frozen small 3D development benchmark. Reproduce with one BLAS thread.
Independent analytic Mie qualification is inherited unmodified; new inversions
are same-discretization noiseless development tests, not continuum recovery.
"""
from pathlib import Path
import os,sys,json,time,resource,traceback
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r7.src.vector_inverse import VectorInverse
METHODS=['ordinary_joint','dual_like_joint','A_material_joint','direct_physical']
OUT=ROOT/'RESULTS/vector_3d'
def dump(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False))
def qualify(model):
    inherited=REPO/'research/research_max_struct_r1/results/backend_qualification/QUALIFICATION.json'
    reference=json.loads(inherited.read_text())
    q=dict(status='COMPLETED',inherited_independent_analytic_Mie=dict(path=str(inherited.relative_to(REPO)),receipt=reference,reused_without_rerun=True,new_cost_LU=0),derivatives={})
    x=np.linspace(.1,1.4,model.N);pred,J,E,lin=model.physical(x);Y=pred*(1.03+.01j)
    # Independent solve in actual CM+RR dipole p coordinates, not chi-density.
    alpha=model.h**3*x/(1-model.s*x);ad=np.repeat(alpha,3)
    B=np.eye(model.M)-ad[:,None]*model.backend.Goff
    p=la.solve(B,ad[:,None]*model.inc,check_finite=False)
    q['CM_RR_equivalence_relative_error']=float(la.norm(p-model.h**3*J)/la.norm(p))
    q['material_variable']='real chi=epsilon_r-1; alpha=Vchi/(1-schi) nonlinear; derivative computed in chi including self'
    rng=np.random.default_rng(715);v=rng.normal(size=model.N);v/=la.norm(v);step=2e-6
    funcs={'physical':lambda z:model.physical(z,Y)}
    for kind in METHODS[:-1]:funcs[kind]=lambda z,k=kind:model.material(z,J,E,k)
    for name,func in funcs.items():
        f,g=func(x);fp=func(x+step*v)[0];fm=func(x-step*v)[0];fd=(fp-fm)/(2*step);analytic=float(g@v)
        q['derivatives'][name]=dict(analytic=analytic,centered_difference=fd,relative_error=float(abs(analytic-fd)/max(abs(fd),abs(analytic),1e-12)),step=step)
    q['new_qualification_counts']=dict(model.counts)
    q['finite_solver_linear_relative_residual']=lin
    q['continuum_certificate']='NOT_CLAIMED; inherited sphere Mie errors are 10.47%/3.66% at h=.04/.03; imaging cube h=.08 is coarser, same-grid development only'
    dump(OUT/'QUALIFICATION.json',q)
    if any(row['relative_error']>2e-5 for row in q['derivatives'].values()):raise RuntimeError('derivative implementation mismatch')
    return q

def targets(model):
    p=model.points
    sphere=2.0*(np.sum(p*p,axis=1)<=.145**2)
    two=1.5*(np.sum((p-np.array([-.105,0,0]))**2,axis=1)<=.105**2)+2.5*(np.sum((p-np.array([.105,.04,.04]))**2,axis=1)<=.105**2)
    # bent bracket: union of three axis-aligned arms, no radial symmetry.
    ij=model.ijk
    shape=((ij[:,0]==1)&(ij[:,1]>=1)&(ij[:,1]<=3)&(ij[:,2]>=1)&(ij[:,2]<=3))|((ij[:,1]==1)&(ij[:,0]>=1)&(ij[:,0]<=3)&(ij[:,2]>=1)&(ij[:,2]<=3))|((ij[:,2]==1)&(ij[:,0]>=1)&(ij[:,0]<=3)&(ij[:,1]>=1)&(ij[:,1]<=3))
    return dict(sphere=sphere,two_spheres=two,non_spherical_bracket=1.8*shape)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    if(OUT/'COMPLETE.json').exists():raise FileExistsError('preserve completed benchmark')
    start=time.perf_counter();cpu=time.process_time();model=VectorInverse()
    config=dict(status='FROZEN_BEFORE_RUN',grid=[5,5,5],N=125,vector_dimension=375,h_lambda=.08,box_chi=[0,4],initial_chi=.25,lambda_smooth=model.lam,lambda_l2=model.l2,rho=model.rho,beta=model.beta,outer_cycles=15,material_inner_maxiter=12,direct_maxiter=180,common_physical_corrector_maxiter=40,source_directions=model.directions.tolist(),source_polarizations=model.polarizations.tolist(),receivers=model.receivers.tolist(),noise='none',target_status='EXPOSED_DEVELOPMENT_SHAPES_NOT_BLIND_HOLDOUT',truth_policy='truth for simulation then scorer only, never passed optimizer/selection',methods=METHODS,threads=1,setup=model.setup)
    dump(OUT/'FROZEN_CONFIG.json',config);qualify(model)
    truths=targets(model);simulation=[]
    for name,truth in truths.items():
        model.reset_counts();Y,_,_,lin=model.physical(truth)
        target=OUT/name;target.mkdir(exist_ok=True)
        np.savez_compressed(target/'truth_for_scoring.npz',chi=truth,points=model.points)
        np.savez_compressed(target/'public_observations.npz',Y=Y,points=model.points,receivers=model.receivers,directions=model.directions,polarizations=model.polarizations)
        simulation.append(dict(target=name,counts=dict(model.counts),linear_relative_residual=lin))
    del truths,truth
    dump(OUT/'SIMULATION_COSTS.json',simulation)
    summaries=[];failures=[]
    for name in ['sphere','two_spheres','non_spherical_bracket']:
        target=OUT/name;Y=np.load(target/'public_observations.npz')['Y']
        for kind in METHODS:
            try:
                arrays,receipt=model.optimize(Y,kind)
                # Seal actual material arrays before scorer has any truth access.
                np.savez_compressed(target/(kind+'_prediction.npz'),**arrays)
                dump(target/(kind+'_receipt.json'),receipt)
                print(json.dumps(dict(target=name,method=kind,wall_s=receipt['wall_s'],residual=receipt['physical_data_residual'],physical_LU=receipt['counts']['physical_LU'])),flush=True)
            except Exception as exc:
                fail=dict(target=name,method=kind,error=repr(exc),traceback=traceback.format_exc(),counts=dict(model.counts));failures.append(fail);dump(target/(kind+'_FAILED.json'),fail)
        # Scoring happens after all four methods have sealed predictions.
        truth=np.load(target/'truth_for_scoring.npz')['chi']
        for kind in METHODS:
            file=target/(kind+'_prediction.npz')
            if not file.exists():continue
            pred=np.load(file);receipt=json.loads((target/(kind+'_receipt.json')).read_text());x=pred['prediction'];raw=pred['raw']
            row=dict(target=name,method=kind,material_NRMSE=float(la.norm(x-truth)/la.norm(truth)),raw_NRMSE=float(la.norm(raw-truth)/la.norm(truth)),material_MAE=float(np.mean(abs(x-truth))),support_IoU=float(np.sum((x>=.75)&(truth>=.75))/max(np.sum((x>=.75)|(truth>=.75)),1)),truth_norm=float(la.norm(truth)),**{k:receipt[k] for k in ['physical_data_residual','physical_phi','pre_corrector_data_residual','wall_s','cpu_s','counts','final_linear_residual']})
            summaries.append(row)
        del truth
    dump(OUT/'SCORES.json',summaries)
    complete=dict(status='COMPLETED' if not failures else 'COMPLETED_WITH_FAILURES',wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,peak_process_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),peak_memory_scope='macOS ru_maxrss bytes; whole sequential benchmark process including forward qualification; per-method peak NOT_MEASURED',setup=model.setup,failures=failures,rows=len(summaries),null_costs=dict(GPU=0,GT_online=0,independent_Mie_rerun=0,mature_vector_CIE_implementation=0,model_specific_parameter_search=0),limitations=['same-grid synthetic/noiseless exposed development set','no continuum convergence claim for imaging','joint controls are vector-adapted alternating baselines, not mature 2D CIE','Dual-like material metric is a local resolvent residual, not an exact historical Dual reproduction','fixed common iteration caps do not mean equal LU, wall time, or optimization convergence','3D isotropic per-voxel unknowns; no tensor material inversion'])
    dump(OUT/'COMPLETE.json',complete)
    print(json.dumps(complete),flush=True)
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
