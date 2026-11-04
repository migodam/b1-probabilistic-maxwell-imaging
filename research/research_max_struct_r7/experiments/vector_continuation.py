"""Preserve initial12; fixed300-iteration common physical continuations.
One additional6^3 sphere with same .4lambda physical box and priors.
Full objective-evaluation iterates are sealed before offline GT scoring.
"""
from pathlib import Path
import os,sys,json,time,resource
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import numpy as np
from scipy import linalg as la
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r7.src.vector_inverse import VectorInverse
from research.research_max_struct_r1.experiments.qualify_backend import mie
METHODS=['ordinary_joint','dual_like_joint','A_material_joint','direct_physical']
OUT=ROOT/'RESULTS/vector_continuation'
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False))

def mie_actual_grids():
    theta=np.linspace(0,np.pi,181);dirs=np.column_stack([np.sin(theta),np.zeros_like(theta),np.cos(theta)])
    _,s2,sigma=mie(.145,2*np.pi,3,np.cos(theta));ref=1j/(2*np.pi)*s2[:,None]*np.column_stack([np.cos(theta),np.zeros_like(theta),-np.sin(theta)])
    rows=[]
    for ns in [5,6]:
        model=VectorInverse(ns,.4/ns);chi=2.*(np.sum(model.points**2,axis=1)<=.145**2)
        Y,J,E,lin=model.physical(chi);p=model.h**3*J.reshape(-1,3,6)
        phase=np.exp(-1j*model.backend.k*dirs@model.points.T);dip=np.einsum('rn,ncp->rcp',phase,p)
        F=model.backend.k**2/(4*np.pi)*(dip-dirs[:,:,None]*np.einsum('rc,rcp->rp',dirs,dip)[:,None,:]);f=F[:,:,0]
        row=dict(nside=ns,total_unknowns=model.N,occupied_cells=int(np.sum(chi>0)),h=model.h,occupied_volume=float(np.sum(chi>0)*model.h**3),analytic_volume=4*np.pi*.145**3/3,mie_far_vector_relative_error=float(la.norm(f-ref)/la.norm(ref)),linear_residual=lin,counts=model.counts,setup=model.setup)
        rows.append(row)
        np.savez_compressed(OUT/('mie_actual_grid_%d.npz'%ns),points=model.points,chi=chi,far=f,mie_far=ref,theta=theta)
    dump(OUT/'MIE_ACTUAL_GRIDS.json',dict(rows=rows,refinement_error_decreased=rows[1]['mie_far_vector_relative_error']<rows[0]['mie_far_vector_relative_error'],gate='REPORT_ACTUAL_ERRORS_NO_CONTINUUM_ACCEPTANCE'))

def continue_one(model,Y,initial,folder,old_receipt):
    model.reset_counts();tick=time.perf_counter();cpu=time.process_time();trials=[];phis=[];data=[];walls=[];cpus=[];lus=[];adj=[]
    def objective(x):
        phi,g=model.physical(x,Y);prior=model.prior(x)[0]
        trials.append(x.copy());phis.append(phi);data.append(max(phi-prior,0.));walls.append(time.perf_counter()-tick);cpus.append(time.process_time()-cpu);lus.append(model.counts['physical_LU']);adj.append(model.counts['physical_adjoint_RHS'])
        return phi,g
    res=minimize(objective,initial,jac=True,method='L-BFGS-B',bounds=[model.box]*model.N,options=dict(maxiter=300,ftol=1e-13,gtol=1e-9,maxls=30))
    # Final point is selected exclusively by public physicalPhi across all trials.
    best=int(np.argmin(phis));final=trials[best].copy();post=model.physical(final,Y,False)
    receipt=dict(status='COMPLETED',counts=model.counts,wall_s=time.perf_counter()-tick,cpu_s=time.process_time()-cpu,optimizer=dict(success=bool(res.success),message=str(res.message),nit=int(res.nit),nfev=int(res.nfev)),common_maxiter=300,physical_phi=float(phis[best]),physical_data_residual=float(np.sqrt(data[best])),final_linear_residual=post[-1],previous_counts=old_receipt['counts'],previous_wall_s=old_receipt['wall_s'],previous_cpu_s=old_receipt['cpu_s'],selection='lowest public physicalPhi trial; no truth available')
    folder.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(folder/'prediction_trace.npz',prediction=final,trials=np.array(trials),physical_phi=phis,data_residual_squared=data,wall_s=walls,cpu_s=cpus,physical_LU=lus,adjoint_RHS=adj,initial=initial)
    dump(folder/'RECEIPT.json',receipt)
    print(json.dumps(dict(folder=str(folder.relative_to(ROOT)),residual=receipt['physical_data_residual'],nit=res.nit,LU=model.counts['physical_LU'],wall=receipt['wall_s'])),flush=True)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    if(OUT/'COMPLETE.json').exists():raise FileExistsError('preserve completed continuation')
    tick=time.perf_counter();cpu=time.process_time()
    dump(OUT/'FROZEN_CONFIG.json',dict(authorization='parent requested common longer physical corrector and one necessary216-DOF sphere after initial12 results',common_maxiter=300,gtol=1e-9,ftol=1e-13,selection='public physicalPhi incumbent',grid216=dict(nside=6,h=.4/6,physical_box_width=.4,targets=['sphere']),preserve_initial12=True,truth_online=False))
    mie_actual_grids()
    # Additional one-grid run, original parameters retained.
    model=VectorInverse(6,.4/6);base=OUT/'grid216_sphere';base.mkdir(exist_ok=True)
    truth=2.*(np.sum(model.points**2,axis=1)<=.145**2);model.reset_counts();Y,_,_,lin=model.physical(truth)
    np.savez_compressed(base/'truth_for_scoring.npz',chi=truth,points=model.points);np.savez_compressed(base/'public_observations.npz',Y=Y,points=model.points)
    dump(base/'SIMULATION_COSTS.json',dict(counts=model.counts,setup=model.setup,linear_residual=lin))
    del truth
    for kind in METHODS:
        arrays,receipt=model.optimize(Y,kind);np.savez_compressed(base/(kind+'_prediction.npz'),**arrays);dump(base/(kind+'_receipt.json'),receipt)
    cases=[(ROOT/'RESULTS/vector_3d'/name,OUT/('grid125_'+name),5,.08) for name in ['sphere','two_spheres','non_spherical_bracket']]+[(base,OUT/'grid216_sphere_continued',6,.4/6)]
    for source,dest,ns,h in cases:
        model=VectorInverse(ns,h);Y=np.load(source/'public_observations.npz')['Y']
        for kind in METHODS:
            arr=np.load(source/(kind+'_prediction.npz'));old=json.loads((source/(kind+'_receipt.json')).read_text());continue_one(model,Y,arr['prediction'],dest/kind,old)
    # All16 continuation prediction traces are now sealed; score only now.
    rows=[]
    for source,dest,ns,h in cases:
        truth=np.load(source/'truth_for_scoring.npz')['chi']
        for kind in METHODS:
            folder=dest/kind;a=np.load(folder/'prediction_trace.npz');r=json.loads((folder/'RECEIPT.json').read_text());old=json.loads((source/(kind+'_receipt.json')).read_text())
            errors=la.norm(a['trials']-truth,axis=1)/la.norm(truth);incumbent=np.minimum.accumulate(a['physical_phi']);idx=np.zeros(len(incumbent),int);best=0
            for i,p in enumerate(a['physical_phi']):
                if p<a['physical_phi'][best]:best=i
                idx[i]=best
            np.savez_compressed(folder/'offline_scored_curve.npz',total_physical_LU=a['physical_LU']+old['counts']['physical_LU'],total_cpu_s=a['cpu_s']+old['cpu_s'],total_wall_s=a['wall_s']+old['wall_s'],incumbent_NRMSE=errors[idx],incumbent_physical_phi=incumbent,incumbent_data_residual=np.sqrt(a['data_residual_squared'][idx]))
            row=dict(case=dest.name,method=kind,initial_NRMSE=float(la.norm(a['initial']-truth)/la.norm(truth)),final_NRMSE=float(la.norm(a['prediction']-truth)/la.norm(truth)),physical_data_residual=r['physical_data_residual'],continuation_LU=r['counts']['physical_LU'],total_physical_LU=r['counts']['physical_LU']+old['counts']['physical_LU'],total_auxiliary_Cholesky=old['counts']['current_normal_factorizations'],continuation_cpu_s=r['cpu_s'],total_cpu_s=r['cpu_s']+old['cpu_s'],nit=r['optimizer']['nit'],status=r['optimizer']['message'])
            rows.append(row)
        del truth
    dump(OUT/'SCORES.json',rows)
    dump(OUT/'COMPLETE.json',dict(status='COMPLETED',wall_s=time.perf_counter()-tick,cpu_s=time.process_time()-cpu,process_peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),continuation_predictions=16,additional_original_grid216_predictions=4,failed_calls=[],cost_scope='includes2 coarse/refined Mie actual-grid forward solves,216base4 reconstructions,16 continuation runs,setup,offline scoring; plots excluded'))
if __name__=='__main__':
    with threadpool_limits(limits=1):main()
