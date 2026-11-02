"""Bounded wiring/stationarity/current-space diagnostics; no GT or inversions."""
from pathlib import Path
from collections import Counter
import json,sys,time,traceback
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_fbg_r3_experiments.src.cie_tsom_baseline import fourier_apcs
from research.research_max_struct_r2.src.cie_global import run_cie,transport_source_metric
OUT=ROOT/'results/diagnostics'

def write(path,obj):
    path.write_text(json.dumps(obj,indent=2,allow_nan=False))

def relative(a,b):return float(la.norm(a-b)/max(la.norm(a),1e-300))

def source_wiring():
    start=time.perf_counter();cpu=time.process_time();models=[]
    rec=dict(status='RUNNING',classification='ONE_NEW_FULL_METRIC_WIRING_DIAGNOSTIC',
             source='src/cie_global.py',iterations=16,prior='same historical real-material prior',
             truth_access='disclosed constant4x4 root only; no blind truth')
    try:
        with np.load(REPO/'research/research_max_struct_r1/data_public/dynamics_operators_4x4.npz',allow_pickle=False)as f:
            ops={k:f[k]for k in ('G','E0','S','grid_shape')}
        root=MaxwellModel(ops);models.append(root)
        Y=root.state(np.full(root.n,.55))['Y'].copy();E=ops['E0'];ns=E.shape[1]
        rng=np.random.default_rng(829);U,_=la.qr(rng.normal(size=(ns,ns))+1j*rng.normal(size=(ns,ns)))
        Od=np.diag(1/np.sum(abs(Y)**2,axis=0));Os=np.diag(1/np.sum(abs(E)**2,axis=0))
        trajectories=[];runs=[]
        for mixed in (False,True):
            public=dict(ops);public['E0']=E@U if mixed else E
            target=Y@U if mixed else Y;model=MaxwellModel(public);models.append(model);trajectory=[]
            cfg=dict(stages=[dict(mf=1,beta_current=1.,beta_material=.2,maxiter=16)],
                dpcs_rtol=1e-3,sourcewise_cg=False,coefficient_tol=0.,miniter=2,
                regularizer=model.regularizer,material_inner_maxiter=60,
                physical_sigma=float(la.norm(Y)),
                data_metric=transport_source_metric(Od,U)if mixed else Od,
                state_metric=transport_source_metric(Os,U)if mixed else Os,
                iteration_callback=lambda st:trajectory.append(st))
            result=run_cie(model,target,np.full(model.n,.1),cfg)
            trajectories.append(trajectory);runs.append(result['execution_receipt'])
        md=[];jd=[];ad=[]
        for a,b in zip(*trajectories):
            md.append(float(la.norm(a['m']-b['m'])))
            jd.append(float(la.norm(a['J']-b['J']@U.conj().T)))
            ad.append(float(la.norm(a['coefficients']-b['coefficients']@U.conj().T)))
        np.savez_compressed(OUT/'SOURCE_WIRING_TRAJECTORIES.npz',unitary=U,
            m_original=np.array([s['m']for s in trajectories[0]]),m_mixed=np.array([s['m']for s in trajectories[1]]),
            J_original=np.array([s['J']for s in trajectories[0]]),J_mixed_back=np.array([s['J']@U.conj().T for s in trajectories[1]]))
        rec.update(status='COMPLETED',material_trajectory_difference=md,current_trajectory_difference=jd,
            coefficient_trajectory_difference=ad,max_material_difference=max(md),max_current_difference=max(jd),
            unitary_defect=float(la.norm(U.conj().T@U-np.eye(ns))),runs=runs,
            interpretation='wiring/covariance diagnostic only; no inverse superiority or newly acquired data')
    except Exception as exc:rec.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc());raise
    finally:
        costs=Counter()
        for m in models:costs.update(m.cost)
        rec.update(wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,physical_model_cost=dict(costs),
            timer_scope='public load,disclosed root,2 fixed16step diagnostics,trajectory serialization; receipt write excluded')
        write(OUT/'SOURCE_WIRING.json',rec)

def projection_setup(public):
    start=time.perf_counter();cpu=time.process_time();cost=Counter();S=public['S']
    Us,sigma,Vh=la.svd(S,full_matrices=False,check_finite=True);cost['receiver_SVD']=1
    rank=int(np.count_nonzero(sigma>1e-3*sigma[0]));Vd=Vh[:rank].conj().T
    bases={}
    for mf in (10,12):
        U,info=fourier_apcs(tuple(public['grid_shape']),mf,Vd,cost)
        Q,R,pivot=la.qr(U,mode='economic',pivoting=True,check_finite=True);cost['APCS_QR']+=1
        diag=abs(np.diag(R));tol=max(U.shape)*np.finfo(float).eps*max(float(diag.max()),1.)
        r=int(np.count_nonzero(diag>tol));Q=Q[:,:r]
        # DPCS and APCS are orthogonal by construction; do not refit SVD modes.
        bases[mf]=(Q,dict(info,retained_QR_rank=r,QR_rank_threshold=float(tol),
                         total_linear_rank=rank+r))
    return Vd,Us[:,:rank],sigma[:rank],bases,dict(cost,wall_seconds=time.perf_counter()-start,
        cpu_seconds=time.process_time()-cpu,DPCS_rank=rank,setup_shared_once=True)

def endpoints():
    start=time.perf_counter();cpu=time.process_time();models=[];rows=[]
    summary=dict(status='RUNNING',classification='POST_PREDICTION_NO_GT_ENDPOINT_DIAGNOSTIC',
        inverse_reruns=0,truth_access='NONE',projection_scope='ordinary scalar-current L2 only, no Hcurl claim',
        stationarity_scope='projected gradient at stored corrected prediction, not a basin certificate',
        auxiliary_scope='stored stage2 J transported to corrected predictionm; not original stage2 defect')
    try:
        manifest=json.loads((ROOT/'data_public/manifest.json').read_text());setup=None
        for record in manifest['records']:
            with np.load(ROOT/'data_public'/record['path'],allow_pickle=False)as f:
                public={k:f[k]for k in ('G','E0','S','Y','grid_shape','sigma','noise_radius')}
            if setup is None:
                Vd,Ud,sig,bases,setup=projection_setup(public)
            Y=public['Y'];JD=Vd@((Ud.conj().T@Y)/sig[:,None]);S=public['S'];G=public['G'];E=public['E0']
            for arm in ('native','global'):
                tick=time.perf_counter();ct=time.process_time();model=MaxwellModel(public);models.append(model)
                with np.load(ROOT/'predictions'/record['observation_id']/(arm+'.npz'),allow_pickle=False)as f:m=f['m']
                with np.load(ROOT/'results/blind'/record['observation_id']/arm/'stage_2.npz',allow_pickle=False)as f:
                    aux=f['J'];stage_m=f['m']
                # One retained physical LU and one original analytic source adjoint.
                state=model.state(m);value,g=model.data_value_gradient(m,Y)
                J=state['j'];pg=m-np.clip(m-g,.02,20.)
                raux=aux-m[:,None]*(E+G@aux);delta=J-aux
                dport=S@delta;physical=state['Y']-Y;aport=S@aux-Y
                row=dict(observation_id=record['observation_id'],arm=arm,status='COMPLETED',
                    physical_objective=float(value),gradient_l2=float(la.norm(g)),gradient_inf=float(la.norm(g,np.inf)),
                    projected_gradient_mapping_l2=float(la.norm(pg)),projected_gradient_mapping_inf=float(la.norm(pg,np.inf)),
                    projected_gradient_step=1.,active_lower_exact=int(np.count_nonzero(m==.02)),active_upper_exact=int(np.count_nonzero(m==20.)),
                    box_violation=max(float(np.max(.02-m)),float(np.max(m-20.)),0.),
                    final_material_minus_original_stage2_material_l2=float(la.norm(m-stage_m)),
                    physical_state_relative_defect=state['state_relative_residual'],physical_residual=float(la.norm(physical)),
                    physical_relative_residual=float(la.norm(physical)/la.norm(Y)),
                    transported_auxiliary_state_defect_l2=float(la.norm(raux)),transported_auxiliary_port_residual=float(la.norm(aport)),
                    physical_current_minus_auxiliary_l2=float(la.norm(delta)),
                    observed_current_residual_amplification=float(la.norm(delta)/max(la.norm(raux),1e-300)),
                    observed_port_residual_amplification=float(la.norm(dport)/max(la.norm(raux),1e-300)),
                    amplification_scope='observed one residual direction, not resolvent operator norm/bound',
                    identity_defect=float(la.norm(delta-m[:,None]*(G@delta)+raux)),projections={})
                for mf,(Q,info)in bases.items():
                    jp=JD+Q@(Q.conj().T@(J-JD))
                    full=Vd@(Vd.conj().T@J)+Q@(Q.conj().T@J)
                    state_defect=jp-m[:,None]*(E+G@jp)
                    row['projections'][str(mf)]=dict(basis=info,
                        fixed_DPCS_affine_relative_L2_error=relative(J,jp),
                        free_DPCS_linear_relative_L2_error=relative(J,full),
                        projected_current_state_defect=float(la.norm(state_defect)),
                        projected_current_port_relative_residual=float(la.norm(S@jp-Y)/la.norm(Y)),
                        fixed_DPCS_coefficient_mismatch_relative=float(la.norm(Vd.conj().T@(J-JD))/max(la.norm(J),1e-300)),
                        interpretation='best L2 projection at final physical material; not a PDE-constrained optimizer or GT projection')
                np.savez_compressed(OUT/(record['observation_id']+'_'+arm+'_physical_endpoint.npz'),m=m,j=J,Q=state['Q'],gradient=g,projected_gradient=pg)
                row.update(wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-ct,physical_model_cost=dict(model.cost))
                rows.append(row);write(OUT/(record['observation_id']+'_'+arm+'.json'),row)
        summary.update(status='COMPLETED',rows=rows,shared_projection_setup=setup,
                       conclusion_policy='inspect projected gradients and optimizer termination; unresolved budget/stationarity is not a true-basin result')
    except Exception as exc:summary.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc(),rows=rows);raise
    finally:
        costs=Counter()
        for m in models:costs.update(m.cost)
        summary.update(wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,physical_model_cost=dict(costs),
                       timer_scope='public+prediction load,one sharedbasis setup,6forward+adjoint endpoints,projections+serialization; summary write excluded')
        write(OUT/'ENDPOINTS.json',summary)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/'SOURCE_WIRING.json').exists()or(OUT/'ENDPOINTS.json').exists():raise FileExistsError('Preserve existing diagnostics; no rerun')
    with threadpool_limits(limits=1):source_wiring();endpoints()

if __name__=='__main__':main()
