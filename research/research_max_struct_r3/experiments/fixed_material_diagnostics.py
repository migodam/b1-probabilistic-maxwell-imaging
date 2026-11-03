"""H1/H2 fixed-material diagnostics on published R2 GT/endpoints; no inversion."""
from pathlib import Path
from collections import Counter
import json,sys,time,traceback
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel

CASE='obs_b20d5a55677d4eca';R2=REPO/'research/research_max_struct_r2'

def dump(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False))

def main():
    out=ROOT/'results/diagnostics';out.mkdir(parents=True,exist_ok=True)
    path=out/'FIXED_MATERIAL_H1_H2.json'
    if path.exists():raise FileExistsError('Preserve diagnostic disposition')
    start=time.perf_counter();cpu=time.process_time();models=[];aux=Counter();rows=[]
    receipt=dict(status='RUNNING',classification='OFFLINE_DIAGNOSTIC_PUBLISHED_R2_GT',case=CASE,
        inverse_runs=0,beta_current=1.,source_metrics='same fixed native wd/ws from observed Y and E',
        truth_access='R2 published m/Y_clean/noise only; no new confirmation GT',
        code='experiments/fixed_material_diagnostics.py',hash_checks=False)
    try:
        with np.load(R2/'data_public'/(CASE+'.npz'),allow_pickle=False)as f:
            public={k:f[k]for k in ('G','E0','S','Y','grid_shape','sigma')}
        with np.load(R2/'sealed_truth'/(CASE+'.npz'),allow_pickle=False)as f:
            gt=f['m'];noise=f['noise'];clean=f['Y_clean']
        anchors={'published_GT':gt};currents={}
        for arm in ('native','global'):
            with np.load(R2/'predictions'/CASE/(arm+'.npz'),allow_pickle=False)as f:anchors[arm]=f['m']
            with np.load(R2/'results/blind'/CASE/arm/'stage_2.npz',allow_pickle=False)as f:currents[arm]=f['J']
        Y=public['Y'];E=public['E0'];G=public['G'];S=public['S'];ns=E.shape[1]
        wd=1/np.sum(abs(Y)**2,axis=0);ws=1/np.sum(abs(E)**2,axis=0)
        Us,sigma,Vh=la.svd(S,full_matrices=False);aux['receiver_SVD']+=1
        rank=int(np.count_nonzero(sigma>1e-3*sigma[0]));D=Vh[:rank];UD=Us[:,:rank];sig=sigma[:rank]
        coeff=(UD.conj().T@Y)/sig[:,None]
        noise_coeff=(UD.conj().T@noise)/sig[:,None]
        SS=S.conj().T@S;SY=S.conj().T@Y;aux['receiver_Gram_builds']+=1
        for tag,m in anchors.items():
            tick=time.perf_counter();ct=time.process_time();model=MaxwellModel(public);models.append(model)
            st=model.state(m);Jp=st['j'];P=1/(1+m);L=st['L'];H=P[:,None]*L;bstate=(P*m)[:,None]*E
            with model._physics.operation('H1_DPCS_TRANSFER',rhs=rank,receiver_transpose_RHS=rank):
                DL=la.lu_solve(st['lu'],D.T,trans=1).T
            A=DL/P[None,:];Gram=A@A.conj().T
            fac=la.cho_factor(Gram,lower=True);aux['H1_source_Gram_Cholesky']+=1
            b=coeff-D@Jp;z=la.cho_solve(fac,b);aux['H1_source_Gram_RHS']+=ns
            floor=np.real(np.sum(b.conj()*z,axis=0));weighted_floor=.5*float(np.sum(ws*floor))
            comparisons={}
            for arm,ja in currents.items():
                defect=P[:,None]*(ja-m[:,None]*(E+G@ja));term=.5*float(np.sum(ws[None,:]*abs(defect)**2))
                comparisons[arm]=dict(weighted_state_term=term,ratio_to_exact_hard_DPCS_floor=term/max(weighted_floor,1e-300),
                    DPCS_constraint_defect=float(la.norm(D@ja-coeff)),
                    semantics='stored stage2 J transported to this anchor m; not original stage2 residual')
            HH=H.conj().T@H;Hb=H.conj().T@bstate;aux['state_Gram_builds']+=1
            Jfull=np.empty_like(E)
            for s in range(ns):
                K=wd[s]*SS+ws[s]*HH;rhs=wd[s]*SY[:,s]+ws[s]*Hb[:,s]
                f=la.cho_factor(K,lower=True);aux['H2_current_Cholesky']+=1
                Jfull[:,s]=la.cho_solve(f,rhs);aux['H2_current_RHS']+=1
            R=H@Jfull-bstate;Rd=S@Jfull-Y
            profile_data=.5*float(np.sum(wd[None,:]*abs(Rd)**2));profile_state=.5*float(np.sum(ws[None,:]*abs(R)**2))
            physical_data=.5*float(np.sum(wd[None,:]*abs(st['Y']-Y)**2));prior=model.reg.evaluate(m)[0]
            normgrad=S.conj().T@(Rd*wd)+H.conj().T@(R*ws)
            diff=Jp-Jfull
            row=dict(anchor=tag,status='COMPLETED',DPCS_rank=rank,H1_per_source_min_state_energy=floor.tolist(),
                H1_weighted_floor=weighted_floor,DPCS_current_mismatch_l2=float(la.norm(b)),
                DPCS_noise_only_mismatch_l2=float(la.norm(noise_coeff)),
                DPCS_nonnoise_material_mismatch_l2=float(la.norm(b-noise_coeff)),
                retained_whitened_measurement_noise_norm=float(la.norm(UD.conj().T@noise)/float(public['sigma'])),
                native_auxiliary_comparisons=comparisons,
                H2_same_metric_exact_profile_data=profile_data,H2_same_metric_exact_profile_state=profile_state,
                H2_same_metric_exact_profile_total=profile_data+profile_state+prior,
                physical_native_metric_data=physical_data,physical_native_metric_total=physical_data+prior,
                physical_original_normalized_Phi=float(.5*la.norm(st['Y']-Y)**2/la.norm(Y)**2+prior),
                physical_relative_residual=float(la.norm(st['Y']-Y)/la.norm(Y)),prior=float(prior),
                profile_to_physical_data_ratio=(profile_data+profile_state)/max(physical_data,1e-300),
                profile_current_normal_equation_gradient_norm=float(la.norm(normgrad)),
                profile_current_to_physical_difference=float(la.norm(diff)),
                profile_current_raw_state_defect=float(la.norm(R/P[:,None])),
                observed_correction_port_amplification=float(la.norm(S@diff)/max(la.norm(R/P[:,None]),1e-300)),
                H1_gram_solve_relative_residual=float(la.norm(Gram@z-b)/max(la.norm(b),1e-300)),
                model_cost=dict(model.cost),wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-ct)
            rows.append(row)
        receipt.update(status='COMPLETED',rows=rows,source_weights=dict(wd=wd.tolist(),ws=ws.tolist()),
            clean_Y_consistency=float(la.norm(Y-clean-noise)),
            interpretation='fixed-material mechanism diagnostic only; exact noisy hard locking conflicts do not establish their dominance in nonlinear failure')
    except Exception as exc:receipt.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc(),rows=rows);raise
    finally:
        cost=Counter()
        for m in models:cost.update(m.cost)
        receipt.update(physical_model_cost=dict(cost),auxiliary_cost=dict(aux),wall_seconds=time.perf_counter()-start,
                       cpu_seconds=time.process_time()-cpu,FFT=0,GPU=0,
                       timer_scope='public/publishedtruth load,SVD/Gram,3physicalanchors,H1/H2fixedlinearwork; ownreceipt write excluded')
        dump(path,receipt)

if __name__=='__main__':
    with threadpool_limits(limits=1):main()
