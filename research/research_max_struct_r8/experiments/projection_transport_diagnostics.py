"""Public state-fixed diagnostics of hard-DPCS and Fourier transport.

Reuses paid NN initial physical states. No forward solve, GT, or parameter fit.
Results are covariates, never a rule for rejecting an initializer.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
import argparse,json,sys,time
import numpy as np
from scipy.linalg import svd
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r8.experiments.run_case import load_public
from research.research_max_struct_r6.src.fft_tsom_ops import NativeFourierAPCS

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cohort',choices=['development','confirmation'],required=True);a=ap.parse_args()
    start,cpu=time.perf_counter(),time.process_time();cache={};counter=Counter();rows=[]
    manifest='warm_development_manifest.json' if a.cohort=='development' else 'warm_all_manifest.json'
    tasks=json.loads((ROOT/'configs'/manifest).read_text())['tasks'];tasks=[t for t in tasks if t['arm']=='a_plus' and t['initialization']=='nn_joint']
    # Public availability qualification only; never reads seal truth metadata.
    for t in tasks:
        r=json.loads((REPO/t['output_dir']/'RECEIPT.json').read_text())
        if r['status']!='SEALED_PUBLIC_PREDICTION':raise RuntimeError('Paid state unavailable '+t['task_id'])
    for t in tasks:
        d=REPO/t['output_dir'];pub=load_public(t['public_files']);key=t['public_files'][0]
        if key not in cache:
            S=pub['S'];u,s,vh=svd(S,full_matrices=False);rank=int(np.count_nonzero(s>.001*s[0]));vd=vh[:rank].conj().T
            U=NativeFourierAPCS((24,24),3,vd,counter);UU=U.apply(np.eye(U.k));pinv=np.linalg.pinv(UU,rcond=1e-12)
            cache[key]=(u[:,:rank],s[:rank],vd,UU,pinv);counter['receiver_SVD']+=1;counter['Fourier_projection_SVD']+=1
        us,s,vd,UU,pinv=cache[key]
        with np.load(d/'initial.npz') as z:m=z['m'].ravel().copy();jp=z['j'].copy();qp=z['Q'].copy()
        G,S,Y,E=pub['G'],pub['S'],pub['Y'],pub['E0'];ns=jp.shape[1];jd=vd@((us.conj().T@Y)/s[:,None]);pd=vd@(vd.conj().T@jp)
        deltaD=jd-pd;deltaN=UU@(pinv@jp)-(jp-pd);zero=jd-jp;truncated=deltaD+deltaN
        norm=float(np.linalg.norm(jp));ws=1/np.sum(abs(E)**2,axis=0);den=np.sum(ws[None,:]*abs(qp)**2,axis=1)
        def act(v):
            counter['Green_columns']+=ns
            return v-m[:,None]*(G@v)
        bD,bN=act(deltaD),act(deltaN)
        record=dict(case_id=t['case_id'],DPCS_rank=len(s),sources=ns,physical_current_norm=norm,
            DPCS_displacement=float(np.linalg.norm(deltaD)/norm),Fourier_displacement=float(np.linalg.norm(deltaN)/norm),
            transported_DPCS_defect=float(np.linalg.norm(bD)/norm),transported_Fourier_defect=float(np.linalg.norm(bN)/norm),
            signed_constitutive_cross_term=float(2*np.vdot(bD,bN).real/norm**2),
            full_APCS_constitutive_defect=float(np.linalg.norm(bD)/norm),
            truncated_constitutive_defect=float(np.linalg.norm(bD+bN)/norm),
            removed_Fourier_constitutive_defect_is_helpful=bool(np.linalg.norm(bD)<np.linalg.norm(bD+bN)),
            fixed_norm='All current/constitutive quantities divided by paid initial physical current Frobenius norm; no GT scaling',
            qualification='State-fixed ordinary LS diagnostic, not actual A+ update or online success predictor')
        for name,dj in [('zero_APCS',zero),('truncated_physical_projection',truncated),('fixed_physical_APCS_anchor',deltaD)]:
            j=jp+dj;q=E+G@j;counter['Green_columns']+=ns
            read=np.clip(np.real(np.sum(ws[None,:]*q.conj()*j,axis=1))/np.sum(ws[None,:]*abs(q)**2,axis=1),.02,20)
            defect=j-m[:,None]*q;linear=np.real(np.sum(ws[None,:]*qp.conj()*defect,axis=1))/den
            step=read-m
            record[name]=dict(current_relative_displacement=float(np.linalg.norm(dj)/norm),constitutive_relative_defect=float(np.linalg.norm(defect)/norm),
                port_change_relative=float(np.linalg.norm(S@dj)/np.linalg.norm(Y)),readout_step_norm=float(np.linalg.norm(step)),
                root_linear_prediction_norm=float(np.linalg.norm(linear)),root_linear_relative_defect=float(np.linalg.norm(step-linear)/max(np.linalg.norm(step),1e-30)))
            counter['receiver_columns']+=ns
        rows.append(record)
    out=ROOT/'RESULTS'/('projection_transport_'+a.cohort);out.mkdir(exist_ok=False)
    (out/'DIAGNOSTICS.json').write_text(json.dumps(dict(rows=rows,truth_access='NONE',selection_effect='NONE; all frozen tasks retained',parameters_fitted=0),indent=2))
    rec=dict(status='PUBLIC_DIAGNOSTIC_COMPLETE',objects=len(rows),created_utc=datetime.now(timezone.utc).isoformat(),cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,LU=0,complex_RHS=0,operation_counts=dict(counter),truth_access='NONE',new_initial_forward=False)
    (out/'RECEIPT.json').write_text(json.dumps(rec,indent=2));print(json.dumps(rec))

if __name__=='__main__':
    with threadpool_limits(limits=1):main()
