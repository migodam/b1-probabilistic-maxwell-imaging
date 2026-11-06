"""Public development diagnostic; primary confirmation truth is never read."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
from collections import Counter
import json,time,sys
import numpy as np
from scipy.linalg import svd
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r7.experiments.score_primary import material_metrics
from research.research_max_struct_r6.src.fft_tsom_ops import NativeFourierAPCS
from research.research_max_struct_r8.experiments.run_case import load_public

def main():
    start,cpu=time.perf_counter(),time.process_time();out=ROOT/'RESULTS/warm_development_analysis';out.mkdir(exist_ok=False)
    cases=json.loads((ROOT/'data_public/warm_development/cases.json').read_text())['records'];tasks=json.loads((ROOT/'configs/warm_development_manifest.json').read_text())['tasks']
    with np.load(ROOT/'data_training/validation.npz') as z:truth=z['m'][:8].reshape(8,576).copy()
    rows=[];projection=[]
    for i,case in enumerate(cases):
        ident=case['case_id'];amp=case['condition']['amplitude'];support=truth[i]>.1+.5*amp
        for t in [t for t in tasks if t['case_id']==ident]:
            d=REPO/t['output_dir'];rec=json.loads((d/'RECEIPT.json').read_text())
            with np.load(d/'prediction.npz') as z:m=z['m'].ravel()
            with np.load(d/'initial.npz') as z:mi=z['m'].ravel()
            met=material_metrics(m,truth[i],support,amp);im=material_metrics(mi,truth[i],support,amp)
            rows.append(dict(case_id=ident,arm=t['arm'],initialization=t.get('initialization','uniform'),cpu_seconds=rec['cpu_seconds'],Phi=rec['physical_objective'],relative_residual=rec['physical_relative_residual'],selected_stage=rec['auxiliary_selected_stage'],initial_nrmse=im['material_nrmse'],**met,corrector_stop=rec['corrector']['stop']))
        # Reconstruct the already specified initial affine projection only.
        d=ROOT/f'RESULTS/warm_development/{ident}__a_plus__nn_joint';pub=load_public(case['public_files'])
        with np.load(d/'initial.npz') as z:m=z['m'].ravel();jp=z['j'].copy();qp=z['Q'].copy()
        with np.load(d/'SOURCE_METRICS.npz') as z:w=z['state_weights'].copy()
        G,S,Y,E=pub['G'],pub['S'],pub['Y'],pub['E0'];u,s,vh=svd(S,full_matrices=False);rank=int(np.count_nonzero(s>.001*s[0]));vd=vh[:rank].conj().T;JD=vd@((u[:,:rank].conj().T@Y)/s[:rank,None]);cost=Counter();U=NativeFourierAPCS((24,24),3,vd,cost);UU=U.apply(np.eye(U.k));a=np.linalg.lstsq(UU,jp-JD,rcond=1e-12)[0]
        pdjp=vd@(vd.conj().T@jp);fixed=JD-pdjp;jfull=JD+jp-pdjp;jnative=JD+UU@a
        base=np.linalg.norm(jp)
        for name,j in [('mf3_projected',jnative),('full_APCS_projected',jfull)]:
            dj=j-jp;q=E+G@j;uconst=j-m[:,None]*q;lin=np.real(np.sum(w[None,:]*qp.conj()*(dj-m[:,None]*(G@dj)),axis=1))/np.sum(w[None,:]*abs(qp)**2,axis=1)
            readout=np.clip(np.real(np.sum(w[None,:]*q.conj()*j,axis=1))/np.sum(w[None,:]*abs(q)**2,axis=1),.02,20.)
            projection.append(dict(case_id=ident,projection=name,DPCS_rank=rank,relative_current_displacement=float(np.linalg.norm(dj)/base),irreducible_DPCS_displacement=float(np.linalg.norm(fixed)/base),complement_representation_displacement=float(np.linalg.norm(dj-fixed)/base),constitutive_defect_relative=float(np.linalg.norm(uconst)/max(np.linalg.norm(j),1e-30)),unregularized_material_readout_step_norm=float(np.linalg.norm(readout-m)),linearized_readout_step_norm=float(np.linalg.norm(lin)),linear_prediction_relative_defect=float(np.linalg.norm((readout-m)-lin)/max(np.linalg.norm(readout-m),1e-30)),scope='STATE_FIXED_DEVELOPMENT_DIAGNOSTIC; ordinary real LS without prior, not an A+ trajectory update; no truth in projection/readout',LU=0,complex_RHS=0))
        np.savez_compressed(out/(ident+'_projection.npz'),m=m,j_physical=jp,j_mf3=jnative,j_full=jfull,JD=JD,DPCS_displacement=fixed)
    (out/'SCORES.json').write_text(json.dumps(dict(rows=rows,scope='OFFLINE_DEVELOPMENT_VALIDATION_ONLY; these eight objects are not independent confirmation'),indent=2));(out/'CURRENT_PROJECTION_DIAGNOSTICS.json').write_text(json.dumps(projection,indent=2))
    (out/'RECEIPT.json').write_text(json.dumps(dict(status='DEVELOPMENT_ANALYSED',cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,LU=0,complex_RHS=0,receiver_SVD=8,Fourier_projection_LS_SVD=8,Green_column_actions=sum(4*c['condition']['sources'] for c in cases),truth_access='TRAINING_VALIDATION_LABELS_ONLY; primary test remains sealed'),indent=2))
    print(json.dumps(dict(rows=len(rows),projection_states=len(projection),output=str(out))))

if __name__=='__main__':main()
