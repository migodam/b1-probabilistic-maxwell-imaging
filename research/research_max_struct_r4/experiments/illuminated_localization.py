"""Formula-driven T-completion refinement using only its illuminated action."""
from pathlib import Path
import sys,time,json,traceback
import numpy as np
from scipy import linalg as la
from collections import Counter
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_fbg_r3_experiments.src.cie_tsom_baseline import _readout
from research.research_max_struct_r2.experiments.run_blind import dump,ALLOW,metric
from research.research_max_struct_r4.src.tmatrix_completion_control import TMatrixCompletionControl
from research.research_max_struct_r4.experiments.cooperative_paths import refine

def run():
    tick,cpu=time.perf_counter(),time.process_time()
    out=ROOT/'results/illuminated_localization';out.mkdir(parents=True,exist_ok=True)
    path=out/'RECEIPT.json'
    if path.exists():raise FileExistsError('Preserve result')
    model=None;aux=Counter()
    receipt=dict(run_id='illuminated_localization',status='RUNNING',stage='DEV_FORMULA_REFINEMENT',truth_access='NONE',
                 case_id='obs_b20d5a55677d4eca',iterations=30,product_cutoff=1e-3,
                 material_projection='classical shared-real CSI row fitting + same public prior; no beta',
                 rationale='minimize illuminated (diag(m)-V)Q rather than Frobenius diag(V); eliminate unused T columns',
                 source_weights='fixed original paper1/||Einc_source||²',
                 corrector='same scaled original corrector1000 as DCTMC control',selection='original Phi only',
                 no_new_method_originality_claim=True)
    dump(path,receipt)
    try:
        with np.load(REPO/f"research/research_max_struct_r2/data_public/{receipt['case_id']}.npz") as f:
            public={k:f[k] for k in f.files if k in ALLOW}
        model=MaxwellModel(public);Y=public['Y'];m=np.full(model.n,float(public['host']))
        projection=TMatrixCompletionControl(model,Y,product_rtol=1e-3)
        weights=1/np.sum(abs(model.E)**2,axis=0)
        bounds=(np.full(model.n,.02),np.full(model.n,20.))
        history=[];materials=[m.copy()];scored=[metric(model,Y,m)]
        for it in range(30):
            state=model.state(m)
            with model._physics.operation('ILLUMINATED_COMPLETION',data_current_projections=1,Green_actions=model.ns):
                encoded=projection.Ua.conj().T@(Y-state['Y'])@projection.Vb
                delta=projection.Va@((projection.mask*encoded)/projection.sa[:,None])@projection.Vb.conj().T
                J=state['j']+delta;Q=model.E+model.G@J
            m,inner=_readout(J,Q,m,None,weights,1.,bounds,lambda z:model.reg.evaluate(z)[:2],True,{},aux)
            row=metric(model,Y,m);scored.append(row);materials.append(m.copy())
            history.append(dict(iteration=it+1,Phi=row['objective'],physical_residual=row['relative_residual'],
                data_compatible_current_residual=float(la.norm(model.S@J-Y)/la.norm(Y)),
                constitutive_defect=float(la.norm(J-m[:,None]*Q)/la.norm(J)),
                current_correction_relative=float(la.norm(delta)/la.norm(state['j'])),**inner))
        best=min(scored,key=lambda x:x['objective'])
        np.savez_compressed(out/'outer_candidates.npz',m=np.asarray(materials))
        dump(out/'HISTORY.json',history)
        final,corr=refine(model,Y,best['m']);dump(out/'CORRECTOR.json',corr)
        selected=min([best,metric(model,Y,final)],key=lambda x:x['objective']);st=model.state(selected['m'])
        np.savez_compressed(out/'prediction.npz',m=selected['m'],Y_prediction=st['Y'])
        receipt.update(status='SEALED_DEVELOPMENT_PREDICTION',outer_selected_Phi=best['objective'],
            outer_selected_residual=best['relative_residual'],outer_selected_iteration=int(np.argmin([x['objective'] for x in scored])),
            final_Phi=selected['objective'],physical_residual=selected['relative_residual'],corrector_stop=corr['stop'],
            retained_complex_coordinates=int(projection.mask.sum()),discarded_data_relative=projection.data_info['discarded_observed_data_relative_norm'])
    except Exception as e:
        receipt.update(status='FAILED',error=repr(e),traceback=traceback.format_exc());raise
    finally:
        receipt.update(physical_cost=dict(model.cost) if model else {},auxiliary_cost=dict(aux),
                       wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-cpu)
        dump(path,receipt)
    print(json.dumps(receipt))
if __name__=='__main__':
    with threadpool_limits(limits=1):run()
