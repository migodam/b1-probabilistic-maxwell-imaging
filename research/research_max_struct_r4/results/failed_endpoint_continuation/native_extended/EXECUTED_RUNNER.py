"""Fixed post-confirmation budget diagnostic, not a modified primary outcome.

One public-failure-selected object, two endpoints, one extra4000 allowance each.
No new auxiliary method, tuning, truth access, or retry after this allowance.
"""
from pathlib import Path
import argparse,json,sys,time,shutil
from collections import Counter
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r4.experiments.cooperative_paths import refine

def main(method):
    case='obs_42cb233e4f224bb3'
    output=ROOT/'results/failed_endpoint_continuation'/method;output.mkdir(parents=True,exist_ok=False)
    config=dict(classification='POST_CONFIRMATION_BUDGET_DIAGNOSTIC',case_id=case,method=method,
        selection='lexicographically first public-data-ineligible case; no material metric selection',
        extra_corrector_maxiter=4000,original_corrector_maxiter=1000,
        scale=1e6,original_gtol=1e-8,ftol=1e-12,material_box=[.02,20.],prior='unchanged',
        scope='Test whether additional original-Maxwell work resolves the failure within a fixed allowance; not a new algorithm')
    (output/'CONFIG.json').write_text(json.dumps(config,indent=2));shutil.copy2(__file__,output/'EXECUTED_RUNNER.py')
    start,cpu=time.perf_counter(),time.process_time();model=None
    receipt=dict(status='RUNNING',configuration=config,truth_access='NONE',failed_attempts=[])
    (output/'RECEIPT.json').write_text(json.dumps(receipt))
    try:
        with np.load(ROOT/'data_public/confirmation/operators.npz',allow_pickle=False) as f:public={k:f[k] for k in f.files}
        with np.load(ROOT/'data_public/confirmation'/(case+'.npz'),allow_pickle=False) as f:public.update({k:f[k] for k in f.files})
        m=np.load(ROOT/'results/confirmation_remote/predictions'/case/(method+'.npz'),allow_pickle=False)['m']
        with threadpool_limits(limits=1):
            model=MaxwellModel(public);value,g=model.data_value_gradient(m,public['Y'])
            receipt['initial']=dict(physical_Phi=value,projected_gradient_inf=float(np.linalg.norm(m-np.clip(m-g,.02,20.),np.inf)))
            final,corr=refine(model,public['Y'],m,maximum=4000)
            history=np.array(corr.pop('history'));value,g=model.data_value_gradient(final,public['Y']);state=model.state(final)
            np.savez_compressed(output/'prediction.npz',m=final,Y_prediction=state['Y'],gradient=g,corrector_history=history)
            receipt.update(status='SEALED_POST_CONFIRMATION_DIAGNOSTIC',corrector=corr,physical_objective=value,
                physical_relative_residual=float(np.linalg.norm(state['Y']-public['Y'])/np.linalg.norm(public['Y'])),
                physical_whitened_residual=float(np.linalg.norm(state['Y']-public['Y'])/float(public['sigma'])),
                projected_gradient_inf=float(np.linalg.norm(final-np.clip(final-g,.02,20.),np.inf)),
                state_relative_residual=float(state['state_relative_residual']),selection='Original Phi, prior unchanged, retain incumbent',
                original_confirmation_predictions_unchanged=True)
    except Exception as exc:
        receipt.update(status='FAILED',error=repr(exc));receipt['failed_attempts'].append(repr(exc));raise
    finally:
        receipt.update(physical_cost=dict(model.cost) if model else {},wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu)
        (output/'RECEIPT.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt),flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--method',choices=['onehop','native_extended'],required=True)
    main(parser.parse_args().method)
