"""Six existing failed endpoints: one exact gradient each, no optimization/GT.

New information sought: raw gradients do not establish box KKT stationarity.
This diagnostic does not modify any frozen prediction or confirmation score.
"""
from pathlib import Path
import json,sys,time
from collections import Counter
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel

def main():
    output=ROOT/'results/failed_endpoint_stationarity';output.mkdir(exist_ok=False)
    start,cpu=time.perf_counter(),time.process_time();rows=[];total=Counter()
    with np.load(ROOT/'data_public/confirmation/operators.npz',allow_pickle=False) as f:ops={k:f[k] for k in f.files}
    with threadpool_limits(limits=1):
        for case in ('obs_42cb233e4f224bb3','obs_63abafa750fc483c'):
            with np.load(ROOT/'data_public/confirmation'/(case+'.npz'),allow_pickle=False) as f:public=dict(ops,**{k:f[k] for k in f.files})
            for method in ('onehop','global','native_extended'):
                m=np.load(ROOT/'results/confirmation_remote/predictions'/case/(method+'.npz'),allow_pickle=False)['m']
                model=MaxwellModel(public);value,g=model.data_value_gradient(m,public['Y'])
                pg=m-np.clip(m-g,.02,20.)
                row=dict(case_id=case,method=method,physical_Phi=value,raw_gradient_inf=float(np.linalg.norm(g,np.inf)),
                    unit_step_projected_gradient_inf=float(np.linalg.norm(pg,np.inf)),
                    near_lower_cells=int(np.sum(m<=.02+1e-8)),near_upper_cells=int(np.sum(m>=20.-1e-8)),
                    physical_cost=dict(model.cost),state_relative_residual=float(model.state(m)['state_relative_residual']))
                rows.append(row);total.update(model.cost)
                np.savez_compressed(output/(case+'__'+method+'.npz'),m=m,gradient=g,projected_gradient=pg)
    record=dict(status='COMPLETED_DIAGNOSTIC_NOT_NEW_PREDICTIONS',truth_access='NONE',
        reason='Raw gradient saved previously is not a box KKT measure; evaluate existing failed endpoints once',
        rows=rows,physical_cost=dict(total),wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,
        objective='same normalized physical Phi and prior',new_optimization_iterations=0,
        caveat='Numerical first-order measure only, not interval-certified KKT/local-minimum or global-uniqueness proof')
    (output/'RECEIPT.json').write_text(json.dumps(record,indent=2));print(json.dumps(record))
if __name__=='__main__':main()
